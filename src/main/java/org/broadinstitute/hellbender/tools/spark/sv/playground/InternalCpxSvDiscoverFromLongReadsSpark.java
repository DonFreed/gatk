package org.broadinstitute.hellbender.tools.spark.sv.playground;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.EnumMap;
import java.util.HashSet;
import java.util.Set;


/**
 * This tool takes a SAM file containing alignments of single-ended long read
 * (be it long read sequencing, or contigs assembled from standard Illumina short reads),
 * searches for reads with complicated split alignments indicating the presence of complex structural variations (cxSV),
 * and outputs a custom file format containing how affected reference segments are rearranged on the provided sample
 * where the long reads are from.
 */
@CommandLineProgramProperties(summary="Parses a SAM file containing long reads alignments, and outputs cxSV rearrangements.",
        oneLineSummary="Parses a long read SAM file, and outputs cxSV rearrangements.",
        usageExample = "InternalCpxSvDiscoverFromLongReadsSpark \\" +
                "-I /path/to/my/dir/longReads.sam \\" +
                "-O /path/to/my/dir/outputDir \\" +
                "-R /path/to/my/reference/reference.2bit --fastaReference /path/to/my/reference/reference.fasta",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class InternalCpxSvDiscoverFromLongReadsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(InternalCpxSvDiscoverFromLongReadsSpark.class);

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
            discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "output directory",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputDir;

    @Argument(doc = "output SAM files",
            shortName = "wSAM",
            fullName = "writeSAM")
    private boolean writeSAMFiles;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public ReadFilter makeReadFilter() {
        return ReadFilterLibrary.MAPPED;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> reads = getReads();
        final SAMFileHeader header = getHeaderForReads();

        // filter alignments and split the gaps
        final JavaRDD<AlignedContig> contigsWithAlignmentsReconstructed =
                InternalFilterLongReadAlignmentsSAMSpark.filterByScore(reads, header, localLogger)
                        .filter(lr -> lr.alignmentIntervals.size()>1).cache();

        // split up the long reads by their possible type of SV
        final EnumMap<RawTypes, JavaRDD<AlignedContig>> splitUpLongReads =
                splitReads(contigsWithAlignmentsReconstructed, localLogger);

        if ( !FileUtils.createDirInBucketToWriteTo(outputDir) )
            throw new GATKException("Could not create directory " + outputDir + " to write results to.");

        if (writeSAMFiles) {
            final Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(header);
            splitUpLongReads.forEach((k, v) -> writeSAM(v, k.name(), reads, headerBroadcast, outputDir, localLogger));
        }

        new InternalVariantDetectorFromLongReadAlignmentsForSimpleInsDel()
                .inferSvAndWriteVCF(splitUpLongReads.get(RawTypes.InsDel), outputDir+"/"+RawTypes.InsDel.name()+".vcf",
                        ctx.broadcast(getReference()), discoverStageArgs.fastaReference, getAuthenticatedGCSOptions(),
                        localLogger);

        new InternalVariantDetectorFromLongReadAlignmentsForSimpleStrandSwitch()
                .inferSvAndWriteVCF(splitUpLongReads.get(RawTypes.Inv), outputDir+"/"+RawTypes.Inv.name()+".vcf",
                        ctx.broadcast(getReference()), discoverStageArgs.fastaReference, getAuthenticatedGCSOptions(),
                        localLogger);
    }

    enum RawTypes {
        MultiConf, Inv, InsDel, SegDupOrMEI, Cpx;
    }

    private static boolean hasEquallyGoodConfigurations(final AlignedContig longread) {
        return longread.hasEquallyGoodConfigurations;
    }

    private static boolean hasOnly2Alignments(final AlignedContig longreadWithOnlyOneConfig) {
        return longreadWithOnlyOneConfig.alignmentIntervals.size() == 2;
    }

    private static boolean isSameChromosomeMapping(final AlignedContig longreadWithOnlyOneConfigAnd2Aln) {
        return longreadWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(0).referenceSpan.getContig()
                .equals(longreadWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(1).referenceSpan.getContig());
    }

    private static boolean isLikelyInvBreakpointOrInsInv(final AlignedContig longreadWithOnlyOneConfigAnd2AlnToSameChr) {
        final AlignmentInterval intervalOne = longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0),
                                intervalTwo = longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(1);
        return intervalOne.forwardStrand ^ intervalTwo.forwardStrand;
    }

    private static boolean isLikelySimpleInsDel(final AlignedContig longreadWithOnlyOneConfigAnd2AlnToSameChr) {
        if (isLikelyInvBreakpointOrInsInv(longreadWithOnlyOneConfigAnd2AlnToSameChr)) return false;

        final AlignmentInterval intervalOne, intervalTwo;
        if (longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0).forwardStrand) {
            intervalOne = longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0);
            intervalTwo = longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(1);
        } else {
            intervalOne = longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(1);
            intervalTwo = longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0);
        }
        return ( intervalOne.referenceSpan.getEnd() <= intervalTwo.referenceSpan.getStart() )
                ||
                ( intervalOne.referenceSpan.getStart() < intervalTwo.referenceSpan.getStart() &&
                        intervalOne.referenceSpan.getEnd() < intervalTwo.referenceSpan.getEnd());
    }

    private static boolean isLikelySegDupOrMEI(final AlignedContig longreadWithOnlyOneConfigAnd2AlnToSameChr) {
        if (isLikelyInvBreakpointOrInsInv(longreadWithOnlyOneConfigAnd2AlnToSameChr)) return false;

        final AlignmentInterval intervalOne, intervalTwo;
        if (longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0).forwardStrand) {
            intervalOne = longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0);
            intervalTwo = longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(1);
        } else {
            intervalOne = longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(1);
            intervalTwo = longreadWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0);
        }
        return intervalOne.referenceSpan.getEnd() > intervalTwo.referenceSpan.getStart() &&
                intervalOne.referenceSpan.getStart() > intervalTwo.referenceSpan.getEnd();
    }

    private static boolean isLikelyCpx(final AlignedContig longreadWithOnlyOneConfig) {

        return !hasOnly2Alignments(longreadWithOnlyOneConfig) || !isSameChromosomeMapping(longreadWithOnlyOneConfig) ||
                !(isLikelyInvBreakpointOrInsInv(longreadWithOnlyOneConfig) || isLikelySegDupOrMEI(longreadWithOnlyOneConfig) ||
                        isLikelySimpleInsDel(longreadWithOnlyOneConfig));
    }

    private static EnumMap<RawTypes, JavaRDD<AlignedContig>> splitReads(final JavaRDD<AlignedContig> longreadsWithAlignmentsReconstructed,
                                                                        final Logger toolLogger) {

        final EnumMap<RawTypes, JavaRDD<AlignedContig>> longreadsByRawTypes = new EnumMap<>(RawTypes.class);

        // long reads with more than 1 best configurations
        longreadsByRawTypes.put(RawTypes.MultiConf,
                longreadsWithAlignmentsReconstructed.filter(InternalCpxSvDiscoverFromLongReadsSpark::hasEquallyGoodConfigurations));


        final JavaRDD<AlignedContig> longreadsWithOnlyOneBestConfig =
                longreadsWithAlignmentsReconstructed.filter(lr -> !hasEquallyGoodConfigurations(lr)).cache();


        // long reads with only 1 best configuration and having only 2 alignments mapped to the same chromosome (after gap split)
        final JavaRDD<AlignedContig> longreadsWithOnlyOneBestConfigAnd2AIToSameChr =
                longreadsWithOnlyOneBestConfig
                        .filter(InternalCpxSvDiscoverFromLongReadsSpark::hasOnly2Alignments)
                        .filter(InternalCpxSvDiscoverFromLongReadsSpark::isSameChromosomeMapping).cache();

        longreadsByRawTypes.put(RawTypes.Inv,
                longreadsWithOnlyOneBestConfigAnd2AIToSameChr.filter(InternalCpxSvDiscoverFromLongReadsSpark::isLikelyInvBreakpointOrInsInv));

        longreadsByRawTypes.put(RawTypes.InsDel,
                longreadsWithOnlyOneBestConfigAnd2AIToSameChr.filter(InternalCpxSvDiscoverFromLongReadsSpark::isLikelySimpleInsDel));

        longreadsByRawTypes.put(RawTypes.SegDupOrMEI,
                longreadsWithOnlyOneBestConfigAnd2AIToSameChr.filter(InternalCpxSvDiscoverFromLongReadsSpark::isLikelySegDupOrMEI));

        // long reads with only 1 best configuration
        longreadsWithOnlyOneBestConfigAnd2AIToSameChr.unpersist();
        longreadsByRawTypes.put(RawTypes.Cpx,
                longreadsWithOnlyOneBestConfig.filter(InternalCpxSvDiscoverFromLongReadsSpark::isLikelyCpx));
        longreadsWithOnlyOneBestConfig.unpersist();

        return longreadsByRawTypes;
    }

    private static void writeSAM(final JavaRDD<AlignedContig> filteredReads, final String rawTypeString,
                                 final JavaRDD<GATKRead> originalLongReads, final Broadcast<SAMFileHeader> headerBroadcast,
                                 final String outputDir, final Logger toolLogger) {

        final Set<String> filteredReadNames = new HashSet<>( filteredReads.map(tig -> tig.contigName).distinct().collect() );
        toolLogger.info(filteredReadNames.size() + " long reads indicating " + rawTypeString);
        final JavaRDD<SAMRecord> splitLongReads = originalLongReads.filter(read -> filteredReadNames.contains(read.getName()))
                .map(read -> read.convertToSAMRecord(headerBroadcast.getValue()));
        FileUtils.writeSAMFile(splitLongReads.collect().iterator(), headerBroadcast.getValue(),
                outputDir+"/"+rawTypeString+".sam", false);
    }
}
