package org.broadinstitute.hellbender.tools.spark.sv.playground;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;


/**
 * This tool takes a SAM file containing alignments of single-ended long read
 * (be it long read sequencing, or contigs assembled from standard Illumina short reads),
 * searches for split alignments indicating the presence of structural variations,
 * and outputs an interpreted annotated single sample VCF.
 */
@CommandLineProgramProperties(summary="Parses a SAM file containing long reads alignments, and outputs an interpreted annotated single sample VCF.",
        oneLineSummary="Parses a long read SAM file, and outputs single sample VCF.",
        usageExample = "InternalDiscoverVariantsFromFilteredContigAlignmentsSAMSpark \\" +
                "-I /path/to/my/dir/longReads.sam -O /path/to/my/dir/structuralVariants.vcf \\" +
                "-R /path/to/my/reference/reference.2bit --fastaReference /path/to/my/reference/reference.fasta",
        programGroup = StructuralVariationSparkProgramGroup.class)
@Deprecated
public final class InternalDiscoverVariantsFromFilteredContigAlignmentsSAMSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(InternalDiscoverVariantsFromFilteredContigAlignmentsSAMSpark.class);

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
            discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "vcf file for interpreted variants", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String vcfOutputFileName;

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

        // filter alignments and split the gaps
        final JavaRDD<AlignedContig> contigsWithAlignmentsReconstructed
                = InternalFilterLongReadAlignmentsSAMSpark.filterByScore(getReads(), getHeaderForReads(), localLogger);

        new ForSimpleInsDel()
                .inferSvAndWriteVCF(contigsWithAlignmentsReconstructed, vcfOutputFileName, ctx.broadcast(getReference()),
                        discoverStageArgs.fastaReference, getAuthenticatedGCSOptions(), localLogger);
    }
}
