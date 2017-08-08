package org.broadinstitute.hellbender.tools.spark.sv.playground;

import avro.shaded.com.google.common.annotations.VisibleForTesting;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKsvVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;


final class ForSimpleStrandSwitch implements VariantDetectorFromLongReadAlignments {

    @SuppressWarnings("unchecked")
    private static final List<String> EMPTY_INSERTION_MAPPINGS = Collections.EMPTY_LIST;

    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    public void inferSvAndWriteVCF(final JavaRDD<AlignedContig> longReads, final String vcfOutputFileName,
                                   final Broadcast<ReferenceMultiSource> broadcastReference, final String fastaReference,
                                   final GCSOptions options, final Logger toolLogger) {

        longReads.cache();
        toolLogger.info(longReads.count() + " chimera indicating either 1) simple strand-switch breakpoints, or 2) inverted duplication.");

        // split between suspected inv dup VS strand-switch breakpoint
        final JavaRDD<VariantContext> simpleStrandSwitchBkpts =
                dealWithSimpleStrandSwitchBkpts(longReads.filter(read -> !isLikelyInvertedDuplication(read)), broadcastReference, toolLogger);
        SVVCFWriter.writeVCF(options, vcfOutputFileName.replace(".vcf", "_simpleSS.vcf"), fastaReference, simpleStrandSwitchBkpts, toolLogger);

        final JavaRDD<VariantContext> invDups =
                dealWithSuspectedInvDup(longReads.filter(ForSimpleStrandSwitch::isLikelyInvertedDuplication), broadcastReference, toolLogger);
        if (invDups != null)
            SVVCFWriter.writeVCF(options, vcfOutputFileName.replace(".vcf", "_invDup.vcf"), fastaReference, simpleStrandSwitchBkpts.union(invDups), toolLogger);

        longReads.unpersist();
    }

    /**
     * @return true iff the overlap of the two AI's made up the CA is strictly greater than half of any of the two AI's ref span.
     */
    @VisibleForTesting
    static boolean isLikelyInvertedDuplication(final AlignedContig longRead) {
        final int[] overlaps = computeOverlaps(longRead);
        return 2 * overlaps[0] > Math.min(longRead.alignmentIntervals.get(0).referenceSpan.size(),
                                          longRead.alignmentIntervals.get(1).referenceSpan.size());
    }

    /**
     * Given a read with only two alignments mapped to the same chromosome, return an array of size 2, with elements
     *   [0] : overlap length on reference
     *   [1] : overlap length on read
     */
    private static int[] computeOverlaps(final AlignedContig longRead) {
        final SimpleInterval refSpanOne = longRead.alignmentIntervals.get(0).referenceSpan,
                             refSpanTwo = longRead.alignmentIntervals.get(1).referenceSpan;

        // dummy number for chr to be used in constructing SVInterval, since input CA has 2 AI & both map to the same chr
        final int dummyChr = 1;
        final SVInterval intOne = new SVInterval(dummyChr, refSpanOne.getStart(), refSpanOne.getEnd() + 1),
                         intTwo = new SVInterval(dummyChr, refSpanTwo.getStart(), refSpanTwo.getEnd() + 1);

        final int overlapOnRef = intOne.overlapLen(intTwo);
        final int overlapOnRead = Math.max(0, 1 + longRead.alignmentIntervals.get(0).endInAssembledContig
                                                - longRead.alignmentIntervals.get(1).startInAssembledContig);

        return new int[]{overlapOnRef, overlapOnRead};
    }

    // =================================================================================================================
    private JavaRDD<VariantContext> dealWithSimpleStrandSwitchBkpts(final JavaRDD<AlignedContig> longReads,
                                                                    final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                    final Logger toolLogger) {
        final JavaPairRDD<ChimericAlignment, byte[]> inversionOrInvertedInsertion =
                longReads
                        .mapToPair(ForSimpleStrandSwitch::convertAlignmentIntervalToChimericAlignment)
                        .filter(Objects::nonNull).cache();

        toolLogger.info(inversionOrInvertedInsertion.count() + " chimera indicating simple strand-switch breakpoints.");

        return inversionOrInvertedInsertion
                        .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2), pair._1))
                        .groupByKey()
                        .mapToPair(noveltyAndEvidence -> inferType(noveltyAndEvidence, broadcastReference.getValue()))
                        .flatMap(noveltyTypeAndEvidence ->
                                AnnotatedVariantProducer
                                        .produceMultipleAnnotatedVcFromNovelAdjacency(noveltyTypeAndEvidence._1,
                                        noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2, broadcastReference));
    }

    /**
     * Taking advantage of the fact that for input read, we know it has only two alignments that map to the same reference
     * chromosome, with strand switch.
     * @return  null if the pair of the alignments are no strong enough to support a strand switch breakpoint
     *                  {@link #splitPairStrongEnoughEvidenceForCA(AlignmentInterval, AlignmentInterval, int, int)},
     *          otherwise a pair {read sequence, chimeric alignment}
     */
    private static Tuple2<ChimericAlignment, byte[]> convertAlignmentIntervalToChimericAlignment
    (final AlignedContig longReadWith2AIMappedToSameChrAndStrandSwitch) {

        final AlignmentInterval intervalOne = longReadWith2AIMappedToSameChrAndStrandSwitch.alignmentIntervals.get(0),
                                intervalTwo = longReadWith2AIMappedToSameChrAndStrandSwitch.alignmentIntervals.get(1);

        if (splitPairStrongEnoughEvidenceForCA(intervalOne, intervalTwo, MORE_RELAXED_ALIGNMENT_MIN_MQ,  MORE_RELAXED_ALIGNMENT_MIN_LENGTH)) {
            return new Tuple2<>(new ChimericAlignment(intervalOne, intervalTwo, EMPTY_INSERTION_MAPPINGS,
                    longReadWith2AIMappedToSameChrAndStrandSwitch.contigName), longReadWith2AIMappedToSameChrAndStrandSwitch.contigSequence);
        } else {
            return null;
        }
    }

    /**
     * Roughly similar to {@link ChimericAlignment#nextAlignmentMayBeNovelInsertion(AlignmentInterval, AlignmentInterval, Integer)}:
     *  1) either alignment may have very low mapping quality (a more relaxed mapping quality threshold);
     *  2) either alignment may consume only a "short" part of the contig, or if assuming that the alignment consumes
     *     roughly the same amount of ref bases and read bases, has isAlignment that is too short
     */
    private static boolean splitPairStrongEnoughEvidenceForCA(final AlignmentInterval intervalOne,
                                                              final AlignmentInterval intervalTwo,
                                                              final int mapQThresholdInclusive,
                                                              final int alignmentLengthThresholdInclusive) {

        if (intervalOne.mapQual < mapQThresholdInclusive || intervalTwo.mapQual < mapQThresholdInclusive)
            return false;

        final int overlap = AlignmentInterval.overlapOnContig(intervalOne, intervalTwo);

        final int x = intervalOne.endInAssembledContig - intervalOne.startInAssembledContig + 1,
                  y = intervalTwo.endInAssembledContig - intervalTwo.startInAssembledContig + 1;

        return Math.min(x - overlap, y - overlap) >= alignmentLengthThresholdInclusive;
    }

    static final class INV55BND extends BreakEndVariantType {

        /**
         * Technically, a strand switch breakpoint should have two VCF records, hence we also have {@code forUpstreamLoc}.
         */
        INV55BND(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc,
                 final ReferenceMultiSource reference) {
            super(getIDString(narl, forUpstreamLoc),
                    Collections.singletonMap(GATKsvVCFConstants.INV55, ""),
                    extractBasesForAltAllele(narl, forUpstreamLoc, reference), true,
                    getTheOtherRefLoc(narl, forUpstreamLoc), true);
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc) {

            return GATKsvVCFConstants.BREAKEND_STR + GATKsvVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    GATKsvVCFConstants.INV55 + GATKsvVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedLeftRefLoc.getContig() + GATKsvVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedLeftRefLoc.getEnd() + GATKsvVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedRightRefLoc.getStart() + GATKsvVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    (forUpstreamLoc ? "1" : "2") ;
        }

        private static byte[] extractBasesForAltAllele(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc,
                                                       final ReferenceMultiSource reference) {
            try {
                final byte[] ref = reference
                        .getReferenceBases(null, forUpstreamLoc ? narl.leftJustifiedLeftRefLoc :
                                                                                narl.leftJustifiedRightRefLoc)
                        .getBases();
                final String ins = narl.complication.getInsertedSequenceForwardStrandRep();
                if (ins.isEmpty()) {
                    return ref;
                } else {
                    return forUpstreamLoc ? ArrayUtils.addAll(ref, ins.getBytes())
                                          : ArrayUtils.addAll(ref, SequenceUtil.reverseComplement(ins).getBytes());
                }
            } catch (final IOException ioex) {
                throw new GATKException("Could not read reference for extracting reference bases.", ioex);
            }
        }

        private static SimpleInterval getTheOtherRefLoc(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc) {
            return forUpstreamLoc ? narl.leftJustifiedRightRefLoc : narl.leftJustifiedLeftRefLoc;
        }
    }

    static final class INV33BND extends BreakEndVariantType {

        /**
         * Technically, a strand switch breakpoint should have two VCF records, hence we also have {@code forUpstreamLoc}.
         */
        INV33BND(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc,
                 final ReferenceMultiSource reference) {
            super(getIDString(narl, forUpstreamLoc),
                    Collections.singletonMap(GATKsvVCFConstants.INV33, ""),
                    extractBasesForAltAllele(narl, forUpstreamLoc, reference), false,
                    getTheOtherRefLoc(narl, forUpstreamLoc), false);
        }

        private static String getIDString(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc) {

            return GATKsvVCFConstants.BREAKEND_STR + GATKsvVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    GATKsvVCFConstants.INV33 + GATKsvVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedLeftRefLoc.getContig() + GATKsvVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedLeftRefLoc.getEnd() + GATKsvVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    narl.leftJustifiedRightRefLoc.getStart() + GATKsvVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                    (forUpstreamLoc ? "1" : "2") ;
        }

        private static byte[] extractBasesForAltAllele(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc,
                                                       final ReferenceMultiSource reference) {
            try {
                final byte[] ref = reference
                        .getReferenceBases(null, forUpstreamLoc ? narl.leftJustifiedLeftRefLoc :
                                                                                narl.leftJustifiedRightRefLoc)
                        .getBases();
                final String ins = narl.complication.getInsertedSequenceForwardStrandRep();
                if (ins.isEmpty()) {
                    return ref;
                } else {
                    return forUpstreamLoc ? ArrayUtils.addAll(ref, ins.getBytes())
                                          : ArrayUtils.addAll(ref, SequenceUtil.reverseComplement(ins).getBytes());
                }
            } catch (final IOException ioex) {
                throw new GATKException("Could not read reference for extracting reference bases.", ioex);
            }
        }

        private static SimpleInterval getTheOtherRefLoc(final NovelAdjacencyReferenceLocations narl, final boolean forUpstreamLoc) {
            return forUpstreamLoc ? narl.leftJustifiedRightRefLoc : narl.leftJustifiedLeftRefLoc;
        }
    }

    private static Tuple2<NovelAdjacencyReferenceLocations, Tuple2<Iterable<SvType>, Iterable<ChimericAlignment>>>
    inferType(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> noveltyAndEvidence,
              final ReferenceMultiSource reference) {

        final NovelAdjacencyReferenceLocations novelAdjacency = noveltyAndEvidence._1;
        final Iterable<ChimericAlignment> chimericAlignments = noveltyAndEvidence._2;
        final BreakEndVariantType bkpt_1, bkpt_2;
        if (novelAdjacency.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE) {
            bkpt_1 = new INV55BND(novelAdjacency, true, reference);
            bkpt_2 = new INV55BND(novelAdjacency, false, reference);
        } else if (novelAdjacency.strandSwitch == StrandSwitch.REVERSE_TO_FORWARD){
            bkpt_1 = new INV33BND(novelAdjacency, true, reference);
            bkpt_2 = new INV33BND(novelAdjacency, false, reference);
        } else {
            throw new GATKException("Wrong type of novel adjacency sent to wrong analysis pathway: no strand-switch being sent to strand-switch path. \n" +
                    Utils.stream(chimericAlignments).map(ChimericAlignment::onErrStringRep).collect(Collectors.toList()));
        }

        return new Tuple2<>(novelAdjacency, new Tuple2<>(Arrays.asList(bkpt_1, bkpt_2), chimericAlignments));
    }






    // =================================================================================================================
    private JavaRDD<VariantContext> dealWithSuspectedInvDup(final JavaRDD<AlignedContig> longReads,
                                                            final Broadcast<ReferenceMultiSource> broadcastReference,
                                                            final Logger toolLogger) {
        final JavaRDD<AlignedContig> invDupSuspects =
                longReads
                        .filter(ForSimpleStrandSwitch::isLikelyInvertedDuplication).cache();

        toolLogger.info(invDupSuspects.count() + " chimera indicating inverted duplication");


        final JavaPairRDD<ChimericAlignment, byte[]> noHomologyOnRead =
                invDupSuspects
                        .filter(lr -> !hasHomology(lr))
                        .mapToPair(ForSimpleStrandSwitch::convertAlignmentIntervalToChimericAlignment)
                        .filter(Objects::nonNull).cache();

        final JavaRDD<AlignedContig> withHomologyOnRead = invDupSuspects.filter(ForSimpleStrandSwitch::hasHomology);

        return null;
    }

    private static boolean hasHomology(final AlignedContig longRead) {
        return longRead.alignmentIntervals.get(0).endInAssembledContig < longRead.alignmentIntervals.get(1).startInAssembledContig;
    }

    private static void forInv55() {

    }

    /**
     * Similar to {@link ChimericAlignment#involvesRefPositionSwitch(AlignmentInterval, AlignmentInterval)} except "equals" case.
     * See which one is more appropriate.
     */
    static boolean involvesRefPositionSwitch(final AlignmentInterval regionWithLowerCoordOnContig,
                                             final AlignmentInterval regionWithHigherCoordOnContig) {

        return regionWithHigherCoordOnContig.referenceSpan.getStart() <= regionWithLowerCoordOnContig.referenceSpan.getStart();
    }
}
