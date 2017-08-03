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


final class InternalVariantDetectorFromLongReadAlignmentsForSimpleStrandSwitch implements InternalVariantDetectorFromLongReadAlignments {

    @SuppressWarnings("unchecked")
    private static final List<String> EMPTY_INSERTION_MAPPINGS = Collections.EMPTY_LIST;

    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    public void inferSvAndWriteVCF(final JavaRDD<AlignedContig> longReads, final String vcfOutputFileName,
                                   final Broadcast<ReferenceMultiSource> broadcastReference, final String fastaReference,
                                   final GCSOptions options, final Logger toolLogger) {

        final JavaPairRDD<ChimericAlignment, byte[]> inversionOrInvertedInsertion =
                longReads
                        .filter(read -> read.alignmentIntervals.size() > 1) // remove contigs who after filtering & gap-split has only one alignment
                        .filter(read -> !isLikelyInvertedDuplication(read)) // split between suspected inv dup VS strand-switch breakpoint
                        .mapToPair(InternalVariantDetectorFromLongReadAlignmentsForSimpleStrandSwitch::convertAlignmentIntervalToChimericAlignment)
                        .filter(Objects::nonNull);

        toolLogger.info(inversionOrInvertedInsertion.count() + " chimera indicating simple strand-switch breakpoints");

        // deal with strand-switch breakpoint
        final JavaRDD<VariantContext> annotatedVariants =
                inversionOrInvertedInsertion
                        .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2), pair._1))
                        .groupByKey()
                        .mapToPair(noveltyAndEvidence -> inferType(noveltyAndEvidence._1, noveltyAndEvidence._2, broadcastReference.getValue()))
                        .flatMap(noveltyTypeAndEvidence ->
                                AnnotatedVariantProducer.produceMultipleAnnotatedVcFromNovelAdjacency(noveltyTypeAndEvidence._1,
                                        noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2, broadcastReference));

        SVVCFWriter.writeVCF(options, vcfOutputFileName, fastaReference, annotatedVariants, toolLogger);

        // deal with suspected inv dup
//        final JavaRDD<AlignedContig> invDupSuspects =
//                longReads
//                        .filter(contig -> contig.alignmentIntervals.size() > 1) // remove contigs who after filtering & gap-split has only one alignment
//                        .filter(InternalVariantDetectorFromLongReadAlignmentsForSimpleStrandSwitch::isLikelyInvertedDuplication);
    }

    /**
     * Taking advantage of the fact that for input read, we know it has only two alignments that map to the same reference
     * chromosome, with strand switch.
     * @return  null if the pair of the alignments are no strong enough to support a strand switch breakpoint,
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
     * @param intervalOne
     * @param intervalTwo
     * @return
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

    /**
     * @return true iff the overlap of the two AI's made up the CA is strictly greater than half of any of the two AI's ref span.
     */
    @VisibleForTesting
    static boolean isLikelyInvertedDuplication(final AlignedContig longRead) {
        final SimpleInterval refSpanOne = longRead.alignmentIntervals.get(0).referenceSpan,
                             refSpanTwo = longRead.alignmentIntervals.get(1).referenceSpan;

        final int dummyChr = 1; // dummy number for chromosome to be used in constructing SVInterval, since we know the input CA has 2 AI and both map to the same chr
        final SVInterval intOne = new SVInterval(dummyChr, refSpanOne.getStart(), refSpanOne.getEnd() + 1),
                         intTwo = new SVInterval(dummyChr, refSpanTwo.getStart(), refSpanTwo.getEnd() + 1);

        final int overlap = intOne.overlapLen(intTwo);
        return 2 * overlap > Math.min(intOne.getLength(), intTwo.getLength());
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

    @SuppressWarnings("unchecked")
    private static Tuple2<NovelAdjacencyReferenceLocations, Tuple2<Iterable<SvType>, Iterable<ChimericAlignment>>>
    inferType(final NovelAdjacencyReferenceLocations novelAdjacency, final Iterable<ChimericAlignment> chimericAlignments,
              final ReferenceMultiSource reference) {

        final BreakEndVariantType bkpt_1, bkpt_2;
        if (novelAdjacency.endConnectionType == NovelAdjacencyReferenceLocations.EndConnectionType.FIVE_TO_FIVE) {
            bkpt_1 = new INV55BND(novelAdjacency, true, reference);
            bkpt_2 = new INV55BND(novelAdjacency, false, reference);
        } else if (novelAdjacency.endConnectionType == NovelAdjacencyReferenceLocations.EndConnectionType.THREE_TO_THREE){
            bkpt_1 = new INV33BND(novelAdjacency, true, reference);
            bkpt_2 = new INV33BND(novelAdjacency, false, reference);
        } else {
            throw new GATKException("Wrong type of novel adjacency sent to wrong analysis pathway: no strand-switch being sent to strand-switch path. \n" +
                    Utils.stream(chimericAlignments).map(ChimericAlignment::onErrStringRep).collect(Collectors.toList()));
        }

        return new Tuple2<>(novelAdjacency, new Tuple2<>(Arrays.asList(bkpt_1, bkpt_2), chimericAlignments));
    }
}
