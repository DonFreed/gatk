package org.broadinstitute.hellbender.tools.spark.sv.playground;

import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark.annotateVariant;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark.inferType;

final class ForSimpleInsDel implements VariantDetectorFromLongReadAlignments {

    public void inferSvAndWriteVCF(final JavaRDD<AlignedContig> longReads, final String vcfOutputFileName,
                                   final Broadcast<ReferenceMultiSource> broadcastReference, final String fastaReference,
                                   final GCSOptions options, final Logger toolLogger){

        // convert to ChimericAlignment, similar to ChimericAlignment.parseOneContig(final AlignedContig, final int), except the head/tail filtering
        final JavaPairRDD<byte[], List<ChimericAlignment>> chimericAlignments =
                longReads
                        .mapToPair(ForSimpleInsDel::convertAlignmentIntervalToChimericAlignment);

        // usual business as in DiscoverVariantsFromContigAlignmentsSAMSpark#discoverVariantsAndWriteVCF()
        final JavaRDD<VariantContext> annotatedVariants =
                chimericAlignments
                        .flatMapToPair(DiscoverVariantsFromContigAlignmentsSAMSpark::discoverNovelAdjacencyFromChimericAlignments)
                        .groupByKey()
                        .mapToPair(noveltyAndEvidence -> inferType(noveltyAndEvidence._1, noveltyAndEvidence._2))
                        .map(noveltyTypeAndEvidence -> annotateVariant(noveltyTypeAndEvidence._1,
                                noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2, broadcastReference));

        SVVCFWriter.writeVCF(options, vcfOutputFileName, fastaReference, annotatedVariants, toolLogger);
    }

    /**
     * Very similar to {@link ChimericAlignment#parseOneContig(AlignedContig, int)}, except the head/tail filtering.
     */
    private static Tuple2<byte[], List<ChimericAlignment>> convertAlignmentIntervalToChimericAlignment
    (final AlignedContig contig) {

        final List<AlignmentInterval> alignmentIntervals = contig.alignmentIntervals;
        final Iterator<AlignmentInterval> iterator = alignmentIntervals.iterator();
        AlignmentInterval current = iterator.next();
        final List<ChimericAlignment> results = new ArrayList<>(alignmentIntervals.size() - 1);
        final List<String> insertionMappings = new ArrayList<>();
        while ( iterator.hasNext() ) {
            final AlignmentInterval next = iterator.next();
            if (nextAlignmentMayBeNovelInsertion(current, next, DEFAULT_MIN_ALIGNMENT_LENGTH)) {
                if (iterator.hasNext()) {
                    insertionMappings.add(next.toPackedString());
                    continue;
                } else {
                    break;
                }
            }
            final boolean isNotSimpleTranslocation = isNotSimpleTranslocation(current, next,
                    determineStrandSwitch(current, next), involvesRefPositionSwitch(current, next));
            if (isNotSimpleTranslocation)
                results.add(new ChimericAlignment(current, next, insertionMappings, contig.contigName));
        }
        return new Tuple2<>(contig.contigSequence,results);
    }
}
