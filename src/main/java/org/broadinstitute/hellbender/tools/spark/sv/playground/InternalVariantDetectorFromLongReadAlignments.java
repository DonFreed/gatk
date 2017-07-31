package org.broadinstitute.hellbender.tools.spark.sv.playground;

import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;

/**
 * A base class for workflow of variant breakpoint detection from split alignments, variant type interpretation,
 * and resolving complication {@link org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakpointComplications}.
 */
interface InternalVariantDetectorFromLongReadAlignments {

    abstract void inferSvAndWriteVCF(final JavaRDD<AlignedContig> longReads, final String vcfOutputFileName,
                                     final Broadcast<ReferenceMultiSource> broadcastReference, final String fastaReference,
                                     final GCSOptions options, final Logger toolLogger);
}
