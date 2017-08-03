package org.broadinstitute.hellbender.tools.spark.sv.playground;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKsvVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Map;

public abstract class BreakEndVariantType extends SvType {
    @Override
    public final String toString() {
        return GATKsvVCFConstants.BREAKEND_STR;
    }

    // see VCF spec 4.2 or 4.3 for BND format ALT allele field for SV
    BreakEndVariantType(final String id, final Map<String, String> typeSpecificExtraAttributes,
                        final byte[] bases, final boolean bracketPointsLeft, final SimpleInterval novelAdjRefLoc,
                        final boolean basesFirst) {
        super(id, constructAltAllele(bases, bracketPointsLeft, novelAdjRefLoc, basesFirst), 0, typeSpecificExtraAttributes);
    }

    private static Allele constructAltAllele(final byte[] bases, final boolean bracketPointsLeft, final SimpleInterval novelAdjRefLoc,
                                             final boolean basesFirst) {
        final String s = bracketPointsLeft ? "]" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getStart() + "]"
                                           : "[" + novelAdjRefLoc.getContig() + ":" + novelAdjRefLoc.getStart() + "[";
        return Allele.create( basesFirst ? new String(bases) + s  : s + new String(bases) );
    }
}
