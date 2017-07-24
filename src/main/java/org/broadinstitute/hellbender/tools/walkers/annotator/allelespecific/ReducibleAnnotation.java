package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.util.List;
import java.util.Map;

/**
 * An interface for annotations that are calculated using raw data across samples, rather than the median (or median of median) of samples values
 *
 * TODO: this interface will eventually contain the actual 'reduce' functionality.
 */
public interface ReducibleAnnotation extends Annotation {
    public abstract String getRawKeyName();

    /**
     * Generate the raw data necessary to calculate the annotation. Raw data is the final endpoint for gVCFs.
     * @param ref the reference context for this annotation
     * @param vc the variant context to annotate
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     */
    public abstract Map<String, Object> annotateRawData(final ReferenceContext ref,
                                                        final VariantContext vc,
                                                        final ReadLikelihoods<Allele> likelihoods);

    /**
     * Combine raw data, typically during the merging of raw data contained in multiple gVCFs as in CombineGVCFs and the
     * preliminary merge for GenotypeGVCFs
     * @param allelesList   The merged allele list across all variants being combined/merged
     * @param listOfRawData The raw data for all the variants being combined/merged
     * @return
     */
    public abstract Map<String, Object> combineRawData(final List<Allele> allelesList, final List <? extends ReducibleAnnotationData> listOfRawData);

    /**
     * Calculate the final annotation value from the raw data
     * @param vc -- contains the final set of alleles, possibly subset by GenotypeGVCFs
     * @param originalVC -- used to get all the alleles for all gVCFs
     * @return
     */
    public abstract Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC);


    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME: there's no way to properly use generics here because of method overriding
    public abstract void calculateRawData(VariantContext vc, final ReadLikelihoods<Allele> likelihoods, ReducibleAnnotationData rawAnnotations);
}
