package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

/**
 * Allele-specific implementation of rank sum test annotations
 */
public abstract class AS_RankSumTest extends RankSumTest implements ReducibleAnnotation {
    private static final Logger logger = Logger.getLogger(AS_RankSumTest.class);
    public static final String SPLIT_DELIM = "\\|"; //String.split takes a regex, so we need to escape the pipe
    public static final String PRINT_DELIM = "|";
    public static final String RAW_DELIM = ",";
    public static final String REDUCED_DELIM = ",";

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        //TODO only raw for now
//        if (AnnotationUtils.walkerRequiresRawData(callingWalker))
            return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getRawKeyName()));
//        else
//            return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                         final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        return annotateRawData(ref, vc, likelihoods);
    }

    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods ) {
        if ( likelihoods == null) {
            return Collections.emptyMap();
        }

        final Map<String, Object> annotations = new HashMap<>();
        final AlleleSpecificAnnotationData<CompressedDataList<Integer>> myData = initializeNewAnnotationData(vc.getAlleles());
        calculateRawData(vc, likelihoods, myData);
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), myData.getAttributeMap());
        if (annotationString == null){
            return Collections.emptyMap();
        }
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    private AlleleSpecificAnnotationData<CompressedDataList<Integer>> initializeNewAnnotationData(final List<Allele> vcAlleles) {
        final Map<Allele, CompressedDataList<Integer>> perAlleleValues = new HashMap<>();
        for (final Allele a : vcAlleles) {
            perAlleleValues.put(a, new CompressedDataList<>());
        }
        final AlleleSpecificAnnotationData<CompressedDataList<Integer>> ret = new AlleleSpecificAnnotationData<>(vcAlleles, perAlleleValues.toString());
        ret.setAttributeMap(perAlleleValues);
        return ret;
    }

    private String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, CompressedDataList<Integer>> perAlleleValues) {
        if (perAlleleValues.values().stream().allMatch(d -> d.isEmpty())){
            return null;
        }
        final StringBuilder annotationString = new StringBuilder();
        for (int i =0; i< vcAlleles.size(); i++) {
            if (i!=0) {
                annotationString.append(PRINT_DELIM);
            }
            final CompressedDataList<Integer> alleleValues = perAlleleValues.get(vcAlleles.get(i));
            annotationString.append(alleleValues);
        }
        return annotationString.toString();
    }

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    @Override
    public void calculateRawData(VariantContext vc, final ReadLikelihoods<Allele> likelihoods, ReducibleAnnotationData myData) {
        if(likelihoods == null) {
            return;
        }

        final int refLoc = vc.getStart();

        final Map<Allele, CompressedDataList<Integer>> perAlleleValues = myData.getAttributeMap();
        for ( final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles() ) {
            if (bestAllele.isInformative() && isUsableRead(bestAllele.read, refLoc)) {
                final OptionalDouble value = getElementForRead(bestAllele.read, refLoc, bestAllele);
                if (value.isPresent() && value.getAsDouble() != INVALID_ELEMENT_FROM_READ && perAlleleValues.containsKey(bestAllele.allele)) {
                    perAlleleValues.get(bestAllele.allele).add((int) value.getAsDouble());
                }
            }
        }
    }

    /**
     *
     * @param vc -- contains the final set of alleles, possibly subset by GenotypeGVCFs
     * @param originalVC -- used to get all the alleles for all gVCFs
     * @return
     */
    public  Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC) {
        if (!vc.hasAttribute(getRawKeyName()))
            return new HashMap<>();

        final String rawRankSumData = vc.getAttributeAsString(getRawKeyName(),null);
        if (rawRankSumData == null)
            return new HashMap<>();

        final Map<String,Object> annotations = new HashMap<>();
        final AlleleSpecificAnnotationData<Histogram> myData = new AlleleSpecificAnnotationData(originalVC.getAlleles(), rawRankSumData);
        parseCombinedDataString(myData);

        final Map<Allele, Double> perAltRankSumResults = calculateReducedData(myData.getAttributeMap(), myData.getRefAllele());
        //shortcut for no ref values
        if (perAltRankSumResults.isEmpty())
            return annotations;
        final String annotationString = makeReducedAnnotationString(vc, perAltRankSumResults);
        annotations.put(getKeyNames().get(0), annotationString);
        return annotations;
    }

    protected void parseCombinedDataString(final ReducibleAnnotationData<Histogram> myData) {
        final String rawDataString = myData.getRawData();
        String rawDataNoBrackets;
        final Map<Allele, Histogram> perAlleleValues = new HashMap<>();
        //Initialize maps
        for (final Allele current : myData.getAlleles()) {
            perAlleleValues.put(current, new Histogram());
        }
        //Map gives back list with []
        if (rawDataString.charAt(0) == '[') {
            rawDataNoBrackets = rawDataString.substring(1, rawDataString.length() - 1);
        }
        else {
            rawDataNoBrackets = rawDataString;
        }
        //rawDataPerAllele is a string representation of the Histogram for each allele in the form value, count, value, count...
        final String[] rawDataPerAllele = rawDataNoBrackets.split(SPLIT_DELIM);
        for (int i=0; i<rawDataPerAllele.length; i++) {
            final String alleleData = rawDataPerAllele[i];
            final Histogram alleleList = perAlleleValues.get(myData.getAlleles().get(i));
            final String[] rawListEntriesAsStringVector = alleleData.split(RAW_DELIM);
            for (int j=0; j<rawListEntriesAsStringVector.length; j+=2) {
                Double value;
                int count;
                if (!rawListEntriesAsStringVector[j].isEmpty()) {
                    value = Double.parseDouble(rawListEntriesAsStringVector[j].trim());
                    if (!rawListEntriesAsStringVector[j + 1].isEmpty()) {
                        count = Integer.parseInt(rawListEntriesAsStringVector[j + 1].trim());
                        if(!value.isNaN())
                            alleleList.add(value,count);
                    }
                }
            }
        }
        myData.setAttributeMap(perAlleleValues);
        myData.validateAllelesList();
    }

    public Map<Allele, Double> calculateReducedData(final Map<Allele, Histogram> perAlleleValues, final Allele ref) {
        final Map<Allele, Double> perAltRankSumResults = new HashMap<>();
        for (final Allele alt : perAlleleValues.keySet()) {
            if (alt.equals(ref, false))
                continue;
            if (perAlleleValues.get(alt) != null)
                perAltRankSumResults.put(alt, perAlleleValues.get(alt).median());
        }
        return perAltRankSumResults;
    }

    public Map<String, Object> combineRawData(final List<Allele> vcAlleles, final List<? extends ReducibleAnnotationData> annotationList) {
        //VC already contains merged alleles from ReferenceConfidenceVariantContextMerger
        final ReducibleAnnotationData combinedData = initializeNewAnnotationData(vcAlleles);

        for (final ReducibleAnnotationData currentValue : annotationList) {
            parseRawDataString(currentValue);
            combineAttributeMap(currentValue, combinedData);

        }
        final Map<String, Object> annotations = new HashMap<>();
        final String annotationString = makeCombinedAnnotationString(vcAlleles, combinedData.getAttributeMap());
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }


    protected void parseRawDataString(final ReducibleAnnotationData<Histogram> myData) {
        final String rawDataString = myData.getRawData();
        String rawDataNoBrackets;
        final Map<Allele, Histogram> perAlleleValues = new HashMap<>();
        //Initialize maps
        for (final Allele current : myData.getAlleles()) {
            perAlleleValues.put(current, new Histogram());
        }
        //Map gives back list with []
        if (rawDataString.charAt(0) == '[') {
            rawDataNoBrackets = rawDataString.substring(1, rawDataString.length() - 1);
        }
        else {
            rawDataNoBrackets = rawDataString;
        }
        //rawDataPerAllele is a per-sample list of the rank sum statistic for each allele
        final String[] rawDataPerAllele = rawDataNoBrackets.split(SPLIT_DELIM);
        for (int i=0; i<rawDataPerAllele.length; i++) {
            final String alleleData = rawDataPerAllele[i];
            final Histogram alleleList = perAlleleValues.get(myData.getAlleles().get(i));
            final String[] rawListEntriesAsStringVector = alleleData.split(",");
            for (int j=0; j<rawListEntriesAsStringVector.length; j++) {
                Double value;
                if (!rawListEntriesAsStringVector[j].isEmpty()) {
                    value = Double.parseDouble(rawListEntriesAsStringVector[j].trim());
                    if(!value.isNaN())
                        alleleList.add(value);
                }
            }
        }
        myData.setAttributeMap(perAlleleValues);
        myData.validateAllelesList();
    }


    protected void combineAttributeMap(final ReducibleAnnotationData<Histogram> toAdd, final ReducibleAnnotationData<Histogram> combined) {
        for (final Allele a : combined.getAlleles()) {
            if (toAdd.hasAttribute(a)) {
                final Histogram alleleData = combined.getAttribute(a);
                if (toAdd.getAttribute(a) != null) {
                    alleleData.add(toAdd.getAttribute(a));
                    combined.putAttribute(a, alleleData);
                }
            }
        }
    }

    private String makeReducedAnnotationString(final VariantContext vc, final Map<Allele,Double> perAltRankSumResults) {
        String annotationString = "";
        for (final Allele a : vc.getAlternateAlleles()) {
            if (!annotationString.isEmpty())
                annotationString += REDUCED_DELIM;
            if (!perAltRankSumResults.containsKey(a))
                logger.warn("ERROR: VC allele not found in annotation alleles -- maybe there was trimming?");
            else
                annotationString += String.format("%.3f", perAltRankSumResults.get(a));
        }
        return annotationString;
    }

    protected String makeCombinedAnnotationString(final List<Allele> vcAlleles, final Map<Allele, Histogram> perAlleleValues) {
        String annotationString = "";
        for (int i = 0; i< vcAlleles.size(); i++) {
            if (vcAlleles.get(i).isReference())
                continue;
            if (i != 0)
                annotationString += PRINT_DELIM;
            final Histogram alleleValue = perAlleleValues.get(vcAlleles.get(i));
            //can be null if there are no ref reads
            if (alleleValue == null)
                continue;
            annotationString += alleleValue.toString();
        }
        return annotationString;
    }
}
