package org.broadinstitute.hellbender.tools.walkers.annotator;


import org.broadinstitute.hellbender.exceptions.GATKException;
import java.util.Arrays;

/**
 * Created by gauthier on 11/1/16.
 *
 * Ported by emeryj 7/24/17
 */
public class Histogram {
    private Double binSize;
    private String precisionFormat;
    private String printDelim;
    final private Double BIN_EPSILON = 0.01;

    private CompressedDataList<Integer> dataList = new CompressedDataList<>();

    public Histogram() {
        this.binSize = 0.1;
        precisionFormat = "%.1f";
    }

    public Histogram(final Double binSize) {
        this.binSize = binSize;
        precisionFormat = "%." + Math.round(-Math.log10(binSize)) + "f";
    }

    public void add(final Double d) {
        long binKey = getBinnedValue(d);
        if (isValidBinKey(binKey))
            dataList.add((int)binKey);
        else
            throw new GATKException("Histogram values are suspiciously extreme.  Failed to add " + d + " to the Histogram.");
    }

    public void add(final Double d, final int count) {
        if (count < 1)
            throw new GATKException("Cannot add non-positive counts to Histogram.");
        long binKey = getBinnedValue(d);
        if (isValidBinKey(binKey))
            dataList.add((int)binKey, count);
        else
            throw new GATKException("Histogram values are suspiciously extreme.  Failed to add " + d + " to the Histogram.");
    }

    public void add(final Histogram h) {
        if (!this.binSize.equals(h.binSize))
            throw new GATKException("Histogram bin sizes are mismatched -- cannot add bin size " + this.binSize + " to " + h.binSize);
        this.dataList.add(h.dataList);
    }

    public Integer get(final Double d) {
        long binKey = getBinnedValue(d);
        if (isValidBinKey(binKey))
            return dataList.getValueCounts().get((int)binKey);
        else
            throw new GATKException("Requested value is suspiciously extreme.  Failed to retrieve " + d + " from the Histogram.");
    }

    /**
     *
     * @return may be null if Histogram is empty
     */

    public Double median() {
        int numItems = 0;
        for(final int count : dataList.valueCounts.values()) {
            numItems += count;
        }
        boolean oddNumberValues = true;
        if(numItems % 2 == 0)
            oddNumberValues = false;
        int medianIndex = (numItems+1)/2;

        int counter = 0;
        Double firstMedian = null;
        for(final Integer key : dataList.valueCounts.keySet()) {
            counter += dataList.valueCounts.get(key);
            if( counter > medianIndex) {
                if (firstMedian == null)
                    return key*binSize;
                else {
                    return (firstMedian+key)/2.0*binSize;
                }
            }
            if( counter == medianIndex) {
                if (oddNumberValues)
                    return key*binSize;
                else {
                    firstMedian = (double) key;
                }
            }
        }
        return null;
    }

    private long getBinnedValue(double d) {
        return Math.round(Math.floor((d+BIN_EPSILON*binSize)/binSize)); //add a little epsilon before division so values exactly on bin boundaries will stay in the same bin
    }

    private boolean isValidBinKey(long binnedValue) {
        return binnedValue <= Integer.MAX_VALUE && binnedValue >= Integer.MIN_VALUE;
    }

    @Override
    public String toString(){
        printDelim = ",";
        String str = "";
        Object[] keys = dataList.valueCounts.keySet().toArray();
        Arrays.sort(keys);
        for (Object i: keys){
            if(!str.isEmpty())
                str+=printDelim;
            str+=(String.format(precisionFormat,(double)(int)i*binSize)+printDelim+dataList.valueCounts.get(i));  //key i needs to be output with specified precision
        }
        return str;
    }
}
