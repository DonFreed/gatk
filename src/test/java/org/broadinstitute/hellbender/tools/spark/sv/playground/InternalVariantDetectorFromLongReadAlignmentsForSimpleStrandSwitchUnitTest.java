package org.broadinstitute.hellbender.tools.spark.sv.playground;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

public class InternalVariantDetectorFromLongReadAlignmentsForSimpleStrandSwitchUnitTest extends BaseTest {

    @Test
    public void test (){

        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(21, 1, 46709983);

        final SAMRecord one =
                ArtificialReadUtils.createArtificialRead(header, "asm024381:tig00001", 20, 6070057, "TGAAATGTATGTGTGAGGTATGCAGTATGTGTGTGAGGTAGTGTGCGATGTGTGTGTAGTGTATGTGGTTGTGTGAGGTATGTGGGGTGTGAGGAATGTATTGTGTATGTGTGATATATACTGATTGTGTGTAAGGGATGTGGAGTGTGTGGCATGTGTGTAAGGTAGGTGTGTGTTGTGTATATGTGAGCTGTATAGTGTCGGGGGGGTGTGAGGTATGTGGTGTATGTTATGTTTGAGATCTAGTGTGTGTGTATGGTGTGTGTGGGAGGTATGTGGGGTGTGTGGTGTGTGGTGTGTATGAGGTATGTAGTGTGAGGTGTGTGATGTGTAGTGTGTGGTGTGGGGTATGTGGTGTATGTGTGAAGTATGTGTTGTGTGATGTGTGGGTGATATTTGGTGCCGTGTGTGTGGTATATGGTGTGTGGTATGAGGTGTGTAGTGTGATATGTGTGGTGTGTAATATGTGGTGTGTGTGTGTGTGTGATATATGGTGTGTGTGGTGTTATGATGTGTGTTGTGAGGTATGTGGTGTCTGTGTGTGATATGTGATTTGGGTGTGAGGTGTGTGTGGTGTGGCGTGTGGTGTGTGTGATGTGATGTGTGTGTGACATGGGGTGGTGCGTGGTGTGGTGTGTGTGGTATGTGGTGGTTGGTGTGTATGTGGTGAGTGAGGGGTGTGTGGTGTGGGTGGTGTGTGTGGTGTGTGTGGTTTGTGGTGTGTGTGGTTTGTGGTGTGTGGTATGTGGTGTGTTGTGTGTGGTTTGTGGTATGGTGTGTGTGGTATGGTTGTGTGTGGTGTGGTGTGTGCTGTGTGTATGGTTTGTGGTGTGTGTGGTGTGT".getBytes(),
                        ArtificialReadUtils.createRandomReadQuals(843), "502M341S").convertToSAMRecord(header);
        one.setMappingQuality(60);
        one.setAttribute("NM", 0);
        one.setAttribute("AS", 502);

        final SAMRecord two =
                ArtificialReadUtils.createArtificialRead(header, "asm024381:tig00001", 20, 43467994, "ACACACCACACACACCACAAACCATACACACAGCACACACCACACCACACACAACCATACCACACACACCATACCACAAACCACACACAACACACCACATACCACACACCACAAACCACACACACCACAAACCACACACACCACACACACCACCCACACCACACACC".getBytes(),
                        ArtificialReadUtils.createRandomReadQuals(167), "167M676H").convertToSAMRecord(header);
        two.setSupplementaryAlignmentFlag(true);
        two.setMappingQuality(60);
        two.setReadNegativeStrandFlag(true);
        two.setAttribute("NM", 0);
        two.setAttribute("AS", 167);

        final AlignmentInterval intervalOne = new AlignmentInterval(one);
        final AlignmentInterval intervalTwo = new AlignmentInterval(two);

        final AlignedContig contig = new AlignedContig("asm024381:tig00001", one.getReadBases(),
                Arrays.asList(intervalOne, intervalTwo), false);

        Assert.assertFalse( InternalVariantDetectorFromLongReadAlignmentsForSimpleStrandSwitch.isLikelyInvertedDuplication(contig) );
    }
}
