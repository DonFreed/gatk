package org.broadinstitute.hellbender.tools.spark.sv.playground;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.UncheckedIOException;
import java.util.Iterator;

class FileUtils {

    static void writeSAMFile(final Iterator<SAMRecord> alignments, final SAMFileHeader header, final String outputName,
                             final boolean preOrdered) {
        try ( final SAMFileWriter writer = createSAMFileWriter(outputName, header, preOrdered) ) {
            alignments.forEachRemaining(writer::addAlignment);
        } catch ( final UncheckedIOException ie) {
            throw new GATKException("Can't write SAM file to the specified location: " + outputName, ie);
        }
    }

    private static SAMFileWriter createSAMFileWriter(final String samFile, final SAMFileHeader header, final boolean preOrdered) {
        final SAMFileWriterFactory factory = new SAMFileWriterFactory();
        final int lastDotIndex = samFile.lastIndexOf('.');
        if (lastDotIndex >= 0) {
            final String extension = samFile.substring(lastDotIndex).toLowerCase();
            if (extension.equals(BamFileIoUtils.BAM_FILE_EXTENSION)) {
                return factory.makeBAMWriter(header, preOrdered, BucketUtils.createFile(samFile));
            } else if (extension.equals(".sam")) {
                return factory.makeSAMWriter(header, preOrdered, BucketUtils.createFile(samFile));
            } else {
                throw new GATKException("unsupported read alignment file name extension (." + extension + ") in requested name: " + samFile);
            }
        } else {
            throw new GATKException("cannot determine the alignment file format from its name: " + samFile);
        }
    }

    static void writeLinesToSingleFile(final Iterator<String> linesToWrite, final String fileName) {
        try ( final OutputStream writer =
                      new BufferedOutputStream(BucketUtils.createFile(fileName)) ) {
            while (linesToWrite.hasNext()) {
                writer.write(linesToWrite.next().getBytes()); writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write "+fileName, ioe);
        }
    }
}
