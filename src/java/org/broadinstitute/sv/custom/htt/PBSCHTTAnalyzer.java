/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2021 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.sv.custom.htt;


import org.broadinstitute.sv.commandline.ArgumentException;
import org.broadinstitute.sv.commandline.CommandLineParser;
import org.broadinstitute.sv.commandline.CommandLineProgram;
import org.broadinstitute.sv.dataset.SAMLocation;
import org.broadinstitute.sv.util.NumericMap;
import org.broadinstitute.sv.util.SequenceUtilities;
import org.broadinstitute.sv.util.io.ErrorCheckingPrintWriter;
import org.broadinstitute.sv.util.sam.SAMUtils;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CloseableIterator;

import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;


/**
 * Custom analysis of HTT STR on Pacbio single cell data.
 * This runs the first stage of the pipeline which implements basic read decoding.
 */
public class PBSCHTTAnalyzer
    extends CommandLineProgram
{
    @Argument(fullName="debug", shortName="debug", required=false,
              doc="Enable debugging (default false)")
    private String mDebugArg = null;
    private boolean mDebug = false;

    @Input(fullName="referenceFile", shortName="R", required=false,
           doc="Reference fasta file")
    private File mReferenceFileArg = null;

    @Input(fullName="inputFile", shortName="I", required=true,
           doc="Input SAM/BAM file(s) or .list file")
    private List<String> mInputFiles = null;

    @Argument(fullName="inputFileIndexCache", shortName="IC", required=false,
              doc="Index cache for accessing remote input files")
    private String mInputFileIndexCache = null;

    @Argument(fullName="scatter", shortName="scatter", required=false,
              doc="Scatter specfication (e.g. k@N for k from 1 to N)")
    private String mScatterArg = null;
    private int mScatterIndex = 0;
    private int mScatterCount = 0;
    private int mScatterReadCount = 0;

    @Output(fullName="statsOutputFile", shortName="statsOutputFile", required=false,
            doc="Output statistics file")
    private File mStatsOutputFile = null;
    private PrintWriter mStatsOutputWriter = null;

    @Output(fullName="detailOutputFile1", shortName="detailOutputFile1", required=false,
            doc="Output detail file")
    private File mDetailOutputFile1 = null;
    private PrintWriter mDetailOutputWriter1 = null;

    @Output(fullName="detailOutputFile2", shortName="detailOutputFile2", required=false,
            doc="Output detail file")
    private File mDetailOutputFile2 = null;
    private PrintWriter mDetailOutputWriter2 = null;

    // For testing
    @Output(fullName="boundaryOutputFile", shortName="boundaryOutputFile", required=false,
            doc="Output boundary file")
    private File mBoundaryOutputFile = null;
    private PrintWriter mBoundaryOutputWriter = null;

    @Argument(fullName="sample", shortName="sample", required=false,
              doc="Sample(s) to extract (or .list file)")
    private List<String> mSampleList = null;

    @Argument(fullName="libraryType", shortName="libraryType", required=false,
              doc="Type of library (default pacbio/10x)")
    private String mLibraryType = "pacbio/10x";


    private static final long PROCESSED_READ_INITIAL_INTERVAL = 1000;
    private static final long PROCESSED_READ_MAX_INTERVAL = 10000;

    private static Logger mLog = Logger.getLogger(PBSCHTTAnalyzer.class);
    private long mProcessedReadCount = 0;
    private long mProcessedReadLogInterval = PROCESSED_READ_INITIAL_INTERVAL;
    private Map<String, SampleStats> mSampleStatsMap = null;
    private String mDebugRead = null; // = "m64425e_230111_174744/2903/ccs";

    // Smith-Waterman parameters
    private boolean mDebugSW = false;
    private Matrix mSWScoringMatrix = null;
    private static final float SW_MATCH = 2f;
    private static final float SW_MISMATCH = -1f;
    private static final float SW_OPEN = 2.5f;
    private static final float SW_EXTEND = 0.5f;

    // We just search for the end of the p7 sequence.
    // This appears to be more reliable because the start is sometimes trimmed from the read.
    // private static final SEQ_P7 = toCharArray("CAAGCAGAAGACGGCATACGAGAT");
    private static final char[] SEQ_P7 = toCharArray("CAAGCAGAAGACGGCATACGAGAT");
    private static final char[] SEQ_P5 = toCharArray("GTGTAGATCTCGGTGGTCGCCGTATCATT");
    // Adapter 1: GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
    private static final char[] SEQ_ADAPTER1_START = toCharArray("GTCTCGTGGG");
    private static final char[] SEQ_ADAPTER1_END = toCharArray("TAAGAGACAG");
    // Adapter 2: AGATCGGAAGAGCGTCGTGTAGCTGTCTCTTATACACATCTGACGCTGCCGACGA
    private static final char[] SEQ_ADAPTER2_START = toCharArray("AGATCGGAAGAG");
    private static final char[] SEQ_ADAPTER2_END = toCharArray("CTGCCGACGA");

    private static final char[] SEQ_EXON1_ANCHOR1 = toCharArray("CCTTCGAGTCCCTCAAGTCCTTCCAG");
    private static final char[] SEQ_EXON1_ANCHOR2_START = toCharArray("CCGCCACCGCCGCCGCCG");
    private static final char[] SEQ_EXON1_ANCHOR2_END = toCharArray("GAGCCGCTGCACCGACC");

    private static final char[] SEQ_REPEAT_LEFT_ANCHOR = toCharArray("CCTCAAGTCCTTCCAG");
    private static final char[] SEQ_CAG = toCharArray("CAG");


    public static void main(String[] args) throws Exception {
        run(new PBSCHTTAnalyzer(), args);
    }

    protected int run() {
        disableJAlignerLogging();

        if (mDebugArg != null) {
            mDebug = CommandLineParser.parseBoolean(mDebugArg);
            if (mDebug) {
                mLog.setLevel(Level.DEBUG);
            }
        }

        if (mScatterArg != null) {
            parseScatterArgument(mScatterArg);
        }

        mSWScoringMatrix = MatrixGenerator.generate(SW_MATCH, SW_MISMATCH);

        mSampleStatsMap = new HashMap<>();
        List<SAMLocation> inputFileList = CommandLineParser.parseSAMLocations(mInputFiles, mInputFileIndexCache, mReferenceFileArg);

        if (inputFileList.isEmpty()) {
            throw new ArgumentException("An empty set of input bam/cram files was specified.");
        }

        List<String> sampleList = mSampleList;
        Set<String> sampleSet = null;
        if (sampleList != null) {
            sampleSet = new HashSet<>(sampleList);
        }

        if (mStatsOutputFile != null) {
            mStatsOutputWriter = new ErrorCheckingPrintWriter(mStatsOutputFile);
            mStatsOutputWriter.println(formatStatsHeader());
        }
        if (mDetailOutputFile1 != null) {
            mDetailOutputWriter1 = new ErrorCheckingPrintWriter(mDetailOutputFile1);
            mDetailOutputWriter1.println(formatDetailHeader());
        }
        if (mDetailOutputFile2 != null) {
            mDetailOutputWriter2 = new ErrorCheckingPrintWriter(mDetailOutputFile2);
            mDetailOutputWriter2.println(formatDetailHeader());
        }
        if (mBoundaryOutputFile != null) {
            mBoundaryOutputWriter = new ErrorCheckingPrintWriter(mBoundaryOutputFile);
            mBoundaryOutputWriter.println(formatBoundaryHeader());
        }

        for (SAMLocation inputFile : inputFileList) {
            info(String.format("Reading input file %s ...", inputFile));
            SamReader reader = inputFile.createSamFileReader(ValidationStringency.SILENT);
            // For now, we assume that the input is localized to the HTT exon.
            CloseableIterator<SAMRecord> iterator = reader.iterator();
            while (iterator.hasNext()) {
                SAMRecord record = iterator.next();

                // We do not immediately stop processing on reads not in the current partition.
                // This is solely so that we can keep an accurate NREADS count while scattering.
                // This could be eliminated without any adverse effects if performance or scalability becomes an issue.
                boolean inScatterPartition = isInScatterPartition(record);
                if (inScatterPartition) {
                    logProcessedRead(record);
                }

                String sample = SAMUtils.getSampleId(record);
                if (sampleSet != null && !sampleSet.contains(sample)) {
                    continue;
                }

                String readId = formatReadId(record);
                SampleStats stats = getSampleStats(sample);
                if (!inScatterPartition) {
                    stats.recordUncountedRead(readId);
                    continue;
                }
                stats.countRead(readId);
                stats.incrementAlignmentCount();
                String seqName = getSequenceName(record);
                if (seqName != null) {
                    stats.incrementSequenceCount(seqName);
                }

                if (!shouldProcessRead(record)) {
                    continue;
                }

                try {
                    processRead(record);
                } catch (RuntimeException exc) {
                    throw new RuntimeException("Error processing read " + formatReadId(record) +
                                               " @" + record.getReferenceName() + ":" + record.getAlignmentStart() +
                                               ": " + exc.getMessage(), exc);
                }
            }
            logScanEnd();
            iterator.close();
            CloserUtil.close(reader);
        }

        if (mStatsOutputWriter != null) {
            emitStatsReport();
        }

        CloserUtil.close(mStatsOutputWriter);
        CloserUtil.close(mDetailOutputWriter1);
        CloserUtil.close(mDetailOutputWriter2);
        CloserUtil.close(mBoundaryOutputWriter);

        info("Processing complete.");
        return 0;
    }

    private void parseScatterArgument(String arg) {
        if (arg == null || arg.length() == 0) {
            throw new IllegalArgumentException("Invalid empty -scatter argument");
        }
        boolean error = false;
        String[] fields = arg.split("@");
        if (fields.length != 2) {
            error = true;
        } else {
            try {
                mScatterIndex = Integer.parseInt(fields[0]);
                mScatterCount = Integer.parseInt(fields[1]);
            } catch (NumberFormatException exc) {
                error = true;
            }
            if (mScatterIndex <= 0 || mScatterCount <= 0 || mScatterIndex > mScatterCount) {
                error = true;
            }
        }
        if (error) {
            throw new ArgumentException("Invalid -scatter argument: " + arg);
        }
        info("Scattering partition " + mScatterIndex + " of " + mScatterCount);
    }

    private boolean isInScatterPartition(SAMRecord record) {
        if (mScatterCount == 0) {
            return true;
        }
        mScatterReadCount++;
        int index = (mScatterReadCount % mScatterCount);
        if (index == 0) {
            index = mScatterCount;
        }
        return (index == mScatterIndex);
    }

    private boolean shouldProcessRead(SAMRecord record) {
        if (record.getReadFailsVendorQualityCheckFlag() ||
            record.getDuplicateReadFlag() ||
            record.isSecondaryAlignment()) {
            return false;
        }
        return true;
    }

    private void processRead(SAMRecord record) {
        boolean onTarget = isOnTarget(record);
        boolean isUnaligned = isUnaligned(record);
        if (!(onTarget || isUnaligned)) {
            return;
        }
        SampleStats stats = getSampleStats(record);
        if (onTarget) {
            stats.incrementOnTargetCount();
        }
        ReadDecoding decoding1 = decodeRead1(record);
        ReadDecoding decoding2 = decodeRead2(decoding1);
        boolean decoded1 = (decoding1 != null && isFullyDecoded(decoding1));
        boolean decoded2 = (decoding2 != null && isFullyDecoded(decoding2));
        if (decoded1) {
            stats.incrementDecodedCount1();
        }
        if (decoded2) {
            stats.incrementDecodedCount2();
        }
        if (decoded1 && decoded2) {
            stats.incrementDecodedCount();
        }
        if (!onTarget && isUnaligned) {
            // Count an unaligned read as on-target only if either the left or right is more than half decoded.
            double decodedFraction1 = getDecodedFieldFraction(decoding1);
            double decodedFraction2 = getDecodedFieldFraction(decoding2);
            if (decodedFraction1 > 0.5 || decodedFraction2 > 0.5) {
                stats.incrementOnTargetCount();
                stats.incrementUnalignedOnTargetCount();
            }
        }
        if (mDetailOutputWriter1 != null) {
            mDetailOutputWriter1.println(formatDetailReportLine(record, decoding1));
        }
        if (mDetailOutputWriter2 != null) {
            mDetailOutputWriter2.println(formatDetailReportLine(record, decoding2));
        }
        if (mBoundaryOutputWriter != null) {
            analyzeRepeatBoundaries(record);
        }
    }

    private boolean isOnTarget(SAMRecord record) {
        return readOverlapsInterval(record, "chr4", 3074826, 3075088);
    }

    private boolean isUnaligned(SAMRecord record) {
        return record.getReadUnmappedFlag();
    }

    private double getDecodedFieldFraction(ReadDecoding decoding) {
        if (decoding == null) {
            return 0;
        }
        int count = getDecodedFieldCount(decoding);
        int denominator = decoding.getPositions().length;
        if (denominator == 0) {
            return 0;
        }
        return (count / (double) denominator);
    }

    private int getDecodedFieldCount(ReadDecoding decoding) {
        if (decoding == null) {
            return 0;
        }
        int count = 0;
        int[] positions = decoding.getPositions();
        for (int i = 0; i < positions.length; i++) {
            if (positions[i] > 0) {
                count++;
            }
        }
        return count;
    }

    private boolean isFullyDecoded(ReadDecoding decoding) {
        int[] positions = decoding.getPositions();
        for (int i = 0; i < positions.length; i++) {
            if (positions[i] <= 0) {
                return false;
            }
        }
        return true;
    }

    private ReadDecoding decodeRead1(SAMRecord record) {
        if (isUnaligned(record)) {
            ReadDecoding forwardDecoding = decodeRead1Internal(record, 'F');
            ReadDecoding reverseDecoding = decodeRead1Internal(record, 'R');
            return getBestDecoding(forwardDecoding, reverseDecoding);
        } else {
            return decodeRead1Internal(record, 'F');
        }
    }

    private ReadDecoding getBestDecoding(ReadDecoding decoding1, ReadDecoding decoding2) {
        if (decoding1 == null) {
            return decoding2;
        } else if (decoding2 == null) {
            return decoding1;
        } else {
            float conf1 = sumConfidences(decoding1.getConfidences());
            float conf2 = sumConfidences(decoding2.getConfidences());
            if (conf2 > conf1) {
                return decoding2;
            } else {
                return decoding1;
            }
        }
    }

    private float sumConfidences(float[] confidences) {
        float sum = 0;
        for (int i = 0; i < confidences.length; i++) {
            if (!Float.isNaN(confidences[i])) {
                sum += confidences[i];
            }
        }
        return sum;
    }

    private ReadDecoding decodeRead1Internal(SAMRecord record, char orientation) {
        boolean debugRead = (mDebugRead != null && mDebugRead.equals(formatReadId(record)));
        if (debugRead) {
            System.out.println("#DBG: decodeRead1(" + formatReadId(record) + ", " + orientation + "):");
        }

        int offset = -1;
        float[] tempConf = new float[1];
        char[] searchBases = null;
        char[] readBases = toCharArray(canonicalizeBases(record.getReadBases()));
        if (orientation == 'R') {
            readBases = SequenceUtilities.reverseComplement(readBases);
        }
        if (debugRead) {
            System.out.println("#DBG: readBases: " + new String(readBases));
        }

        int p7_start = 0;
        int p7_end = 0;
        float p7_conf = Float.NaN;
        searchBases = extractReadField(readBases, 1, StrictMath.min(50, readBases.length));
        offset = searchBasesSW(searchBases, SEQ_P7, 0, -1, tempConf);
        if (offset >= 0) {
            p7_start = 1;
            p7_end = offset + SEQ_P7.length;
            p7_conf = tempConf[0];
        }
        char[] p7_bases = extractReadField(readBases, p7_start, p7_end);

        int a1_start = 0;
        int a1_end = 0;
        float a1_conf = Float.NaN;
        int searchOffset = p7_end;
        searchBases = extractReadField(readBases, p7_end + 1, StrictMath.min(p7_end + 100, readBases.length));
        offset = searchBasesSW(searchBases, SEQ_ADAPTER1_START, 0, -1, tempConf);
        if (offset >= 0) {
            a1_start = searchOffset + offset + 1;
            a1_conf = tempConf[0];
        }
        offset = searchBasesSW(searchBases, SEQ_ADAPTER1_END, 0, -1, tempConf);
        if (offset >= 0) {
            a1_end = searchOffset + offset + SEQ_ADAPTER1_END.length;
            a1_conf = StrictMath.min(a1_conf, tempConf[0]);
        } else {
            a1_conf = Float.NaN;
        }
        char[] a1_bases = extractReadField(readBases, a1_start, a1_end);

        int i7_start = 0;
        int i7_end = 0;
        float i7_conf = Float.NaN;
        if (p7_end > 0 && a1_start > 0 && a1_start > p7_end + 1) {
            i7_start = p7_end + 1;
            i7_end = a1_start - 1;
            i7_conf = StrictMath.min(p7_conf, a1_conf);
        }
        char[] i7_bases = extractReadField(readBases, i7_start, i7_end);

        searchOffset = StrictMath.max(0, readBases.length - 200);
        searchBases = extractReadField(readBases, searchOffset + 1, readBases.length);

        int p5_start = 0;
        int p5_end = 0;
        float p5_conf = Float.NaN;
        offset = searchBasesSW(searchBases, SEQ_P5, 0, -1, tempConf);
        if (offset >= 0) {
            p5_start = searchOffset + offset + 1;
            p5_end = readBases.length;
            p5_conf = tempConf[0];
        }
        char[] p5_bases = extractReadField(readBases, p5_start, p5_end);

        int a2_start = 0;
        int a2_end = 0;
        float a2_conf = Float.NaN;
        offset = searchBasesSW(searchBases, SEQ_ADAPTER2_START, 0, -1, tempConf);
        if (offset >= 0) {
            a2_start = searchOffset + offset + 1;
            a2_conf = tempConf[0];
        }
        offset = searchBasesSW(searchBases, SEQ_ADAPTER2_END, 0, -1, tempConf);
        if (offset >= 0) {
            a2_end = searchOffset + offset + SEQ_ADAPTER2_END.length;
            a2_conf = StrictMath.min(a2_conf, tempConf[0]);
        } else {
            a2_conf = Float.NaN;
        }
        char[] a2_bases = extractReadField(readBases, a2_start, a2_end);

        int i5_start = 0;
        int i5_end = 0;
        float i5_conf = Float.NaN;
        if (a2_end > 0 && p5_start > 0 && p5_start > a2_end + 1) {
            i5_start = a2_end + 1;
            i5_end = p5_start - 1;
            i5_conf = StrictMath.min(p5_conf, a2_conf);
        }
        char[] i5_bases = extractReadField(readBases, i5_start, i5_end);

        int umicbc_start = 0;
        int umicbc_end = 0;
        float umicbc_conf = Float.NaN;
        if (a2_start > 0) {
            int umicbcMinLength = 28;
            int umicbcScanLength = 35;
            if (mLibraryType.equals("pacbio/pipseq")) {
                int umicbcLength = 54; // chemistry v4 requires >= 54
                umicbc_start = StrictMath.max(1, a2_start - umicbcLength);
                umicbc_end = a2_start - 1;
                umicbc_conf = 5.0f;
            } else {
                int scanStart = StrictMath.max(1, a2_start - umicbcScanLength);
                int scanEnd = StrictMath.max(1, a2_start - umicbcMinLength);
                for (umicbc_start = scanStart; umicbc_start < scanEnd; umicbc_start++) {
                    if (readBases[umicbc_start - 1] != 'A') {
                        break;
                    }
                }
                umicbc_end = a2_start - 1;
                // Similar to other high confidence values
                umicbc_conf = 5.0f;
            }
        }
        char[] umicbc_bases = extractReadField(readBases, umicbc_start, umicbc_end);

        int read_start = 0;
        int read_end = 0;
        float read_conf = Float.NaN;
        if (a1_end > 0) {
            read_start = a1_end + 1;
        }
        if (umicbc_start > 0) {
            read_end = umicbc_start - 1;
        }
        if (read_start > 0 && read_end > 0) {
            read_conf = StrictMath.min(a1_conf, umicbc_conf);
        }
        char[] read_bases = elideLongField(extractReadField(readBases, read_start, read_end));

        int fieldCount = 8;
        int[] positions = new int[2*fieldCount];
        float[] confs = new float[fieldCount];
        char[][] subFields = new char[fieldCount][];
        int idx = 0;
        positions[2*idx] = p7_start;
        positions[2*idx+1] = p7_end;
        confs[idx] = p7_conf;
        subFields[idx++] = p7_bases;
        positions[2*idx] = i7_start;
        positions[2*idx+1] = i7_end;
        confs[idx] = i7_conf;
        subFields[idx++] = i7_bases;
        positions[2*idx] = a1_start;
        positions[2*idx+1] = a1_end;
        confs[idx] = a1_conf;
        subFields[idx++] = a1_bases;
        positions[2*idx] = read_start;
        positions[2*idx+1] = read_end;
        confs[idx] = read_conf;
        subFields[idx++] = read_bases;
        positions[2*idx] = umicbc_start;
        positions[2*idx+1] = umicbc_end;
        confs[idx] = umicbc_conf;
        subFields[idx++] = umicbc_bases;
        positions[2*idx] = a2_start;
        positions[2*idx+1] = a2_end;
        confs[idx] = a2_conf;
        subFields[idx++] = a2_bases;
        positions[2*idx] = i5_start;
        positions[2*idx+1] = i5_end;
        confs[idx] = i5_conf;
        subFields[idx++] = i5_bases;
        positions[2*idx] = p5_start;
        positions[2*idx+1] = p5_end;
        confs[idx] = p5_conf;
        subFields[idx++] = p5_bases;

        return new ReadDecoding(record, orientation, positions, confs, subFields);
    }

    private ReadDecoding decodeRead2(ReadDecoding input) {
        SAMRecord record = input.getRecord();
        int readStart = input.getPositions()[6];
        int readEnd = input.getPositions()[7];
        if (readStart == 0 || readEnd == 0 || readEnd < readStart) {
            return null;
        }

        boolean debugRead = (mDebugRead != null && mDebugRead.equals(formatReadId(record)));
        if (debugRead) {
            System.out.println("#DBG: decodeRead2(" + formatReadId(record) + "):");
        }

        int offset = -1;
        float[] tempConf = new float[1];

        char orientation = input.getOrientation();
        char[] readBases = toCharArray(canonicalizeBases(record.getReadBases()));
        if (orientation == 'R') {
            readBases = SequenceUtilities.reverseComplement(readBases);
        }
        char[] coreBases = extractReadField(readBases, readStart, readEnd);
        char[] searchBases = null;

        if (debugRead) {
            System.out.println("#DBG: coreBases: " + formatBases(coreBases));
        }

        // This is the left end of exon 1.
        int reg1_start = 0;
        int reg1_end = 0;
        float reg1_conf = Float.NaN;
        searchBases = extractReadField(coreBases, 1, StrictMath.min(100, coreBases.length));
        offset = searchBasesSW(searchBases, SEQ_EXON1_ANCHOR1, 0, -1, tempConf);
        if (offset >= 0) {
            reg1_start = readStart + offset;
            reg1_end = readStart + offset + SEQ_EXON1_ANCHOR1.length - 3 - 1;
            reg1_conf = tempConf[0];
        }
        char[] reg1_bases = extractReadField(readBases, reg1_start, reg1_end);

        if (debugRead) {
            System.out.println("#DBG: reg1: " + formatBases(reg1_bases));
        }

        // This is the right end of exon 1.
        int reg2_start = 0;
        int reg2_end = 0;
        float reg2_conf = Float.NaN;
        offset = searchBasesSW(coreBases, SEQ_EXON1_ANCHOR2_START, 0, -1, tempConf);
        if (offset >= 0) {
            int priorOffset = offset;
            offset = searchBackwardsForCodon(coreBases, offset, SEQ_CAG, 50);
            if (debugRead) {
                System.out.println("#DBG: search backwards " + formatReadId(record) + " " + priorOffset + " -> " + offset);
            }
            reg2_start = readStart + offset;
            reg2_conf = tempConf[0];
        }
        offset = searchBasesSW(coreBases, SEQ_EXON1_ANCHOR2_END, 0, -1, tempConf);
        if (offset >= 0) {
            reg2_end = readStart + offset + SEQ_EXON1_ANCHOR2_END.length - 1;
            reg2_conf = StrictMath.min(tempConf[0], reg2_conf);
        }
        char[] reg2_bases = extractReadField(readBases, reg2_start, reg2_end);

        if (debugRead) {
            System.out.println("#DBG: reg2: " + formatBases(reg2_bases));
        }

        int repeat_start = 0;
        int repeat_end = 0;
        float repeat_conf = Float.NaN;
        if (reg1_end > 0 && reg2_start > 0 && reg2_start > reg1_end + 1) {
            repeat_start = reg1_end + 1;
            repeat_end = reg2_start - 1;
            repeat_conf = StrictMath.min(reg1_conf, reg2_conf);
        }
        char[] repeat_bases = extractReadField(readBases, repeat_start, repeat_end);

        // This is the transcript beyond the end of exon 1.
        int reg3_start = 0;
        int reg3_end = 0;
        float reg3_conf = Float.NaN;
        if (reg2_end > 0 && reg2_end < readEnd) {
            reg3_start = reg2_end + 1;
            reg3_end = readEnd;
            reg3_conf = reg2_conf;
        }
        char[] reg3_bases = extractReadField(readBases, reg3_start, reg3_end);

        int fieldCount = 4;
        int[] positions = new int[2*fieldCount];
        float[] confs = new float[fieldCount];
        char[][] subFields = new char[fieldCount][];
        int idx = 0;
        positions[2*idx] = reg1_start;
        positions[2*idx+1] = reg1_end;
        confs[idx] = reg1_conf;
        subFields[idx++] = reg1_bases;
        positions[2*idx] = repeat_start;
        positions[2*idx+1] = repeat_end;
        confs[idx] = repeat_conf;
        subFields[idx++] = repeat_bases;
        positions[2*idx] = reg2_start;
        positions[2*idx+1] = reg2_end;
        confs[idx] = reg2_conf;
        subFields[idx++] = reg2_bases;
        positions[2*idx] = reg3_start;
        positions[2*idx+1] = reg3_end;
        confs[idx] = reg3_conf;
        subFields[idx++] = reg3_bases;

        return new ReadDecoding(record, orientation, positions, confs, subFields);
    }

    private int searchBackwardsForCodon(char[] bases, int offset, char[] codon, int limit) {
        final int codonLength = codon.length;
        final int deltaLimit = StrictMath.min(limit, offset - codonLength);
        for (int delta = 0; delta <= deltaLimit; delta++) {
            int testOffset = offset - delta - codonLength;
            if (basesMatch(bases, testOffset, codon, 0, codonLength)) {
                return (offset-delta);
            }
        }
        return offset;
    }

    private boolean basesMatch(char[] bases1, int offset1, char[] bases2, int offset2, int length) {
        for (int i = 0; i < length; i++) {
            if (bases1[offset1+i] != bases2[offset2+i]) {
                return false;
            }
        }
        return true;
    }

    // This method takes read positions in 1-based coordinates.
    // It also explicitly does error checking on the positions so that callers do not need to.
    private char[] extractReadField(char[] bases, int start, int end) {
        if (bases == null || start <= 0 || start > bases.length || end <= 0 || end > bases.length || end < start) {
            return null;
        }
        int resultLength = end - start + 1;
        char[] result = new char[resultLength];
        System.arraycopy(bases, start - 1, result, 0, resultLength);
        return result;
    }

    private char[] elideLongField(char[] bases) {
        if (bases == null || bases.length <= 20) {
            return bases;
        }
        int resultLength = StrictMath.min(23, bases.length);
        char[] result = new char[resultLength];
        Arrays.fill(result, '.');
        System.arraycopy(bases, 0, result, 0, 10);
        System.arraycopy(bases, bases.length - 10, result, result.length - 10, 10);
        return result;
    }

    private void emitStatsReport() {
        if (mStatsOutputWriter == null) {
            return;
        }
        List<String> sampleList = new ArrayList<>(mSampleStatsMap.keySet());
        Collections.sort(sampleList);
        for (String sample : sampleList) {
            mStatsOutputWriter.println(formatStatsReportLine(mSampleStatsMap.get(sample)));
        }
    }

    private String formatStatsHeader() {
        StringBuilder builder = new StringBuilder();
        addField(builder, "SAMPLE");
        addField(builder, "NREADS");
        addField(builder, "NRECORDS");
        addField(builder, "NALIGNED");
        addField(builder, "NONTARGET");
        addField(builder, "NUNALIGNEDONTARGET");
        addField(builder, "NDECODE1");
        addField(builder, "NDECODE2");
        addField(builder, "NDECODE");
        return builder.toString();
    }

    private String formatStatsReportLine(SampleStats stats) {
        int nReads = stats.getReadCount();
        int nRecords = stats.getAlignmentCount();
        Integer nUnaligned = stats.getSequenceCountMap().get("*");
        int nAligned = nRecords - (nUnaligned == null ? 0 : nUnaligned);
        StringBuilder builder = new StringBuilder();
        addField(builder, stats.getSample());
        addField(builder, nReads);
        addField(builder, nRecords);
        addField(builder, nAligned);
        addField(builder, stats.getOnTargetCount());
        addField(builder, stats.getUnalignedOnTargetCount());
        addField(builder, stats.getDecodedCount1());
        addField(builder, stats.getDecodedCount2());
        addField(builder, stats.getDecodedCount());
        return builder.toString();
    }

    private String formatDetailHeader() {
        StringBuilder builder = new StringBuilder();
        addField(builder, "SAMPLE");
        addField(builder, "READNAME");
        addField(builder, "READLEN");
        addField(builder, "ALIGNLEN");
        addField(builder, "ALIGNSTART");
        addField(builder, "ALIGNEND");
        addField(builder, "ORIENTATION");
        addField(builder, "POSITIONS");
        addField(builder, "CONFS");
        addField(builder, "FIELDS");
        return builder.toString();
    }

    // TBD: Rename to formatDecodingReportLine ?
    private String formatDetailReportLine(SAMRecord record, ReadDecoding decoding) {
        String sample = SAMUtils.getSampleId(record);
        String readName = formatReadId(record);
        int readLength = getReadLength(record);
        int alignLength = getAlignmentLength(record);
        int alignStart = getAlignmentStart(record);
        int alignEnd = getAlignmentEnd(record);

        StringBuilder builder = new StringBuilder();
        addField(builder, formatString(sample));
        addField(builder, formatString(readName));
        addField(builder, readLength);
        addField(builder, alignLength);
        addField(builder, formatPosition(alignStart));
        addField(builder, formatPosition(alignEnd));
        if (decoding == null) {
            addField(builder, "NA");
            addField(builder, "NA");
            addField(builder, "NA");
            addField(builder, "NA");
        } else {
            addField(builder, formatOrientation(decoding.getOrientation()));
            addField(builder, formatPositions(decoding.getPositions()));
            addField(builder, formatConfidences(decoding.getConfidences()));
            addField(builder, formatReadFields(decoding.getSubFields()));
        }
        return builder.toString();
    }

    private String formatBoundaryHeader() {
        StringBuilder builder = new StringBuilder();
        addField(builder, "SAMPLE");
        addField(builder, "READNAME");
        addField(builder, "READLEN");
        addField(builder, "LEFTPOS");
        addField(builder, "RIGHTPOS");
        addField(builder, "LEFTCONF");
        addField(builder, "RIGHTCONF");
        return builder.toString();
    }

    private String formatBoundaryReportLine(SAMRecord record, int leftPos, int rightPos, float leftConf, float rightConf) {
        StringBuilder builder = new StringBuilder();
        addField(builder, SAMUtils.getSampleId(record));
        addField(builder, formatReadId(record));
        addField(builder, getReadLength(record));
        addField(builder, leftPos);
        addField(builder, rightPos);
        addField(builder, formatDouble("%1.1f", leftConf));
        addField(builder, formatDouble("%1.1f", rightConf));
        return builder.toString();
    }

    private String formatPositions(int[] positions) {
        if (positions == null || positions.length == 0) {
            return "NA";
        }
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < positions.length; i++) {
            if (i > 0) {
                builder.append(',');
            }
            builder.append(positions[i]);
        }
        return builder.toString();
    }

    private String formatPosition(int position) {
        if (position <= 0) {
            return "NA";
        }
        return Integer.toString(position);
    }

    private String formatConfidences(float[] confidences) {
        if (confidences == null || confidences.length == 0) {
            return "NA";
        }
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < confidences.length; i++) {
            if (i > 0) {
                builder.append(',');
            }
            builder.append(formatDouble("%1.1f", (double) confidences[i]));
        }
        return builder.toString();
    }

    private String formatReadFields(char[][] fields) {
        if (fields == null || fields.length == 0) {
            return "NA";
        }
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < fields.length; i++) {
            if (i > 0) {
                builder.append(',');
            }
            builder.append(formatReadField(fields[i]));
        }
        return builder.toString();
    }

    private String formatReadField(char[] field) {
        if (field == null || field.length == 0) {
            return "NA";
        }
        return new String(field);
    }

    private boolean readOverlapsInterval(SAMRecord record, String seqName, int start, int end) {
        if (!record.getReferenceName().equals(seqName)) {
            return false;
        }
        int pos1 = record.getAlignmentStart();
        int pos2 = record.getAlignmentEnd();
        return (pos1 <= end && pos2 >= start);
    }

    // This method analyzes the read boundaries from the alignments
    // with an extra search if necessary for the left boundary.
    // This is just used for comparison with full read decoding.
    private void analyzeRepeatBoundaries(SAMRecord record) {

        float leftConf = Float.NaN;
        float rightConf = Float.NaN;
        float defaultConfidence = 5.0f;

        int leftPos = getLeftBoundaryPosition(record);
        int rightPos = getRightBoundaryPosition(record);

        if (leftPos > 0) {
            leftConf = defaultConfidence;
        }
        if (rightPos > 0) {
            rightConf = defaultConfidence;
        }
        if (leftPos <= 0 && rightPos > 0) {
            float[] tempConf = new float[1];
            leftPos = searchForLeftBoundaryPosition(record, rightPos, tempConf);
            if (leftPos > 0) {
                leftConf = tempConf[0];
            }
        }
        if (mBoundaryOutputWriter != null) {
            mBoundaryOutputWriter.println(formatBoundaryReportLine(record, leftPos, rightPos, leftConf, rightConf));
        }
    }

    // Run a more complex heuristic scan for the left boundary if the aligner truncated the read.
    private int searchForLeftBoundaryPosition(SAMRecord record, int rightPos, float[] confidence) {
        // We use as the left anchor the bases to the left of the first CAG plus one CAG unit.
        int anchorOffset = SEQ_REPEAT_LEFT_ANCHOR.length - 3;
        char[] readBases = toCharArray(canonicalizeBases(record.getReadBases()));
        int endOffset = StrictMath.max(0, rightPos - 1 - anchorOffset);
        int offset = searchBasesSW(readBases, SEQ_REPEAT_LEFT_ANCHOR, 0, endOffset, confidence);
        int leftPos = offset + 1;
        if (leftPos > 0) {
            leftPos += anchorOffset;
        }
        return leftPos;
    }

    // Always search in forward direction
    // Allow consecutive best scores (return midpoint going left on ties)
    private int searchBasesSW(char[] bases, char[] query, int startOffset, int endOffset, float[] confidenceOut) {
        if (bases == null) {
            return -1;
        }
        final String queryString = new String(query);
        final int queryLength = query.length;
        final int lastOffset = StrictMath.min(bases.length - queryLength, (endOffset >= 0) ? endOffset : bases.length);
        final char[] buffer = new char[queryLength];
        int bestOffset = -1;
        int bestOffsetRun = 0;
        float bestScore = -1;
        float bestScore2 = -1;
        if (mDebugSW) {
            System.out.println("#DBG: SW start Q=" + new String(query) + " T=" + new String(bases));
            System.out.println("#DBG: SW scan offsets " + startOffset + " - " + lastOffset);
        }
        for (int offset = startOffset; offset <= lastOffset; offset++) {
            System.arraycopy(bases, offset, buffer, 0, queryLength);
            float score = scoreAlignmentSW(queryString, new String(buffer));
            if (score > bestScore) {
                bestOffset = offset;
                bestOffsetRun = 1;
                bestScore2 = bestScore;
                bestScore = score;
            } else if (score == bestScore) {
                if (offset == bestOffset + bestOffsetRun) {
                    bestOffsetRun++;
                } else {
                    bestOffsetRun = 1;
                    bestScore2 = bestScore;
                    bestScore = score;
                }
            } else if (score > bestScore2) {
                bestScore2 = score;
            }
            if (mDebugSW) {
                System.out.println("#DBG SW @" + offset + " score=" + score +
                                   " best=" + bestScore + " @" + bestOffset +
                                   " run=" + bestOffsetRun + " best2=" + bestScore2);
            }
        }
        float conf = bestScore - bestScore2;
        if (conf <= 0) {
            bestOffset = -1;
        }
        if (confidenceOut != null && bestScore2 >= 0) {
            confidenceOut[0] = conf;
        }
        /***
        We no longer return the middle of the best run, but rather the leftmost position.
        if (bestOffset >= 0) {
            // Use the middle of the best run, ties go to the left.
            bestOffset += bestOffsetRun/2;
        }
        ***/
        return bestOffset;
    }

    private float scoreAlignmentSW(String bases1, String bases2) {
        // For some unknown reason, this is sometimes not sticky ...
        disableJAlignerLogging();

        Sequence seq1 = new Sequence(bases1);
        Sequence seq2 = new Sequence(bases2);
        Alignment align = null;
        try {
            align = SmithWatermanGotoh.align(seq1, seq2, mSWScoringMatrix, SW_OPEN, SW_EXTEND);
        } catch (Exception e) {
            e.printStackTrace();
            return Float.NaN;
        }
        float score = align.getScore();
        int initialGapLength = align.getStart1();
        if (initialGapLength > 0) {
            score -= SW_OPEN + (initialGapLength * SW_EXTEND);
        }
        int terminalGapLength = bases1.length() - computeNonGapLength(align.getSequence1()) - initialGapLength;
        if (terminalGapLength > 0) {
            score -= SW_OPEN + (terminalGapLength * SW_EXTEND);
        }
        return score;
    }

    private int computeNonGapLength(char[] seq) {
        int length = 0;
        final int seqlen = seq.length;
        for (int i = 0; i < seqlen; i++) {
            if (!isGap(seq[i])) {
                length++;
            }
        }
        return length;
    }

    private boolean isGap(char value) {
        return (value == '-');
    }

    // Crude first stab. Should probably require more matching bases (for short reads).
    // For long reads, we perhaps do not want to require the aligned base to match - needs investigation.
    // Decided to make the boundary surround the CAG repeat exactly (even though first base of next codon on the right is also a C - CCG).
    // Also including the synonymous CAA repeat that often appears near the right end.
    private int getLeftBoundaryPosition(SAMRecord record) {
        int offset = getMatchingReadOffset(record, "chr4", 3074876, 'C');
        if (offset < 0) {
            return 0;
        }
        return (offset + 1 + 1);
    }

    private int getRightBoundaryPosition(SAMRecord record) {
        int offset = getMatchingReadOffset(record, "chr4", 3074940, 'C');
        if (offset < 0) {
            return 0;
        }
        return (offset + 1 - 1);
    }

    private int getMatchingReadOffset(SAMRecord record, String targetContig, int targetPos, char targetBase) {
        if (record.getReadUnmappedFlag() || !record.getContig().equals(targetContig)) {
            return -1;
        }
        if (!(record.getAlignmentStart() <= targetPos && record.getAlignmentEnd() >= targetPos)) {
            return -1;
        }
        byte[] bases = record.getReadBases();
        Cigar cigar = record.getCigar();
        int refPos = record.getAlignmentStart();
        int readOffset = 0;
        for (CigarElement ce : cigar.getCigarElements()) {
            final CigarOperator op = ce.getOperator();
            int opLength = ce.getLength();
            final boolean opConsumesRead = op.consumesReadBases();
            final boolean opConsumesReference = op.consumesReferenceBases();
            for (int i = 0; i < opLength; i++) {
                if (refPos == targetPos && opConsumesRead && opConsumesReference) {
                    char readBase = (char) canonicalizeBase(bases[readOffset]);
                    if (readBase == targetBase) {
                        return readOffset;
                    } else {
                        return -1;
                    }
                }
                if (opConsumesRead) {
                    readOffset++;
                }
                if (opConsumesReference) {
                    refPos++;
                }
            }
        }
        return -1;
    }

    private String getSequenceName(SAMRecord record) {
        // We take advantage of the fact that this returns "*" for unaligned reads.
        String seqName = record.getReferenceName();
        return seqName;
    }

    private int getReadLength(SAMRecord record) {
        return record.getReadLength();
    }

    private int getAlignmentLength(SAMRecord record) {
        if (isUnaligned(record)) {
            return 0;
        }
        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();
        return (alignmentEnd - alignmentStart + 1);
    }

    private int getAlignmentStart(SAMRecord record) {
        if (isUnaligned(record)) {
            return 0;
        }
        return record.getAlignmentStart();
    }

    private int getAlignmentEnd(SAMRecord record) {
        if (isUnaligned(record)) {
            return 0;
        }
        return record.getAlignmentEnd();
    }

    private static char[] toCharArray(String s) {
        if (s == null) {
            return null;
        }
        return s.toCharArray();
    }

    private char[] toCharArray(byte[] bytes) {
        if (bytes == null) {
            return null;
        }
        final int length = bytes.length;
        final char[] result = new char[length];
        for (int i = 0; i < length; i++) {
            result[i] = (char) (bytes[i] & 0xFF);
        }
        return result;
    }

    private String toBaseString(byte[] bases) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < bases.length; i++) {
            builder.append(toBaseChar(bases[i]));
        }
        return builder.toString();
    }

    private char toBaseChar(byte base) {
        return ((char) (base & 0xFF));
    }

    private byte[] canonicalizeBases(String bases) {
        int length = bases.length();
        byte[] result = new byte[length];
        for (int i = 0; i < length; i++) {
            result[i] = canonicalizeBase(bases.charAt(i));
        }
        return result;
    }

    private byte[] canonicalizeBases(byte[] bases) {
        int length = bases.length;
        byte[] result = new byte[length];
        for (int i = 0; i < length; i++) {
            result[i] = canonicalizeBase(bases[i]);
        }
        return result;
    }

    private byte canonicalizeBase(char base) {
        return (byte) (Character.toUpperCase(base) & 0xFF);
    }

    private byte canonicalizeBase(byte base) {
        return canonicalizeBase((char) (base & 0xFF));
    }

    private SampleStats getSampleStats(SAMRecord record) {
        String sample = SAMUtils.getSampleId(record);
        if (sample == null) {
            return null;
        }
        return getSampleStats(sample);
    }

    private SampleStats getSampleStats(String sample) {
        SampleStats stats = mSampleStatsMap.get(sample);
        if (stats == null) {
            stats = new SampleStats(sample);
            mSampleStatsMap.put(sample, stats);
        }
        return stats;
    }

    private void logScanEnd() {
        info(String.format("Processed %d read%s.", mProcessedReadCount, (mProcessedReadCount == 1) ? "" : "s"));
        mProcessedReadCount = 0;
        mProcessedReadLogInterval = PROCESSED_READ_INITIAL_INTERVAL;
    }

    private void logProcessedRead(SAMRecord record) {
        mProcessedReadCount++;
        if ((mProcessedReadCount % mProcessedReadLogInterval) == 0) {
            String seqName = formatSequenceName(record);
            int pos = record.getAlignmentStart();
            info(String.format("Processed %d reads @%s:%d ...", mProcessedReadCount, seqName, pos));
            updateProcessedReadLogInterval();
        }
    }

    private void updateProcessedReadLogInterval() {
        // We log at 10k intervals to start, then switch to 100k at 1M reads and 1M at 10M reads.
        if (mProcessedReadCount == 10 * mProcessedReadLogInterval && mProcessedReadLogInterval < PROCESSED_READ_MAX_INTERVAL) {
            mProcessedReadLogInterval = mProcessedReadLogInterval * 10;
        }
    }

    private String formatReadId(SAMRecord record) {
        String readId = record.getReadName();
        if (record.getReadPairedFlag()) {
            if (record.getFirstOfPairFlag()) {
                readId = readId + "/1";
            } else {
                readId = readId + "/2";
            }
        }
        return readId;
    }

    private String formatSequenceName(SAMRecord record) {
        String seqName = record.getContig();
        if (seqName == null) {
            // Unaligned read
            seqName = "*";
        }
        return seqName;
    }

    private String formatOrientation(char orientation) {
        if (orientation == 0) {
            return "NA";
        }
        return Character.toString(orientation);
    }

    private String formatDouble(String format, double value) {
        if (Double.isNaN(value)) {
            return "NA";
        }
        return String.format(format, value);
    }

    private String formatBases(char[] value) {
        if (value == null) {
            return null;
        }
        return new String(value);
    }

    private String formatString(Object value) {
        if (value == null) {
            return "NA";
        } else {
            return value.toString();
        }
    }

    private void addField(StringBuilder builder, Object value) {
        if (builder.length() > 0) {
            builder.append('\t');
        }
        builder.append((value == null) ? "NA" : value.toString());
    }

    private void debug(String message) {
        mLog.debug(message);
    }

    private void info(String message) {
        mLog.info(message);
    }

    private void warn(String message) {
        mLog.warn(message);
    }

    /**
     * JAligner logs by default.
     * This method disables the logging.
     */
    private static void disableJAlignerLogging() {
        // Disable JAligner logging, which is on by default.
        // System.out.println("#DBG: disabling JAligner logging ...");
        java.util.logging.Logger.getLogger("jaligner").setLevel(java.util.logging.Level.OFF);
    }

    private static class ReadDecoding
    {
        private SAMRecord mRecord;
        private char mOrientation;
        private int[] mPositions;
        private float[] mConfidences;
        private char[][] mSubFields;

        ReadDecoding(SAMRecord record, char orientation, int[] positions, float[] confidences, char[][] subFields) {
            mRecord = record;
            mOrientation = orientation;
            mPositions = positions;
            mConfidences = confidences;
            mSubFields = subFields;
        }

        public SAMRecord getRecord() { return mRecord; }
        public char getOrientation() { return mOrientation; }
        public int[] getPositions() { return mPositions; }
        public float[] getConfidences() { return mConfidences; }
        public char[][] getSubFields() { return mSubFields; }
    }

    private static class SampleStats
    {
        private String mSample = null;
        private int mReadCount = 0;
        private int mAlignmentCount = 0;
        private int mOnTargetCount = 0;
        private int mUnalignedOnTargetCount = 0;
        private int mDecodedCount = 0;
        private int mDecodedCount1 = 0;
        private int mDecodedCount2 = 0;
        private Set<String> mReadIdSet = null;
        private Map<String, Integer> mSequenceCountMap = null;

        SampleStats(String sample) {
            mSample = sample;
            mReadIdSet = new HashSet<>();
            mSequenceCountMap = new HashMap<>();
        }

        public String getSample() { return mSample; }
        public int getReadCount() { return mReadCount; }
        public int getAlignmentCount() { return mAlignmentCount; }
        public int getOnTargetCount() { return mOnTargetCount; }
        public int getUnalignedOnTargetCount() { return mUnalignedOnTargetCount; }
        public int getDecodedCount() { return mDecodedCount; }
        public int getDecodedCount1() { return mDecodedCount1; }
        public int getDecodedCount2() { return mDecodedCount2; }
        public Map<String, Integer> getSequenceCountMap() { return mSequenceCountMap; }

        public void countRead(String readId) {
            if (mReadIdSet.add(readId)) {
                mReadCount++;
            }
        }

        public void recordUncountedRead(String readId) {
            // Record that we observed a read, but do not count it.
            mReadIdSet.add(readId);
        }

        public void incrementSequenceCount(String seqName) {
            NumericMap.add(mSequenceCountMap, seqName, 1);
        }

        public void incrementAlignmentCount() { mAlignmentCount++; }
        public void incrementOnTargetCount() { mOnTargetCount++; }
        public void incrementUnalignedOnTargetCount() { mUnalignedOnTargetCount++; }
        public void incrementDecodedCount() { mDecodedCount++; }
        public void incrementDecodedCount1() { mDecodedCount1++; }
        public void incrementDecodedCount2() { mDecodedCount2++; }
    }
}
