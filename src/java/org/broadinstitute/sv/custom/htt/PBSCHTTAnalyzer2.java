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
import org.broadinstitute.sv.util.GenomeInterval;
import org.broadinstitute.sv.util.NumericMap;
import org.broadinstitute.sv.util.Unchecked;
import org.broadinstitute.sv.util.fasta.IndexedFastaFile;
import org.broadinstitute.sv.util.io.ErrorCheckingPrintWriter;
import org.broadinstitute.sv.util.io.TextFileLineIterator;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;

import htsjdk.samtools.util.CloserUtil;

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
 * Separate utility to parse PBSCHTTAnalyzer output files and do further downstream analysis.
 */
public class PBSCHTTAnalyzer2
    extends CommandLineProgram
{
    @Argument(fullName="debug", shortName="debug", required=false,
              doc="Enable debugging (default false)")
    private String mDebugArg = null;
    private boolean mDebug = false;

    @Argument(fullName="verbose", shortName="verbose", required=false,
              doc="Generate verbose output (default false)")
    private String mVerboseArg = null;
    private boolean mVerbose = false;

    @Input(fullName="referenceFile", shortName="R", required=false,
           doc="Reference fasta file")
    private File mReferenceFileArg = null;
    private IndexedFastaFile mReferenceFile = null;

    @Input(fullName="detailFile1", shortName="detailFile1", required=false,
           doc="Input detail file from PBSCHTTAnalyzer")
    private File mInputDetailFile1 = null;

    @Input(fullName="detailFile2", shortName="detailFile2", required=false,
           doc="Input detail2 file from PBSCHTTAnalyzer")
    private File mInputDetailFile2 = null;

    @Input(fullName="indexFile", shortName="indexFile", required=true,
           doc="Input file of reaction indexes")
    private File mReactionIndexFile = null;

    @Input(fullName="cbcFile", shortName="cbcFile", required=false,
           doc="Input files of known cell barcodes")
    private List<File> mInputCBCFileList = null;

    // For testing new index matching logic
    @Input(fullName="indexMatchMode", shortName="indexMatchMode", required=false,
           doc="Index matching mode for testing")
    private Integer mIndexMatchMode = null;

    @Output(fullName="outputFile1", shortName="outputFile1", required=false,
            doc="Output report file")
    private File mOutputFile1 = null;
    private ErrorCheckingPrintWriter mOutputWriter1 = null;

    @Output(fullName="outputFile2", shortName="outputFile2", required=false,
            doc="Output report file")
    private File mOutputFile2 = null;
    private ErrorCheckingPrintWriter mOutputWriter2 = null;

    @Output(fullName="outputFile3", shortName="outputFile3", required=false,
            doc="Output report file")
    private File mOutputFile3 = null;
    private ErrorCheckingPrintWriter mOutputWriter3 = null;

    // Smith-Waterman parameters
    private boolean mDebugSW = false;
    private Matrix mSWScoringMatrix = null;
    private static final float SW_MATCH = 2f;
    private static final float SW_MISMATCH = -1f;
    private static final float SW_OPEN = 2.5f;
    private static final float SW_EXTEND = 0.5f;

    // Categorization score thresholds, based on SW parameters above, set from empirical distributions
    private static final float POLYA_STRONG_THRESHOLD = 15.0f;
    private static final float INTRON1_STRONG_THRESHOLD = 20.0f;
    private static final float INTRON1_WEAK_THRESHOLD = 10.0f;
    private static final float INTRON1_NEG_THRESHOLD = 5.0f;
    private static final float EXON2_STRONG_THRESHOLD = 15.0f;
    private static final float EXON2_WEAK_THRESHOLD = 5.0f;
    private static final float EXON2_NEG_THRESHOLD = 0.0f;
    private static final float EXON3_STRONG_THRESHOLD = 15.0f;
    private static final float EXON3_WEAK_THRESHOLD = 15.0f;
    private static final float EXON3_NEG_THRESHOLD = 10.0f;

    // Constants guiding analysis of read through transcripts
    // The intron interval is currently hard-wired for hg38.
    private static String HTT_INTRON1_INTERVAL = "chr4:3075089-3086938";
    private static final int READ_THROUGH_QUERY_LENGTH = 20;
    private static final int READ_THROUGH_SEARCH_RADIUS = 200;

    private String[] mReactionIds = null;
    private String[][] mReactionIndexes = null;
    private String[][] mReactionIndexesTranspose = null;
    private Set<String>[] mKnownCellBarcodes = null;

    private GenomeInterval mIntron1Interval = null;
    private char[] mIntron1Bases = null;

    // For debugging (to be removed)
    private final boolean DEBUG_REPUNITS = false;
    private Map<String, Integer> mRepUnitGapMap = new HashMap<>();

    private static Logger mLog = Logger.getLogger(PBSCHTTAnalyzer2.class);


    public static void main(String[] args) throws Exception {
        run(new PBSCHTTAnalyzer2(), args);
    }

    protected int run() {
        disableJAlignerLogging();

        if (mDebugArg != null) {
            mDebug = CommandLineParser.parseBoolean(mDebugArg);
            if (mDebug) {
                mLog.setLevel(Level.DEBUG);
            }
        }
        if (mVerboseArg != null) {
            mVerbose = CommandLineParser.parseBoolean(mVerboseArg);
        }

        if (mReferenceFileArg != null) {
            mReferenceFile = IndexedFastaFile.open(mReferenceFileArg);
        }

        mSWScoringMatrix = MatrixGenerator.generate(SW_MATCH, SW_MISMATCH);
        info("Reading reaction info from " + mReactionIndexFile);
        parseReactionIndexFile(mReactionIndexFile);

        if (mInputCBCFileList != null) {
            int cbcTotal = 0;
            int nInputFiles = mInputCBCFileList.size();
            mKnownCellBarcodes = Unchecked.<Set<String>[]>cast(new Set[nInputFiles]);
            for (int i = 0; i < nInputFiles; i++) {
                mKnownCellBarcodes[i] = new HashSet<>(parseCellBarcodes(mInputCBCFileList.get(i)));
                cbcTotal += mKnownCellBarcodes[i].size();
            }
            info("Read " + cbcTotal + " cell barcodes from " + nInputFiles + " input files.");
        }

        if (mOutputFile1 != null) {
            mOutputWriter1 = new ErrorCheckingPrintWriter(mOutputFile1);
            mOutputWriter1.println(formatOutputHeader1());
        }
        if (mOutputFile2 != null) {
            mOutputWriter2 = new ErrorCheckingPrintWriter(mOutputFile2);
            mOutputWriter2.println(formatOutputHeader2());
        }
        if (mOutputFile3 != null) {
            if (mReferenceFile == null) {
                warn("No reference file supplied: Disabling analysis of read-through transcripts");
            } else {
                mIntron1Interval = GenomeInterval.parse(HTT_INTRON1_INTERVAL);
                mIntron1Bases = toCharArray(mReferenceFile.getSequence(mIntron1Interval));
                if (mIntron1Bases == null) {
                    throw new RuntimeException("Failed to read sequence for interval " + mIntron1Interval + " from " + mReferenceFileArg);
                }
                mOutputWriter3 = new ErrorCheckingPrintWriter(mOutputFile3);
                mOutputWriter3.println(formatOutputHeader3());
            }
        }

        if (mInputDetailFile1 != null) {
            processDetailFile1(mInputDetailFile1);
        }
        if (mInputDetailFile2 != null) {
            processDetailFile2(mInputDetailFile2);
        }

        if (DEBUG_REPUNITS) {
            for (Map.Entry<String, Integer> entry : mRepUnitGapMap.entrySet()) {
                System.out.println(String.format("#DBG: GAPMAP\t%s\t%d", entry.getKey(), entry.getValue()));
            }
        }

        CloserUtil.close(mOutputWriter1);
        CloserUtil.close(mOutputWriter2);
        CloserUtil.close(mOutputWriter3);
        return 0;
    }

    private void parseReactionIndexFile(File inputFile) {
        List<String> reactionIdList = new ArrayList<>();
        List<String[]> reactionIndexList = new ArrayList<>();
        TextFileLineIterator lineIterator = new TextFileLineIterator(inputFile);
        while (lineIterator.hasNext()) {
            String line = lineIterator.next();
            String[] fields = line.split("\t");
            if (fields.length != 3) {
                throw new RuntimeException("Invalid input line (wrong number of fields) in file " + inputFile + ": " + line);
            }
            if (lineIterator.getLineCount() == 1) {
                if (!fields[0].equals("RXNID")) {
                    throw new RuntimeException("Invalid header line in file " + inputFile + ": " + line);
                }
                continue;
            }
            String rxnId = fields[0];
            String index1 = fields[1];
            String index2 = fields[2];
            String[] indexes = new String[2];
            indexes[0] = index1;
            indexes[1] = index2;
            reactionIdList.add(rxnId);
            reactionIndexList.add(indexes);
        }
        int rxnCount = reactionIdList.size();
        mReactionIds = new String[rxnCount];
        mReactionIndexes = new String[rxnCount][];
        for (int i = 0; i < rxnCount; i++) {
            mReactionIds[i] = reactionIdList.get(i);
            mReactionIndexes[i] = reactionIndexList.get(i);
        }
        mReactionIndexesTranspose = transpose(mReactionIndexes);
    }

    private List<String> parseCellBarcodes(File inputFile) {
        List<String> cbcList = new ArrayList<>();
        TextFileLineIterator lineIterator = new TextFileLineIterator(inputFile);
        while (lineIterator.hasNext()) {
            String line = lineIterator.next();
            if (lineIterator.getLineCount() == 1 && line.length() > 0 && line.startsWith("#")) {
                continue;
            }
            String cbc = line.trim();
            cbcList.add(cbc);
        }
        return cbcList;
    }

    private void processDetailFile1(File inputFile) {
        TextFileLineIterator lineIterator = new TextFileLineIterator(inputFile);
        while (lineIterator.hasNext()) {
            String line = lineIterator.next();
            String[] fields = line.split("\t");
            if (fields.length != 10) {
                throw new RuntimeException("Invalid input line (wrong number of fields) in file " + inputFile + ": " + line);
            }
            if (lineIterator.getLineCount() == 1 && fields.length > 0 && fields[0].equals("SAMPLE")) {
                continue;
            }
            String sample = fields[0];
            String readName = fields[1];
            int readLength = Integer.parseInt(fields[2]);
            int alignedLength = Integer.parseInt(fields[3]);
            int alignmentStart = parsePosition(fields[4]);
            int alignmentEnd = parsePosition(fields[5]);
            char orientation = parseOrientation(fields[6]);
            int[] positions = parsePositionVector(fields[7]);
            // This field currently unused.
            // float[] confs = parseFloatVector(fields[8]);
            String[] readFields = parseStringVector(fields[9]);

            if (readFields.length != 8) {
                throw new RuntimeException("Invalid input line (wrong number of read subfields) in file " + inputFile + ": " + line);
            }
            String indexSeq1 = readFields[1];
            String indexSeq2 = readFields[6];
            String[] ecIndexes = errorCorrectIndexes(indexSeq1, indexSeq2);
            int rxnNumber = getReactionNumber(ecIndexes);
            String rxnId = null;
            if (rxnNumber > 0) {
                rxnId = mReactionIds[rxnNumber - 1];
            }

            String umicbc = readFields[4];
            String cbc = null;
            String umi = null;
            String eccbc = null;
            if (umicbc != null && umicbc.length() == 28) {
                umi = umicbc.substring(0,12);
                cbc = umicbc.substring(12);
                int cbcset = 0;
                if (mKnownCellBarcodes != null && rxnNumber > 0) {
                    if (matchBarcodes(cbc, mKnownCellBarcodes[rxnNumber - 1])) {
                        eccbc = cbc;
                    }
                }
            }
            if (mOutputWriter1 != null) {
                mOutputWriter1.println(formatOutputRecord1(readName, rxnNumber, rxnId,
                                                           indexSeq1, indexSeq2, ecIndexes[0], ecIndexes[1],
                                                           cbc, eccbc, umi, umicbc));
            }
        }
        lineIterator.close();
    }

    private int getReactionNumber(String[] ecIndexes) {
        if (ecIndexes == null) {
            return 0;
        }
        int rxn = 0;
        for (int i = 0; i < mReactionIndexes.length; i++) {
            if (Arrays.equals(ecIndexes, mReactionIndexes[i])) {
                if (rxn > 0) {
                    // Multiple indexes match
                    rxn = 0;
                    break;
                }
                rxn = i+1;
            }
        }
        return rxn;
    }

    private void processDetailFile2(File inputFile) {
        TextFileLineIterator lineIterator = new TextFileLineIterator(inputFile);
        while (lineIterator.hasNext()) {
            String line = lineIterator.next();
            String[] fields = line.split("\t");
            if (fields.length != 10) {
                throw new RuntimeException("Invalid input line (wrong number of fields) in file " + inputFile + ": " + line);
            }
            if (lineIterator.getLineCount() == 1 && fields.length > 0 && fields[0].equals("SAMPLE")) {
                continue;
            }
            String sample = fields[0];
            String readName = fields[1];
            int readLength = Integer.parseInt(fields[2]);
            int alignedLength = Integer.parseInt(fields[3]);
            int alignmentStart = parsePosition(fields[4]);
            int alignmentEnd = parsePosition(fields[5]);
            char orientation = parseOrientation(fields[6]);
            int[] positions = parsePositionVector(fields[7]);
            // This field currently unused.
            // float[] confs = parseFloatVector(fields[8]);
            String[] readFields = parseStringVector(fields[9]);

            float repLength = Float.NaN;
            float repPurity = Float.NaN;
            boolean repHasCAA = false;
            int repUnits = -1;
            String category = null;

            int polyACount = 0;
            int polyACount1 = 0;
            int tailLength = 0;
            float exon2Score = Float.NaN;
            float exon3Score = Float.NaN;
            float intron1Score = Float.NaN;

            if (readFields != null) {
                if (readFields.length != 4) {
                    throw new RuntimeException("Invalid input line (wrong number of read subfields) in file " + inputFile + ": " + line);
                }
                String repeat = readFields[1];
                if (repeat != null) {
                    int repCAASuffixLength = getCAASuffixLength(repeat);
                    if (repCAASuffixLength > 0) {
                        repHasCAA = true;
                        repeat = repeat.substring(0, repeat.length() - repCAASuffixLength);
                    }
                    repLength = repeat.length() / 3.0f;
                    repPurity = computeRepeatPurity(repeat);
                    repUnits = computeRepeatUnits(repeat);
                }
                // Moved inline from function to be able to emit other debugging/QC info (exon2Score, etc.)
                String exon1Tail = readFields[3];
                String exon2 = "AAAGAAAGAACTTTCA";
                String exon3 = "AAATTCTCCAGAATTT";
                String intron1 = "GTGAGTTTGGGCCCGC";
                if (exon1Tail != null) {
                    polyACount = getPolyALength(exon1Tail, 0);
                    polyACount1 = getPolyALength(exon1Tail, 1);
                    tailLength = exon1Tail.length();
                    String overhang = exon1Tail.substring(0, StrictMath.min(16, exon1Tail.length()));
                    exon2Score = scoreAlignmentSW(overhang, exon2);
                    exon3Score = scoreAlignmentSW(overhang, exon3);
                    intron1Score = scoreAlignmentSW(overhang, intron1);
                    if (polyACount1 >= POLYA_STRONG_THRESHOLD) {
                        category = "POLYA";
                    } else if (intron1Score >= INTRON1_STRONG_THRESHOLD ||
                               (intron1Score >= INTRON1_WEAK_THRESHOLD && exon2Score < EXON2_NEG_THRESHOLD && exon3Score < EXON3_NEG_THRESHOLD)) {
                        category = "INTRON1";
                    } else if (exon2Score >= EXON2_STRONG_THRESHOLD ||
                               (exon2Score >= EXON2_WEAK_THRESHOLD && exon3Score < EXON3_NEG_THRESHOLD && intron1Score < INTRON1_NEG_THRESHOLD)) {
                        category = "EXON2";
                    } else if (exon3Score >= EXON3_STRONG_THRESHOLD ||
                               (exon3Score >= EXON3_WEAK_THRESHOLD && exon2Score < EXON2_NEG_THRESHOLD && intron1Score < INTRON1_NEG_THRESHOLD)) {
                        category = "EXON3";
                    } else {
                        category = "OTHER";
                    }
                }
                if (exon1Tail != null && category != null && category.equals("INTRON1")) {
                    analyzeReadThrough(readName, readLength, exon1Tail);
                }
            }
            if (mOutputWriter2 != null) {
                mOutputWriter2.println(formatOutputRecord2(readName, readLength, alignmentStart, alignmentEnd,
                                                           repLength, repPurity, repUnits, repHasCAA,
                                                           category, polyACount, polyACount1, tailLength,
                                                           exon2Score, exon3Score, intron1Score));
            }
        }
        lineIterator.close();
    }

    private void analyzeReadThrough(String readName, int readLength, String exon1Tail) {

        if (mIntron1Bases == null || mOutputWriter3 == null) {
            return;
        }

        int queryLength = READ_THROUGH_QUERY_LENGTH;
        int searchRadius = READ_THROUGH_SEARCH_RADIUS;
        String trimmedTail = trimPolyA(exon1Tail);
        int tailLength = trimmedTail.length();
        int polyALength = exon1Tail.length() - tailLength;
        char[] readBases = toCharArray(trimmedTail.substring(StrictMath.max(0, tailLength - queryLength)));
        char[] refBases = mIntron1Bases;
        int intronStartPos = mIntron1Interval.getStart();

        int startOffset = StrictMath.max(0, tailLength - searchRadius);
        int endOffset = tailLength + searchRadius;
        float[] confOut = new float[1];
        int offset = searchBasesSW(refBases, readBases, startOffset, endOffset, confOut);
        int endPos = (offset < 0) ? 0 : (intronStartPos + offset + readBases.length - 1);
        float endScore = confOut[0];

        if (mOutputWriter3 != null) {
            String outReadBases = trimmedTail.substring(StrictMath.max(0, tailLength - 16));
            String outRefBases = null;
            if (endPos > 0) {
                outRefBases = mReferenceFile.getSequence(mIntron1Interval.getSequenceName(), StrictMath.max(1, endPos - 16 + 1), endPos);
            }
            mOutputWriter3.println(formatOutputRecord3(readName, readLength, tailLength, polyALength, endPos, endScore, outReadBases, outRefBases));
        }
    }

    private String trimPolyA(String sequence) {
        int length = sequence.length();
        for (int i = 0; i < length; i++) {
            char ch = sequence.charAt(length-i-1);
            if (!(ch == 'A' || ch == 'a')) {
                if (i == 0) {
                    return sequence;
                } else {
                    return sequence.substring(0, length-i);
                }
            }
        }
        return "";
    }

    private float computeRepeatPurity(String sequence) {
        if (sequence == null || sequence.equals("")) {
            return Float.NaN;
        }
        int numerator = 0;
        int length = sequence.length();
        int index = 0;
        while (index < length) {
            if (index <= length-3) {
                String codon = sequence.substring(index, index+3);
                if (isRepeatCodon(codon)) {
                    numerator += 3;
                    index += 3;
                    continue;
                }
            }
            if (index == 0) {
                // Check for prefixes (give benefit of the doubt) but only test for CAG as prefix
                if (length >= 2 && sequence.substring(0,2).equals("GC")) {
                    numerator += 1;
                    index += 1;
                    continue;
                }
                if (length >= 3 && sequence.substring(0,3).equals("AGC")) {
                    numerator += 2;
                    index += 2;
                    continue;
                }
            }
            if (length - index < 3) {
                // Check for suffixes (give benefit of the doubt)
                int tailLength = length - index;
                String suffix = sequence.substring(index, length);
                if (tailLength == 1 && isRepeatCodon(suffix + "AG")) {
                    numerator += tailLength;
                } else if (tailLength == 2 && isRepeatCodon(suffix + "G")) {
                    numerator += tailLength;
                }
                // After processing tail, we are done
                break;
            }
            index++;
        }
        float purity = numerator / (float) length;
        // System.out.println("#DBG: pur = " + String.format("%1.5f", purity) + " num = " + numerator + " den = " + length + " seq= " + sequence);
        return purity;
    }

    private boolean isRepeatCodon(String seq) {
        return (seq.equals("CAG") || seq.equals("CAA"));
    }

    private int getCAASuffixLength(String sequence) {
        if (sequence.endsWith("CAACAG")) {
            return 6;
        }
        return 0;
    }

    private int computeRepeatUnits(String sequence) {
        if (sequence == null || sequence.equals("")) {
            return 0;
        }
        int repUnits = 0;
        int index = 0;
        int lastIndex = 0;
        final int length = sequence.length();
        while (index < length) {
            if (index <= length-3) {
                String codon = sequence.substring(index, index+3);
                if (codon.equals("CAG")) {
                    if (index > lastIndex && DEBUG_REPUNITS) {
                        String gap = sequence.substring(lastIndex, index);
                        NumericMap.add(mRepUnitGapMap, gap, 1);
                    }
                    repUnits += 1;
                    index += 3;
                    lastIndex = index;
                    continue;
                }
            }
            index++;
        }
        if (index > lastIndex && DEBUG_REPUNITS) {
            String gap = sequence.substring(lastIndex, index);
            NumericMap.add(mRepUnitGapMap, gap, 1);
        }
        return repUnits;
    }

    private int getPolyALength(String sequence, int maxMismatches) {
        int count = 0;
        int mismatches = 0;
        final int length = sequence.length();
        for (int i = 0; i < length; i++) {
            char base = sequence.charAt(i);
            if (base != 'A') {
                mismatches++;
                if (mismatches > maxMismatches) {
                    break;
                }
            }
            count++;
        }
        return count;
    }

    private String[] errorCorrectIndexes(String index1, String index2) {
        if (mIndexMatchMode == null || mIndexMatchMode.intValue() == 1) {
            // Default to the old way of error correcting and matching indexes.
            return errorCorrectIndexesV1(index1, index2);
        } else if (mIndexMatchMode.intValue() == 2) {
            return errorCorrectIndexesV2(index1, index2);
        } else {
            throw new ArgumentException("Unrecognized index match mode: " + mIndexMatchMode);
        }
    }

    private String[] errorCorrectIndexesV1(String index1, String index2) {
        int which = matchBothIndexesExact(index1, index2, mReactionIndexesTranspose[0], mReactionIndexesTranspose[1]);
        if (which >= 0) {
            return mReactionIndexes[which];
        }

        int which1 = matchIndexesOriginal(index1, mReactionIndexesTranspose[0], 0.5f);
        int which2 = matchIndexesOriginal(index2, mReactionIndexesTranspose[1], 0.5f);
        if (which1 >= 0 && which2 >= 0 && which1 != which2) {
            which1 = -1;
            which2 = -1;
        } else if (which1 < 0 && which2 >= 0) {
            which1 = which2;
        } else if (which2 < 0 && which1 >= 0) {
            which2 = which1;
        }
        if (which1 != which2) {
            throw new RuntimeException("Internal error: indexes do not match");
        }
        if (which1 < 0) {
            return new String[2];
        }
        return mReactionIndexes[which1];
    }

    private String[] errorCorrectIndexesV2(String index1, String index2) {
        int which = matchBothIndexesExact(index1, index2, mReactionIndexesTranspose[0], mReactionIndexesTranspose[1]);
        if (which >= 0) {
            return mReactionIndexes[which];
        }
        String match1 = matchIndexes(index1, mReactionIndexesTranspose[0], 0.5f);
        String match2 = matchIndexes(index2, mReactionIndexesTranspose[1], 0.5f);
        String[] result = { match1, match2 };
        return(result);
    }

    private String[] errorCorrectIndexesV2Obsolete(String index1, String index2) {
        int which = matchBothIndexesExact(index1, index2, mReactionIndexesTranspose[0], mReactionIndexesTranspose[1]);
        if (which >= 0) {
            return mReactionIndexes[which];
        }
        String ecIndex1 = errorCorrectIndex(index1, mReactionIndexesTranspose[0], 0.5f);
        String ecIndex2 = errorCorrectIndex(index2, mReactionIndexesTranspose[1], 0.5f);
        String[] ecIndexes = new String[2];
        ecIndexes[0] = ecIndex1;
        ecIndexes[1] = ecIndex2;
        return(ecIndexes);
    }

    private String errorCorrectIndex(String query, String[] targets, float threshold) {
        if (query == null || query.length() == 0) {
            return null;
        }
        // Quick check for an exact match
        final int targetCount = targets.length;
        for (int i = 0; i < targetCount; i++) {
            if (targets[i].equals(query)) {
                return targets[i];
            }
        }

        int bestIndex = -1;
        float bestScore = Float.NaN;
        float bestScore2 = Float.NaN;
        for (int i = 0; i < targetCount; i++) {
            if (bestIndex >= 0 && targets[i].equals(targets[bestIndex])) {
                // Ignore redundant targets
                continue;
            }
            float score = scoreAlignmentSW(query, targets[i]);
            // System.out.println("#DBG: errorCorrectIndex " + query + "," + i + ": score " + score);
            if (Float.isNaN(score)) {
                continue;
            }
            if (bestIndex < 0 || score >= bestScore) {
                bestIndex = i;
                bestScore2 = bestScore;
                bestScore = score;
            } else if (score > bestScore2) {
                bestScore2 = score;
            }
        }
        if (!Float.isNaN(bestScore2)) {
            float conf = bestScore - bestScore2;
            if (conf < threshold) {
                bestIndex = -1;
            }
        }
        if (bestIndex < 0) {
            return null;
        }
        return targets[bestIndex];
    }

    private String[][] transpose(String[][] input) {
        if (input == null) {
            return null;
        }
        int dim1 = input.length;
        if (dim1 == 0) {
            throw new RuntimeException("Cannot transpose one dimensional matrix");
        }
        int dim2 = input[0].length;
        for (int i = 0; i < dim1; i++) {
            if (input[i].length != dim2) {
                throw new RuntimeException("Cannot transpose ragged matrix");
            }
        }
        String[][] result = new String[dim2][];
        for (int i = 0; i < dim2; i++) {
            result[i] = new String[dim1];
            for (int j = 0; j < dim1; j++) {
                result[i][j] = input[j][i];
            }
        }
        return result;
    }

    private boolean matchBarcodes(String query, Set<String> targets) {
        return targets.contains(query);
    }

    private int matchBarcodes(String query, String[] targets) {
        // For now, exact match only on barcodes.
        if (query == null || query.length() == 0) {
            return -1;
        }
        // Quick check for exact matches (assumes targets are unique)
        final int targetCount = targets.length;
        for (int i = 0; i < targetCount; i++) {
            if (targets[i].equals(query)) {
                return i;
            }
        }
        return -1;
    }

    private int matchBothIndexesExact(String query1, String query2, String[] targets1, String[] targets2) {
        final int targetCount = targets1.length;
        if (targetCount != targets2.length) {
            throw new RuntimeException("Internal error: mismatched target arrays");
        }
        if (query1 == null || query1.length() == 0 || query2 == null || query2.length() == 0) {
            return -1;
        }
        int bestIndex = -1;
        for (int i = 0; i < targetCount; i++) {
            if (targets1[i].equals(query1) && targets2[i].equals(query2)) {
                if (bestIndex >= 0) {
                    return -1;
                }
                bestIndex = i;
            }
        }
        return bestIndex;
    }

    private String matchIndexes(String query, String[] targets, float threshold) {
        if (query == null || query.length() == 0) {
            return null;
        }
        // Quick check for exact matches
        final int targetCount = targets.length;
        for (int i = 0; i < targetCount; i++) {
            if (targets[i].equals(query)) {
                return(query);
            }
        }
        String bestMatch = null;
        float bestScore = Float.NaN;
        float bestScore2 = Float.NaN;
        for (int i = 0; i < targetCount; i++) {
            if (bestMatch != null && targets[i].equals(bestMatch)) {
                continue;
            }
            float score = scoreAlignmentSW(query, targets[i]);
            // System.out.println("#DBG: matchIndexes " + query + "," + i + ": score " + score);
            if (Float.isNaN(score)) {
                continue;
            }
            if (bestMatch == null || score >= bestScore) {
                bestMatch = targets[i];
                bestScore2 = bestScore;
                bestScore = score;
            } else if (score > bestScore2) {
                bestScore2 = score;
            }
        }
        if (!Float.isNaN(bestScore2)) {
            float conf = bestScore - bestScore2;
            if (conf < threshold) {
                bestMatch = null;
            }
        }
        // System.out.println("#DBG: matchIndexes " + query + ": SW " + bestMatch);
        return bestMatch;
    }

    private int matchIndexesOriginal(String query, String[] targets, float threshold) {
        if (query == null || query.length() == 0) {
            return -1;
        }
        // Quick check for exact matches
        // If targets are non-unique, return an ambiguous match.
        final int targetCount = targets.length;
        int bestIndex = -1;
        for (int i = 0; i < targetCount; i++) {
            if (targets[i].equals(query)) {
                // System.out.println("#DBG: matchIndexes " + query + ": exact " + i);
                if (bestIndex >= 0) {
                    return -1;
                }
                bestIndex = i;
            }
        }
        if (bestIndex >= 0) {
            return bestIndex;
        }
        float bestScore = Float.NaN;
        float bestScore2 = Float.NaN;
        for (int i = 0; i < targetCount; i++) {
            float score = scoreAlignmentSW(query, targets[i]);
            // System.out.println("#DBG: matchIndexes " + query + "," + i + ": score " + score);
            if (Float.isNaN(score)) {
                continue;
            }
            if (bestIndex < 0 || score >= bestScore) {
                bestIndex = i;
                bestScore2 = bestScore;
                bestScore = score;
            } else if (score > bestScore2) {
                bestScore2 = score;
            }
        }
        if (!Float.isNaN(bestScore2)) {
            float conf = bestScore - bestScore2;
            if (conf < threshold) {
                bestIndex = -1;
            }
        }
        // System.out.println("#DBG: matchIndexes " + query + ": SW " + bestIndex);
        return bestIndex;
    }

    // Always search in forward direction
    // Allow consecutive best scores (return midpoint going left on ties)
    private int searchBasesSW(char[] bases, char[] query, int startOffset, int endOffset, float[] confidenceOut) {
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
        if (bestOffset >= 0) {
            // Use the middle of the best run, ties go to the left.
            bestOffset += bestOffsetRun/2;
        }
        return bestOffset;
    }

    private float scoreAlignmentSW(String bases1, String bases2) {
        // For some unknown reason, this is Sometimes not sticky ...
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

    private int[] parsePositionVector(String text) {
        if (text.equals("NA")) {
            return null;
        }
        String[] fields = text.split(",");
        int[] result = new int[fields.length];
        for (int i = 0; i < fields.length; i++) {
            result[i] = parsePosition(fields[i]);
        }
        return result;
    }

    private int parsePosition(String text) {
        if (text == null || text.equals("") || text.equals("NA")) {
            return -1;
        }
        return Integer.parseInt(text);
    }

    private char parseOrientation(String text) {
        if (text == null || text.equals("") || text.equals("NA")) {
            return 0;
        }
        if (text.length() != 1) {
            throw new RuntimeException("Invalid orientation value: " + text);
        }
        return text.charAt(0);
    }

    private String[] parseStringVector(String text) {
        if (text.equals("NA")) {
            return null;
        }
        String[] fields = text.split(",");
        String[] result = new String[fields.length];
        for (int i = 0; i < fields.length; i++) {
            result[i] = fields[i];
            if (fields[i] != null && fields[i].equals("NA")) {
                result[i] = null;
            }
        }
        return result;
    }

    private String formatOutputHeader1() {
        StringBuilder builder = new StringBuilder();
        addField(builder, "READNAME");
        addField(builder, "RXN");
        addField(builder, "RXNID");
        addField(builder, "INDEX1");
        addField(builder, "INDEX2");
        addField(builder, "ECINDEX1");
        addField(builder, "ECINDEX2");
        addField(builder, "CBC");
        addField(builder, "ECCBC");
        addField(builder, "UMI");
        addField(builder, "RAWUMI");
        return builder.toString();
    }

    private String formatOutputRecord1(String readName, int rxnNumber, String rxnId,
                                       String index1, String index2, String ecIndex1, String ecIndex2,
                                       String cbc, String eccbc, String umi, String umicbc) {
        StringBuilder builder = new StringBuilder();
        addField(builder, formatString(readName));
        addField(builder, rxnNumber);
        addField(builder, formatString(rxnId));
        addField(builder, formatString(index1));
        addField(builder, formatString(index2));
        addField(builder, formatString(ecIndex1));
        addField(builder, formatString(ecIndex2));
        addField(builder, formatString(cbc));
        addField(builder, formatString(eccbc));
        addField(builder, formatString(umi));
        addField(builder, formatString(umicbc));
        return builder.toString();
    }

    private String formatOutputHeader2() {
        StringBuilder builder = new StringBuilder();
        addField(builder, "READNAME");
        addField(builder, "READLENGTH");
        addField(builder, "ALIGNSTART");
        addField(builder, "ALIGNEND");
        addField(builder, "REPUNITS");
        addField(builder, "REPHASCAA");
        addField(builder, "REPLENGTH");
        addField(builder, "REPPURITY");
        addField(builder, "CATEGORY");
        addField(builder, "POLYA0");
        addField(builder, "POLYA1");
        addField(builder, "TAILLEN");
        addField(builder, "EXON2");
        addField(builder, "EXON3");
        addField(builder, "INTRON1");
        return builder.toString();
    }

    private String formatOutputRecord2(String readName, int readLength, int alignStart, int alignEnd,
                                       float repLength, float repPurity, int repUnits, boolean repHasCAA,
                                       String category, int polyACount, int polyACount1, int tailLength,
                                       float exon2Score, float exon3Score, float intron1Score) {
        StringBuilder builder = new StringBuilder();
        addField(builder, formatString(readName));
        addField(builder, readLength);
        addField(builder, formatPosition(alignStart));
        addField(builder, formatPosition(alignEnd));
        addField(builder, repUnits);
        addField(builder, (repHasCAA ? "T" : "F"));
        addField(builder, formatDouble("%1.1f", repLength));
        addField(builder, formatDouble("%1.3f", repPurity));
        addField(builder, formatString(category));
        addField(builder, polyACount);
        addField(builder, polyACount1);
        addField(builder, tailLength);
        addField(builder, formatDouble("%1.1f", exon2Score));
        addField(builder, formatDouble("%1.1f", exon3Score));
        addField(builder, formatDouble("%1.1f", intron1Score));
        return builder.toString();
    }

    private String formatOutputHeader3() {
        StringBuilder builder = new StringBuilder();
        addField(builder, "READNAME");
        addField(builder, "READLENGTH");
        addField(builder, "TAILLEN");
        addField(builder, "POLYALEN");
        addField(builder, "ENDPOS");
        addField(builder, "ENDSCORE");
        addField(builder, "READBASES");
        addField(builder, "REFBASES");
        return builder.toString();
    }

    private String formatOutputRecord3(String readName, int readLength, int tailLength, int polyALength,
                                       int endPos, float endScore, String readBases, String refBases) {
        StringBuilder builder = new StringBuilder();
        addField(builder, formatString(readName));
        addField(builder, readLength);
        addField(builder, tailLength);
        addField(builder, polyALength);
        addField(builder, endPos);
        addField(builder, formatDouble("%1.1f", endScore));
        addField(builder, formatString(readBases));
        addField(builder, formatString(refBases));
        return builder.toString();
    }

    private String formatPosition(int position) {
        if (position <= 0) {
            return "NA";
        }
        return Integer.toString(position);
    }

    private String formatDouble(String format, double value) {
        if (Double.isNaN(value)) {
            return "NA";
        }
        return String.format(format, value);
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

    private static char[] toCharArray(String s) {
        if (s == null) {
            return null;
        }
        return s.toCharArray();
    }

    private void debug(String message) {
        mLog.debug(message);
    }

    private void info(String message) {
        if (mVerbose || mDebug) {
            mLog.info(message);
        }
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
}
