package Project2_github.Nucleotide_Sequence_Alignment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

public class SN {
    private static final int MATCH = 1;
    private static final int MISMATCH = -1;
    private static final int GAP_OPEN = -10;
    private static final float GAP_EXTEND = -0.5f;

    public static void main(String[] args) {
        String file = "sequence.fasta";

        try {
            Map<String, String> sequences = parseFasta(file);
            // Check if there are exactly two sequences in the input file
            if (sequences.size() != 2) {
                System.err.println("The input FASTA file must contain exactly two sequences.");
                System.exit(1);
            }

            // Retrieve the two sequences from the parsed sequences map
            String[] ids = sequences.keySet().toArray(new String[0]);
            String seq1 = sequences.get(ids[0]);
            String seq2 = sequences.get(ids[1]);

            // Align the sequences and store the result
            AlignmentResult result = alignSequences(seq1, seq2);

            // Print the alignment result
            System.out.println(">" + ids[0] + " vs " + ids[1]);
            System.out.println(result.sequence1);
            System.out.println(result.sequence2);
            System.out.println("Total penalty: " + result.penalty);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // Method to parse a FASTA file containing two sequences
    private static Map<String, String> parseFasta(String filename) throws IOException {
        Map<String, String> sequences = new LinkedHashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String line;
            String currentID = null;
            StringBuilder currentSeq = new StringBuilder();
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (currentID != null) {
                        sequences.put(currentID, currentSeq.toString());
                        currentSeq.setLength(0);
                    }
                    currentID = line.substring(1);
                } else {
                    currentSeq.append(line);
                }
            }
            if (currentID != null) {
                sequences.put(currentID, currentSeq.toString());
            }
        }
        return sequences;
    }

    // Alignment result, including the aligned sequences and the penalty score
    private static class AlignmentResult {
        String sequence1;
        String sequence2;
        float penalty;

        // Constructor for the AlignmentResult class
        AlignmentResult(String sequence1, String sequence2, float penalty) {
            this.sequence1 = sequence1;
            this.sequence2 = sequence2;
            this.penalty = penalty;
        }
    }

    // Method to align two sequences using dynamic programming
    private static AlignmentResult alignSequences(String seq1, String seq2) {
        int m = seq1.length();
        int n = seq2.length();
        // Initialize the dynamic programming, gap1, and gap2 matrices
        float[][] dp = new float[m + 1][n + 1];
        float[][] gap1 = new float[m + 1][n + 1];
        float[][] gap2 = new float[m + 1][n + 1];
        // Initialize the first column of the dp and gap1 matrices
        for (int i = 1; i <= m; i++) {
            gap1[i][0] = GAP_OPEN + i * GAP_EXTEND;
            dp[i][0] = gap1[i][0];
        }
        // Initialize the first row of the dp and gap2 matrices
        for (int j = 1; j <= n; j++) {
            gap2[0][j] = GAP_OPEN + j * GAP_EXTEND;
            dp[0][j] = gap2[0][j];
        }
        // Calculate the alignment scores using dynamic programming
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                // Check if the characters match or mismatch, and update the score accordingly
                int match = seq1.charAt(i - 1) == seq2.charAt(j - 1) ? MATCH : MISMATCH;
                dp[i][j] = dp[i - 1][j - 1] + match;
                // Update gap1 and gap2 matrices using the previously calculated values and
                // penalties
                gap1[i][j] = Math.max(dp[i - 1][j] + GAP_OPEN, gap1[i - 1][j] + GAP_EXTEND);
                gap2[i][j] = Math.max(dp[i][j - 1] + GAP_OPEN, gap2[i][j - 1] + GAP_EXTEND);
                // Update the current cell of the dp matrix with the maximum score from the
                // three matrices
                dp[i][j] = Math.max(dp[i][j], Math.max(gap1[i][j], gap2[i][j]));
            }
        }
        // Create StringBuilders to store the aligned sequences
        StringBuilder alignedSeq1 = new StringBuilder();
        StringBuilder alignedSeq2 = new StringBuilder();
        // Traceback through the dp matrix to build the aligned sequences
        int i = m, j = n;
        while (i > 0 || j > 0) {
            // If we're moving diagonally, append the corresponding characters to the
            // aligned sequences
            if (i > 0 && j > 0
                    && dp[i][j] == dp[i - 1][j - 1] + (seq1.charAt(i - 1) == seq2.charAt(j - 1) ? MATCH : MISMATCH)) {
                alignedSeq1.append(seq1.charAt(i - 1));
                alignedSeq2.append(seq2.charAt(j - 1));
                i--;
                j--;
                // If we're moving vertically, append the character from seq1 to alignedSeq1 and
                // a gap to alignedSeq2
            } else if (i > 0 && dp[i][j] == gap1[i][j]) {
                alignedSeq1.append(seq1.charAt(i - 1));
                alignedSeq2.append('_');
                i--;
                // If we're moving horizontally, append a gap to alignedSeq1 and the character
                // from seq2 to alignedSeq2
            } else if (j > 0) {
                alignedSeq1.append('_');
                alignedSeq2.append(seq2.charAt(j - 1));
                j--;
            }
        }
        // Reverse the aligned sequences to obtain the correct order
        alignedSeq1.reverse();
        alignedSeq2.reverse();
        // Return a new AlignmentResult object containing the aligned sequences and the
        // final penalty score
        return new AlignmentResult(alignedSeq1.toString(), alignedSeq2.toString(), dp[m][n]);
    }
}    
