import java.util.*;
import java.io.*;

public class SN {
    private static final float MATCH = 1;
    private static final float MISMATCH = -4;
    private static final float GAP_OPEN = -10;
    private static final float GAP_EXTEND = -0.5f;

    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in); 
        System.out.println("enter the name of the fasta file you wish to use"); 
        String file = scan.nextLine(); 


        try {
            Map<String, String> sequences = parseFasta(file);

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

    // printing array
    public static void pa(float[][] arr) {
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr[i].length; j++) {
                System.out.printf("%-15s", arr[i][j]);
            }
            System.out.println();
        }
    }

    // add val to whole row
    public static void rr(float[][] arr, int rowIndex, float value) {
        for (int i = 1; i < arr[rowIndex].length; i++) {
            if (rowIndex != 0) {
                return;
            }
            arr[rowIndex][i] = value;
        }
    }

    // add val to whole column
    public static void rc(float[][] arr, int colIndex, float value) {
        for (int i = 1; i < arr.length; i++) {
            if (colIndex != 0) {
                return;
            }
            arr[i][colIndex] = value;
        }
    }

    // Method to align two sequences using dynamic programming

    private static AlignmentResult alignSequences(String seq1, String seq2) {
        float ni = Float.NEGATIVE_INFINITY;
        int m = seq1.length();
        int n = seq2.length();
        // Initialize the dynamic programming, gap1, and gap2 matrices
        float[][] dp = new float[m + 1][n + 1];
        float[][] gap1 = new float[m + 1][n + 1];
        float[][] gap2 = new float[m + 1][n + 1];
        // Initialize the first column of the dp and gap1 matrices
        rr(dp, 0, ni);
        rc(dp, 0, ni);
        rr(gap1, 0, ni);
        rc(gap2, 0, ni);
        // in the first column at the ith row will start w gap penalty starting with row
        // 1

        for (int i = 1; i <= m; i++) {
            gap1[i][0] = GAP_OPEN + (i * GAP_EXTEND);
        }
        // Initialize the first row of the dp and gap2 matrices
        for (int j = 1; j <= n; j++) {
            gap2[0][j] = GAP_OPEN + (j * GAP_EXTEND);
        }

        // Calculate the alignment scores using dynamic programming

        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                // Check if the characters match or mismatch, and update the score accordingly
                float match = seq1.charAt(i - 1) == seq2.charAt(j - 1) ? MATCH : MISMATCH;

                // dp[i][j] = dp[i - 1][j - 1] + match;
                // Update gap1 and gap2 matrices using the previously calculated values and
                // penalties
                gap1[i][j] = Math.max(GAP_OPEN + GAP_EXTEND + dp[i - 1][j],
                        Math.max(GAP_EXTEND + gap1[i - 1][j], GAP_OPEN + GAP_EXTEND + gap2[i - 1][j]));
                gap2[i][j] = Math.max(GAP_OPEN + GAP_EXTEND + dp[i][j - 1],
                        Math.max(GAP_OPEN + GAP_EXTEND + gap1[i][j - 1], GAP_EXTEND + gap2[i][j - 1]));
                // Update the current cell of the dp matrix with the maximum score from the
                // three matrices

                dp[i][j] = Math.max(dp[i - 1][j - 1], Math.max(gap1[i][j], gap2[i][j])) + match;
            }
        }
        // Create StringBuilders to store the aligned sequences
        StringBuilder alignedSeq1 = new StringBuilder();
        StringBuilder alignedSeq2 = new StringBuilder();
        // Traceback through the dp matrix to build the aligned sequences

        int i = m, j = n;

        while (i > 0 || j > 0) {
            float maxVal = Math.max(dp[i][j], Math.max(gap1[i][j], gap2[i][j]));

            if (i > 0 && j > 0
                    && maxVal == dp[i - 1][j - 1] + (seq1.charAt(i - 1) == seq2.charAt(j - 1) ? MATCH : MISMATCH)) {
                alignedSeq1.append(seq1.charAt(i - 1));
                alignedSeq2.append(seq2.charAt(j - 1));
                i--;
                j--;
            } else if (i > 0 && maxVal == gap1[i][j]) {
                alignedSeq1.append(seq1.charAt(i - 1));
                alignedSeq2.append('_');
                i--;
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