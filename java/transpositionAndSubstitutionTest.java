import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Scanner;
import org.jblas.util.Permutations;
import org.w3c.dom.css.Rect;

public class transpositionAndSubstitutionTest {
    public static void main(String[] args) throws Exception {

        String cipherDir = "../data/simple_cipher/simple_408";
        String trigramDir = "../data/trigram.txt";
        String bigramDir = "../data/bigram";
        double[][][] trigram = util.readTrigram(trigramDir);
        double[][] bigram = util.readBigram(bigramDir);

        int nrows = 24;
        int ncols = 17;

        Scanner reader = new Scanner(new BufferedReader(new FileReader(cipherDir)));
//        int[][] plain = new int[nrows][ncols];
//        for (int i = 0; i < nrows; i++) {
//            for (int j = 0; j < ncols; j++) {
//                plain[i][j] = reader.nextInt();
//            }
//        }
        int[] plain = new int[nrows*ncols];
        for (int i = 0; i < nrows * ncols; i++) {
            plain[i] = reader.nextInt();
        }

        int[][] cipher = new int[nrows][ncols];
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                cipher[i][j] = reader.nextInt();
            }
        }
        int[] perm = Permutations.randomPermutation(17);
        int[][] permutedCipher = RectangularCipher.transpose(cipher, perm);

        RectangularCipher cipherSolver = new RectangularCipher(permutedCipher, bigram, trigram);
        System.out.println("plain text score = "+cipherSolver.logProbArray(plain));
        cipherSolver.iterativeDecipher(10);
    }
}
