
/**
 * A Rectangular Cipher is a cipher text that fits into a rectangle. Example are the Zodiac 408 and 340
 * This class wraps around a 2d array, and supports column transposition.
 */

import org.apache.commons.math3.util.FastMath;
import org.jblas.util.Permutations;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.*;
import java.util.Random;

public class RectangularCipher {
    public int[][] oldmatrix;
    public int[][] matrix;
    int nrows;
    int ncols;
    long rng = 2018;

    FastTHMM hmmSolver;
    TravelingSalesmanSolver tsSolver = new TravelingSalesmanSolver();

    double[][] bigram;
    double[][][] trigram;

    public RectangularCipher(int[][] arr, double[][] bigram, double[][][] trigram){
        oldmatrix = arr.clone();
        nrows = oldmatrix.length;
        for (int i = 0; i < arr.length; i++) {
            assert(arr[0].length == arr[i].length);
        }
        ncols = arr[0].length;
        matrix = arr.clone();

        this.bigram = bigram;
        this.trigram =trigram;
    }

    public int[][] transpose(int[] cindices){
        assert(cindices.length == ncols);
        int[][] newarray = new int[nrows][ncols];
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                newarray[i][j] = matrix[i][cindices[j]];
            }
        }
        return newarray;
    }

    public int[] getCol(int c){
        int[] res = new int[nrows];
        for (int i = 0; i < nrows; i++) {
            res[i] = matrix[i][c];
        }
        return res;
    }

    public void setCol(int c, int[] val){
        assert(val.length == nrows);
        for (int i = 0; i < nrows; i++) {
            matrix[i][c] = val[i];
        }
    }

    public int[] flatten(){
        int[] res = new int[nrows * ncols];
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                res[i * ncols + j] = matrix[i][j];
            }
        }
        return res;
    }

    public void shapen(int[] arr){
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                matrix[i][j] = arr[i * ncols + j];
            }
        }
    }

    public double getCol2ColDistance(int i, int j){
        double logP = 0;
        for (int k = 0; k < nrows; k++) {
            logP += FastMath.log(bigram[matrix[k][i]][matrix[k][j]]);
        }
        return logP;
    }

    public double[][] getCol2ColDistanceMatrix(){
        double[][] distance = new double[ncols][ncols];
        for (int i = 0; i < ncols; i++) {
            for (int j = 0; j < ncols; j++) {
                distance[i][j] = getCol2ColDistance(i,j);
            }
        }
        return distance;
    }

    public double deSubstitution(int nHidden, int nObs){
        int[] seq = flatten();
        hmmSolver = new FastTHMM(trigram, nHidden, nObs, seq);
        hmmSolver.train(200, false,rng, false);
        int[] proposal = hmmSolver.viterbi();
        shapen(proposal);
        return hmmSolver.logProbFromAlpha();
    }

    public double deTransposition(){
        double[][] distance = getCol2ColDistanceMatrix();
        tsSolver.solve(ncols, distance);
        int[] tour = tsSolver.getTour();
        matrix = transpose(tour);
        return tsSolver.getScore();
    }

    public double[] iterativeDecipher(int iter){
        LinkedList<Double> logScore = new LinkedList<>();
        double currScore = Double.NEGATIVE_INFINITY;
         logScore.add(currScore);
        for (int i = 0; i < iter; i++) {
            System.out.println("Iteration "+i);
            oldmatrix = util.deepClone(matrix);
            System.out.print("Solving substitution...");
            deSubstitution(26, 26);
            currScore = logProb();
            System.out.println("Score = "+currScore);
            System.out.println(Arrays.deepToString(matrix));
//            if ( logScore.getLast() >= currScore){
//                matrix = oldmatrix;
//                break;
//            }
             logScore.add(currScore);
            System.out.print("Solving transposition...");
            deTransposition();
            currScore = logProb();
            System.out.println("Score = "+currScore);
            System.out.println(Arrays.deepToString(matrix));
//            if ( logScore.getLast() >= currScore){
//                matrix = oldmatrix;
//                break;
//            }
             logScore.add(currScore);
        }
        double[] arrayScores =  logScore.stream().mapToDouble(i->i).toArray();

        return arrayScores;
    }

    public static int[][] readFromFile(String cipherDir, int nrows, int ncols) throws FileNotFoundException {
        Scanner reader = new Scanner(new BufferedReader(new FileReader(cipherDir)));
        int[][] plain = new int[nrows][ncols];
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                plain[i][j] = reader.nextInt();
            }
        }
        int[][] cipher = new int[nrows][ncols];
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                plain[i][j] = reader.nextInt();
            }
        }
        return cipher;
    }

    public static int[] permutation(int range){
        int[] perm = Permutations.randomPermutation(range);
        return perm;
    }

    public static int[] permuteKey(String method, long seed) {
        int[] dict = new int[26];
        if (method.equals("random")) {
            Random rng = new Random(seed);
            LinkedList<Integer> temp = new LinkedList<>();
            for (int i = 0; i < 26; i++) {
                temp.push(i);
            }
            for (int i = 25; i >= 0; i--) {
                double r = rng.nextDouble();
                dict[i] = temp.remove((int) (r * (i + 1)));
            }
        } else if (method.equals("caesar")) {
            for (int i = 0; i < 26; i++) {
                dict[i] = (i + 3) % 26;
            }
        }
        return dict;
    }

    public double logProb(){
        int[] seq = flatten();
        double lp = 0;
        for (int i = 2; i < seq.length; i++) {
            lp += FastMath.log(trigram[seq[i-2]][seq[i-1]][seq[i]]);
        }
        return lp;
    }

    public double logProbArray(int[] seq){
        double lp = 0;
        for (int i = 2; i < seq.length; i++) {
            lp += FastMath.log(trigram[seq[i-2]][seq[i-1]][seq[i]]);
        }
        return lp;
    }
}
