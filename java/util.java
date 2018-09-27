import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by drproduck on 9/13/17.
 */
public class util {

    //to extract characters from the Brown corpus

    public static int[] textProcess(String path, boolean nospace, int n_skip, int max) throws IOException {
        BufferedReader rd = new BufferedReader(new FileReader(path));
        for (int i = 0; i < n_skip; i++) {
            rd.read();
        }
        int c = rd.read();
        List<Integer> list = new ArrayList<>();
        int iter = 0;
        boolean added = false;
        while (c != -1 && iter < max) {
            added = false;
            if (Character.isAlphabetic(c)) {
                if (Character.isUpperCase(c)) {
                    list.add(c - 65);
                }
                else list.add(c - 97);
                added = true;
            } else if (!nospace && Character.isSpaceChar(c)) {
                list.add(26);
                added = true;
            }
            c = rd.read();
            if (added) iter ++;
        }
        int[] res = new int[list.size()];
        for (int i = 0; i < res.length; i++) {
            res[i] = list.get(i);
        }
        return res;
    }

    // display cipher in a more decorous way

    public static String display2(double[][] array) {
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < array[0].length; j++) {
            sb.append("--------");
        }
        sb.append("\n");
        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[i].length; j++) {
                sb.append(String.format(" %.3f |", array[i][j]));
            }
            sb.append("\n");
            for (int j = 0; j < array[0].length; j++) {
                sb.append("--------");
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    // digram frequency extracted from the brown corpus, default is next 1000000 characters (skip the first 1000000)

    public static double[][] digram() throws IOException {
        int[] chars = util.textProcess("/home/drproduck/Documents/hmm/src/brown_nolines.txt", true, 1000000, 1000000);
        System.out.println(chars.length);
        double[][] digraphFreq = new double[26][26];
        int a = chars[0];
        int b;
        for (int i = 1; i < chars.length; i++) {
            b = chars[i];
            digraphFreq[a][b] += 1;
            a = b;
        }
        double rowSum;
        for (int i = 0; i < 26; i++) {
            rowSum = 0;
            for (int j = 0; j < 26; j++) {
                digraphFreq[i][j] += 5;
                rowSum += digraphFreq[i][j];
            }
            for (int j = 0; j < 26; j++) {
                digraphFreq[i][j] /= rowSum;
            }
        }
        return digraphFreq;
    }

    public static double[][][] trigram(int len, int skip) throws Exception {
        int[] chars = util.textProcess("/home/drproduck/Documents/hmm/src/brown_nolines.txt", true, skip, len);
        double[][][] trigramFreq = new double[26][26][26];
        int a = chars[0];
        int b = chars[1];
        int c;
        for (int i = 2; i < chars.length; i++) {
            c = chars[i];
            trigramFreq[a][b][c] +=1;
            a = b;
            b = c;
        }
        double s;
        for (int i = 0; i < 26; i++) {
            for (int j = 0; j < 26; j++) {
                s = 0;
                for (int k = 0; k < 26; k++) {
                    trigramFreq[i][j][k] += 1;
                    s += trigramFreq[i][j][k];
                }
                for (int k = 0; k < 26; k++) {
                    trigramFreq[i][j][k] /= s;
                }
            }
        }
        return trigramFreq;
    }

    public static double[] stationary_pi() throws Exception {
        int[] chars = util.textProcess("/home/drproduck/Documents/hmm/src/brown_nolines.txt", true, 1000000, 1000000);
        double[] pi = new double[26];
        double sum = 0.0;
        for (int i = 0; i < chars.length; i++) {
            pi[chars[i]] += 1.0;
            sum += 1.0;
        }
        for (int i = 0; i < pi.length; i++) {
            pi[i] /= sum;
        }
        double s = 0;
        for (int i = 0; i < pi.length; i++) {
            if (pi[i] == 0) {
                throw new Exception("a probability is 0");
            }
            s += pi[i];
        }
        if (Math.abs(s - 1) > 0.0001) throw new Exception(String.format("not normalized, s = %f\n", s));
        return pi;
    }

    //randomize a key for for simple substitution

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

    //simple substitution key score

    public static double keyScore(int[] dict, double[][] B) {
        double score = 0.0;
        int[] keys = new int[26];
        for (int i = 0; i < B.length; i++) {
            Pair<Double, Integer> x = arrayMax(B[i]);
            keys[i] = x.idx;
        }
        for (int i = 0; i < 26; i++) {
            if (keys[i] == dict[i]) {
                score ++;
            }
        }
        return score / 26;
    }

    //convenient method

    public static Pair<Double, Integer> arrayMax(double[] a) {
        double max = Double.NEGATIVE_INFINITY;
        int idx = -1;
        Pair<Double, Integer> res = new Pair<>();
        for (int i = 0; i < a.length; i++) {
            if (max < a[i]) {
                max = a[i];
                idx = i;
            }
        }
        res.val = max;
        res.idx = idx;
        return res;
    }

    static class Pair<A, B> {
        A val;
        B idx;
    }

    /// convenient method

    public static double[][] deepClone(double[][] a) {;
        int d1 = a.length;
        int d2 = a[0].length;
        double[][] res = new double[d1][d2];
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                res[i][j] = a[i][j];
            }
        }
        return res;
    }

    public static int[][] deepClone(int[][] a) {
        int d1 = a.length;
        int d2 = a[0].length;
        int[][] res = new int[d1][d2];
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                res[i][j] = a[i][j];
            }
        }
        return res;
    }

    //read file Z408.txt to integer array.

    public static int[] read408(String dir) throws FileNotFoundException {
        Scanner rd = new Scanner(new BufferedReader(new FileReader(dir)));
        List<Integer> l = new ArrayList();
        while (rd.hasNext()) {
            l.add(rd.nextInt()-1);
        }
        int[] res = new int[l.size()];
        for (int i = 0; i < l.size(); i++) {
            res[i] = l.get(i);
        }
        return res;
    }

    //read a plain text of 408 I found on the web

    public static int[] plain(String dir) throws IOException {
        BufferedReader rd = new BufferedReader(new FileReader(dir));
        int c = rd.read();
        List<Integer> l = new ArrayList<>();
        while (c != -1) {
            if (Character.isAlphabetic(c)) {
                if (Character.isUpperCase(c)) {
                    l.add(c - 65);
                }
            }
            c = rd.read();
        }
        int[] res = new int[l.size()];
        for (int i = 0; i < l.size(); i++) {
            res[i] = l.get(i);
        }
        return res;
    }

    //make the true key for homophonic substitution for scoring, not to be confused with making a test key. some checks are added
    //since some cipher character in 408 match with more than 1 key, the array of sets is needed

    public static Set<Integer>[] makeDict() throws Exception {
        int[] cipher = read408("408.txt");
        int[] plain = plain("plain.txt");
//        System.out.println(Arrays.toString(cipher));
//        System.out.println(Arrays.toString(plain));
        Set<Integer>[] dict = new Set[53];
        for (int i = 0; i < 53; i++) {
            dict[i] = new HashSet<Integer>();
        }
        if (cipher.length != 408 || plain.length != 408) throw new Exception("has to be 408 characters!");
        for (int i = 0; i < 408; i++) {
            dict[cipher[i]].add(plain[i]);
        }
        return dict;
    }

    public static void main(String[] args) throws Exception {
        Set<Integer>[] dict = makeDict();
        System.out.println(Arrays.deepToString(dict));

    }

    //scoring a 408 key. Here I chose to divide the score by 53: the number of cipher characters. It makes sense right?

    public static double homophonicScore(Set<Integer>[] dict, double[][] B) {
        int d1 = B.length;
        int d2 = B[0].length;
        double score = 0;
        for (int i = 0; i < d2; i++) {
            double max = Double.NEGATIVE_INFINITY;
            int idx = -1;
            for (int j = 0; j < d1; j++) {
                if (B[j][i] > max) {
                    max = B[j][i];
                    idx = j;
                }
            }
            if (dict[i].contains(idx)) score ++;
        }
        return score / 53.0;
    }

    public static double[] readPi(String dir) throws Exception {
        Scanner rd = new Scanner(new BufferedReader(new FileReader(dir)));
        double[] pi = new double[26];
        for (int i = 0; i < 26; i++) {
            pi[i] = rd.nextDouble();
        }
        return pi;
    }

    public static double[][] readBigram(String dir) throws Exception {
        Scanner rd = new Scanner(new BufferedReader(new FileReader(dir)));
        double[][] ret = new double[26][26];
        for (int i = 0; i < 26; i++) {
            for (int j = 0; j < 26; j++) {
                ret[i][j] = rd.nextDouble();
            }
        }
        return ret;
    }

    public static double[][][] readTrigram(String dir) throws Exception {
        Scanner rd = new Scanner(new BufferedReader(new FileReader(dir)));
        double[][][] ret = new double[26][26][26];
        for (int i = 0; i < 26; i++) {
            for (int j = 0; j < 26; j++) {
                for (int k = 0; k < 26; k++) {
                    ret[i][j][k] = rd.nextDouble();
                }
            }
        }
        return ret;
    }
}
