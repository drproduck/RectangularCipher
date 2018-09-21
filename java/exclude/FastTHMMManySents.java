import org.apache.commons.math3.util.FastMath;

import java.util.Random;

/**
 * Trigram exclude.HMM. a similar normalization is used instead of log-scale and the ensuing logSumExp.
 */

public class FastTHMMManySents {

    double[][][] A;
    double[][] B;
    int nt;
    int nw;
    double[][][] ap, bt;
    double[] c;
    double[][][] digamma;
    double[][][][] trigamma;
    double[][] gamma;
    int[][] sentences;
    double logProb;
    double oldLogProb;
    boolean firstRestart;
    static double eps = 0.1; // smoothing constant

    public FastTHMMManySents(double[][][] transition, int nTag, int nWord, int[][] sentences) {
        nt = nTag;
        nw = nWord;
        A = transition;
        B = new double[nt][nw];
        ap = new double[len - 1][nt][nt];
        bt = new double[len - 1][nt][nt];
        c = new double[len - 1];
        digamma = new double[len][nt][nt];
        trigamma = new double[len][nt][nt][nt];
        gamma = new double[len][nt];
        this.sentences = sentences;

        firstRestart = true;
    }

//    public getSentencesStatistics(int seq){
//
//    }

    private void alpha(int[] seq) {
        int len = seq.length;
        c[0] = 0;
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < nt; j++) {
                ap[0][i][j] = B[i][seq[0]] * B[j][seq[1]];
                c[0] += ap[0][i][j];
            }
        }
        c[0] = 1 / c[0];
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < nt; j++) {
                ap[0][i][j] *= c[0];
            }
        }
        for (int t = 1; t < len - 1; t++) {
            c[t] = 0;
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    ap[t][i][j] = 0;
                    for (int k = 0; k < nt; k++) {
                        ap[t][i][j] += ap[t - 1][k][i] * A[k][i][j] * B[j][seq[t + 1]];
                    }
                    c[t] += ap[t][i][j];
                }
            }
            c[t] = 1 / c[t];
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    ap[t][i][j] *= c[t]; // p(O0,O1,x0,x1) / p(O0,O1) = p(x0,x1|O0,O1)
                }
            }
        }
    }

    public double logProbFromAlpha() {
        logProb = 0;
        for (int i = 0; i < len - 1; i++) {
            logProb -= FastMath.log(c[i]);
        }
        return logProb;
    }

    public void beta(int[] seq) {
        int len = seq.length;
        // nominal code
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < nt; j++) {
                bt[len - 2][i][j] = 1;
            }
        }
        for (int t = len - 3; t >= 0; t--) {
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    bt[t][i][j] = 0;
                    for (int k = 0; k < nt; k++) {
                        bt[t][i][j] += bt[t + 1][j][k] * A[i][j][k] * B[k][seq[t + 2]];
                    }
                    bt[t][i][j] *= c[t + 1]; // p(O2|x0,x1) / p(O2|O0,O1)
                }
            }
        }
    }

    public void digamma() {
        for (int t = 0; t < len - 1; t++) {
            for (int i = 0; i < nt; i++) {
                gamma[t][i] = 0;
                for (int j = 0; j < nt; j++) {
                    digamma[t][i][j] = ap[t][i][j] * bt[t][i][j]; // p(O0,01,O2,x0,x1) / p(O0,O1,O2) = p(x0,x1|O)
                    gamma[t][i] += digamma[t][i][j]; // p(x0|O) = \sum_{x1} p(x0,x1|O)
                }
            }
        }
        for (int j = 0; j < nt; j++) {
            gamma[len - 1][j] = 0;
            for (int i = 0; i < nt; i++) {
                gamma[len - 1][j] += digamma[len - 2][i][j]; // p(x2|O) = \sum_{x1} p(x1,x2|O)
            }
        }
    }

    public void trigamma() {
        for (int t = 0; t < len - 2; t++) {
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    for (int k = 0; k < nt; k++) {
                        trigamma[t][i][j][k] = ap[t][i][j] * A[i][j][k] * B[k][seq[t + 2]] * bt[t + 1][j][k] * c[t + 1];
                        // (p(O0,O1,x0,x1)/p(O0,O1)) * p(x2|x0,x1) * p(O2|x2) * (p(O3|x1,x2)/p(O3|O1,O2)) / p(O2|O0,O1)
                    }
                }
            }
        }
    }

    public void reEstimate(boolean trainA) {
        if (firstRestart) {
            alpha();
            logProbFromAlpha();
            oldLogProb = logProb;
            firstRestart = false;
        }
        beta();
        digamma();
        if (trainA) {
            trigamma();
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    double denom = 0;
                    for (int k = 0; k < nt; k++) {
                        A[i][j][k] = 0;
                        for (int t = 0; t < len - 2; t++) {
                            A[i][j][k] += trigamma[t][i][j][k];
                        }
                        denom += A[i][j][k];
                    }
                    for (int k = 0; k < nt; k++) {
                        A[i][j][k] /= denom;
                    }
                }
            }
        }
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < nw; j++) {
                double denom = eps * nw;
                double numer = 0;
                for (int t = 0; t < len; t++) {
                    if (seq[t] == j) {
                        numer += gamma[t][i];
                    }
                    denom += gamma[t][i];
                }
                numer += eps;
                B[i][j] = numer / denom;
            }
        }
        alpha();
        logProbFromAlpha();
    }

    public void train(int iter, boolean trainA, long seed, boolean verbose) {
        init(seed, trainA);
        long start;
        long stop;
        for (int it = 0; it < iter; it++) {
            if (verbose) System.out.printf("iteration: %d\n", it);
            start = System.nanoTime();
            reEstimate(trainA);
            stop = System.nanoTime();
            if (verbose) {
                System.out.println(logProb);
                System.out.printf("time: %f\n", (stop - start) / 1e9);
            }
            if (it != 0 && logProb <= oldLogProb) {
                System.out.printf("something may be wrong: logprob = %f, oldlogprob = %f\n", logProb, oldLogProb);
                break;
            }
        }
    }

    public void init(long seed, boolean initA) {
        Random rng = new Random(seed);
        double s;
        for (int i = 0; i < nt; i++) {
            s = 0;
            for (int j = 0; j < nw; j++) {
                B[i][j] = rng.nextDouble();
                s += B[i][j];
            }
            for (int j = 0; j < nw; j++) {
                B[i][j] /= s;
            }
        }
        if (initA) {
            A = new double[nt][nt][nt];
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    s = 0;
                    for (int k = 0; k < nt; k++) {
                        A[i][j][k] = rng.nextDouble();
                        s += A[i][j][k];

                    }
                    for (int k = 0; k < nt; k++) {
                        A[i][j][k] /= s;
                    }
                }
            }
        }
    }

    public int[] viterbi() {
        double c = 0;

        int[][][] dp = new int[len - 1][nt][nt];
        double[][] head = new double[nt][nt];
        double[][] temp = new double[nt][nt];

        // i -> j
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < nt; j++) {
                head[i][j] = ap[0][i][j];
            }
        }

        for (int t = 0; t < len - 2; t++) {
            c = 0;
            // i -> (j -> k), fix j and k, search over i
            for (int k = 0; k < nt; k++) {
                for (int j = 0; j < nt; j++) {
                    double x = 0;
                    double max = Double.NEGATIVE_INFINITY;
                    int argmax = -1;
                    for (int i = 0; i < nt; i++) {
                        x = head[i][j] * A[i][j][k];
                        if (x > max) {
                            max = x;
                            argmax = i;
                        }
                    }
                    dp[t][j][k] = argmax;
                    temp[j][k] = max * Math.pow(B[k][seq[t + 2]], 3);
                    c += temp[j][k];
                }
            }
            //scale and move back to head
            c = 1 / c;
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    head[i][j] = c * temp[i][j];
                }
            }
        }

        //get the last state of the best DP sequence
        int argmaxi = -1;
        int argmaxj = -1;
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < nt; j++) {
                if (head[i][j] > max) {
                    max = head[i][j];
                    argmaxi = i;
                    argmaxj = j;
                }
            }
        }
        //trace back
        int[] resn = new int[len];
        resn[len - 2] = argmaxi;
        resn[len - 1] = argmaxj;
        for (int t = len - 3; t >= 0; t--) {
            argmaxi = dp[t][argmaxi][argmaxj];
            argmaxj = argmaxi;
            resn[t] = argmaxi;
        }
        char[] resc = new char[len];
        for (int i = 0; i < len; i++) {
            resc[i] = (char) (resn[i] + 65);
        }
        return resn;
    }

    public static void main(String[] args) throws Exception {
        int[] pt = util.plain("/home/drproduck/Documents/exclude.HMM/data/408plaincleaned");
        int[] ct = util.read408("/home/drproduck/Documents/exclude.HMM/data/408ciphercleaned");
        int nTag = 26;
        int nWord = 54;
        int len = 408;

        double[][][] trigram = util.readTrigram("/home/drproduck/Documents/exclude.HMM/data/trigram.txt");
        long start = System.nanoTime();
        FastTHMMManySents FHMM = new FastTHMMManySents(trigram, nTag, nWord, ct);
        long stop = System.nanoTime();
        FHMM.train(200, false, -3076155353333121539L, true);
//        FHMM.train(200, false, 8781939572407739913L, true);
//        FHMM.train(200, false, 1209845257843231593L, true);
//        FHMM.train(200, false, 3738420990656387694L, true);

        stop = System.nanoTime();
        System.out.printf("training time: %f\n", (stop - start) / 1e9);

        /*
        Use Viterbi decode
         */
//        int[] sol = FHMM.viterbi();
//        double acc = 0;
//        for (int i = 0; i < sol.length; i++) {
//            if (sol[i] == pt[i]) {
//                acc += 1;
//            }
//        }
//        System.out.printf("accuracy = %f\n", acc / len);
//        StringBuilder sb = new StringBuilder();
//        for (int i = 0; i < sol.length; i++) {
//            sb.append(Character.toString((char) (sol[i] + 65)));
//        }
//        System.out.println(sb.toString());
//        System.out.println(util.display2(FHMM.B));
//        BufferedWriter out = new BufferedWriter(new FileWriter("B_1236"));
//        for (int i = 0; i < 26; i++) {
//            for (int j = 0; j < 54; j++) {
//                out.write(String.format("%f ", FHMM.B[i][j]));
//            }
//            out.newLine();
//        }
//        out.close();

        /*
        Use highest-by-column decode
         */
        int[] argmax = new int[FHMM.nw];
                for (int k = 0; k < FHMM.nw; k++) {
                    double max = FHMM.B[0][k];
                    int arg = 0;
                    for (int l = 1; l < FHMM.nt; l++) {
                        if (max < FHMM.B[l][k]) {
                            max = FHMM.B[l][k];
                            arg = l;
                        }
                    }
                    argmax[k] = arg;
                }
                int[] sol_num = new int[len];
                double acc = 0;
                for (int k = 0; k < len; k++) {
                    sol_num[k] = argmax[ct[k]];
                    if (sol_num[k] == pt[k]) {
                        acc++;
                    }
                }
                acc /= len;
        System.out.println(acc);
    }
}
