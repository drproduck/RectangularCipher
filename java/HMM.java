import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.Random;

/**
 * Created by drproduck on 9/3/17.
 */
public class HMM {

    public double[][][] digamma;
    public double[][] A, B, ap, bt, gamma;
    public double[] c, pi;
    public int[] seq;

    public int len, nt, nw;
    public double oldLogProb;
    public double logProb;
    public boolean firstRestart;
    static double eps = 0.1;

    public HMM(int n_states, int n_obs, int[] observations) {
        this(n_states, n_obs, observations.length);
        seq = observations;
    }

    public HMM(int n_states, int n_obs, int time) {
        len = time;
        this.nt = n_states;
        this.nw = n_obs;
        A = new double[nt][nt];
        B = new double[nt][nw];
        pi = new double[nt];
        ap = new double[len][nt];
        bt = new double[len][nt];
        c = new double[len];
        gamma = new double[len][nt];
        digamma = new double[len][nt][nt];
        firstRestart = true;
    }

    public void alpha() {
        c[0] = 0;
        for (int i = 0; i < nt; i++) {
            ap[0][i] = pi[i] * B[i][seq[0]];
            c[0] += ap[0][i];
        }

        c[0] = 1 / c[0];
        for (int i = 0; i < nt; i++) {
            ap[0][i] = c[0] * ap[0][i];
        }

        for (int t = 1; t < len; t++) {
            c[t] = 0;
            for (int i = 0; i < nt; i++) {
                ap[t][i] = 0;
                for (int j = 0; j < nt; j++) {
                    ap[t][i] += ap[t - 1][j] * A[j][i];
                }
                ap[t][i] = ap[t][i] * B[i][seq[t]]; // p(h_t = i, O_t | O_1^{t-1})
                c[t] += ap[t][i]; // p(O_t | O_1^{t-1}) probability of observation at time t given previous observations
            }
            //scale
            c[t] = 1 / c[t];
            for (int i = 0; i < nt; i++) {
                ap[t][i] *= c[t]; // p(h_t = i | O_1^t) probability of state i at time t given observations up to time t
            }
        }
    }

    public int[] viterbi() {
        double c = 0;

        int[][] dp = new int[len][nt];
        double[] head = new double[nt];
        double[] temp = new double[nt];

        for (int i = 0; i < nt; i++) {
            head[i] = pi[i] * B[i][seq[0]];
            c += ap[0][i];
        }

        c = 1 / c;
        for (int i = 0; i < nt; i++) {
            head[i] *= c;
        }

        for (int t = 1; t < len; t++) {
            c = 0;
            for (int i = 0; i < nt; i++) {
                double x = 0;
                double max = Double.NEGATIVE_INFINITY;
                int lst = -1;
                for (int j = 0; j < nt; j++) {
                    x = head[j] * A[j][i];
                    if (x > max) {
                        max = x;
                        lst = j;
                    }
                }
                dp[t][i] = lst;

                temp[i] = max * B[i][seq[t]];
                c += temp[i];
            }
            //scale
            c = 1 / c;
            for (int i = 0; i < nt; i++) {
                temp[i] *= c;
            }

            //copy temp back to head
            for (int i = 0; i < nt; i++) {
                head[i] = temp[i];
            }
        }

        //get the last state of the best DP sequence
        int bs = -1;
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < nt; i++) {
            if (head[i] > max) {
                max = head[i];
                bs = i;
            }
        }
        //trace back
        int[] res = new int[len];
        res[len - 1] = bs;
        for (int t = len - 1; t > 0; t--) {
            bs = dp[t][bs];
            res[t - 1] = bs;
        }

        return res;
    }

    // best sequence in exclude.HMM sense
    public int[] besthmm() {
        alpha();
        beta();
        gamma();
        int[] res = new int[len];
        for (int t = 0; t < len; t++) {
            double max = Double.NEGATIVE_INFINITY;
            int p = -1;
            for (int i = 0; i < nt; i++) {
                if (gamma[t][i] > max) {
                    max = gamma[t][i];
                    p = i;
                }
            }
            res[t] = p;
        }
        return res;
    }

    public void beta() {
        // Let βT −1(i) = 1, scaled by cT −1
        for (int i = 0; i < nt; i++) {
            bt[len - 1][i] = c[len - 1];
        }
        for (int t = len - 2; t >= 0; t--) {
            for (int i = 0; i < nt; i++) {
                bt[t][i] = 0;
                for (int j = 0; j < nt; j++) {
                    bt[t][i] = bt[t][i] + A[i][j] * B[j][seq[t + 1]] * bt[t + 1][j];
                }
                // scale βt(i) with same scale factor as αt(i)
                bt[t][i] = c[t] * bt[t][i];
            }
        }
    }

    public void gamma() {
        for (int t = 0; t < len - 1; t++) {
            double denom = 0;
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    denom = denom + ap[t][i] * A[i][j] * B[j][seq[t + 1]] * bt[t + 1][j];
                }
            }
            for (int i = 0; i < nt; i++) {
                gamma[t][i] = 0;
                for (int j = 0; j < nt; j++) {
                    digamma[t][i][j] = (ap[t][i] * A[i][j] * B[j][seq[t + 1]] * bt[t + 1][j]) / denom;
                    gamma[t][i] = gamma[t][i] + digamma[t][i][j];
                }
            }
        }

        // special case
        double denom = 0;
        for (int i = 0; i < nt; i++) {
            denom = denom + ap[len - 1][i];
        }
        for (int i = 0; i < nt; i++) {
            gamma[len - 1][i] = ap[len - 1][i] / denom;
        }
    }

    //reEstimate the model only after alpha, beta and gamma passes
    public void reEstimate(boolean trainA) {
        if (firstRestart) {
            alpha();
            logProbFromAlpha();
            oldLogProb = logProb;
            firstRestart = false;
        }
        beta();
        gamma();
        for (int i = 0; i < nt; i++) {
            pi[i] = gamma[0][i];
        }

        if (trainA) {
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    double numer = 0;
                    double denom = 0;
                    for (int t = 0; t < len - 1; t++) {
                        numer = numer + digamma[t][i][j];
                        denom = denom + gamma[t][i];
                    }
                    A[i][j] = numer / denom;
                }
            }
        }

        // inefficient
//        for (int i = 0; i < nt; i++) {
//            for (int j = 0; j < nw; j++) {
//                double numer = 0;
//                double denom = 0;
//                for (int t = 0; t < len; t++) {
//                    if (seq[t] == j) {
//                        numer = numer + gamma[t][i];
//                    }
//                    denom = denom + gamma[t][i];
//                }
//                B[i][j] = numer / denom;
//            }
//        }

        // add smoothing constant
        double[][] BB = new double[nt][nw];
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < nw; j++) {
                BB[i][j] = eps;
            }
        }
        for (int i = 0; i < nt; i++) {
            double denom = nw * eps;
            for (int t = 0; t < len; t++) {
                BB[i][seq[t]] += gamma[t][i];
                denom += gamma[t][i];
            }
            for (int j = 0; j < nw; j++) {
                BB[i][j] /= denom;
            }
        }
        B = BB;

        alpha();
        logProbFromAlpha();
    }

    // compute logProb
    public void logProbFromAlpha() {
        logProb = 0;
        for (int i = 0; i < len; i++) {
            logProb = logProb + Math.log(c[i]);
        }
        logProb = -logProb;
    }

    //compute probability of seq
    public double score(int[] obs) {
        seq = obs;
        alpha();
        double p = 1;
        for (int i = 0; i < c.length; i++) {
            p *= c[i];
        }
        return 1.0 / p;
    }

    public void train(int iter, long seed, boolean trainA, boolean verbose) {
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

    public double jointProb(int[] states) throws Exception {
        if (states.length != len) throw new Exception("number of states must correspond to number of observations");
        double p = 1;
        p *= pi[states[0]] * B[states[0]][seq[0]];
        for (int t = 1; t < len; t++) {
            p *= A[states[t - 1]][states[t]] * B[states[t]][seq[t]];
        }
        return p;
    }

    // randomly initialize pi, A and B so that the rows sum to 1

    public void init(long seed, boolean initA) {
        Random rng = new Random(seed);
        double rem = 1.0;
        for (int i = 0; i < nt - 1; i++) {
            pi[i] = rng.nextDouble() * rem;
            rem = rem - pi[i];
        }
        pi[nt - 1] = rem;

        if (initA) {
            for (int i = 0; i < nt; i++) {
                rem = 1.0;
                for (int j = 0; j < nt - 1; j++) {
                    A[i][j] = rng.nextDouble() * rem;
                    rem = rem - A[i][j];
                }
                A[i][nt - 1] = rem;
            }
        }

        for (int i = 0; i < nt; i++) {
            rem = 1.0;
            for (int j = 0; j < nw - 1; j++) {
                B[i][j] = rng.nextDouble() * rem;
                rem = rem - B[i][j];
            }
            B[i][nw - 1] = rem;
        }
    }


    public static void main(String[] args) throws Exception {
        HMM test = new HMM(2, 3, 3);
        test.seq = new int[]{1, 0, 2};
        test.A = new double[][]{{0.7, 0.3}, {0.4, 0.6}};
        test.B = new double[][]{{0.1, 0.4, 0.5}, {0.7, 0.2, 0.1}};
        test.pi = new double[]{0.0, 1.0};
        test.alpha();

        System.out.println(Arrays.toString(test.viterbi()));
        System.out.println(Arrays.toString(test.besthmm()));
        System.out.println(test.sc);

        int[] txt_obs = util.textProcess("/home/drproduck/Documents/hmm/src/brown_nolines.txt", false, 0, 100000);
        System.out.println("text processed successfully");
        HMM trainer = new HMM(26, 27, txt_obs);
        System.out.println("declare successfully");
        //System.out.println(Arrays.deepToString(trainer.A));
        //System.out.println(Arrays.deepToString(trainer.B));
        trainer.train(100, 5678, true, true);
        BufferedWriter outputWriter = new BufferedWriter(new FileWriter("output26.txt"));
        outputWriter.write(Arrays.deepToString(trainer.A));
        outputWriter.newLine();
        outputWriter.write(Arrays.deepToString(trainer.B));
        outputWriter.newLine();
        outputWriter.write(Arrays.toString(trainer.pi));
        outputWriter.newLine();
        outputWriter.flush();
        outputWriter.close();

    }

    double sc = 0;
}


