import org.jblas.DoubleMatrix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LPSolution;
import scpsolver.problems.LPWizard;
import scpsolver.problems.LPWizardConstraint;


public class TravelingSalesmanSolver {
    public int n;
    public LPSolution lastSolution = null;
    private LinearProgramSolver solver  = SolverFactory.newDefault();

    private DoubleMatrix addDummyNode(double[][] arr){
        DoubleMatrix A = DoubleMatrix.zeros(this.n, this.n);
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                A.put(i,j,arr[i][j]);
            }
        }
        return A;
    }

    public double solve(int n, double[][] arr){
        n = n + 1;
        this.n = n;

        DoubleMatrix A = addDummyNode(arr);

        String var;
        LPWizard lp = new LPWizard();

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                var = "x"+" "+i+" "+j;
                lp.plus(var,A.get(i,j));
                lp.setBoolean(var);
            }
        }

        LPWizardConstraint c;

        // constraint for each row
        for (int i = 0; i < n; i++) {
            c = lp.addConstraint("row sum "+i,1,"=");
            for (int j = 0; j < n; j++) {
                var = "x"+" "+i+" "+j;
                c.plus(var, 1);
            }
        }

        // constraint for each column
        for (int i = 0; i < n; i++) {
            c = lp.addConstraint("column sum "+i,1,"=");
            for (int j = 0; j < n; j++) {
                var = "x"+" "+j+" "+i;
                c.plus(var,1);
            }
        }

        addSubtourConstraint(lp);

        lp.setMinProblem(false);

        lastSolution = lp.solve(solver);
        return lastSolution.getObjectiveValue();
    }

    private void addSubtourConstraint(LPWizard lp){
        // bounds for dummy variables

        for (int i = 1; i < n; i++) {
            String ui = "u"+" "+i;
            lp.addConstraint("lower bound for u"+i,2,"<=").plus(ui,1);
            lp.addConstraint("upper bound for u"+i,n,">=").plus(ui,1);
            lp.setInteger(ui);
        }

        // dummy variable constraint to eliminate subtour
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                String ui = "u"+" "+i;
                String uj = "u"+" "+j;
                String var = "x"+" "+i+" "+j;
                lp.addConstraint("dummy pair"+i+j,n-2,">=").plus(ui,1).plus(uj,-1).plus(var,n-1);
            }
        }

    }

    public int[] getTour(){
        Map<Integer, Integer> dict = new HashMap<>();

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (lastSolution.getInteger("x"+" "+i+" "+j) == 1){
                    dict.put(i, j);
                }
            }
        }
        int[] tour = new int[n-1];
        int start = n-1;
        for (int i = 0; i < n-1; i++) {
            tour[i] = dict.get(start);
            start = tour[i];
        }

        return tour;
    }

    public double getScore() {
        return lastSolution.getObjectiveValue();
    }

    public static void main(String[] args) throws IOException {
        double t0 = System.nanoTime();
        int n = 17;

        double[][] arr = new double[n][n];

        File file = new File("distance");

        Scanner br = new Scanner(new BufferedReader(new FileReader(file)));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                arr[i][j] = br.nextDouble();
            }
        }

        TravelingSalesmanSolver solver = new TravelingSalesmanSolver();
        System.out.println(solver.solve(17, arr));
        System.out.println(Arrays.toString(solver.getTour()));
        System.out.println(System.nanoTime() - t0);

        System.out.println(solver.solve(3, new double[][]{{-1,-2,-3},{-3,-4,-5},{-5,-6,-7}}));
        System.out.println(Arrays.toString(solver.getTour()));
        System.out.println(System.nanoTime() - t0);


    }

}
