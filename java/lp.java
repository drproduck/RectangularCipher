import scpsolver.constraints.LinearBiggerThanEqualsConstraint;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

import java.util.Arrays;

public class lp {
    public static void main(String[] args) {
        double t0 = System.nanoTime();
        LinearProgram lp = new LinearProgram(new double[]{5.0, 10.0});
        lp.addConstraint(new LinearBiggerThanEqualsConstraint(new double[]{3.0, 1.0}, 8.0, "c1"));
        lp.addConstraint(new LinearBiggerThanEqualsConstraint(new double[]{0.0, 4.0}, 4.0, "c2"));
        lp.addConstraint(new LinearSmallerThanEqualsConstraint(new double[]{2.0, 0.0}, 2.0, "c3"));
        lp.setMinProblem(true);
        LinearProgramSolver solver = SolverFactory.newDefault();
        double[] sol = solver.solve(lp);
        System.out.println(Arrays.toString(sol));
        System.out.println(System.nanoTime() - t0);
    }
}
