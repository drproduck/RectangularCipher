

public class TransitionKernel {

    public TransitionKernel(int nHidden, int nObserved, int[] plainSeq, int[] cipherSeq){
        int[][] cache = new int[nHidden][nObserved];
        assert(plainSeq.length == cipherSeq.length);
    }

}
