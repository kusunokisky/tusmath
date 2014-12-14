package sjtools;
/**
 * 数値計算に必要な値をまとめるクラス
 * @author takeru
 *
 */
class NumericalData {
	/**
	 * 係数行列
	 */
	private double[][] a;
	/**
	 * 初期ベクトル
	 */
	private double[] initX;
	/**
	 * 右辺項
	 */
	private double[] b;
	/**
	 * 許容誤差
	 */
	private double eps;
	/**
	 * 最大反復回数
	 */
	private int maxN;
	/**
	 * 収束判定条件を満たすかどうか
	 */
	private boolean convergenceCriterion;
	/**
	 * 使用する収束判定条件
	 */
	private ConvergenceCriterion con;
	/**
	 * 使用するノルム
	 */
	private Norm norm;
	/**
	 * 反復計算処理において反復回数を保持する
	 * 収束しない場合は -1 がはいる
	 */
	private int count = 0;

	public NumericalData(double[][] a,double[] b) {
		this.a = a;
		this.b = b;
	}
	double[][] getA() {
		return a;
	}
	/**
	 * 係数行列(Ax=bのA)を設定する
	 * @param a
	 */
	void setA(double[][] a) {
		this.a = a;
	}
	double[] getInitX() {
		return initX;
	}
	/**
	 * 初期ベクトルを設定する
	 * @param initX
	 */
	void setInitX(double[] initX) {
		this.initX = initX;
	}
	double[] getB() {
		return b;
	}
	/**
	 * 右辺項(Ax=bのb)を設定する
	 * @param b
	 */
	void setB(double[] b) {
		this.b = b;
	}
	double getEps() {
		return eps;
	}
	/**
	 * 許容誤差を設定する
	 * @param eps
	 */
	void setEps(double eps) {
		this.eps = eps;
	}
	int getMaxN() {
		return maxN;
	}
	/**
	 * 最大反復回数を設定する
	 * @param maxN
	 */
	void setMaxN(int maxN) {
		this.maxN = maxN;
	}
	ConvergenceCriterion getCon() {
		return con;
	}
	/**
	 * 収束判定条件(誤差,残差,相対誤差,相対残差)を設定する
	 * @param con
	 */
	void setCon(ConvergenceCriterion con) {
		this.con = con;
	}
	Norm getNorm() {
		return norm;
	}
	/**
	 * 使用するノルムの種類を設定する
	 * @param norm
	 */
	void setNorm(Norm norm) {
		this.norm = norm;
	}
	int getCount() {
		return count;
	}
	/**
	 * 反復回数を設定する
	 * @param count
	 */
	void setCount(int count) {
		this.count = count;
	}
	/**
	 * 収束判定条件を満たすかを返す
	 * @return
	 */
	boolean fullfilConvergence() {
		return convergenceCriterion;
	}
	/**
	 * 収束したか設定する
	 * @param convergence
	 */
	void setConvergence(boolean convergenceCriterion) {
		this.convergenceCriterion = convergenceCriterion;
	}


}
