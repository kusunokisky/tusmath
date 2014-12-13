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
	 * 使用する収束判定
	 */
	private Convergence con;
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
	void setA(double[][] a) {
		this.a = a;
	}
	double[] getInitX() {
		return initX;
	}
	void setInitX(double[] initX) {
		this.initX = initX;
	}
	double[] getB() {
		return b;
	}
	void setB(double[] b) {
		this.b = b;
	}
	double getEps() {
		return eps;
	}
	void setEps(double eps) {
		this.eps = eps;
	}
	int getMaxN() {
		return maxN;
	}
	void setMaxN(int maxN) {
		this.maxN = maxN;
	}
	Convergence getCon() {
		return con;
	}
	void setCon(Convergence con) {
		this.con = con;
	}
	Norm getNorm() {
		return norm;
	}
	void setNorm(Norm norm) {
		this.norm = norm;
	}
	int getCount() {
		return count;
	}
	void setCount(int count) {
		this.count = count;
	}

}
