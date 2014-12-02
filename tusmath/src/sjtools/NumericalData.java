package sjtools;

import sjtools.CalcNumerical.Convergence;
import sjtools.CalcNumerical.Norm;

class NumericalData {
	private double[][] a;
	private double[] initX;
	private double[] b;
	private double eps;
	private int maxN;
	private Convergence con;
	private Norm norm;
	
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

}
