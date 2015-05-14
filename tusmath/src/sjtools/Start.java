package sjtools;

public class Start {
	public static void main(String[] args) {
	double[][] a = {{2,3,5,7},{3,5,7,11},{5,7,11,13},{7,11,13,17}};
	double[][] a_ =CalcNumerical.tridiagonal(a);
	CalcTool.printMat(a_);
	}
}

