package sjtools;

import sjtools.CalcNumerical.Convergence;
import sjtools.CalcNumerical.Norm;


public class Start {
	public static void main(String[] args) {
		double[][] a = {{6,1,1,1,0},
						{1,7,1,1,1},
						{0,1,8,1,1},
						{0,0,1,9,1},
						{0,0,0,0,10}};
		double[] b = {9,11,11,11,11};
		double[] x = {0,0,0,0,0};
		double[][] a_ = {{4,-2,0},
						 {-2,4,4},
						 {4,-2,4}};
		double[] b_ = {1,1,1};
		double[] x_ = {10,10,0};
		NumericalData sordata = new NumericalData(a, b);
		sordata.setEps(1.0E-10);
		sordata.setMaxN(100);
		sordata.setInitX(x);
		sordata.setCon(Convergence.RELATIVEERROR);
		sordata.setNorm(Norm.INFINITY);
		CalcTool.printVec(CalcNumerical.jacobi(sordata));
		CalcTool.printVec(CalcNumerical.gaussSeidel(sordata));
		System.out.println("SOR法");
		CalcTool.printVec(CalcNumerical.sor(sordata,1.0));
		
	}
	
}
