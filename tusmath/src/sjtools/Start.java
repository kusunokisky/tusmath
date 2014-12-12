package sjtools;

import sjtools.CalcNumerical.Convergence;
import sjtools.CalcNumerical.Norm;



public class Start {
	public static void main(String[] args) {
		int n = 100;
		double d[][] = new double[n][n];
		for(int i = 0;i < n;i++){
			for(int j = 0;j < n; j++){
				if(i == j){
					d[i][j] = i+1;
				}
			}
		}
		System.out.println("[1.1]");
		System.out.println("||D||1 = " + Calc.matNorm1(d));
		System.out.println("||D||2 = " + 100);
		System.out.println("||D||∞ = " + Calc.matNormInf(d));
		
		double[] y = new double[n];
		for(int i = 0;i < n;i++){
			y[i] = i + 1;
		}
		System.out.println("[1.2]");
		System.out.println( Calc.vecNormP(y, 3) );
		
		n = 5;
		double[][] a = new double[n][n];
		double[] b = new double[n];
		for(int i = 0;i < n;i++){
			for(int j = 0;j < n;j++){
				if(i == j){
					a[i][j] = i+1;
				}
				else{
					a[i][j] = 1;
				}
			}
			b[i] = i+5;
		}
		a[0][0] = 1.0E-4;
		b[0] = 4.0001;
		double[] x = new double[n];
		for(int i = 0;i < n;i++){
			x[i] = 1;
		}
		System.out.println("[2.1]");
		CalcTool.printVec(x);
		System.out.println("[2.2]");
		System.out.println( Calc.vecNorm2( Calc.subVec(x, CalcNumerical.gauss(a, b)) ) );
		System.out.println("[2.3]");
		System.out.println( Calc.vecNorm2( Calc.subVec(x, CalcNumerical.partialPivotGauss(a, b)) ) );
		n = 12;
		double[][] h = CalcTool.createHilbert(12);
		double[] hx = new double[n];
		double[] hb = new double[n];
		for(int i = 0; i < n;i++){
			for(int j = 0; j < n;j++){
				hb[i] += h[i][j];
			}
			hx[i] = 1;
		}
		NumericalData hdata = new NumericalData(h,hb);
		hdata.setNorm(Norm.INFINITY);
		
		System.out.println("[3.1]");
		System.out.println(Calc.conditionNumber(hdata.getA(),hdata.getNorm()));
		System.out.println("[3.2]");
		double[] hx_ = CalcNumerical.partialPivotGauss(hdata.getA(), hdata.getB());
		hdata.setCon(Convergence.RELATIVERESIDUAL);
		System.out.println(CalcTool.calcConvergence(hdata, hb, Calc.matVec(h, hx_)));
		System.out.println("[3.3]");
		hdata.setCon(Convergence.RELATIVEERROR);
		System.out.println(CalcTool.calcConvergence(hdata, hx, hx_));
		System.out.println("[3.4");
		System.out.println("理論的なことは省略");
		
		System.out.println("[4]");
		System.out.println("プログラムを参照");
		n = 200;
		double[][] ta1 = new double[n][n];
		double[] tb1 = new double[n];
		for(int i = 0; i < n;i++){
			for(int j = 0; j < n;j++){
				if(i == j){
					ta1[i][j] = -1;
				}else if(i == j-1){
					ta1[i][j] = -2;
				}else if(i == j+1){
					ta1[i][j] = 1;
				}
			}
			tb1[i] = 1;
		}
		double[] tx1 = new double[n];
		NumericalData t1data = new NumericalData(ta1, tb1);
		t1data.setInitX(tx1);
		t1data.setEps(1.0E-8);
		t1data.setMaxN(500);
		t1data.setCon(Convergence.RELATIVERESIDUAL);
		t1data.setNorm(Norm.INFINITY);
		System.out.println("[5.1]");
		for(int i = 1;i < 20 ; i++){
			CalcNumerical.sor(t1data, i / 10.0);
			System.out.println("ω = " + (i / 10.0) + " ; 反復回数 : " + t1data.getCount());
		}
		n = 3;
		double[][] ta2 = new double[n][n];
		double[] tb2 = new double[n];
		for(int i = 0; i < n;i++){
			for(int j = 0; j < n;j++){
				if(i == j){
					ta2[i][j] = -1;
				}else if(i == j-1){
					ta2[i][j] = -2;
				}else if(i == j+1){
					ta2[i][j] = 1;
				}
			}
			tb2[i] = 1;
		}
		double[] tx2 = new double[n];
		NumericalData t2data = new NumericalData(ta2, tb2);
		t2data.setInitX(tx2);
		t2data.setEps(1.0E-8);
		t2data.setMaxN(500);
		t2data.setCon(Convergence.RELATIVERESIDUAL);
		t2data.setNorm(Norm.INFINITY);
		System.out.println("[5.2]");
		CalcNumerical.jacobi(t2data);
		System.out.println(t2data.getCount());
		System.out.println("理論的なのことは省略");
		
		n = 100;
		double[][] la = new double[n][n];
		for(int i = 0; i < n;i++){
			for(int j = 0; j < n;j++){
				if(i == j){
					la[i][j] = 2;
				}else if(i == j-1){
					la[i][j] = 1;
				}else if(i == j+1){
					la[i][j] = 1;
				}
			}
		}
		double lsum = 1;
		la = CalcNumerical.luDecomposition(la);
		for(int i = 0; i < n;i++){
			lsum *= la[i][i];
		}
		System.out.println("[6]");
		System.out.println(lsum);

	}
	
}
