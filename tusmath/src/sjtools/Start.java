package sjtools;



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
		NumericalData sordata = new NumericalData(a_, b_);
		double[][] ta = new double[20][20];
		int g = 100;
		for(int i = 0;i < ta.length;i++){
			for(int j = 0; j < ta[i].length;j++){
				if(i == j){
					ta[i][j] = 2;
				}else if((i-1) == j){
					ta[i][j] = g;
				}else if((i+1) == j){
					ta[i][j] = 1;
				}else{
					ta[i][j] = 0;
				}
			}
		}
		CalcTool.printMat(CalcNumerical.inverse(ta));
		System.out.println(Calc.matNorm1(ta));
		System.out.println(Calc.matNorm1(ta)*Calc.matNorm1(CalcNumerical.inverse(ta)));
	}
	
}
