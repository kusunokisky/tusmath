package sjtools;

class Example {
	/*=================================================================================
	 * x : 誤差	x~ = x + Δx : 近似値
	 * Ax=b の真の解 x	AΔx =Δb
	 *Δx = A^-1 * Δb ⇒ ||Δx|| ≦||A^-1|| ||b||,||A|| ||x|| ≧||b||
	 *||Δx||/||x|| ≦||A|| ||A^-1|| ||Δb||/||b|| = κ(A)||Δb|| / ||b||
	 *=================================================================================
	 *||I|| = 1 ||A|| < 1
	 *⇒1/(1+||A||) ≦ ||(I±A)^-1||≦1/(1-||A||)
	 *---------------------------------------------------------------------------------
	 *X : Aの近似逆行列 R : I-XA
	 *||R||<1⇒||A^-1||≦||X|| / (1-||R||)
	 *=================================================================================
	 *(A+ΔA)(x+Δx) = (b + Δb),||ΔA|| ||A^-1|| < 1
	 *⇒||Δx||/||x||≦κ(A)(||ΔA||/||A|| + ||Δb||/||b||) / (1 - κ(A)||ΔA||/||A||)
	 *proof
	 *
	 *=================================================================================
	 *D : 対角行列  E : 狭義上三角行列  F : 狭義下三角行列
	 *A = D + E + Fに分離し
	 *Jacobi法
	 *x^(m+1) = -D^-1(E+F)x^m + D^-1 b
	 *Gauss-Seidel法
	 *x^(m+1) = D^-1(b-Ex^(m+1)-Fx^(m))
	 *SOR法
	 *x^(m+1) = (1-ω)x^m + ωx_^(m+1) (x_ はGauss-Seidel法と同様に求める) 
	 *---------------------------------------------------------------------------------
	 *jacobi法
	 *T = -D^-1(E+F) : 反復式
	 *ρ(T) > 1 ⇒ jacobi法は収束しない
	 *Gauss-Seidel法
	 *T = -(D+E)^-1 F : 反復式
	 *ρ(T) > 1 ⇒ Gauss-Seidel法は収束しない ρ(T):スペクトル半径 = max{|λ|}
	 *SOR法
	 *L = D^-1 E | U = D^-1 F
	 *T = (I+ωL)^-1{(1-ω)I - ωU} : 反復式
	 *ρ(T) ≧ |ω-1| が成立
	 *proof
	 *φ(λ) = det(λI - T)とする
	 *(I+ωL)は正則でdet()(I+ωL) = 1
	 *φ(λ) = det(I+ωL)・det(λI - T) = det{(λ + ω - 1)I + λωL + ωU}
	 *λ = 0とするとTの固有値λ_1,λ_2,....,λ_nに対して
	 *(-1)^nλ_1,λ_2,....,λ_n = det{(ω-1)I + ωU} = (ω-1)^n
	 *ρ(T) = max|λn|≧|ω-1|となる
	 *=================================================================================
	 *
	 *
	 */
	static void bTest(){
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
		hdata.setCon(ConvergenceCriterion.RELATIVERESIDUAL);
		System.out.println(CalcTool.calcConvergence(hdata, hx_, Calc.matVec(h, hx_)));
		System.out.println("[3.3]");
		hdata.setCon(ConvergenceCriterion.RELATIVEERROR);
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
		t1data.setCon(ConvergenceCriterion.RELATIVEERROR);
		t1data.setNorm(Norm.INFINITY);
		System.out.println("[5.1]");
		for(int i = 1;i < 20 ; i++){
			CalcNumerical.sor(t1data, i / 10.0);
			System.out.println("ω = " + (i / 10.0) + " ; 反復回数 : " +t1data.getCount());
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
		t2data.setCon(ConvergenceCriterion.RELATIVERESIDUAL);
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
	static void H24ATest(){
		int n = 20;
		double[][] a = new double[n][n];
		for(int i = 0; i < n;i++){
			for(int j = 0; j < n;j++){
				if(i==j){
					a[i][j] = 10;
				}else{
					a[i][j] = Math.sqrt(Math.pow(i+1, 2) + Math.pow(j+1,2)) / 5;
				}
			}
		}
		System.out.println(Calc.matNorm1(a));
		System.out.println(Calc.matNormInf(a));
		System.out.println(Calc.matNormFrobenius(a));
		double [][] b = CalcNumerical.luDecomposition(a);
		double sum =1;
		for(int i = 0; i < n;i++){
			sum *= b[i][i];
		}
		System.out.println(sum);
		double[] t = {1,1,1,1,1};
		System.out.println(Calc.vecNorm2(t));
		double[][] a5 = { {3,-2,0},
						  {-2,3,-2},
						  {0,-2,3} };
		double[] b5 = {1,1,1};
		System.out.println(Calc.vecNorm1(Calc.subVec(b5, Calc.matVec(a5, CalcNumerical.partialPivotGauss(a5, b5)))));
		NumericalData data  = new NumericalData(a5, b5);
		double[] initX = {1,-1,1};
		data.setInitX(initX);
		data.setEps(1.0E-7);
		data.setMaxN(200);
		data.setNorm(Norm.ONE);
		data.setCon(ConvergenceCriterion.RESIDUAL);
		CalcNumerical.jacobi(data);
		System.out.println(data.fullfilConvergence());
		double[] x_ = CalcNumerical.gaussSeidel(data);
		CalcTool.printVec(x_);
		System.out.println(data.getCount());
		System.out.println(data.fullfilConvergence());
		System.out.println(Calc.vecNorm1(Calc.subVec(b5, Calc.matVec(a5, x_))));
		double[][] h = new double[10][10];
		for(int i =0 ;i < 10 ;i++){
			for(int j = 0;j < 10;j++){
				h[i][j] = Math.sqrt(Math.pow(i+1, 2) + Math.pow(j+1, 2) - 1);
			}
		}
		CalcTool.printMat(h);
		double[][] luh = CalcNumerical.luDecomposition(h);
		sum = 1;
		for(int i = 0;i < 10;i++){
			sum*= luh[i][i];
		}
		System.out.println(sum);
		
	}
	static void test(){
		double x[] = new double[49];
		for(int i = 0; i < 49 ; i++){
			x[i] = Math.sqrt(i + 1);
		}
		//System.out.println(Calc.vecNorm2(x));
		int n = 100;
		double[][] a = new double[n][n];
		for(int i = 0; i < n;i++){
			for(int j = 0; j < n;j++){
					a[i][j] = i + j +2;
			}
		}
		//CalcTool.printMat(a);
		//System.out.println(Calc.matNormFrobenius(a));
		n = 6;
		double[][] h = CalcTool.createHilbert(n);
		//System.out.println(Calc.conditionNumber(h, Norm.INFINITY));
		double[] hb = new double[n]; 
		double[] trueX = new double[n];
		for(int i = 0; i < n ; i++){
		 for(int j = 0; j < n;j++){
			 hb[i] += h[i][j];
		 }
		 trueX[i] = 1;
		}
		//CalcTool.printVec(CalcNumerical.gauss(h, hb));
		//CalcTool.printVec(CalcNumerical.partialPivotGauss(h, hb));
		//CalcTool.printVec(CalcNumerical.completePivotGauss(h, hb));
		double[] nearX = CalcNumerical.partialPivotGauss(h,hb);
		NumericalData hdata = new NumericalData(h, hb);
		hdata.setCon(ConvergenceCriterion.RELATIVEERROR);
		hdata.setNorm(Norm.INFINITY);
		//System.out.println(CalcTool.calcConvergence(hdata, trueX, nearX));
		System.out.println(Calc.vecNormInf(Calc.subVec(trueX, nearX)) / Calc.vecNormInf(trueX));
		hb[0] += 0.001 + hb[0];
		nearX = CalcNumerical.partialPivotGauss(h,hb);
		System.out.println(Calc.vecNormInf(Calc.subVec(trueX, nearX)) / Calc.vecNormInf(trueX));
		System.out.println("////////////////////////3//////////////////////");
		n=3;
		double[][] ta = {{1,-2,2},{-1,1,-1},{-2,-2,1}};
		double[] tb = {1,1,1};
		NumericalData tdata = new NumericalData(ta, tb);
		double[] initX = { 0,0,0};
		tdata.setInitX(initX);
		tdata.setCon(ConvergenceCriterion.RELATIVEERROR);
		tdata.setNorm(Norm.INFINITY);
		tdata.setEps(1.0E-12);
		tdata.setMaxN(50);
		CalcTool.printVec(CalcNumerical.jacobi(tdata));
		System.out.println(tdata.fullfilConvergence());
		System.out.println(tdata.getCount());
		
		double[][] pa = {{0,-1,2},{-4,5,6},{8,9,10}};
		double[] pb = {3,7,11};
		NumericalData pdata = new NumericalData(pa, pb);
		pdata.setInitX(initX);
		pdata.setCon(ConvergenceCriterion.RELATIVEERROR);
		pdata.setNorm(Norm.TWO);
		pdata.setEps(1.0E-8);
		pdata.setMaxN(50);
		CalcNumerical.jacobi(pdata);
		System.out.println(pdata.fullfilConvergence());
		n = 100;
		double[][] ta1 = new double[n][n];
		double[] tb1 = new double[n];
		double[] initXp = new double[n];
		int ka = 9;
		for(int i = 0; i < n;i++){
			for(int j = 0; j < n;j++){
				if(i == j){
					ta1[i][j] = ka;
				}else if(i == j-1){
					ta1[i][j] = 4;
				}else if(i == j+1){
					ta1[i][j] = 4;
				}
			}
			tb1[i] = 1;
			initXp[i] = 0;
		}
		NumericalData kdata = new NumericalData(ta1, tb1);
		kdata.setInitX(initXp);
		kdata.setCon(ConvergenceCriterion.RELATIVEERROR);
		kdata.setNorm(Norm.INFINITY);
		kdata.setEps(1.0E-8);
		kdata.setMaxN(500);
		for(int i = 1;i < 20 ; i++){
			CalcNumerical.sor(kdata, i / 10.0);
			System.out.println("ω = " + (i / 10.0) + " ; 反復回数 : " +kdata.getCount());
		}
		double[][] la = new double[n][n];
		for(int i = 0; i < n;i++){
			for(int j = 0; j < n;j++){
				if(i == j){
					la[i][j] = 4;
				}else if(i == j-1){
					la[i][j] = 2;
				}else if(i == j+1){
					la[i][j] = 2;
				}else if(i == j-2){
					la[i][j] = 1;
				}else if(i == j+2){
					la[i][j] = 1;
				}
			}
		}
		double[][] ila = CalcNumerical.inverse(la);
		double sum = 0;
		CalcTool.printMat(ila);
		for(int i= n-1; 0 < i;i--){
			sum += ila[i][i];
		}
		System.out.println(sum);
		double[][] te = {{1,0,0},{-1,1,0},{-2,-2,1}};
		double[][] tr = {{0,-2,2},{0,0,-1},{0,0,0}};
		CalcTool.printMat(CalcNumerical.inverse(te));
		CalcTool.printMat(Calc.multipleMat(CalcNumerical.inverse(te) , tr));
	}
	
}
