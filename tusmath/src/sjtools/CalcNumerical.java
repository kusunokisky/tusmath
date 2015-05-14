package sjtools;

import java.util.Arrays;

public class CalcNumerical{

	protected CalcNumerical() {
		throw new UnsupportedOperationException();
	}
	/**
	 * ホーナー法を用いて多項式の値を計算
	 * @param x
	 * @param coefficient x^n ～ x^0の係数
	 * @return xを代入した値
	 */
	static double hornerRule(double x,double[] coefficient){
		double result = 0;
		for(int i = coefficient.length-1 ; 0 < i ; i-- ){
			result = result * x + coefficient[i];
		}
		return result;
	}
	/**
	 * 二分法による収束判定
	 * @param a 初期点1
	 * @param b 初期点2
	 * @param eps 許容誤差
	 * @param f 関数
	 * */
	 static double bisectionMethod(double a,double b,double eps,Function f){
		double a_ = a;
		double b_ = b;
		double c = 0;
		int count = 0;
		while(Math.abs(b_ - a_) / 2.0 >= eps){


				c  = (a_ + b_) / 2.0;
				double dep = f.fx(a_) * f.fx(c);

				if(dep > 0){
					a_ = c;
				}else if(dep <0){
					b_ = c;
				}else{
					break;
				}
				count++;

		}
		System.out.println("n = " + count);
		return c;
	}
	/**
	 * ニュートン法の相対誤差による収束判定
	 * <p>
	 * step1   : 初期値 x0,許容誤差 eps,最大反復回数 Nを与える<br>
	 * </p>
	 * <p>
	 * Step2   : k=0,1,...,Nについて以下を繰り返す<br>
	 * Step2-1 : Xk+1 = Xk - f(Xk)/f`(Xk)を計算する<br>
	 * Step2-2 : |(Xk+1 - Xk) / Xk+1| < eps ならばStep3へ
	 * </p>
	 * <p>
	 * Step3   : Xk+1を解とする
	 * </p>
	 * @param x 初期値
	 * @param eps 許容誤差
	 * @param n 最大反復回数
	 * @param f 関数
	 * @return xk+1 解
	 *
	 *
	 * */
	static double relativeErrorNewton(double x,double eps,int n,Function f){
		double xk = x;
		double xk1 = 0;
		//DecimalFormat df1 = new DecimalFormat("##.########################################");
		for(int i = 0;i < n;i++){
			xk1 = xk - f.fx(xk)/f.dfx(xk);
			//System.out.println(df1.format( Math.abs( (xk1 -xk) / xk1) ) + ",");
			if(Math.abs((xk1 -xk)/xk1) < eps){
				System.out.println("n = " + (i+1));
				break;
			}

			xk = xk1;
			//System.out.println(xk);
		}
		return xk1;
	}
	/**
	 * ニュートン法の残差による収束判定
	 * <p>
	 * step1   : 初期値 x0,許容誤差 eps,最大反復回数 Nを与える<br>
	 * </p>
	 * <p>
	 * Step2   : k=0,1,...,Nについて以下を繰り返す<br>
	 * Step2-1 : |f(Xk)| < eps ならばStep3へ<br>
	 * Step2-2 : Xk+1 = Xk - f(Xk)/f`(Xk)を計算する
	 * </p>
	 *
	 * <p>
	 * Step3   : Xkを解とする
	 * </p>
	 *
	 * @param x 初期値
	 * @param eps 許容誤差
	 * @param n 最大反復回数
	 * @param f 関数
	 * @return xk 解
	 * */
	static double residualNewton(double x,double eps,int n,Function f){
		double xk =  x;
		for(int i = 0;i <=  n;i++){
			if(Math.abs(f.fx(xk)) < eps){
				break;
			}else{
				xk = xk - f.fx(xk)/f.dfx(xk);
			}
		}
		return xk;
	}
	static double relativeErrorParallelChord(double x,double eps,int n,Function f){
		double xk = x;
		double xk1 =0;
		for(int i = 0;i < n;i++){
			xk1 = xk -f.fx(xk)/f.dfx(x);
			if(Math.abs((xk1 -xk)/xk1) < eps){
				System.out.println("n = " + (i+1));
				break;
			}
			xk = xk1;
		}
		return xk1;
	}
	static double residualParallelChord(double x,double eps,int n,Function f){
		double xk = x;
		for(int i = 0;i <= n;i++){
			if(Math.abs(f.fx(xk)) < eps){
				System.out.println("n = " + i );
				break;
			}else{
				xk = xk - f.fx(xk)/f.dfx(x);
			}
		}

		return xk;
	}
	/**
	 * セカント法の相対誤差による収束判定
	 * @param x0 初期点1
	 * @param x1 初期点2
	 * @param eps 許容誤差
	 * @param n 最大反復回数
	 * @param f 関数
	 * @return 解
	 */
	static double relativeErrorSecant(double x0,double x1,double eps,int n,Function f){
		double xk = x0;
		double xk1 = x1;
		for(int i = 0;i < n;i++){
			double xk1_ = xk1;
			xk1 = xk - ( f.fx(xk)*(xk1-xk) / ( f.fx(xk1) - f.fx(xk) ) );
			xk = xk1_;
			if(( Math.abs( (xk1-xk)/xk1) ) < eps) {
				System.out.println("n = " + (i+1));
				break;
			}
		}
		return xk1;
	}
	/**
	 * セカント法の残差による収束判定
	 * @param x0 初期点1
	 * @param x1 初期点2
	 * @param eps 許容誤差
	 * @param n 最大反復回数
	 * @param f 関数
	 * @return 解
	 */
	static double residualSecant(double x0,double x1,double eps,int n,Function f){
		double xk = x0;
		double xk1 = x1;
		double temp = 0.0;
		int count = 0;
		while(Math.abs(f.fx(xk) ) >= eps){
			temp = xk - f.fx(xk)*(xk-xk1)/(f.fx(xk)-f.fx(xk1));
			xk = xk1;
			xk1 = temp;
			count++;
			if(count == n){
				break;
			}
		}
		//System.out.println("n = " +  count);
		return xk;
	}
	/**
	 * 部分ピポッド選択付ガウスの消去法を使いAx=bの解ベクトルを求める
	 * @param a 方程式の係数行列
	 * @param b 方程式の右辺項
	 * @return 解ベクトル
	 */
	static double[] partialPivotGauss(double a[][],double[] b){
		double[][] a_ = CalcTool.arrayCopy(a);
		double[] b_ = Arrays.copyOf(b, b.length);
		//前進消去
		for(int k = 0;k < a_.length;k++){
			//Pivot選択
			int pivot = k;
			double max = Math.abs(a_[k][k]);
			for(int m = k;m < a.length-1;m++){
				double dep = Math.abs(a_[m+1][k]);
				if(max < dep){
					max = dep;
					pivot = m+1;
				}
			}
			if(k != pivot){
				CalcTool.transformationG(a_, k, pivot);
				CalcTool.swap(b_, k, pivot);
			}
			for(int i = k+1;i < a_.length;i++){
				double t = a_[i][k]/a_[k][k];
				for(int j = k+1 ;j < a_[i].length;j++){
					a_[i][j] -=  t * a_[k][j];
				}
				b_[i] -= t * b_[k];
			}
		}
		//後退代入
		for(int i = b_.length-1; i >= 0 ;i--){
			for(int j = i+1; j < b_.length ;j++){
				b_[i] -= a_[i][j]*b_[j];
			}
			b_[i] /= a_[i][i];
		}
		return b_;
	}
	/**
	 * 完全ピポッド選択付ガウスの消去法を使いAx=bの解ベクトルを求める
	 * @param a 方程式の係数行列
	 * @param b 方程式の右辺項
	 * @return 解ベクトル : 1次元配列
	 */
	static double[] completePivotGauss(double a[][],double[] b){
		double[][] a_ = CalcTool.arrayCopy(a);
		double[] b_ = Arrays.copyOf(b, b.length);
		int[] base = CalcTool.createBaseData(b.length);
		//前進消去法
		for(int k = 0;k < a_.length;k++){
			//ピポットの初期化
			int rowPivot = k;
			int columnPivot = k;
			//最大値の設定
			double max = Math.abs(a_[k][k]);
			//作業対象内からの最大要素の検索
			for(int i = k; i < a_.length;i++){
				for(int j = k;j<a_[i].length;j++){
					double dep = Math.abs(a_[i][j]);
					if(max < dep){
						max = dep;
						rowPivot = i;
						columnPivot = j;
					}
				}
			}
			//System.out.println("i" + rowPivot + "j" + columnPivot);
			if(k != rowPivot){
				CalcTool.transformationG(a_, k, rowPivot);
				CalcTool.swap(b_, k, rowPivot);
			}
			if(k != columnPivot){
				CalcTool.transformationGC(a_, k, columnPivot);
				CalcTool.swap(base, k, columnPivot);
			}
			for(int i = k+1;i < a_.length;i++){
				double t = a_[i][k]/a_[k][k];
				for(int j = k+1 ;j < a_[i].length;j++){
					a_[i][j] -=  t * a_[k][j];
				}
				b_[i] -= t * b_[k];
			}
			CalcTool.printVec(b_);
		}
		CalcTool.printMat(a_);
		//後退代入
		for(int i = b_.length-1; i >= 0 ;i--){
			for(int j = i+1; j < b_.length ;j++){
				b_[i] -= a_[i][j]*b_[j];
			}
			b_[i] /= a_[i][i];
		}
		//解ベクトルの解の入れ替えを行う
		//CalcTool.baseBubbleSort(b_, base);
		System.out.println(base[0]);
		System.out.println(base[1]);
		System.out.println(base[2]);
		System.out.println(base[3]);
		System.out.println(base[4]);
		System.out.println();
		CalcTool.baseQuickSort(b_, base, 0, b_.length-1);
		return b_;
	}
	/**
	 * ガウスの消去法を使い解を求めた配列を返す
	 *
	 * @param a 方程式の係数配列
	 * @param b 方程式の数値
	 * @return 解の配列
	 */
	static double[] gauss(double a[][],double b[]){
		//配列のコピーを作成する
		double a_[][] = CalcTool.arrayCopy(a);
		double b_[] = Arrays.copyOf(b, b.length);
		//前進消去法
		for(int k = 0;k < a_.length-1;k++){
			for(int i = k+1;i < a_.length;i++){
				double t = a_[i][k]/a_[k][k];
				for(int j = k+1 ;j < a_[i].length;j++){
					a_[i][j] -=  t * a_[k][j];
				}
				b_[i] -= t * b_[k];
			}
		}
		//後退代入
		for(int i = b_.length-1; i >= 0 ;i--){
			for(int j = i+1; j < b_.length ;j++){
				b_[i] -= a_[i][j]*b_[j];
			}
			b_[i] /= a_[i][i];
		}
		return b_;
	}
	/**
	 * LU分解された行列を利用して解を求める
	 *
	 * @param lu Lu分解された行列
	 * @param b 右辺の項
	 * @return 求められた方程式の解
	 */
	static double[] solveLU(double[][] lu,double[]b){
		double y[] = Arrays.copyOf(b, b.length);
		/*`int plus = 0;
		int sub = 0;
		int multi = 0;
		int wari = 0;
		*/
		for(int i = 1; i < lu.length ;i++){
			for(int j = 0; j < i ;j++){
				//System.out.println("i : " + i+"j : " + j);
				 y[i] -= lu[i][j]*y[j];
				// multi++;
				// sub++;
			}
		}
		for(int i = y.length-1; i >= 0 ;i--){
			for(int j = i+1; j < y.length ;j++){
				//System.out.println("i : " + i+"j : " + j);
				y[i] -= lu[i][j]*y[j];
				//multi++;
				//sub++;
			}
			y[i] /= lu[i][i];
			//System.out.println("i = " + i);
			//wari++;
			//System.out.println("/ = " + wari);
		}
		//System.out.println("解を求めるときの演算回数");
		//System.out.println(plus +" : " + sub+" : " + multi+" : " + wari);
		return y;
	}
	/**
	 * LU分解を実行する
	 * @param a LU分解したい行列
	 * @return L ＼ U という形に分解された行列
	 */
	static double[][] luDecomposition(double[][] a){
		double a_[][] = CalcTool.arrayCopy(a);
		for(int k = 0; k < a_.length-1 ;k++){
			for(int i = k+1; i < a_.length;i++){
				a_[i][k] = a_[i][k] / a_[k][k];
				for(int j = k + 1;j < a_.length;j++){
					a_[i][j] -=  a_[i][k] * a_[k][j];
				}
			}
		}
		return a_;
	}
	/**
	 * LU分解を利用してAx=bの解であるxのベクトルを求める
	 * @param a 方程式の係数行列
	 * @param b 方程式の右辺項
	 * @return 解ベクトル
	 */
	static double[] luDecomposition(double[][] a, double[] b){
		return solveLU(luDecomposition(a),b);
	}
	/**
	 * LU分解された行列を利用して逆行列を求める
	 * @param lu lu分解された行列
	 * @return  逆行列
	 */
	private static double[][] luInverse(double lu[][]){
		double[][] e  = CalcTool.unitMat(lu.length);
		for(int i = 0;i < lu.length;i++){

			for(int k =0;k < lu.length;k++){
				for(int j = 0; j < k ;j++){
					e[k][i] -= lu[k][j]*e[j][i];
				}
			}
			for(int k = lu.length-1; k >= 0; k--){
				for(int j = k+1; j<lu.length; j++){
					e[k][i] -= lu[k][j]*e[j][i];
				}
				e[k][i] /= lu[k][k];
			}

		}
		return e;
	}
	/**
	 * 正則な行列の逆行列を求める
	 * @param a 正則な行列
	 * @return 逆行列
	 */
	static double[][] inverse(double[][] a){
		return luInverse(luDecomposition(a));
	}
	/**
	 * ヤコビ法により解ベクトルを求める
	 * A = D+E+F (D:対角行列 E:狭義上三角行列 F:狭義下三角行列)
	 * T = D^-1(E+F) : 反復行列
	 *  * @return 解ベクトル
	 */
	static double[] jacobi(NumericalData data){
		double[][] a = CalcTool.arrayCopy(data.getA());
		double[] b = Arrays.copyOf(data.getB(), data.getB().length);
		double[] dep = Arrays.copyOf(data.getInitX(), data.getInitX().length);
		double[] result = new double[dep.length];
		double sum;
		int count = 0;
		while(true){
			for(int i = 0;i < result.length;i++){
				sum = 0;
				for(int j = 0;j < result.length;j++){
					if(i != j){
						sum += a[i][j]*dep[j];
					}
				}
				result[i] = (b[i] - sum)/a[i][i];
			}
			count++;
			//収束判定
			if(CalcTool.calcConvergence(data, result, dep) < data.getEps()){
				data.setConvergence(true);
				data.setCount(count);
				break;
			}
			//反復判定
			if( data.getMaxN() < count ){
				data.setConvergence(false);
				data.setCount(-1);
				return new double[0];
			}
			dep = Arrays.copyOf(result, result.length);
		}
		return result;
	}
	/**
	 * ガウス・ザイデル法により解ベクトルを求める
	 * @param data 数値計算用データ
	 * @return 解ベクトル
	 */
	static double[] gaussSeidel(NumericalData data){
		double[][] a = CalcTool.arrayCopy(data.getA());
		double[] b = Arrays.copyOf(data.getB(), data.getB().length);
		double[] dep = Arrays.copyOf(data.getInitX(), data.getInitX().length);
		double[] result = new double[dep.length];
		double sum;
		int count = 0;
		while(true){
			for(int i = 0;i < result.length;i++){
				sum = 0;

				for(int j = 0;j < result.length;j++){
					if(j < i){
						sum += a[i][j] * result[j];
					}else if(j > i){
						sum += a[i][j] * dep[j];
					}
				}
				result[i] = (b[i] - sum)/a[i][i];
			}
			count++;
			//収束判定
			if(CalcTool.calcConvergence(data, result, dep) < data.getEps()){
				data.setConvergence(true);
				data.setCount(count);
				break;
			}
			//反復判定
			if( data.getMaxN() < count ){
				data.setConvergence(false);
				data.setCount(-1);
				return new double[0];
			}
			dep = Arrays.copyOf(result, result.length);
		}
		return result;
	}
	/** SOR法を用いて解ベクトルを求める
	 * @param data 数値計算用データ
	 * @param accelerate 加速パラメータ
	 * @return 解ベクトル
	 */
	static double[] sor(NumericalData data,double accelerate){
		double[][] a = CalcTool.arrayCopy(data.getA());
		double[] b = Arrays.copyOf(data.getB(), data.getB().length);
		double[] dep = Arrays.copyOf(data.getInitX(), data.getInitX().length);
		double[] result = new double[dep.length];
		double sum;
		int count = 0;
		while(true){
			for(int i = 0;i < result.length;i++){
				sum = 0;

				for(int j = 0;j < result.length;j++){
					if(j < i){
						sum += a[i][j] * result[j];
					}else if(j > i){
						sum += a[i][j]*dep[j];
					}
				}
				result[i] = (b[i] - sum)/a[i][i];
				result[i] = (1 - accelerate) * dep[i] + (accelerate * result[i]);
			}
			count++;
			//収束判定
			if(CalcTool.calcConvergence(data, result, dep) < data.getEps()){
				data.setConvergence(true);
				data.setCount(count);
				break;
			}
			//反復判定
			if( data.getMaxN() < count ){
				data.setConvergence(false);
				data.setCount(-1);
				//空配列をかえす
				return new double[0];
			}
			dep = Arrays.copyOf(result, result.length);

		}
		data.setCount(count);
		return result;
	}
	/**
	 * 修正コレスキー分解された行列を利用して解を求める
	 * @param ldl
	 * @param b
	 * @return
	 */
	static double[] solveCholesky(double[][] ldl,double[] b){
		double[][] ldl_ = CalcTool.arrayCopy(ldl);
		double y[] = Arrays.copyOf(b, b.length);
		double x[] = new double[y.length];
		/*int plus = 0;
		int sub = 0;
		int multi = 0;
		int wari = 0;*/
		for(int i = 0; i < ldl.length ;i++){
			for(int j = 0; j < i ;j++){
				 y[i] -= ldl_[i][j]*y[j];
				 //multi++;
				 //sub++;
			}
		}
		for(int i = ldl.length-1;i >=0;i--){
			x[i] = y[i]*ldl[i][i];
			//multi++;
			for(int j = i+1;j<ldl.length;j++){
				x[i] -= ldl[j][i]*x[j];
				//multi++;
				//sub++;
			}
		}
		//System.out.println(plus +" : " + sub+" : " + multi+" : " + wari);
		return x;
		//System.out.println(plus +" : " + sub+" : " + multi+" : " + wari);
		//System.out.println(plus+sub+multi+wari);
		//return y;
	}
	/**
	 * 修正コレスキー分解を利用してAx=bの解であるxのベクトルを求める
	 * @param a 方程式の係数行列
	 * @param b 方程式の右辺項
	 * @return 解ベクトル
	 */
	static double[] choleskyDecomposition(double[][] a, double[] b){
		return solveCholesky(choleskyDecomposition(a),b);
	}
	/**
	 * 修正コレスキー分解をする
	 * @param a
	 * @return
	 */
	static double[][] choleskyDecomposition(double[][] a){
		double[][] a_ = CalcTool.arrayCopy(a);
		double[] w= new double[a_.length];
		for (int j = 0; j < a_.length; j++) {
			for (int i = 0; i < j ; i++) {
				w[i] = a_[j][i];
				for (int k = 0; k < i ; k++) {
					w[i] -= a_[i][k] * w[k];
					}
				a_[j][i] = w[i] * a_[i][i];
				}
			double t = a_[j][j];
			for (int k = 0; k < j; k++) {
				t -= a_[j][k] * w[k];
			}
			a_[j][j] = 1.0 / t;
		}
		return a_;
	}
	static double[] power(NumericalData data){
		double[] old_x = data.getInitX();
		double[] new_x = new double[old_x.length];
		old_x = CalcTool.normalized(old_x, Norm.TWO);
		new_x = Calc.matVec(data.getA(), old_x);
		int max = 0;
		for(int i = 1; i< old_x.length;i++){
			if(old_x[max] < old_x[i]){
				max = i;
			}
		}
		double old_eigen = new_x[max] / old_x[max];
		old_x = CalcTool.normalized(new_x, Norm.TWO);
		max = 0;
		double new_eigen = 0;
		int maxN = data.getMaxN();
		for(int i = 1;i < maxN;i++){
			new_x = Calc.matVec(data.getA(), old_x);
			for(int j = 1; j < old_x.length;j++){
				if(old_x[max] < old_x[j]){
					max = j;
				}
			}
			new_eigen = new_x[max] / old_x[max];
			if(Math.abs(new_eigen - old_eigen)/Math.abs(new_eigen) < data.getEps()){
				break;

			}
			old_x = CalcTool.normalized(new_x, Norm.TWO);
			old_eigen = new_eigen;
			max = 0;
		}
		System.out.println(new_eigen);
		return new_x;
	}
	/**
	 * 実対称行列を3重対角行列へと変換する
	 * @param a 実対称行列
	 * @return 3重対角行列
	 */
	static double[][] tridiagonalization(double a[][]){


		return null;
	}
	static int sign(double x){
		//
		if(x < 0){
			return -1;
		}else{
			return 1;
		}
	}
	static double[][] tridiagonal(double[][] a){
		double[][] a_ = CalcTool.arrayCopy(a);
		int n = a.length-2;
		int t = a.length-1;

		double[] u_ = new double[t];
		double[][] p_ = new double[t][t];
		double[][] pa = new double[t][t];
		double[][] pap = new double[t][t];
		for(int k = 0; k < n;k++,t--){
			//変更行列の位置
			int r = k+1;
			//ノルム合わせの値S
			double s = 0;
			for(int i = r;i<a.length;i++){
				s += Math.pow(a[i][k],2);
			}
			s = Math.sqrt(s);
			//副対角成分の値
			double na = -sign(a_[r][k])*s;
			//ベクトルuを定める
			for(int i = 0;i < t;i++){
				if(i == 0){
					u_[i] = sign(a_[r][k])*Math.sqrt((1+Math.abs(a_[r][k])/s)/2);
				}else{
					u_[i] = a[r+i][k]/(2*s*Math.abs(u_[0]));
				}
			}
			//ハウスホルダー小行列を定める
			for(int i = 0;i < t;i++){
				for(int j = 0;j<t;j++){
					if(j==i){
							p_[i][j] = 1.0-2*u_[i]*u_[j];
					}else{
							p_[i][j] = -2*u_[i]*u_[j];
					}
				}
			}
			//PAPの小行列を作成
			//P*Aの積
			for(int i = 0;i < t;i++ ){
				for(int j = 0;j<t;j++){
					double sum = p_[i][0]*a_[r][j+r];;
					for(int l =1;l < t;l++){
						sum += p_[i][l]*a_[l+r][j+r];
					}
					pa[i][j] = sum;
				}
			}
			System.out.println("pa");
			CalcTool.printMat(pa);
			//PA * Pの積
			for(int i = 0;i < t;i++ ){
				for(int j = 0;j<t;j++){
					double sum = 0;
					for(int l =0;l < t;l++){
						sum += pa[i][l]*p_[l][j];
					}
					pap[i][j] = sum;
				}
			}
			CalcTool.printMat(pap);
			//行列ＡにPAP小行列を差し込んでいく
			for(int i = 0;i < a.length;i++){
				for(int j = 0; j < a.length;j++){
					if(i == k){
						if(j==r){
							a_[i][j] = na;
						}else if(r < j){
							a_[i][j] = 0;
						}
					}else if(i==r){
						if(j==k){
							a_[i][j] = na;
						}else if(r <= j){
							a_[i][j] = pap[i-r][j-r];
						}
					}else if(r < i){
						if(j==k){
							a_[i][j] = 0;
						}else if(r <= j){
							a_[i][j] = pap[i-r][j-r];
						}
					}
				}
			}
		}
		return a_;
	}
}

