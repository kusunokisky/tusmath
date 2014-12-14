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
		}
		//後退代入
		for(int i = b_.length-1; i >= 0 ;i--){
			for(int j = i+1; j < b_.length ;j++){
				b_[i] -= a_[i][j]*b_[j];
			}
			b_[i] /= a_[i][i];
		}
		//解ベクトルの解の入れ替えを行う
		//CalcTool.baseBubbleSort(b_, base);
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
		for(int k = 0;k < a_.length;k++){
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

		for(int i = 0; i < lu.length ;i++){
			for(int j = 0; j < i ;j++){
				 y[i] -= lu[i][j]*y[j];
			}
		}
		for(int i = y.length-1; i >= 0 ;i--){
			for(int j = i+1; j < y.length ;j++){
				y[i] -= lu[i][j]*y[j];
			}
			y[i] /= lu[i][i];
		}
		return y;
	}
	/**
	 * LU分解を実行する
	 * @param a LU分解したい行列
	 * @return L ＼ U という形に分解された行列
	 */
	static double[][] luDecomposition(double[][] a){
		double a_[][] = CalcTool.arrayCopy(a);
		for(int k = 0; k < a_.length ;k++){
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
}
