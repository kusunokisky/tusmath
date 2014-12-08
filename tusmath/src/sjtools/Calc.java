package sjtools;

import sjtools.CalcNumerical.Norm;


/**
 * 配列の演算処理など処理を行う
 * @author takeru
 */
public class Calc {

	protected Calc() {
		throw new UnsupportedOperationException();
	}
	/**
	 * ベクトルxをスカラーc倍する
	 * @param c スカラー
	 * @param x ベクトル行列(1次元配列)
	 * @return 	スカラー倍されたベクトル行列<br/>
	 * 			第二引数がnullの場合はnull
	 */
	static double[] scalarMultiple(double c,double x[]){
		if(x == null)return null;
		double x_[] = x;
		for(int i = 0;i < x_.length ; i++){
			x_[i] = c * x_[i];
		}
		return x_;
	}
	/**
	 * ベクトル同士の加算(x+y)を計算する
	 * @param x ベクトル行列(1次元配列)
	 * @param y ベクトル行列(1次元配列)
	 * @return 	加算されたベクトル<br/>
	 * 			第一引数もしくは第二引数がnullの場合null
	 */
	static double[] addVec(double x[],double y[]){
		if(x == null || y== null)return null;
		if(x.length != y.length)throw new UnsupportedOperationException("第一引数と第二引数の配列の長さが異なるため");
		double result[] = new double[x.length];
		for(int i = 0;i<result.length;i++){
			result[i] = x[i] + y[i];
		}
		return result;
	}
	/**
	 *  ベクトル同士の減算(x-y)を計算する
	 * @param x ベクトル行列(1次元配列)
	 * @param y ベクトル行列(1次元配列)
	 * @return 	減算されたベクトル<br/>
	 * 			第一引数もしくは第二引数がnullの場合null
	 */
	static double[] subVec(double x[],double y[]){
		if(x == null || y== null)return null;
		if(x.length != y.length)throw new UnsupportedOperationException("第一引数と第二引数の配列の長さが異なるため");
		double result[] = new double[x.length];
		for(int i = 0;i<result.length;i++){
			result[i] = x[i] - y[i];
		}
		return result;
	}
	/**
	 *  ベクトル同士の内積(x, y)を計算する
	 * @param x ベクトル行列(1次元配列)
	 * @param y ベクトル行列(1次元配列)
	 * @return 	内積の値
	 * @throws UnsupportedOperationException 配列の長さが異なるとき<br/>
	 * 										 
	 */
	static double innProb(double x[],double y[]){
		if(x.length != y.length)throw new UnsupportedOperationException();
		double result = 0;
		for(int i = 0;i<x.length;i++){
			result += x[i] * y[i];
		}
		return result;
	}
	/**
	 * 行列aとベクトルxの積を計算する
	 * @param a
	 * @param x
	 * @return
	 */
	static double[] matVec(double a[][] , double x[]){
		for(int i = 0;i < a.length;i++)if(a[i].length != x.length)throw new UnsupportedOperationException();
		double result[] = new double[a.length];
		for(int i  = 0;i <a.length;i++){
			for(int j = 0;j < x.length;j++){
				result[i] += a[i][j] * x[j] ;
			}
		}
		return result;
	}
	/**
	 * 行列aとベクトルx,bに対して,b-Ax(残差)を計算する
	 * @param a
	 * @param x
	 * @param b
	 * @return
	 */
	static double[] residual(double a[][],double x[] , double b[]){
		return subVec(b, matVec(a, x) );
	}
	/**
	 * 行列同士の加算(a+b)を計算する
	 * @param a
	 * @param b
	 * @return
	 */
	static double[][] addMat(double a[][],double b[][]){
		for(int i = 0;i < a.length;i++)if(a[i].length != b.length)throw new UnsupportedOperationException();
		double result[][] = new double[a.length][a[0].length];
		for(int i = 0; i < result.length;i++){
			for(int j = 0; j < result[i].length;j++){
				result[i][j] = a[i][j] + b[i][j];
			}
		}
		return result;
	}
	/**
	 * 行列同士の積(ab)を計算する
	 * @param a
	 * @param b
	 * @return
	 */
	static double[][] multipleMat(double a[][],double b[][]){
		for(int i = 0;i < a.length;i++)if(a[i].length != b.length)throw new UnsupportedOperationException();
		double result[][] = new double[a.length][b[0].length];
		for(int i = 0; i < a.length;i++){
			for(int j = 0; j < b[0].length ; j++){
				for(int k = 0; k < a[0].length;k++){
					result[i][j] += a[i][k] * b[k][j];
				}
			}
		}
		return result;
	}
	/**
	 * ベクトルのノルムを計算する
	 * @param x ベクトル
	 * @param type ノルムの種類
	 * @return
	 */
	static double vecNorm(double x[],Norm type){
		switch (type) {
		case ONE:
			return vecNorm1(x);
		case TWO:
			return vecNorm2(x);
		case INFINITY:
			return vecNormInf(x);

		default:
			throw new IllegalArgumentException();
		}
	}
	/**
	 * ベクトルの１ノルムを計算する
	 * @param x
	 * @return
	 */
	static double vecNorm1(double x[]){
		double result = 0;
		for(int i = 0; i < x.length;i++){
			result += Math.abs(x[i]);
		}
		return result;
	}
	/**
	 * ベクトルの２ノルムを計算する
	 * @param x
	 * @return
	 */
	static double vecNorm2(double x[]){
		double result = 0;
		for(int i = 0;i < x.length;i++){
			result += Math.pow(Math.abs(x[i]),2);
		}
		return Math.sqrt(result);
	}
	static double vecNormP(double x[] ,int p){
		double result = 0;
		for(int i = 0;i < x.length;i++){
			result += Math.pow(Math.abs(x[i]), p);
		}
		return Math.pow(result, 1.0 / p);
	}
	/**
	 * ベクトルの∞ノルムを計算する
	 * @param x
	 * @return
	 */
	static double vecNormInf(double x[]){
		double result = x[0];
		for(int i = 1 ; i < x.length;i++){
			double dep = Math.abs(x[i]);
			if( result < dep){
				result = dep;
			}
		}
		return result;
	}

	/**
	 * 最大公約数を求める
	 * @param a
	 * @param b
	 * @return 最大公約数
	 */
	static int gcd(int a,int b){
		int result = 0;
		if(a > b){
			int r = 0;
			while(true){
				r = a % b;

				if(r == 0){
					result = b;
					break;
				}
				a = b;
				b = r;
			}
		}else if(b > a){
			int r = 0;
			while(true){
				r = b % a;
				if(r == 0){
					result = a;
					break;
				}
				b = a;
				a = r;
			}
		}else{
			result = a;
		}
		return result;
	}

	/**
	 *最大列和
	 */
	static double matNorm1(double x[][]){
		//最初の列を追加
		double result = 0;
		for(int i = 0;i < x.length;i++){
			result += Math.abs(x[i][0]);
		}
		//残りの列を判定していく
		for(int i = 1; i < x[0].length;i++){
			double sum = 0;
			for(int j = 0;j < x.length;j++){
				sum += Math.abs(x[j][i]);
			}
			if(result < sum){
				result = sum;
			}
		}
		return result;
	}
	/**
	 * 最大行和
	 *
	 * @param result 最大行和
	 * @param x n*nの正方行列 n>2
	 */
	static double matNormInf(double x[][]){
		//最初の行を追加
		double result = 0;
		for(int i = 0; i < x[0].length;i++){
			result += Math.abs(x[0][i]);
		}
		//残りの行を判定していく
		for(int i = 1; i < x.length;i++){
			double sum = 0;
			for(int j = 0; j < x[i].length;j++){
				sum += Math.abs(x[i][j]);
			}
			if(result < sum){
				result = sum;
			}
		}
		return result;
	}
	static double conditionNumber(double[][] a,Norm norm){
		switch (norm) {
		case ONE:
			return Calc.matNorm1(a) * Calc.matNorm1( CalcNumerical.inverse(a) );
		case INFINITY:
			return Calc.matNormInf(a) * Calc.matNormInf( CalcNumerical.inverse(a) );
		default:
			throw new UnsupportedOperationException();
		}

	}
	/**
	 * 相対誤差を返す
	 * @param x 真値
	 * @param x_ 近似値
	 * @return 相対誤差
	 */
	static double relativeError(double x , double x_){
		return Math.abs( (x_ - x) / x ) ;
	}
	/**
	 * 絶対誤差を返す
	 * @param x 真値
	 * @param x_ 近似値
	 * @return 絶対誤差
	 */
	static double absoluteError(double x , double x_){
		return Math.abs(x - x_);
	}
}
