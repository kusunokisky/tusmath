package sjtools;

import java.util.Arrays;
/**
 * 配列を使い演算を行う上で必要な処理を行う
 * @author takeru
 */

public class CalcTool {
	private CalcTool() {}
	/**
	 * ベクトル(1次元配列)をコンソール出力する
	 * @param x
	 */
	static void printVec(double x[]){
		if(x == null){
			System.out.println("不正な値です");
			return;
		}
		for(int i = 0; i < x.length;i++){
			System.out.print(x[i] + " ");
		}
		System.out.println();
	}
	/**
	 * 行列(2次元配列)をコンソール出力する
	 * @param a
	 */
	static void printMat(double a[][]){
		if(a == null){
			System.out.println("不正な値です");
			return;
		}
		for(int i = 0; i < a.length;i++){
			for(int j = 0; j < a[i].length;j++){
				System.out.print(a[i][j] + " ");
			}
			System.out.println();
		}
		System.out.println();
	}
	/**
	 * n行n列の単位行列を作成する
	 * @param n 大きさ
	 * @return n×nの単位行列
	 */
	static double[][] unitMat(int n){
		double[][] unit_a = new double[n][n];
		for(int i = 0;i < unit_a.length;i++){
			for(int j = 0; j < unit_a[i].length;j++){
				if(i == j){
					unit_a[i][j] = 1;
				}else{
					unit_a[i][j] = 0;
				}
			}
		}
		return unit_a;
	}
	/**
	 * 逆行列を返す
	 * @param a 正則な行列
	 * @return aに対する逆行列
	 */
	static double[][] inverse(double a[][]){
		//行列のコピー
		double[][] a_ = arrayCopy(a);
		//単位行列の作成
		double[][] inv_a = unitMat(a_.length);
		//掃き出し法
		for(int i = 0 ; i < a_.length;i++){
			//0ではない先頭要素を選択する
			if(a_[i][i] == 0){
				for(int k = i+1 ;k < a_.length;k++ ){
					if(a_[k][i] != 0){
						CalcTool.transformationG(inv_a, i, k);
						CalcTool.transformationG(a_, i, k);
						break;
					}
				}
			}
			//先頭要素以外を0にする
			for(int k = 0;k < a_.length;k++){
				if(k!=i){
					CalcTool.transformationH(inv_a, k, i, -a_[k][i]/a_[i][i]);
					CalcTool.transformationH(a_, k, i, -a_[k][i]/a_[i][i]);
				}
			}
		}
		//元の要素を正規化する
		for(int i = 0; i < a_.length ; i++){
			CalcTool.transformationF(inv_a, i, 1/a_[i][i]);
			CalcTool.transformationF(a_, i, 1/a_[i][i]);
		}
		return inv_a;
	}
	/**
	 * 転置した行列を返す
	 * @param a 転置される行列
	 * @return 転置した行列
	 */
	static double[][] transpose(double[][] a){
		double[][] trans_a = new double[a[0].length][a.length];
		for(int i = 0; i < a.length;i++){
			for(int j = 0; j < a[i].length;j++){
				trans_a[i][j] = a[j][i];
			}
		}
		return trans_a;
	}
	/**
	 * 数値を分数で表示する
	 *
	 * @param x
	 */
	@Deprecated
	static void printFraction(double x){
		//数値を文字列に変換
		int ix = (int)x;
		double smallX = x - ix;
		String s = String.valueOf(smallX);
		double t = Math.pow(10, s.length() - 1);
		int nume = (int)( smallX * Math.pow( 10, s.length()-1 ) );
		int deco = (int)t;
		while(true){
			int gcd = Calc.gcd(nume,deco);
			if(gcd == 1)break;

			nume /= gcd;
			deco /= gcd;
		}
		System.out.println(nume + "/" + deco);
	}
	/**
	 * 正方行列の複製を作成する
	 *
	 * @param a コピーしたい2次元配列
	 * @return コピーされた2次元配列
	 */
	static  double[][] arrayCopy(double a[][]){
		double a_[][] = new double[a.length][];
		for(int i = 0; i < a.length ;i++){
			a_[i] = Arrays.copyOf(a[i], a[i].length);
		}
		return a_;
	}
	static double[][] scalarMat(double[][] x , double c){
		double[][] x_ = CalcTool.arrayCopy(x);
		for(int i = 0;i < x_.length;i++){
			transformationF(x_, i, c);
		}
		return x_;
	}
	static double[][] copyMat(double a[][],int n){
		int t = a.length-n;
		double[][] a_ = new double[t][t];
		for(int k = 0;k<t;k++){
			for(int l = 0;l<t;l++){
				a_[k][l] = a[k+n][l+n];
			}
		}
		return a_;
	}
	/**
	 * i行目をc倍する
	 * @param x 行列
	 * @param i c倍される行
	 * @param c 倍率
	 * @return 作業が施された行列
	 */
	static double[][] transformationF(double[][] x,int i, double c){
		for(int k = 0;k < x[i].length;k++){
			x[i][k] *= c;
		}
		return x;
	}
	/**
	 * i行目をj行目を入れ替える
	 * @param x 行列
	 * @param i 入れ替えられる行(1,2,,,,)
	 * @param j 入れ替える行(1,2,,,,)
	 */
	static void transformationG(double[][] x,int i, int j){
		if(x[i].length != x[j].length)throw new UnsupportedOperationException(i +"行目と" + j + "行目の配列の長さが異なる");
		for(int k = 0;k < x[i].length;k++){
			double dep = x[i][k];
			x[i][k] = x[j][k];
			x[j][k] = dep;
		}
	}
	/**
	 * i列目とj列目を入れ替える
	 * @param x 元となるデータ
	 * @param i 交換したい列
	 * @param j 交換される列
	 */
	static void transformationGC(double[][] x,int i ,int j){
		for(int k = 0 ;k < x.length ;k++){
			double dep = x[k][i];
			x[k][i] = x[k][j];
			x[k][j] = dep;
		}
	}
	/**
	 * i行目にj行目のc倍を加える
	 * @param x 行列
	 * @param i 加えられる行(1,2,,,,)
	 * @param j 加える行(1,2,,,,)
	 * @param c 倍率
	 * @return 作業を施された行列
	 */
	static double[][] transformationH(double[][] x,int i, int j ,double c){
		double x_[][] = x;
		for(int k = 0;k < x[i].length;k++){
			x_[i][k] += c * x_[j][k];
		}
		return x_;
	}
	/**
	 * 配列の要素を任意の位置と交換する
	 * @param x 元となるデータ
	 * @param i 交換したい位置
	 * @param j 交換される位置
	 */
	static void swap(double x[],int i,int j){
		double dep = x[i];
		x[i] = x[j];
		x[j] = dep;
	}
	/**
	 * 配列の要素を任意の位置と交換する
	 * @param x 元となるデータ
	 * @param i 交換したい位置
	 * @param j 交換される位置
	 */
	static void swap(int x[],int i,int j){
		int dep = x[i];
		x[i] = x[j];
		x[j] = dep;
	}
	/**
	 * 初期基準データを作成する
	 * {0,1,2,..}
	 * @param length データの長さ
	 * @return base 基準データ
	 */
	static int[] createBaseData(int length){
		int[] base = new int[length];
		for(int i = 0;i < length;i++){
			base[i] = i;
		}
		return base;
	}
	/**
	 * バブルソートを利用し,基準データを元にソートを行う
	 * @param x ソートしたいデータ
	 * @param base 変更された基準データ
	 */
	static void baseBubbleSort(double[] x,int[] base){
		for(int i = 0; i < base.length;i++){
			for(int j = base.length - 1; j > i;j--){
				if(base[j-1] > base[j]){
					swap(x, j-1, j);
					swap(base, j-1, j);
				}
			}
		}
	}
	/**
	 * クイックソートを利用し,基準データを元にソートを行う
	 * @param x ソートしたいデータ
	 * @param base 変更された基準データ
	 * @param left
	 * @param right
	 */
	public static void baseQuickSort(double[] x,int[] base, int left, int right) {
		int curleft = left;
		int curright = right;
		int pivot = base[(curleft + curright) / 2];

		do {
			while (base[curleft] < pivot) {
				curleft++;
			}

			while (base[curright] > pivot) {
				curright--;
			}
			if (curleft <= curright) {
				swap (base, curleft++, curright--);
				swap (x, curleft-1, curright+1);
			}
		} while (curleft <= curright);

		if (left < curright) {
			baseQuickSort(x, base, left, curright);
		}

		if (curleft < right) {
			baseQuickSort(x, base, curleft, right);
		}
	}
	/**
	 * 収束判定の値を返す
	 * 誤差 : ||真値-近似値||
	 * 残差 : ||b-Ax||(x:真値)
	 * @param data
	 * @param true_x
	 * @param near_x
	 * @return 収束判定の値
	 * @exception UnsupportedOperationException 使用できない動作
	 */
	static double calcConvergence(NumericalData data,double[] true_x,double[] near_x){
		switch (data.getCon()) {
		case ERROR:
			return Calc.vecNorm(Calc.subVec(true_x, near_x), data.getNorm());
		case RESIDUAL:
			return Calc.vecNorm(Calc.residual(data.getA(), true_x, data.getB()), data.getNorm()) ;
		case RELATIVEERROR:
			return Calc.vecNormInf(Calc.subVec(true_x, near_x)) / Calc.vecNormInf(true_x);
			//return Calc.vecNorm(Calc.subVec(true_x, near_x), data.getNorm()) / Calc.vecNorm(true_x, data.getNorm());

		case RELATIVERESIDUAL:
			return Calc.vecNorm(Calc.residual(data.getA(), true_x, data.getB()), data.getNorm())
				/ Calc.vecNorm(data.getB(), data.getNorm());
		default:
			throw new UnsupportedOperationException();
		}
	}
	/**
	 * ヒルバート行列を作成する
	 * @param n
	 * @return
	 */
	static double[][] createHilbert(int n){
		double[][] hilbert = new double[n][n];
		for(int i = 0;i < n ;i++){
			for(int j = 0;j < n;j++){
				hilbert[i][j] = 1.0/(i+j+1.0);
			}
		}
		return hilbert;
	}

}