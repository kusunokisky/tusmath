package sjtools;


class ImaginaryCalc extends Calc {
	private ImaginaryCalc() {
		throw new UnsupportedOperationException();
	}
	/**
	 * 複素数の絶対値を返す
	 * @param real 実部
	 * @param imaginary 虚部
	 * @return 絶対値
	 */
	static double abs(double real,double imaginary){
		return Math.sqrt(Math.pow(real, 2) + Math.pow(imaginary,2));
	}
	/**
	 * 虚数を含んだベクトルの１ノルムを計算する
	 * @param data 複素数
	 * @return
	 */
	static double vecNorm1(ImaginaryData[] data){
		double result = 0;
		for(int i = 0; i < data.length;i++){
			result += abs(data[i].getReal(),data[i].getImaginary());
		}
		return result;
	}
	/**
	 * 虚数を含んだベクトルの２ノルムを計算する
	 * @param data 複素数
	 * @return
	 */
	static double vecNorm2(ImaginaryData[] data){
		double result = 0;
		for(int i = 0; i < data.length ; i++){
			result += Math.pow( Math.pow(data[i].getReal(), 2) + Math.pow(data[i].getImaginary(), 2) ,2);
		}
		return Math.sqrt(result);
	}
	/**
	 * 虚数を含んだベクトルのpノルムを計算する
	 * @param x ベクトル(1次元配列)
	 * @param p 数値
	 * @return 第一引数のpノルム
	 */
	static double vecNormP(ImaginaryData[] data ,int p){
		double result = 0;
		for(int i = 0;i < data.length;i++){
			result += Math.pow(abs(data[i].getReal(),data[i].getImaginary()), p);
		}
		return Math.pow(result, 1.0 / p);
	}
	/**
	 * 虚数を含んだベクトルの∞ノルムを計算する
	 * @param x
	 * @return
	 */
	static double vecNormInf(ImaginaryData[] data){
		double result = abs(data[0].getReal(),data[0].getImaginary());
		for(int i = 1 ; i < data.length;i++){
			double dep = abs(data[i].getReal(),data[i].getImaginary());
			if( result < dep){
				result = dep;
			}
		}
		return result;
	}
}
