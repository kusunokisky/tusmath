package sjtools;


class ImaginaryCalc extends Calc {
	private ImaginaryCalc() {
		throw new UnsupportedOperationException();
	}
	/**
	 * 複素共役を返す
	 * @param x 複素数
	 * @return 複素共役
	 */
	static ImaginaryData conjugate(ImaginaryData x){ 
		return new ImaginaryData(x.getReal(), -x.getImaginary());
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
	 * 複素数の加算を行う
	 * @param x 複素数
	 * @param y 複素数
	 * @return 加算された複素数
	 */
	static ImaginaryData add(ImaginaryData x,ImaginaryData y){
		return new ImaginaryData(x.getReal() + y.getReal(), x.getImaginary() + y.getImaginary());
	}
	/**
	 * 複素数の減算を行う
	 * @param x 複素数
	 * @param y 複素数
	 * @return 減算された複素数
	 */
	static ImaginaryData subtract(ImaginaryData x,ImaginaryData y){
		return new ImaginaryData(x.getReal() - y.getReal(), x.getReal() - y.getReal());
	}
	/**
	 * 複素数の乗算を行う
	 * (a+bi) * (c + di) = (ac-bd)+(ad+bc)i
	 */
	static ImaginaryData multiple(ImaginaryData x,ImaginaryData y){
		return new ImaginaryData(x.getReal()*y.getReal() - x.getImaginary()*y.getImaginary(),
								x.getReal()*y.getImaginary() + x.getImaginary()*y.getReal());
	}
	/**
	 * 複素数の除算を行う
	 * @param x 分子とする複素数
	 * @param y 分母とする複素数
	 * @return y / x
	 */
	static ImaginaryData division(ImaginaryData x,ImaginaryData y){
		double a = y.getReal();
		double b = y.getImaginary();
		double denom = Math.pow(a, 2) + Math.pow(b, 2);
		return new ImaginaryData( (a*x.getReal() + b*x.getImaginary() ) / denom,
								  (a*x.getImaginary() - b*x.getReal() ) / denom);
	}
	/**
	 * 複素数ベクトルの加算を計算する
	 * @param x 複素1次元配列
	 * @param y 複素1次元配列
	 * @return 複素1次元配列
	 */
	static ImaginaryData[] addVec(ImaginaryData[] x,ImaginaryData[] y){
		if(x.length != y.length)throw new UnsupportedOperationException();
		ImaginaryData[] result = new ImaginaryData[x.length];
		for(int i = 0;i < result.length;i++){
			result[i] = add(x[i],y[i]);
		}
		return result;
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
