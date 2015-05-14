package ks3;
/**
 * 計算数学３に対応した数値計算アルゴリズム
 *
 */
class CalcNumerical {

	CalcNumerical() {
		// TODO 自動生成されたコンストラクター・スタブ
	}
	/**
	 * ラグランジュ関数を用いて係数を算出し、Pn(x)の値を求める
	 */
	static double[] newtonCoefficient(double startPoint,double endPoint,int n){
		double[] coe = new double[ n+1 ];
		double[] ff = new double[ n+1 ];
		double[] x  = new double[ n+1 ];
		for(int i = 0;i <= n ; i++){		
			x[i] = startPoint + (endPoint - startPoint) / (double)n * (double)i;
			ff[i] = x[i] * x[i] * x[i] * x[i];
			System.out.printf("x[%2d] = %10.7f      f[%2d] = %10.7f\n",i,x[i],i,ff[i]);
		}
		for(int i = 0;i<= n; i++){	
			coe[i] = ff[i];
			for(int j = 0;j < i;j++){
				coe[i] -= coe[j];
				coe[i] /= x[i] - x[j];
			}
			System.out.printf("c[%2d] = %10.7f\n",i,coe[i]);
		}
		double p = coe[n];
		double t = (x[n-1] + x[n]) * 0.5;//初期値
		for(int i = n-1;i >= 0;i--){
			p *= t - x[i];
			p += coe[i];
		}
		System.out.printf("P_%2d[%10.7f] = %10.7f\n" ,n,t,p);
		System.out.println(t*t*t*t);
		return coe;
				
	}
	/**
	 * エルミート補間関数を用いて係数を算出し、Pn(x)の値を求める
	 */
	static double[] hermitianCoefficient(double startPoint,double endPoint,int n){
		double[] coe = new double[ n+1 ];
		double[] ff = new double[ n+1 ];
		double[] d = new double[n+1];
		double[] gg = new double[n+1];
		double[] x  = new double[ n+1 ];
		for(int i = 0;i <= n ; i++){		
			x[i] = startPoint + (endPoint - startPoint) / (double)n * (double)i;
			ff[i] = x[i] * x[i] * x[i] * x[i];
			gg[i] = 4 * x[i] * x[i] * x[i];
			System.out.printf("x[%2d] = %10.7f      f[%2d] = %10.7f      f'[%2d] = %10.7f\n",i,x[i],i,ff[i],i,gg[i]);
		}
		for(int i = 0;i<= n; i++){	
			coe[i] = ff[i];
			d[i] = gg[i];
			for(int j = 0;j < i;j++){
				coe[i] -= coe[j];
				coe[i] /= x[i ]- x[j];
				d[i] -= coe[j];
				d[i] /= x[i ]- x[j];
				coe[i] -= d[j];
				coe[i] /= x[i ]- x[j];
				d[i] -= coe[j];
				d[i] /= x[i ]- x[j];
			}
			System.out.printf("c[%2d] = %10.7f    d[%2d] = %10.7f\n",i,coe[i],i,d[i]);
		}
		double p = 0;//coe[n]
		double t = (x[n-1] + x[n]) * 0.5;//初期値
		for(int i = n-1;i >= 0;i--){
			p *= t - x[i];
			p += d[i];
			p *= t - x[i];
			p += coe[i];
		}
		System.out.printf("P_%2d[%10.7f] = %10.7f\n" ,n,t,p);
		System.out.println(t*t*t*t);
		return coe;
				
	}

}
