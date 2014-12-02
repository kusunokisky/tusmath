package sjtools;

import java.util.HashMap;

class Change {
	/**
	 * 第一引数に対して36進数までの文字から数値変換用のマップを追加する
	 * @param map 対象のマップ
	 * @return 36進数までの文字から数値変換が追加されたマップ
	 */
	private static HashMap<String,String> addMap(HashMap<String,String> map){
		String[] s = {"A","B","C","D","E","F","G","H","I","J","K","L","M",
				"N","O","P","Q","R","S","T","U","V","W","X","Y","Z"};
		String i_ ;
		for(int i = 0;i < 36;i++){
			i_ = String.valueOf(i);
			if(i >= 10){
				map.put( s[i-10] , i_ );
				map.put( i_ ,  s[i-10]);
			}else{
				map.put(i_,i_);
			}
		}
		return map;
	}

	/**
	 *	N進数をM進数に変換する 
	 * @param num N進数表記された文字列
	 * @param nBase 変換前の基数
	 * @param mBase 変換後の基数
	 * @return M進数に変換された文字列
	 */

	static String nToM(String num,int nBase,int mBase){
		return Change.tenToN(Change.nToTen(num, nBase), mBase);		

	}
	/**
	 * n進数を10進数に変換する
	 * @param num 変換対象の文字列
	 * @param base 数列の基数
	 * @param map 文字を数値に変換させるためのマップ
	 * @return 10進数に変換された数
	 */
	private static double nToTen(String num ,int base,HashMap<String,String> map) throws NumberFormatException{
		if(num == null){
			throw new NumberFormatException("numが" + num + "です");
		}
		if(base < Character.MIN_RADIX){
			throw new NumberFormatException("baseが" + Character.MIN_RADIX + "よりも小さいです");
		}
		if(Character.MAX_RADIX < base){
			throw new NumberFormatException("baseが" + Character.MAX_RADIX + "よりも大きいです" );
		}
		//整数と小数に分ける
		String[] s = num.split("\\.",0);

		/*								*
		 *		 整数部分の処理			*
		 *								*/

		//正負の判定
		boolean pn = true;								//true==正 false==負
		if(s[0].indexOf("-") != -1){
			s[0] = s[0].replaceAll("-", "");			//-の文字を削除
			pn = false;									//負の判定
		}
		//10進数に変換
		double iResult = calcNnum(s[0], base, map);

		/*								*
		 * 		 小数部分の処理			*
		 * 								*/

		double dResult = 0;								//変換後の格納
		//小数部分の存在確認
		if(s.length == 2){
			//10進数に変換
			dResult = calcDnum(s[1], base, map);
		}
		//整数部分と小数部分を結合
		double result = iResult + dResult;
		//正負の設定
		if(!pn){
			result *= -1;
		}
		return result;
	}
	/**
	 * n進数を10進数に変換する
	 * @param num 変換対象の文字列
	 * @param base 数列の基数
	 * @return 10進数に変換された数
	 */
	static double nToTen(String num ,int base){
		HashMap<String, String> map = new HashMap<String, String>();
		addMap(map);
		return nToTen(num, base, map);
	}
	/**
	 * 10進数をn進数に変換する
	 * @param n 10進数表記の数値
	 * @param base 変換先の基数
	 * @param map 36進数までの文字から数値変換が追加されたマップ
	 * @return n進数に変換された文字列
	 */
	private static String tenToN(double n,int base ,HashMap<String,String> map){
		//正負の判定
		double n_ = n;
		boolean pn = true;
		if(n_  < 0){
			n_ *= -1;
			pn = false;
		}
		//整数と小数に分ける
		int in = (int)n_;
		double dn = n_ - in;
		//変換後の格納
		int[] nin = new int[10000];		//現状では1000桁とおいてあるが本来ならば対数をとって求める必要がある
		int[] ndn = new int[10000];		//
		//整数部分の処理
		//n進数へ変換
		int it = 0;
		while(in != 0&&it < 10000){
			nin[it] = in % base;
			in = in / base;
			it++;
		}
		//数値を文字に変換
		String[] nis = new String[it];
		for(int i = 0;i < nis.length;i++){
			nis[i] = map.get( String.valueOf( nin[i] ) );
		}
		//文字の結合及び反転
		String nResult = convineString( reverseString( nis ) );

		//小数部分の処理
		//n進数へ変換
		int dt = 0;
		while(dn != 0&&dt<10000){
			ndn[dt] = (int)(base*dn);
			dn = base * dn - ndn[dt];
			dt++;
		}
		String[] nds = new String[dt];
		for(int i = 0; i < nds.length;i++){
			nds[i] = map.get( String.valueOf( ndn[i] ) );
		}
		String dResult = convineString(nds);
		//正負の訂正
		if(!pn){
			nResult = convineString("-", nResult);
		}
		//結果の出力
		
		//両方が空の場合
		if(dResult.isEmpty() && nResult.isEmpty()){
			return null;
		}//整数が空の場合
		else if(nResult.isEmpty()){
			return "0." + dResult;
		}//小数が空の場合
		else if(dResult.isEmpty()){
			return nResult;
		}else{
			return nResult + "." + dResult;
		}
	}
	/**
	 * 10進数をn進数に変換する
	 * @param n 10進数表記の数値
	 * @param base 変換先の基数
	 * @return n進数に変換された文字列
	 */
	static String tenToN(double n,int base){
		HashMap<String, String> map = new HashMap<String, String>();
		addMap(map);
		return tenToN(n, base, map);
	}
	/**
	 * n進数の整数を10進数に変換
	 * @param s 変換対象の文字列
	 * @param base 変換対象の基数
	 * @param map 文字を数値に変換させるためのマップ
	 * @return 10進数に変換された数
	 */
	static private double calcNnum(String s,int base,HashMap<String, String> map){
		double sum = 0;
		String[] s_ = sepalateString(s);
		for(int i = 0;i < s_.length;i++){
			if(map.get(s_[i]) == null){
				throw new NumberFormatException("s[" + i + "]が" + s_[i] + "です");
			}
			sum = sum * base + Double.parseDouble( map.get( s_[i] ) );
		}
		return sum;
	}
	/**
	 * n進数の小数を10進数に変換
	 * @param s 変換対象の文字列
	 * @param base 変換対象の基数
	 * @param map 文字を数値に変換させるためのマップ
	 * @return 10進数に変換された数
	 */
	static private double calcDnum(String s,int base ,HashMap<String, String> map){
		double sum = 0;
		String[] s_ = sepalateString(s);
		for(int i = s_.length-1; 0 <= i;i--){
			sum += Double.parseDouble( map.get( s_[i] ) ) * Math.pow(base, -(i+1));
		}
		return sum;
	}
	/**
	 * 文字列を１文字ずつに分割
	 * @param s 分割したい文字列
	 * @return 1文字ずつに分割された文字列配列
	 */
	 static private String[] sepalateString(String s){
	 
		char c[] = s.toCharArray();
		String[] s_ = new String[c.length];
		for(int i = 0; i < s.length();i++){
			s_[i] = String.valueOf(c[i]);
		}
		return s_;
	}
	/**
	 * 文字列を反転させる
	 * @param s
	 * @return 反転された文字列
	 * @throws IllegalArgumentException
	 */
	static String reverseString(String s)throws IllegalArgumentException{
		if(s == null){
			throw new IllegalArgumentException("引数が null です");
		}
		char[] c = s.toCharArray();
		String[] s_ = new String[c.length];
		for(int i = s.length();0 < i;i--){
			s_[i-1] = String.valueOf(c[s.length() - i]);
		}
		return convineString(s_);
	}
	/**
	 * 文字列を反転させる
	 * @param s
	 * @return 反転された文字列
	 * @throws IllegalArgumentException
	 */
	static private String[] reverseString(String[] s)throws IllegalArgumentException{
		if(s == null){
			throw new IllegalArgumentException("引数が null です");
		}
		String[] s_ = new String[s.length];
		for(int i = s.length;0 < i;i--){
			if(s[s.length-i] == null){
				throw new IllegalArgumentException("配列第"+(i-1) + "が null です" );
			}
			s_[i-1] = s[s.length-i];
		}
		return s_;
	}
	/**
	 * 文字列配列を文字列に変換（結合）する
	 * @param s
	 * @return 結合された文字列
	 * @throws IllegalArgumentException
	 */
	static private String convineString(String[] s)throws IllegalArgumentException{
		StringBuffer buf = new StringBuffer();
		for(int i = 0;i < s.length;i++){
			if(s[i] == null){
				throw new IllegalArgumentException("配列"+ i + "がnullです");
			}
			buf.append(s[i]);
		}
		return buf.toString();
	}
	/**
	 * 2つの文字を結合
	 * @param s1 前半文字
	 * @param s2 後半文字
	 * @return 結合された文字
	 */
	static private String convineString(String s1,String s2){
		StringBuffer buf = new StringBuffer();
		buf.append(s1);
		buf.append(s2);
		return buf.toString();
	}
}
