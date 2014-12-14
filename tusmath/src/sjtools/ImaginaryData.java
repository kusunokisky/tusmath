package sjtools;

class ImaginaryData {
	private double real;
	private double imaginary;
	
	ImaginaryData(double real,double imaginary) {
		this.real = real;
		this.imaginary = imaginary;
	}

	double getReal() {
		return real;
	}

	void setReal(double real) {
		this.real = real;
	}

	double getImaginary() {
		return imaginary;
	}

	void setImaginary(double imaginary) {
		this.imaginary = imaginary;
	}
}
