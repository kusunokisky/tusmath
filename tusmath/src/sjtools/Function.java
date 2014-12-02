package sjtools;

abstract class Function {
	double coefficient[];
	Function(double coefficient[]) {
		this.coefficient = coefficient;
	}
	abstract double fx(double x);
	abstract double dfx(double x);
}
