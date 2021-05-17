#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

double integral(const double x) {
	return 1.0 / (1.0 + (25.0 * (x * x)));
	//Have defined a function to return the f(x) of the integral provided in our question
}

const double I_exact = (2.0 / 5.0) * double(atan(5.0));
//This is the exact value of the definite integral when integrated, to be used when comparing the exact answer of the integration with the methods computed numerically

double inner_product(const vector<double>& u, const vector<double>& v);

double Trapezium_rule(const double a, const double b, const double N, double func(const double)) {

	const double delta_x = (b - a) / double(N);
	//Declaring a variable that computes delta x. To be used when working out the weights
	//cout << "delta x has value of " << delta_x << endl;

	vector<double> x(N + 1);
	//Declaring a vector x of size N+1 to hold entries of double numbers for the values of xi
	for (int i = 0; i < x.size(); i++) {
		x[i] = double(a) + double(i) * (double(b) - double(a)) / double(N);
		//I am computing x[i] using the equidistant nodes formula mentioned in lectures { xi = a + i*(b-a)/N, i = 0,1,....,N
	}

	vector<double> w(N + 1);
	//Declaring a vector w of size N+1 to hold entries of double numbers for the weights wi
	for (int i = 0; i < w.size(); i++) {
		if (i == 0 || i == N) {
			w[i] = delta_x / 2.0;
		}
		else {
			w[i] = delta_x;
		}
		//I have made an if and else statement to compute the weights wi for when i = 0 or i = N, and 1 <= i <= N-1
	}

	//cout << "i" << "\t" << setw(15) << left << "x" << setw(15) << left << " w " << setw(18) << left << "\tfx" << endl;

	vector<double> fx(N + 1.0);
	for (int i = 0; i < fx.size(); i++) {
		fx[i] = func(x[i]);
		//This vector computes the value of fx[i] by inputting the values of x[i] into the function
		//cout << i << "\t" << setprecision(8) << setw(15) << left << x[i] << "\t" << setw(15) << left << w[i] << "\t" << setw(15) << left << fx[i] << endl;
	}

	const double I_trapezium = inner_product(w, fx);
	//To compute the integral I using the trapezium rule, I use the function defined in Question 2a and compute the inner product of vectors wi and fi

	cout << "The value of the integral computed using the composite trapezium rule is " << setprecision(12) << I_trapezium << endl;

	//cout << "The exact value of the integral is " << setprecision(12) << I_exact << endl;

	cout << "The difference between the numerical integral and the exact integral is " << "";

	return I_trapezium - I_exact;
	//The function returns the difference between the Integral computed using the trapezium rule and the exact integral
}


int main() {
	const double a = -1.0;
	//Lower bound of the definite integral
	const double b = 1.0;
	//Upper bound of the definite integral
	const double N = 63.0;
	//Total number of equidistant points 

	cout << Trapezium_rule(a, b, N, integral);
	//The trapezium rule accepts arguments: lower bound of the integral, upper bound of the integral, total number of points and the f(x) of the integral being computed 
	cout << endl;

}

double inner_product(const vector<double>& u, const vector<double>& v) {
	double sum1;
	sum1 = 0.0;

	if (u.size() == v.size())
		for (int i = 0; i < u.size(); i++) {
			sum1 += double(u[i]) * double(v[i]);
		}
	else {
		cout << "size of vectors aren't equal, therefore, cannot compute." << endl;
	}
	return sum1;
}
//This is the function inner_product. This codes functionality was demonstrated in question 2a
