#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

double integral(const double x);

double inner_product(const vector<double>& u, const vector<double>& v);

const double I_exact = (2.0 / 5.0) * double(atan(5.0));

double f_prime_minus(const double a, const double b) {
	//To use the composite Hermite integration rule we need to compute [f'(a) - f'(b)] for any given function, therefore, I have defined a function to do exactly this
	vector<double> f_prime_vec(2);
	//I have defined a vector f_prime_vec to store the values of f’(a) and f’(b)
	double minus_primes = 0;
	for (int i = 0; i < f_prime_vec.size(); i++) {
		//The derivative of f’(x) = -50x / (1+25x^2)^2
		f_prime_vec[0] = (-50.0 * double(a) / (pow(1.0 + 25.0 * double((a * a)), 2.0)));
		//The derivative of f(x) at x = a is f’(a) = -50(a) / (1+25(a)^2)^2
		f_prime_vec[1] = (-50.0 * double(b) / (pow(1.0 + 25.0 * double((b * b)), 2.0)));
		//The derivative of f(x) at x = b is f’(b) = -50(b) / (1+25(b)^2)^2
		minus_primes = f_prime_vec[0] - f_prime_vec[1];
		//Finally, I compute the value of f’(a) – f’(b) by subtracting f_prime_vec[0] and f_prime_vec[1]
	}
	return minus_primes;
}

double Hermite_integ(const double a, const double b, const double N, double func(const double), double diff_func(const double, const double)) {

	const double delta_x = (double(b) - double(a)) / double(N);
	//cout << "delta x has value of " << delta_x << endl;

	vector<double> x(N + 1);
	for (int i = 0; i < x.size(); i++) {
		x[i] = double(a) + double(i) * (double(b) - double(a)) / double(N);
	}

	vector<double> w(N + 1.0);
	for (int i = 0; i < w.size(); i++) {
		if (i == 0 || i == N) {
			w[i] = double(delta_x) / double(2.0);
		}
		else {
			w[i] = double(delta_x);
		}
	}

	//cout << "i" << "\t" << setw(15) << left << "x" << setw(15) << left << " w " << setw(18) << left << "\tfx" << endl;

	vector<double> fx(N + 1);
	for (int i = 0; i < fx.size(); i++) {
		fx[i] = func(x[i]);
		//cout << i << "\t" << setprecision(8) << setw(15) << left << x[i] << "\t" << setw(15) << left << w[i] << "\t" << setw(15) << left << fx[i] << endl;
	}

	const double I_hermite = inner_product(w, fx) + ((delta_x * delta_x) / 12.0) * double(diff_func(a, b));
	//To compute the integral I using the Hermite integration rule, I use the function defined in Question 2a and compute the inner product of vectors wi and fi and add delta_x^2/12*[f’(a) – f’(b)]

	cout << "The value of the integral computed using the composite Hermite integration rule is " << setprecision(12) << I_hermite << endl;

	//cout << "The exact value of the integral is " << setprecision(12) << I_exact << endl;

	cout << "The difference between the numerical integral and the exact integral is " << "";

	return I_hermite - I_exact;
	//The function returns the difference between the Integral computed using the Hermite integration rule and the exact integral
}



int main() {

	const double a = -1.0;
	//Lower bound of the definite integral
	const double b = 1.0;
	//Upper bound of the definite integral
	const double N = 63.0;
	//Total number of equidistant points

	cout << Hermite_integ(a, b, N, integral, f_prime_minus);
	//The Hermite integration rule accepts arguments: lower bound of the integral, upper bound of the integral, total number of points, the f(x) of the integral being computed and the difference between f’(a) and f’(b)
	cout << endl;
}

double integral(const double x) {
	return 1.0 / (1.0 + (25.0 * (x * x)));
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
