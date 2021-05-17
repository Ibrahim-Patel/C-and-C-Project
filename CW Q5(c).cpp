#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

double integral(const double x);

const double I_exact = (2.0 / 5.0) * atan(5);

double inner_product(const vector<double>& u, const vector<double>& v);

double Clenshaw_curtis(const double a, const double b, const double N, double func(const double)) {

	const double delta_x = (b - a) / double(N);
	//cout << "delta x has value of " << delta_x << endl;

	vector<double> theta(N + 1);
	//Declaring a vector theta of size N+1 to hold entries of double numbers for the values of theta_i
	for (int i = 0; i < theta.size(); i++) {
		const double pi = 3.14159265358979;
		theta[i] = (double(i) * pi) / double(N);
		//This for loop computes the value of theta at each i using the formula provided in the question
	}

	vector<double> x(N + 1);
	for (int i = 0; i < x.size(); i++) {
		x[i] = double(-cos(theta[i]));
	}

	vector<double> w(N + 1);
	//Declaring a vector ‘w’ of size N+1 to hold entries of double numbers for the values of the weights wi

	for (int i = 0; i < w.size(); i++) {
		if (i == 0 || i == N) {
			w[i] = 1.0 / (double(N) * double(N));
		}
		else {
			double sum = 0.0;
			for (int k = 1; k <= ((double(N) - 1.0) / 2.0); k++) {
				sum += (2.0 * cos(2.0 * double(k) * theta[i])) / (4.0 * (double(k) * double(k)) - 1.0);
				//To compute the weights wi for 1 <= i <= N – 1 I first compute the value of the sum (2cos(2k(thetai))/4(k^2)-1 from k = 1 to (N – 1)/2
				w[i] = (2.0 / double(N) * (1.0 - sum));
				//After computing the sum, I work out the value of wi by doing 2/N*(1 – sum)
				//I have made an if and else statement to compute the weights wi for when i = 0 or i = N, and 1 <= i <= N-1
			}
		}
	}

	//cout << "i" << "\t" << setw(15) << left  << "theta" << setw(15) << left << "x" << setw(15) << left << "  w " << setw(18) << left << "\tfx" << endl;


	vector<double> fx(N + 1);
	for (int i = 0; i < fx.size(); i++) {
		fx[i] = func(x[i]);
		//cout << i << "\t" << setprecision(8) << setw(15) << left << theta[i] << setw(15) << left << x[i] << "\t" << setw(15) << left << w[i] << "\t" << setw(15) << left << fx[i] << endl;

	}

	const double I_clenshawcurtis = inner_product(w, fx);
	//To compute the integral I using the Clenshaw-Curtis quadrature rule, I use the function defined in Question 2a and compute the inner product of vectors wi and fi


	cout << "The value of the integral computed using the Clenshaw-Curtis quadrature rule is " << setprecision(12) << I_clenshawcurtis << endl;

	//cout << "The exact value of the integral is " << setprecision(12) << I_exact << endl;

	cout << "The difference between the numerical integral and the exact integral is " << "";

	return I_clenshawcurtis - I_exact;
	//The function returns the difference between the Integral computed using the Clenshaw-Curtis quadrature rule and the exact integral

}


int main() {
	const double a = -1.0;
	const double b = 1.0;
	const double N = 63.0;

	cout << Clenshaw_curtis(a, b, N, integral);
	//The Clenshaw-Curtis quadrature rule accepts arguments: lower bound of the integral, upper bound of the integral, total number of points and the f(x) of the integral being computed
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

double integral(const double x) {
	return 1.0 / (1.0 + (25.0 * (x * x)));
}