#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
using namespace std;

double integral(const double x);

const double I_exact = (2.0 / 5.0) * atan(5);

double max_c(const double a, const double b, double func(const double)) {
	//I have defined this function as such, to return the largest value of f(x). 
	//Rather than computing the value analytically I have adopted this method because it can be used to calculate the maximum value any function can have

	const double arbitrary_num = 1000.0;
	//Declaring a variable arbitrary_num to declare the size of my vector as well as to compute the maximum f(x) using arbitrary_num of points (1001)

	vector<double> x(arbitrary_num);
	for (int i = 0; i < x.size(); i++) {
		x[i] = a + i * (b - a) / double(arbitrary_num);
		//Calculating points for x[i] using the equidistant method
	}

	vector<double> fx(arbitrary_num);
	for (int i = 0; i < fx.size(); i++) {
		fx[i] = func(x[i]);
		//This will compute all of the values of the function f(x) by inputting the values of xi 
	}
	return *max_element(fx.begin(), fx.end());
	//This returns the maximum value of my vector fx and hence the maximum value of any function
	//I have used the asterisk (*) before the max_element() function as this returns an iterator.
	//https://stackoverflow.com/questions/9874802/how-can-i-get-the-max-or-min-value-in-a-vector
}
double Monte_carlo(const double a, const double b, double func(const double)) {

	double area_rec = max_c(a, b, func) * (b - a);
	//Declared the variable area_rec to compute the area of the rectangle between a <= x <= b and 0 <= fx <= c

	const int seed = 93;
	//I am setting my seed as 93 to initialise the sequence of uniform random numbers 
	mt19937_64 rand(seed);
	//I am using the Mersenne Twister algorithm to generate my uniform random numbers because it is better than the linear congruential method
	//mt19937_64 is the object that will generate the uniform random numbers
	//I have declared my variable to be called rand and have supplied my seed to it
	//Combining the two will generate uniform random numbers
	uniform_real_distribution<double> X(a, b);
	//I am generating uniform random numbers for X between the values of a and b
	uniform_real_distribution<double> W(0, max_c(a, b, func));
	//I am generating uniform random numbers for W between 0 and the maximum of f(x) (a <= x <= b)

	const int Monte_carlo_N = 10000;
	//The number of times I want my for loop to run and hence also the amount of uniform random variables I want to generate
	int M = 0;
	for (int i = 0; i < Monte_carlo_N; i++) {
		const double x = X(rand);
		//This will generate uniform random numbers for x between a and b
		const double w = W(rand);
		//This will generate uniform random numbers for w between 0 and c_max
		//Therefore, these combined will generate points randomly within the region specified (in a rectangle where -1 < x < 1 and 0 < y < 1)

		if (w <= func(x))
			//The if statement will ensure that all of the points that are under the function curve will be accounted for by the counter M as M will increase in value by 1 if and only if W <= f(X) 
			M++;
	}

	const double I_MC = (double(M) * area_rec) / double(Monte_carlo_N);

	cout << setprecision(16) << "The value of the integral computed using the hit and miss Monte Carlo method with N = 10000 is " << I_MC << endl;

	//cout << "The exact value of the integral is " << setprecision(12) << I_exact << endl;

	cout << "The difference between the numerical integral and the exact integral is " << "";

	return I_MC - I_exact;
}
int main() {
	const double a = -1.0;
	const double b = 1.0;
	const double N = 63.0;

	cout << (setprecision(16)) << Monte_carlo(a, b, integral) << endl;

	//cout << max_c(a, b, integral);
}

double integral(const double x) {
	return 1.0 / (1.0 + (25.0 * (x * x)));
}
