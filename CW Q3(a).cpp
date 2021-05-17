#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

double function_exp(const double x) {
	//I have declared a function for f(x) which returns e^-(x^2) to test my programs functionality 
	return exp(-(x * x));
}

double f_prime_exp(const double x) {
	// I have declared a function for f’(x) which returns -2x * e^-(x^2) to test my programs functionality
	return -2.0 * double(x) * double(exp(-(x * x)));
}

vector<double> error(const int N, double func(const double), double diff_func(const double)) {
	//I have created a vector function to return the error vector, this is so that I can use this function for Question 3b
	//I have also defined it in this manner so that 'vector<double> error' can be used for any function such as e^-x^2 (as given in the question)
	//Also, because of the layout of my function I can output any vector x, fx, f'x numeric, f'x analytic and error vector easily for any function and any number of points (N+1)
	//From all the methods I experimented with I found this method to be the neatest and most efficient

	double delta_x = 2.0 / double(N);
	//The value for delta x, given in the question

	vector<double> x(N + 1);
	for (int i = 0; i < x.size(); i++) {
		x[i] = (2.0 * double(i) - double(N)) / double(N);
		//cout << "Values of xi for each element in the vector \t" << i << "\t" << x[i] << endl;
	}
	//The values of xi are assumed to be located at the grid-points (2i-N) / N, for each i the value is stored in the i’th entry of the vector x
	//cout << endl;

	vector<double> fx(N + 1);
	for (int i = 0; i < fx.size(); i++) {
		fx[i] = func(x[i]);
		//cout << "Values of fx for each element in the vector \t" << i << "\t" << fx[i] << endl;
	}
	//This vector stores the values of f(xi) by inputting the values xi, computed in the previous vector, into the argument func 
	//cout << endl;

	vector<double> f_prime_numeric(N + 1);
	for (int i = 0; i < f_prime_numeric.size(); i++) {
		if (i == 0) {
			f_prime_numeric[0] = ((-3.0 * fx[0]) + (4.0 * fx[1]) - (fx[2])) / (2.0 * delta_x);
		}

		else if (i > 0 && i <= (N - 1)) {
			f_prime_numeric[i] = (fx[i + 1] - fx[i - 1]) / (2.0 * delta_x);
		}

		else {
			f_prime_numeric[N] = ((fx[N - 2]) - (4.0 * fx[N - 1]) + (3.0 * fx[N])) / (2.0 * delta_x);
		}
		//cout << "Values of f'x numeric for each element in the vector \t" << i << "\t" << f_prime_numeric[i] << endl;
	}
	//This vector stores the numerical values of f'(x) using the finite differences method covered in the interpolating lecture
	//cout << endl;

	vector<double> f_prime_analytic(N + 1);
	for (int i = 0; i < f_prime_analytic.size(); i++) {
		f_prime_analytic[i] = diff_func(x[i]);
		//cout << "Values of f'x analytic for each element in the vector \t" << i << "\t" << f_prime_analytic[i] << endl;
	}
	//This vector stores the analytic values of f'(x) by inputting the values of xi, into the differential function argument 
//cout << endl;

	vector<double> error_vector(N + 1);
	for (int i = 0; i < error_vector.size(); i++) {
		error_vector[i] = f_prime_numeric[i] - f_prime_analytic[i];
		cout << setprecision(12) << "Error values for each element in the vector \t " << i << "\t" << error_vector[i] << endl;
	}
	//This vector stores the error values, which is f'numerical(xi) - f'analytical(xi)
	return error_vector;
	//Returns the error vector computed in the previous step
	//This is the only method I found that will allow me to output vectors, as the void function is not recognised as a vector when trying to input the error vector into the weighted norm object defined for question 3(b)
}


int main() {
	
	//This is the vector I have defined to output the answers for Question 3(a) in the neatest manner possible
	//I have also done it as such so that if I wanted to output the vectors within the function for more points, then I can easily do so by changing the limit within the for loop
	//Simply run this vector for however many values of N and for whatever function you would like to test it with
	vector<double> Question_3a(1);
	for (int j = 16; j < 18; j *= 2) {
		error(j - 1, function_exp, f_prime_exp);
		cout << endl;
	}

}
