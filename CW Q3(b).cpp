#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

double function_exp(const double x) {
	return exp(-(x * x));
}

double f_prime_exp(const double x) {
	return -2.0 * x * exp(-(x * x));
}

vector<double> error(const int N, double func(const double), double diff_func(const double)) {

	double delta_x = 2.0 / double(N);

	vector<double> x(N + 1);
	for (int i = 0; i < x.size(); i++) {
		x[i] = (2.0 * double(i) - double(N)) / double(N);
	}

	vector<double> fx(N + 1);
	for (double i = 0; i < fx.size(); i++) {
		fx[i] = func(x[i]);
	}

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
	}

	vector<double> f_prime_analytic(N + 1);
	for (int i = 0; i < f_prime_analytic.size(); i++) {
		f_prime_analytic[i] = diff_func(x[i]);
	}

	vector<double> error_vector(N + 1);
	for (int i = 0; i < error_vector.size(); i++) {
		error_vector[i] = f_prime_numeric[i] - f_prime_analytic[i];
	}
	//This vector stores the error values, which is f'numerical(xi) - f'analytical(xi)
	return error_vector;
	//Returns the error vector computed in the previous step
}

//This is the class from Question 2(b) to work out the weighted norm of two vectors
//I have simply renamed the object to Mean_error and declared the object before int main() and defined it under int main() as this functions functionality was already demonstrated in the previous question
class Mean_error {
private:
	int m;
public:
	Mean_error() : m(1) {}
	Mean_error(int a) : m(a) {}
	double operator()(const vector<double> u, const
		vector<double> v) const;
};

int main() {

	//This is the vector I have defined to output the answers for Question 3(b) in the neatest manner possible
	//I have also done it as such so that if I wanted to output the mean error for more points (N) then I can easily do so by changing the limit within the for loop

	vector<double> Question_3b(5);
	//This for loop will allow me to output the mean error <e> for N + 1 = 16, 32, 64, 128, 256.
	for (int j = 16; j < 257; j *= 2) {
		const Mean_error N_mean_error(1);
		//In order to work out the mean error, the value for m in weighted norm m will be 1
		double error_N = N_mean_error(error(j - 1, function_exp, f_prime_exp), error(j - 1, function_exp, f_prime_exp));
		//The mean error is to be computed for vectors ei
		error_N /= ((double(j) - 1.0) + 1.0);
		//Once the sum is computed this value has to be divided by N + 1
		cout << "The mean error for N + 1 = " << j << " is " << setprecision(12) << error_N << endl;
		cout << "N^2<e> is " << setprecision(12) << double(j - 1.0) * double(j - 1.0) * error_N << endl;
		cout << endl;
	}

	cout << "As we can see from the values that have been outputted; as the value of N increases the mean error decreases, therefore as <e> decreases 1/N^2 also decreases proportionally to <e> " << endl;
	cout << " As we can see from the values of N^2<e> that have been output the value of N^2<e> is approximately constant " << endl;

}

double Mean_error::operator()(const vector<double> u, const
	vector<double> v) const {
	double l = 0.0;
	double sum2 = 0.0;

	if (u.size() == v.size()) {
		for (int i = 0; i < u.size(); i++) {
			sum2 += pow(abs(double(u[i]) * double(v[i])), double(m) / 2.0);
		}
		l = pow(sum2, 1.0 / double(m));
	}
	else {
		cout << "size of vectors aren't equal, therefore, cannot compute." << endl;
	}
	return l;
};

