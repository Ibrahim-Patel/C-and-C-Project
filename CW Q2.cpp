#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#define eps 0.000001
//Defined eps as such so that the compiler will replace references to these constants with the defined value at compile time 
//https://www.arduino.cc/reference/en/language/structure/further-syntax/define/
using namespace std;

double inner_product(const vector<double>& u, const vector<double>& v) {
	//The function inner_product takes u and v as arguments and u and v are passed as references so the compiler doesn’t need to copy the vectors.

	double sum1;
	sum1 = 0.0;

	if (u.size() == v.size())
		//I have made an if statement to ensure that the sizes of the vectors are the same, as the dot product of two vectors can't be computed if the sizes of the vectors are different

		for (int i = 0; i < u.size(); i++) {
			//As the size of vectors u and v in our example are 3 the conditions for sum are held as we are counting from 0 to 3 not inclusive
			//cout << u[i] << " * " << v[i] << " = " << u[i] * v[i] << endl;
			sum1 += double(u[i]) * double(v[i]);
			//At each iteration of the for loop the value of the sum is increased by the value of u[i] * v[i]. Therefore, on the final for loop the value of the sum will be equal to the inner product
			//cout << sum1 << endl;
		}
	else {
		cout << "size of vectors aren't equal, therefore, cannot compute." << endl;
	}
	return sum1;
}

class Weighted_norm {
private:
	int m;
	//I have declared my member variable m to be private so that it can only be accessed via member functions and not individually 
public:
	Weighted_norm() : m(1) {}
	//As the m’th root of the weighted norm expression can't be zero I have set the value of m to be 1 when a value for m has not been specified
	Weighted_norm(int a) : m(a) {}
	//Set this constructor as such so that the weighted norm accepts ‘a’ as a parameter for the value of m

	double operator()(const vector<double> u, const
		vector<double> v) const {
		//This member function takes two vectors as arguments and as the function won't change the member data, the member function can be declared as a const function
		double l = 0.0;
		double sum2 = 0.0;

		if (u.size() == v.size()) {
			//I have made an if statement to ensure that the sizes of the vectors are the same, as the dot product can't be computed if the sizes of the vectors are different
			for (int i = 0; i < u.size(); i++) {
				sum2 += pow(abs(double(u[i]) * double(v[i])), double(m) / 2.0);
				//At each iteration of the for loop the value of the sum is increased by the value of |u[i] * v[i]| ^ m/2
				//As I am doing a division of m/2, both these numbers are integers, therefore, they must be converted into doubles
			}
			l = pow(sum2, 1.0 / double(m));
			//I have declared the variable ‘l’ to calculate the m’th root of the sum computed in the previous for loop. As I am doing a division of 1/m, both these numbers are integers, therefore, they must be converted into doubles
		}
		else {
			cout << "size of vectors aren't equal, therefore, cannot compute." << endl;
		}
		return l;
	}
};


int main() {
	vector<double> u = { 2, 7, 2 };
	vector<double> v = { 3, 1, 4 };
	cout << "The inner product of vectors u and v is " << inner_product(u, v) << endl;

	const int m1 = 1;
	//Value we are assigning our member variable m to be
	//const Weighted_norm m1_weighted_norm;
	//If we didn’t declare m, then by default the value of m would be 1 as I set the constructor as such

	const Weighted_norm m1_weighted_norm(m1);
	//Before we use the object, Weighted_norm, we first have to declare an instance of the object, which in our case is m1_weighted_norm
	//Accessing our weighted norm object and using our parameter m1 as the input for m 
	double l1 = m1_weighted_norm(u, v);
	//Assigning vectors u and v to be the arguments in the member function
	cout << "The weighted norm for l1(u,v) is " << l1 << endl;

	const int m2 = 2;
	const Weighted_norm m2_weighted_norm(m2);
	double l2 = m2_weighted_norm(u, v);
	cout << "The weighted norm for l2(u,v) is " << l2 << endl;

	const double error = (l2 * l2) - inner_product(u, v);
	//l2 is a double variable and therefore doesn’t need to be converted

	cout << "The quantity l2(u,v)^2 equals the inner product of u*v as the error value is " << error << endl;

	/*if (abs((l2 * l2) - inner_product(u, v)) == 0) {
		//I have defined my own error value to check whether the two expressions are roughly equal as we can't use == because there'll be a very slight difference in the value between them
		cout << "the quantity l2(u,v)^2 equals " << l2 * l2 << ", thus equalling the inner product of u*v" << endl;
	}
	else
		cout << "the quantity l2(u,v)^2 does not equal the inner product of u*v" << endl;*/
}
