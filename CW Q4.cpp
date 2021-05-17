#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <tuple>
using namespace std;

tuple<double, double> step_rk4(const double x, const double z, const double h, const double step,
    double z_prime(const double, const double, const double),
    double h_prime(const double, const double, const double)) {
    //This function is declared to be a tuple so that there can exist two return values, which are accessed by their positioning in the tuple
    //http://www.cplusplus.com/reference/tuple/

    const double k1 = step * z_prime(z, h, x);
    const double l1 = step * h_prime(z, h, x);
    const double k2 = step * z_prime(z + k1 / 2.0, h + l1 / 2.0, x + step / 2.0);
    const double l2 = step * h_prime(z + k1 / 2.0, h + l1 / 2.0, x + step / 2.0);
    const double k3 = step * z_prime(z + k2 / 2.0, h + l2 / 2.0, x + step / 2.0);
    const double l3 = step * h_prime(z + k2 / 2.0, h + l2 / 2.0, x + step / 2.0);
    const double k4 = step * z_prime(z + k3, h + l3, x + step);
    const double l4 = step * h_prime(z + k3, h + l3, x + step);

    return { z + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0, h + (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0 };
    //This code multiplies each k and l by 'step' so that the code is neater
    //The first element in the tuple increments the value of z in each time step by (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    //The second element in the tuple increments the value of h in each time step by (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0
    //This method of computing RK4 is defined in lectures and each k1,k2,k3,k4 are 4 slope estimates between xi and xi+1 and the value of z and h are incremented by the average of the 4 slope estimates in the manner specified by the RK4 method to give a more accurate approximation to the curve
    //This is the general method for computing the next step of the fourth order Runge-Kutta method in a system of ODE's and I have left it in this manner so that no matter what system of ODE's need to be computed they can easily be computed by simply defining them and inputting them into the step-RK4 function
}

double z_prime(const double z, const double h, const double x) {
    if (x == 0) {
        //As the equation is singular at x = 0, we use L'hopital's rule to work out the value of (2/x)h'(x) at x = 0. Therefore, after differentiating both the numerator and denominator we get 2h''(x)/1 and after rearranging our second order ODE we calculate the value of h''(0) to be equal to -1/3 
        return -1.0 / 3.0;
    }
    else {
        return ((-2.0 / x) * z) - h;
    }
}

double h_prime(const double z, const double h, const double x) {
    return z;
}

class Error_norm {
private:
    int m;
public:
    Error_norm() : m(1) {}
    Error_norm(int a) : m(a) {}
    double operator()(const vector<double> u, const
        vector<double> v) const;
};

int main() {

    const double x0 = 0.0;
    //Initial value of x
    const double xN = 3.14159265358979;
    //Final value of x
    const int N = 100;
    //Number of points
    const double step = (xN - x0) / double(N);

    const double z0 = 0.0;
    //Initial condition for z
    vector<double> z(N + 2, z0);
    //Declared a double vector z which has an initial element z0 and a size of N+2

    const double h0 = 1.0;
    //Initial condition for h
    vector<double> h(N + 2, h0);
    //Declared a double vector h which has an initial element h0 and a size of N+2

    vector<double> x(N + 2, x0);
    //Declared a double vector x which has an initial element x0 and a size of N+2

    for (int i = 0; i <= N; i++) {
        auto [znext, hnext] = step_rk4(x[i], z[i], h[i], step, z_prime, h_prime);
        //Given our initial value of x, we compute the value of z[i + 1] and h[i + 1] in the next time step of the RK4 method, defined above
        z[i + 1] = znext;
        h[i + 1] = hnext;
        //Our initial value of z0 and h0 were defined initially, therefore the RK4 method computes the values of the next z and h
        x[i + 1] = (double(i) + 1.0) * (xN - x0) / double(N);
        //Have defined the x[i + 1] in this manner because this will compute x by multiplying i+1 with the equidistant formula so even if the function were to be changed, the elements of the vector x would not need to be changed 
        //Once one time step has been computed and we have worked out the value of z[i + 1] and h[i + 1] we then update our value of x so that the RK4 method can work out the values of the next z and h
    }

    vector<double> h_exact(N + 2, 1.0);
    //As sin(x)/x at x = 0 is indeterminate we use L'Hopital's rule to compute its value at x = 0. Therefore, by differentiating both the numerator and denominator we get cos(x)/1 and at x = 0 this has a value of one. 
    //For this reason I have declared my vector with double variables to have values of 1.0 and these values will be updated with my for loop as required
    for (int i = 1; i <= N; i++) {
        //I have defined my for loop to be from 1 to N as the initial value of h_exact is 1.0 and already defined
        h_exact[i] = sin(x[i]) / x[i];
        //sin(x) is a double function by default and I have declared my vectors to be double so for this reason I don’t need to convert them to doubles in the fraction
        //cout << "h_exact " << i << "\t" << h_exact[i] << endl;
    }

    cout << "i" << "\t" << setw(23) << left << "x" << setw(20) << left << " h " << setw(20) << left << "\te" << endl;

    vector<double> e(N + 2);
    for (int i = 0; i <= N; i++) {
        e[i] = h[i] - h_exact[i];
        if (i % 10 == 0) {
            //prints values of x, h(x), e(x) for t = 0, 10, .., 100 by setting an if statement to only return those values of i which when divided by 10 leaves a remainder of 0
            cout << i << "\t" << setprecision(10) << setw(20) << left << x[i] << "\t" << setw(20) << left << h[i] << "\t" << setw(20) << left << e[i] << endl;
        }
    }
    cout << endl;

    Error_norm error1(1);
    double error_norm_l1 = error1(e, e);
    cout << "error norm l1 = " << error_norm_l1 << endl;
    cout << endl;

    Error_norm error2(2);
    double error_norm_l2 = error2(e, e);
    cout << "error norm l2 = " << error_norm_l2 << endl;
    //Have explained the uses of classes in Question 2

}

double Error_norm::operator()(const vector<double> u, const
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

