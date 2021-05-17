#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <tuple>
using namespace std;

tuple<double, double> step_rk2(const double t, const double q, const double p, const double step,
    double q_prime(const double, const double, const double),
    double p_prime(const double, const double, const double)) {
    //This function is declared to be a tuple so that there can exist two return values which are accessed by their positioning in the tuple
    //http://www.cplusplus.com/reference/tuple/

    const double k1 = step * q_prime(q, p, t);
    const double l1 = step * p_prime(q, p, t);
    const double k2 = step * q_prime(q + k1 / 2.0, p + l1 / 2.0, t + step / 2.0);
    const double l2 = step * p_prime(q + k1 / 2.0, p + l1 / 2.0, t + step / 2.0);

    return { q + k2, p + l2 };
    //This code multiplies each k and l by 'step' so that the code is neater
    //The first element in the tuple increments the value of q in each time step by k2
    //The second element in the tuple increments the value of p in each time step by l2
    //This method of computing Runge-Kutta 2 midpoint method is defined in lectures and we approximate p and q at the midpoint of the interval to give us the slope k2 and l2, we then use this slope k2 and l2 to compute the next value of q and p   
    //This is the general method for computing the next step of the Runge-Kutta 2 midpoint method in a system of ODE's and I have left it in this manner so that no matter what system of ODE's need to be computed they can easily be computed by simply defining them and inputting them into the step-RK2 function
}

double q_prime(const double q, const double p, const double t) {
    return p;
}

double p_prime(const double q, const double p, const double t) {
    return -q;
}

int main() {
    const double t0 = 0.0;
    //Initial time
    const double tN = 100.0;
    //Final time
    const int N = 1000;
    //Number of points is N+1; therefore, the number of intervals is N
    const double step = (tN - t0) / double(N);
    //Time step is 0.1, therefore, minus the difference between the initial time and final time and divide by N

    const double q0 = 0;
    //Initial condition for q0
    vector<double> q(N + 2, q0);
    //Declared a vector q which has an initial element q0 and a size of N+2

    const double p0 = sqrt(2);
    //Initial condition for p0
    vector<double> p(N + 2, p0);
    //Declared a vector p which has an initial element p0 and a size of N+2

    vector<double> t(N + 2, t0);
    //Declared a vector t which has an initial element t0 and a size of N+2

    for (int i = 0; i <= N; i++) {
        auto [qnext, pnext] = step_rk2(t[i], q[i], p[i], step, q_prime, p_prime);
        q[i + 1] = qnext;
        p[i + 1] = pnext;
        //Our initial value of p0 and q0 were defined initially, therefore, the RK2 midpoint method computes the values of the next q and p
        t[i + 1] = (double(i) + 1.0) * (double(tN) - double(t0)) / double(N);
        //Have defined the t[i + 1] in this manner because this will compute t by multiplying i+1 with the equidistant formula so even if the function were to be changed this would not need to be changed 
        //Once one time step has been computed and we have worked out the value of q[i + 1] and p[i + 1] we then update our value of t so that the RK2 midpoint method can work out the next values of q and p
    }

    vector<double> E(N + 2);
    for (int i = 0; i <= N; i++) {
        E[i] = (1.0 / 2.0) * ((p[i] * p[i]) + (q[i] * q[i]));
    }

    vector<double> e(N + 2);
    for (int i = 0; i <= N; i++) {
        e[i] = E[i] - E[0];
    }

    cout << setw(7) << left << "t" << setw(18) << left << "Position q(t)" << setw(18) << left << "Momentum p(t) " << setw(18) << left << "Energy E(t)" << setw(18) << left << "Difference e(t)" << setw(18) << endl;
    cout << setprecision(7) << setw(7) << left << 0 << setw(18) << left << q[0] << setw(18) << left << p[0] << setw(18) << left << E[0] << setw(18) << left << e[0] << setw(18) << endl;

    for (int j = 10; j <= N; j *= 10) {
        cout << setprecision(12) << setw(7) << left << t[j] << setw(18) << left << q[j] << setw(18) << left << p[j] << setw(18) << left << E[j] << setw(18) << left << e[j] << setw(18) << endl;
    }
}
