#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double rofx(const double x) {
    //I have specified the return type of my function rofx to be double before defining it

    double eps = pow(double(10), double(-12));
    //Initialising a tolerance parameter of 10^-12 using the pow function located in the <cmath> library 
    //The value of epsilon is an extremely small number, therefore, I will take this number to be equivalent to zero and use it accordingly in my for loop

    double r_n = double(0.8) * x;
    //Initialised the variable r_n to be my initial guess of [r(n)]  

    double r_n1 = x - double(log(abs(r_n - 1.0)));
    //Initialised the variable r_n1 to be equal to x – ln|r(n) – 1|    

    double iterations = 0.0;
    while (double(abs(r_n1 - r_n)) >= eps) {
        //The while loop will keep running until the condition (abs(r_n1 - r_n) >= eps) is met
        //When the distance between rn+1 and rn is less than 10^-12, we can assume that the value of r has converged as the distance between r(n) and r(n+1) is less than the tolerance level

        r_n = r_n1;
        //After one iteration the value of r(n) will become r(n+1) and after the second iteration will be r(n+2) and so on

        r_n1 = x - double(log(abs((r_n - 1))));
        //After one iteration the value of r(n+1) will become r(n+2) and after the second iteration will be r(n+3) and so on

        iterations++;
        //Initialised a variable ‘iterations’ before my while loop which increments its value by 1 each time the while loop runs to count how many times the while loop has run
    }

    cout << "The function rofx converged in " << iterations << " iterations. " << endl;

    const double error = (r_n1 + double(log(abs(r_n1 - 1)))) - x;

    //Initialised a variable called error to compute the difference between x and r(n+1) + ln|rn – 1|

    cout << "The equation x = r + ln|r - 1| is satisfied when substituting rn+1 into the equation as we get " << error << endl;


    /* if (double((r_n1 + log(abs((r_n1 - 1)))) - x) <= eps)
         //I have adopted this method as the digits of accuracy won't be precise, therefore, I believe the method I have adopted for checking whether r can be substituted into the Eq. is effective as it accounts for minor discrepancies
         cout << "The equation x = r + ln|r - 1| is satisfied when substituting r in as we get " << r_n1 + log(abs((r_n1 - 1))) << " which is the same value as x" << endl;
     else
         cout << "The equation x = r + ln|r - 1| is not satisfied when substituting r into the equation" << endl;*/

    return r_n1;
}


int main() {
    double x;
    cout << "Input a coordinate value for x which is strictly greater than 2 " << endl;
    cin >> x;
    cout << "The final value of r is: " << setprecision(16) << rofx(x) << endl;
}
