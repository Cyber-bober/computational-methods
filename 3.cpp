#include <bits/stdc++.h>
using namespace std;

double lagrange(double xp, double x[], double f[], int n) {
    double L = 0;
    for (int i = 0; i < n; i++) {
        double term = f[i];
        for (int j = 0; j < n; j++)
            if (j != i)
                term *= (xp - x[j]) / (x[i] - x[j]);
        L += term;
    }
    return L;
}

int main() {
    double x[]  = {1, 2, 3, 4};
    double f[]  = {3, 10, 29, 66};

    double x1[] = {1.5, 2.5, 3.5, 4.5};
    double f1[] = {5.375, 17.625, 44.875, 93.125};

    int n = 4;

    double xp = 2.5;

    cout << "По основным точкам:      L(" << xp << ") = " << lagrange(xp, x,  f,  n) << endl;
    cout << "По промежуточным точкам: L(" << xp << ") = " << lagrange(xp, x1, f1, n) << endl;

    return 0;
}
