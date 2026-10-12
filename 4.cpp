#include <bits/stdc++.h>
using namespace std;

// Таблица разделённых разностей
vector<vector<double>> divided_diff(double x[], double f[], int n) {
    vector<vector<double>> t(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) t[i][0] = f[i];
    for (int k = 1; k < n; k++)
        for (int i = 0; i + k < n; i++)
            t[i][k] = (t[i+1][k-1] - t[i][k-1]) / (x[i+k] - x[i]);
    return t;
}

// Значение многочлена Ньютона в точке xp
double newton(double xp, double x[], vector<vector<double>>& t, int n) {
    double res = t[0][0], prod = 1;
    for (int k = 1; k < n; k++) {
        prod *= (xp - x[k-1]);
        res += t[0][k] * prod;
    }
    return res;
}

void print_table(double x[], vector<vector<double>>& t, int n) {
    cout << "Таблица ТРР:\n";
    for (int i = 0; i < n; i++) {
        cout << setw(8) << x[i];
        for (int k = 0; k < n - i; k++)
            cout << setw(12) << t[i][k];
        cout << "\n";
    }
}

int main() {
    double x[]  = {1, 2, 3, 4};
    double f[]  = {3, 10, 29, 66};
    double x1[] = {1.5, 2.5, 3.5, 4.5};
    double f1[] = {5.375, 17.625, 44.875, 93.125};
    int n = 4;

    auto t  = divided_diff(x,  f,  n);
    auto t1 = divided_diff(x1, f1, n);

    cout << "Основные точки\n";
    print_table(x, t, n);
    cout << "Среднеузловые значения:\n";
    for (int i = 0; i + 1 < n; i++) {
        double xm = (x[i] + x[i+1]) / 2;
        cout << "  N(" << xm << ") = " << newton(xm, x, t, n) << "\n";
    }

    cout << "\nПромежуточные точки\n";
    print_table(x1, t1, n);
    cout << "Среднеузловые значения:\n";
    for (int i = 0; i + 1 < n; i++) {
        double xm = (x1[i] + x1[i+1]) / 2;
        cout << "  N(" << xm << ") = " << newton(xm, x1, t1, n) << "\n";
    }
}
