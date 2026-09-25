#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int main() {
    double x[] = {1, 2, 3, 4};
    double f[] = {3, 10, 29, 66};
    int n = 4;

    cout << "=== Исходная таблица f (2x4) ===" << endl;
    cout << "x: ";
    for (int i = 0; i < n; i++) cout << setw(8) << x[i];
    cout << endl;
    cout << "f: ";
    for (int i = 0; i < n; i++) cout << setw(8) << f[i];
    cout << endl << endl;

    // 2) Построение РМК СЛАУ (расширенная матрица)
    // P3(x) = a3*x^3 + a2*x^2 + a1*x + a0
    // Система: для каждого xi: a3*xi^3 + a2*xi^2 + a1*xi + a0 = fi
    double A[4][5]; // 4x5 расширенная матрица
    for (int i = 0; i < n; i++) {
        A[i][0] = pow(x[i], 3);
        A[i][1] = pow(x[i], 2);
        A[i][2] = x[i];
        A[i][3] = 1.0;
        A[i][4] = f[i];
    }

    cout << "=== РМК СЛАУ (4x5) ===" << endl;
    cout << "a3      a2      a1      a0      |  f" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 5; j++) {
            cout << setw(8) << A[i][j];
            if (j == 3) cout << "  |";
        }
        cout << endl;
    }
    cout << endl;

    // 3) Решение СЛАУ методом Гаусса
    // Прямой ход
    for (int k = 0; k < n; k++) {
        // Поиск главного элемента
        int maxRow = k;
        for (int i = k + 1; i < n; i++)
            if (fabs(A[i][k]) > fabs(A[maxRow][k])) maxRow = i;
        swap(A[k], A[maxRow]);

        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j <= n; j++)
                A[i][j] -= factor * A[k][j];
        }
    }

    // Обратный ход
    double a[4]; // a[0]=a3, a[1]=a2, a[2]=a1, a[3]=a0
    for (int i = n - 1; i >= 0; i--) {
        a[i] = A[i][n];
        for (int j = i + 1; j < n; j++)
            a[i] -= A[i][j] * a[j];
        a[i] /= A[i][i];
    }

    cout << "=== Коэффициенты ===" << endl;
    cout << "a3 = " << a[0] << endl;
    cout << "a2 = " << a[1] << endl;
    cout << "a1 = " << a[2] << endl;
    cout << "a0 = " << a[3] << endl << endl;

    // 4) Таблица P3(x) 2x7
    // Точки: x0, (x0+x1)/2, x1, (x1+x2)/2, x2, (x2+x3)/2, x3
    double xp[7];
    xp[0] = x[0];
    xp[1] = (x[0] + x[1]) / 2.0;
    xp[2] = x[1];
    xp[3] = (x[1] + x[2]) / 2.0;
    xp[4] = x[2];
    xp[5] = (x[2] + x[3]) / 2.0;
    xp[6] = x[3];

    cout << "=== Таблица P3(x) (2x7) ===" << endl;
    cout << "x:  ";
    for (int i = 0; i < 7; i++) cout << setw(8) << xp[i];
    cout << endl;
    cout << "P3: ";
    for (int i = 0; i < 7; i++) {
        double val = a[0]*pow(xp[i],3) + a[1]*pow(xp[i],2) + a[2]*xp[i] + a[3];
        cout << setw(8) << val;
    }
    cout << endl;

    return 0;
}
