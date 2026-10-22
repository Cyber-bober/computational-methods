#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

int main() {
    const int n = 5;
    double v = 2.0;
    double eps = 1e-2;  // 10^(-2)

    // Матрица A: диагональ = v+3*i, остальные = диагональ * eps
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    vector<double> b(n);

    for (int i = 0; i < n; i++) {
        double diag = v + 3.0 * i;
        b[i] = diag;
        for (int j = 0; j < n; j++) {
            if (i == j)
                A[i][j] = diag;
            else
                A[i][j] = diag * eps;
        }
    }

    // Вывод исходной матрицы
    cout << "Исходная матрица A и вектор b:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            cout << setw(10) << fixed << setprecision(4) << A[i][j];
        cout << " | " << b[i] << endl;
    }
    cout << endl;

    // Прямой ход Гаусса с выводом промежуточных матриц
    for (int k = 0; k < n; k++) {
        // Выбор главного элемента
        int maxRow = k;
        for (int i = k + 1; i < n; i++)
            if (fabs(A[i][k]) > fabs(A[maxRow][k]))
                maxRow = i;
        swap(A[k], A[maxRow]);
        swap(b[k], b[maxRow]);

        // Исключение
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++)
                A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }

        cout << "После шага k=" << k + 1 << ":" << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                cout << setw(10) << fixed << setprecision(6) << A[i][j];
            cout << " | " << b[i] << endl;
        }
        cout << endl;
    }

    // Обратный ход
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }

    // Вывод решения
    cout << "Решение x:" << endl;
    for (int i = 0; i < n; i++)
        cout << "x[" << i + 1 << "] = " << fixed << setprecision(6) << x[i] << endl;

    // Проверка
    cout << "\nПроверка A*x (должно совпадать с b = [2,5,8,11,14]):" << endl;
    // Восстановим исходную матрицу для проверки
    for (int i = 0; i < n; i++) {
        double diag = v + 3.0 * i;
        double sum = 0;
        for (int j = 0; j < n; j++) {
            double aij = (i == j) ? diag : diag * eps;
            sum += aij * x[j];
        }
        cout << "Строка " << i + 1 << ": " << sum
             << " (должно быть " << diag << ")" << endl;
    }

    return 0;
}
