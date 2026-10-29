#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

// Структура результата метода Гаусса
struct GaussResult {
    vector<double> x;       // Решение системы
    double determinant;     // Определитель матрицы A
    bool is_singular;       // Флаг вырожденности
    int swaps;              // Количество перестановок строк
};

/**
 * Решение СЛАУ Ax=b методом Гаусса с выбором главного элемента
 * и вычислением определителя.
 * 
 * @param A Матрица коэффициентов (будет изменена в процессе)
 * @param b Вектор правой части (будет изменен)
 * @param eps Точность для проверки главного элемента (например, 1e-9)
 */
GaussResult solveGauss(vector<vector<double>> A, vector<double> b, double eps) {
    int n = A.size();
    GaussResult res;
    res.swaps = 0;
    res.is_singular = false;
    res.determinant = 1.0;

    // прямой ход
    for (int k = 0; k < n; ++k) {
        // 1. Выбор главного элемента в столбце k (строки k..n-1)
        int max_row = k;
        double max_val = fabs(A[k][k]);
        
        for (int i = k + 1; i < n; ++i) {
            if (fabs(A[i][k]) > max_val) {
                max_val = fabs(A[i][k]);
                max_row = i;
            }
        }

        // 2. Проверка условия |a_kk| > eps
        if (max_val < eps) {
            res.is_singular = true;
            res.determinant = 0.0;
            return res;
        }

        // 3. Перестановка строк, если необходимо
        if (max_row != k) {
            swap(A[k], A[max_row]);
            swap(b[k], b[max_row]);
            res.swaps++;
        }

        // 4. Исключение переменных
        for (int i = k + 1; i < n; ++i) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // вычисление определителя
    // Формула: Delta_A = (-1)^nu * a_11 * a_22^(1) * ... * a_nn^(n-1)
    for (int i = 0; i < n; ++i) {
        res.determinant *= A[i][i];
    }
    if (res.swaps % 2 != 0) {
        res.determinant = -res.determinant;
    }

    // обратный ход
    res.x.resize(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= A[i][j] * res.x[j];
        }
        res.x[i] = sum / A[i][i];
    }

    return res;
}

int main() {
    const int n = 5;
    double v = 2.0;
    double eps = 1e-9;

    // Генерация матрицы A и вектора b
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    vector<double> b(n);

    cout << "Формирование тестовой системы Ax=b:\n";
    for (int i = 0; i < n; ++i) {
        double diag = v + 3.0 * i;
        b[i] = diag;
        for (int j = 0; j < n; ++j) {
            if (i == j)
                A[i][j] = diag;
            else
                A[i][j] = diag * 1e-2;
        }
    }

    // Вывод исходной системы
    cout << fixed << setprecision(4);
    cout << "Матрица A:\n";
    for (const auto& row : A) {
        for (double val : row) cout << setw(10) << val;
        cout << "\n";
    }
    cout << "\nВектор b: ";
    for (double val : b) cout << val << " ";
    cout << "\n\n";

    // Решение
    GaussResult result = solveGauss(A, b, eps);

    if (result.is_singular) {
        cout << "Главный элемент меньше eps (" << eps << ").\n";
    } else {
        cout << "Определитель det(A) = " << setprecision(10) << result.determinant << "\n";
        cout << "Количество перестановок строк (nu) = " << result.swaps << "\n\n";
        
        cout << "Решение системы x:\n";
        for (int i = 0; i < n; ++i) {
            cout << "x[" << i+1 << "] = " << setprecision(8) << result.x[i] << "\n";
        }
        // Простая проверка подстановкой
        cout << "\nПроверка (невязка ||Ax - b||):\n";
        double error = 0.0;
        // Восстанавливаем исходную матрицу для проверки
        for (int i = 0; i < n; ++i) {
            double diag = v + 3.0 * i;
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                double aij = (i == j) ? diag : diag * 1e-2;
                sum += aij * result.x[j];
            }
            error += fabs(sum - b[i]);
        }
        cout << "Суммарная абсолютная невязка: " << scientific << error << "\n";
    }
    return 0;
}
