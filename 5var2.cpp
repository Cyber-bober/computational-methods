#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Структура для хранения коэффициентов сплайна
struct SplineCoeff {
    double a, b, c, d;
};

int main() {
    // === Входные данные (изменяемые)
    vector<double> x = {1, 2, 3, 4};       // Узлы
    vector<double> f = {3, 10, 29, 66};    // Значения функции
    int n = x.size() - 1;                  // Количество сплайнов (интервалов)
    
    if (n < 1) return 1;

    int N = 4 * n; // Общее число неизвестных (12 для n=3)
    vector<vector<double>> A(N, vector<double>(N + 1, 0.0)); // Расширенная матрица (A|B)

    // Заполнение матрицы
    for (int i = 0; i < n; ++i) {
        double h = x[i+1] - x[i];
        int base = 4 * i; // Смещение для текущей группы (a_i, b_i, c_i, d_i)

        // 1. Условия интерполяции (левый конец)
        A[2*i][base] = 1.0;             // a_i = f_i
        A[2*i][N] = f[i];

        // 2. Условия интерполяции (правый конец)
        // a_i + b_i*h + c_i*h^2 + d_i*h^3 = f_{i+1}
        A[2*i+1][base]     = 1.0;
        A[2*i+1][base + 1] = h;
        A[2*i+1][base + 2] = h * h;
        A[2*i+1][base + 3] = h * h * h;
        A[2*i+1][N]        = f[i+1];
    }

    // 3. Условия гладкости (стыковка производных)
    for (int i = 0; i < n - 1; ++i) {
        double h_curr = x[i+1] - x[i];
        double h_next = x[i+2] - x[i+1];
        
        int base_curr = 4 * i;
        int base_next = 4 * (i + 1);
        int row_smooth = 2 * n + i; // Строки после условий интерполяции

        // S'_i(x_{i+1}) = S'_{i+1}(x_{i+1})
        // b_i + 2*c_i*h_curr + 3*d_i*h_curr^2 = b_{i+1}
        A[row_smooth][base_curr + 1] = 1.0;
        A[row_smooth][base_curr + 2] = 2.0 * h_curr;
        A[row_smooth][base_curr + 3] = 3.0 * h_curr * h_curr;
        A[row_smooth][base_next + 1] = -1.0;

        // S''_i(x_{i+1}) = S''_{i+1}(x_{i+1})
        // 2*c_i + 6*d_i*h_curr = 2*c_{i+1}
        int row_smooth2 = 2 * n + (n - 1) + i;
        A[row_smooth2][base_curr + 2] = 2.0;
        A[row_smooth2][base_curr + 3] = 6.0 * h_curr;
        A[row_smooth2][base_next + 2] = -2.0;
    }

    // 4. Граничные условия (естественный сплайн: S''=0 на краях)
    int row_bc_left = 4 * n - 2;
    int row_bc_right = 4 * n - 1;
    
    // Левая граница: 2*c_0 = 0
    A[row_bc_left][2] = 2.0; 
    
    // Правая граница: 2*c_n + 6*d_n*h_n = 0
    double h_last = x[n] - x[n-1];
    A[row_bc_right][4*(n-1) + 2] = 2.0;
    A[row_bc_right][4*(n-1) + 3] = 6.0 * h_last;

    // Метод Гаусса
    for (int col = 0; col < N; ++col) {
        // Поиск главного элемента
        int max_row = col;
        for (int k = col + 1; k < N; ++k)
            if (fabs(A[k][col]) > fabs(A[max_row][col])) max_row = k;
        swap(A[col], A[max_row]);

        // Приведение к треугольному виду
        for (int k = col + 1; k < N; ++k) {
            double factor = A[k][col] / A[col][col];
            for (int j = col; j <= N; ++j)
                A[k][j] -= factor * A[col][j];
        }
    }

    // Обратный ход
    vector<double> X(N);
    for (int i = N - 1; i >= 0; --i) {
        double sum = A[i][N];
        for (int j = i + 1; j < N; ++j)
            sum -= A[i][j] * X[j];
        X[i] = sum / A[i][i];
    }

    cout << fixed << setprecision(4);
    cout << "Коэффициенты сплайнов:\n";
    for (int i = 0; i < n; ++i) {
        cout << "S" << i+1 << "(x) на [" << x[i] << ", " << x[i+1] << "]:\n";
        cout << "  a = " << X[4*i]     << "\n";
        cout << "  b = " << X[4*i + 1] << "\n";
        cout << "  c = " << X[4*i + 2] << "\n";
        cout << "  d = " << X[4*i + 3] << "\n\n";
    }

    return 0;
}
