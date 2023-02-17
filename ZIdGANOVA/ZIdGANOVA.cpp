#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

vector<double> gaussian_elimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        // Поиск максимального элемента в столбце i
        double max_el = abs(A[i][i]);
        int max_row = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > max_el) {
                max_el = abs(A[j][i]);
                max_row = j;
            }
        }

        // Перестановка строки с максимальным элементом на i-ую строку
        if (max_row != i) {
            swap(A[i], A[max_row]);
            swap(b[i], b[max_row]);
        }

        // Приведение матрицы к треугольному виду
        for (int j = i + 1; j < n; j++) {
            double f = A[j][i] / A[i][i];
            b[j] -= f * b[i];
            for (int k = i; k < n; k++) {
                A[j][k] -= f * A[i][k];
            }
        }
    }

    // Обратный ход
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}

void print_matrix_and_b(vector<vector<double>> A, vector<double> b) {
    cout << "Получившаяся матрица" << endl;
    string tabulation = "";
    for (int i = 0; i < A.size(); i++) {
        tabulation += "\t";
    }
    cout << "A" << tabulation << "| " << "b" << endl;
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A.size(); j++) {
            cout << A[i][j] << "\t";
        }
        cout << "| " << b[i] << endl;
    }
}

void print_matrix(vector<vector<double>> A) {
    int n = A.size();
    cout << "Matrix A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << A[i][j] << " ";
        }
        cout << endl;
    }
}

double determinant(vector<vector<double>> A) {
    int n = A.size();
    double det = 1.0;

    for (int i = 0; i < n; i++) {
        int max_row = i;

        // Ищем строку с максимальным элементом в текущем столбце
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[max_row][i])) {
                max_row = j;
            }
        }

        // Если нашли строку с максимальным элементом, меняем текущую строку с нею местами
        if (max_row != i) {
            swap(A[i], A[max_row]);
            det *= -1;
        }

        // Приводим матрицу к треугольному виду
        for (int j = i + 1; j < n; j++) {
            double c = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= c * A[i][k];
            }
        }

        // Умножаем определитель на элемент главной диагонали
        det *= A[i][i];
    }

    return det;
}

unsigned menu() {
    cout << "Лабораторный проект №1. Вариант 4. Разложение на основе гауссова исключения по столбцам" << endl;
    cout << "Авторы: Нуштаев Никита Петрович, Яин Максим Викторович, ПРИ-О-21/1" << endl;
    cout << "1) Ввод данных" << endl;
    cout << "2) Факторизация" << endl;
    cout << "3) Решение СЛАУ" << endl;
    cout << "4) Определитель" << endl;
    cout << "5) Обращение через AX=E(элементарные преобразования)" << endl;
    cout << "6) Эксперимент 1" << endl;
    cout << "7) Эксперимент 2" << endl;
    cout << "8) Эксперимент 3" << endl;
    cout << "9) Выход" << endl;

    unsigned menuUk;
    cin >> menuUk;
    if (menuUk < 11 && menuUk > 0) {
        return menuUk;
    }
    else
    {
        return 0;
    }
}

void matrix_inverse(vector<vector<double>> A, vector<vector<double>>& A_inv) {
    int n = A.size();
    vector<vector<double>> E(n, vector<double>(n, 0));
    A_inv = vector<vector<double>>(n, vector<double>(n, 0));

    // Создаем единичную матрицу
    for (int i = 0; i < n; i++) {
        E[i][i] = 1;
    }

    // Прямой ход метода Гаусса
    for (int k = 0; k < n; k++) {
        int imax = k;
        double max_val = abs(A[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (abs(A[i][k]) > max_val) {
                max_val = abs(A[i][k]);
                imax = i;
            }
        }

        if (A[imax][k] == 0) {
            cout << "Ошибка матрица сингулярна" << endl;
            return;
        }

        if (imax != k) {
            swap(A[k], A[imax]);
            swap(E[k], E[imax]);
        }

        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; j++) {
                A[i][j] = A[i][j] - factor * A[k][j];
            }
            A[i][k] = 0;
            for (int j = 0; j < n; j++) {
                E[i][j] = E[i][j] - factor * E[k][j];
            }
        }
    }

    // Обратный ход метода Гаусса
    for (int k = n - 1; k >= 0; k--) {
        for (int j = 0; j < n; j++) {
            A_inv[k][j] = E[k][j] / A[k][k];
        }
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < n; j++) {
                A_inv[i][j] = A_inv[i][j] - A[k][j] * A_inv[k][j] / A[k][k];
            }
        }
    }
}

void faktorization_matrix() {}

void ExperimentOne() {
    system("cls");
    cout << "Эксперимент 1" << endl;
    system("pause");
}

void ExperimentTwo() {
    system("cls");
    cout << "Эксперимент 2" << endl;
    system("pause");
}

void ExperimentThree() {
    system("cls");
    cout << "Эксперимент 3" << endl;
    system("pause");
}

int main() {
    vector<double> x;
    vector<vector<double>> A_inv;
    double det;
    setlocale(LC_ALL, "ru");
    while (true) {
        unsigned sch;
        system("cls");
        sch = menu();
        if (sch == 0) {
            cout << "Вы ввели что-то не то" << endl;
        }
        else {
            if (sch == 1) {
                system("cls");
                int n;
                cout << "Введите количество уравнений: ";
                cin >> n;
                vector<vector<double>> A(n, vector<double>(n));
                vector<double> b(n);
                cout << "Введите коэффициенты уравнений: " << endl;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        cin >> A[i][j];
                    }
                    cin >> b[i];
                }
                print_matrix_and_b(A, b);
                x = gaussian_elimination(A, b);
                det = determinant(A);
                matrix_inverse(A, A_inv);
                system("pause");
            }
            if (sch == 2) {
                faktorization_matrix();
            }
            if (sch == 3) {
                system("cls");
                cout << "решение: ";
                for (int i = 0; i < x.size(); i++) {
                    cout << x[i] << " ";
                }
                cout << endl;
                system("pause");
            }
            if (sch == 4) {
                system("cls");
                cout << "Определитель матрицы - " << det << endl;
                system("pause");
            }
            if (sch == 5) {
                system("cls"); 
                print_matrix(A_inv);
                system("pause");
            }
            if (sch == 6) {
                ExperimentOne();
            }
            if (sch == 7) {
                ExperimentTwo();
            }
            if (sch == 8) {
                ExperimentThree();
            }
            if (sch == 9) {
                return 0;
            }
        }
    }
    return 0;
}

