import numpy as np
from scipy.sparse import random

def print_menu():
    print("Меню:")
    print("1. Генерация матрицы P")
    print("2. Разложение Холесского")
    print("3. Решение системы линейных уравнений")
    print("4. Решение системы линейных уравнений с заданной матрицей P")
    print("5. Эксперимент 1: Количество арифметических операций")
    print("6. Эксперимент 2: Решение СЛАУ с заполненной матрицей P")
    print("7. Эксперимент 3: Решение СЛАУ с разреженной матрицей P")
    print("0. Выход")

def print_A_and_b(A, b):
    print("Матрица A:")
    print(A)
    print("Вектор b:")
    print(b)

def input_A_and_b():
    print("1) автоматически 2)вручную")
    choise = int(input())
    n = int(input("Введите размерность матрицы: "))
    A = np.zeros((n, n))
    b = np.zeros(n)
    if choise == 1:
        print("Введите элементы матрицы A:")
        for i in range(n):
            for j in range(n):
                A[i][j] = float(input(f"A[{i}][{j}]: "))

        print("Введите элементы вектора b:")
        for i in range(n):
            b[i] = float(input(f"b[{i}]: "))
    if choise == 2:
        A = np.random.rand(n, n)
        b = np.random.rand(n)
    else:
        print("ошибка ввода матрица заполнена автоматически")
        A = np.random.rand(n, n)
        b = np.random.rand(n)
    return A, b

def generate_matrix_P(n):
    return np.random.rand(n, n)

def cholesky_decomposition(A):
    n = len(A)
    L = np.zeros_like(A)

    for i in range(n):
        for j in range(i+1):
            if i == j:
                sum_term = np.sum(L[i][:j]**2)
                L[i][j] = np.sqrt(A[i][i] - sum_term)
            else:
                sum_term = np.sum(L[i][:j] * L[j][:j])
                L[i][j] = (A[i][j] - sum_term) / L[j][j]

    return L

def solve_linear_equations(A, b):
    n = len(A)
    x = np.zeros(n)

    for k in range(n):
        for j in range(k):
            b[k] -= A[k][j] * x[j]
        x[k] = b[k] / A[k][k]

    return x

def solve_linear_equations_with_P(P, A, b):
    try:
        n = len(A)
        A_modified = P @ A
        x = np.zeros(n)

        for k in range(n):
            for j in range(k):
                b[k] -= A_modified[k][j] * x[j]
            x[k] = b[k] / A_modified[k][k]

        return x
    except Exception as e:
        print("Ошибка", e)
    return x

def experiment_1(A,b):
    n = A.shape[0]
    print_A_and_b(A,b)

    # Решение системы линейных уравнений Ax = b
    x = np.linalg.solve(A, b)
    print("x - ")
    print(x)
    # Оценка количества операций
    num_operations = n ** 3 / 3 + n ** 2 / 2
    print("Количество арифметических операций:", num_operations)

def experiment_2(A,b):
    n = A.shape[0]
    print_A_and_b(A, b)
    P = np.random.rand(n, n)
    print("p - ")
    print(P)
    # Решение системы линейных уравнений Ax = b с матрицей P
    A_modified = P @ A
    x = np.linalg.solve(A_modified, b)

    print("Решение системы линейных уравнений с заполненной матрицей P:")
    print(x)

def experiment_3(A,b):
    n = A.shape[0]
    print_A_and_b(A,b)
    # Генерация разреженной матрицы P
    density = float(input("Введите плотность разреженной матрицы P (от 0 до 1): "))
    P = random(n, n, density=density).toarray()

    # Решение системы линейных уравнений Ax = b с разреженной матрицей P
    A_modified = P @ A
    x = np.linalg.solve(A_modified, b)

    print("Решение системы линейных уравнений с разреженной матрицей P:")
    print(x)

def main():
    A,b = input_A_and_b()
    while True:
        print_menu()
        choice = input("Выберите пункт из меню (0-7): ")

        if choice == '0':
            print("Программа завершена.")
            break
        elif choice == '1':
            n = int(input("Введите размерность матрицы P: "))
            P = generate_matrix_P(n)
            print("Матрица P:")
            print(P)
        elif choice == '2':
            # Запрос исходной матрицы A
            L = cholesky_decomposition(A)
            print("Разложение Холесского матрицы A:")
            print(L)
        elif choice == '3':
            # Запрос исходной матрицы A и вектора b
            x = solve_linear_equations(A, b)
            print("Решение системы линейных уравнений Ax = b:")
            print(x)
        elif choice == '4':
            # Запрос матрицы P, исходной матрицы A и вектора b
            n = int(input("Введите размерность матрицы P: "))
            P = generate_matrix_P(n)
            print("Матрица P:")
            print(P)
            x = solve_linear_equations_with_P(P, A, b)
            print("Решение системы линейных уравнений с матрицей P:")
            print(x)
        elif choice == '5':
            experiment_1()
        elif choice == '6':
            experiment_2()
        elif choice == '7':
            experiment_3(A,b)
        else:
            print("Некорректный ввод. Пожалуйста, выберите пункт из меню.")

if __name__ == "__main__":
    main()
