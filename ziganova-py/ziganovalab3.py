import numpy as np

def display_menu():
    print("Меню:")
    print("1. Ввод матрицы A и вектора b")
    print("2. Факторизация матрицы")
    print("3. Решение СЛАУ")
    print("4. Вычисление определителя матрицы")
    print("5. Обращение матрицы")
    print("6. Эксперимент1")
    print("7. Эксперимент2")
    print("8. Эксперимент3")
    print("9. Эксперимент4")
    print("0. Выход")

def enter_matrix():
    rows = int(input("Введите количество строк матрицы: "))
    cols = int(input("Введите количество столбцов матрицы: "))
    matrix = np.zeros((rows, cols))

    print("Введите элементы матрицы:")
    for i in range(rows):
        for j in range(cols):
            matrix[i, j] = float(input(f"A[{i+1}, {j+1}]: "))

    vector = np.zeros(rows)
    print("Введите элементы вектора b:")
    for i in range(rows):
        vector[i] = float(input(f"b[{i+1}]: "))

    return matrix, vector

def matrix_factorization(matrix):
    n = matrix.shape[0]
    L = np.eye(n)  # Единичная матрица L
    U = np.copy(matrix)  # Копия исходной матрицы

    for k in range(n - 1):
        for i in range(k + 1, n):
            if U[k, k] == 0:
                raise ValueError("Матрица не может быть факторизована с использованием LU-разложения.")
            L[i, k] = U[i, k] / U[k, k]
            for j in range(k, n):
                U[i, j] -= L[i, k] * U[k, j]

    return L, U

def solve_linear_equations(matrix, vector):
    L, U = matrix_factorization(matrix)

    # Решение системы Ly = b
    n = matrix.shape[0]
    y = np.zeros(n)
    y[0] = vector[0]
    for i in range(1, n):
        sum_val = np.dot(L[i, :i], y[:i])
        y[i] = vector[i] - sum_val

    # Решение системы Ux = y
    x = np.zeros(n)
    x[n - 1] = y[n - 1] / U[n - 1, n - 1]
    for i in range(n - 2, -1, -1):
        sum_val = np.dot(U[i, i + 1:], x[i + 1:])
        x[i] = (y[i] - sum_val) / U[i, i]

    return x

def compute_determinant(matrix):
    L, U = matrix_factorization(matrix)

    # Определитель равен произведению диагональных элементов матрицы U
    determinant = np.prod(np.diagonal(U))

    return determinant

def rotate(matrix, i, j, c, s):
    # Применение вращения Гивенса к матрице matrix
    n = matrix.shape[1]
    for k in range(n):
        temp = c * matrix[i, k] + s * matrix[j, k]
        matrix[j, k] = -s * matrix[i, k] + c * matrix[j, k]
        matrix[i, k] = temp

def invert_matrix(matrix):
    n = matrix.shape[0]
    identity = np.eye(n)
    augmented_matrix = np.hstack((matrix, identity))
    for i in range(n):
        for j in range(i + 1, n):
            if augmented_matrix[j, i] != 0:
                r = np.sqrt(augmented_matrix[i, i]**2 + augmented_matrix[j, i]**2)
                c = augmented_matrix[i, i] / r
                s = augmented_matrix[j, i] / r
                rotate(augmented_matrix, i, j, c, s)

    for i in range(n - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            if augmented_matrix[j, i] != 0:
                r = np.sqrt(augmented_matrix[i, i]**2 + augmented_matrix[j, i]**2)
                c = augmented_matrix[i, i] / r
                s = augmented_matrix[j, i] / r
                rotate(augmented_matrix, i, j, c, s)

    inverted_matrix = augmented_matrix[:, n:]

    return inverted_matrix

def experiment1():
    # Реализуйте эксперимент 1
    pass

def experiment2():
    # Реализуйте эксперимент 2
    pass

def experiment3():
    # Реализуйте эксперимент 3
    pass

def experiment4():
    # Реализуйте эксперимент 4
    pass

# Основной код программы
matrix = None
vector = None

while True:
    display_menu()
    choice = input("Выберите пункт меню: ")

    if choice == '1':
        matrix, vector = enter_matrix()
    elif choice == '2':
        if matrix is not None:
            factorized_matrix = matrix_factorization(matrix)
            print("Факторизованная матрица:")
            print(factorized_matrix)
        else:
            print("Ошибка: матрица не была введена.")
    elif choice == '3':
        if matrix is not None and vector is not None:
            solution_vector = solve_linear_equations(matrix, vector)
            print("Вектор-решение:")
            print(solution_vector)
        else:
            print("Ошибка: матрица или вектор не были введены.")
    elif choice == '4':
        if matrix is not None:
            determinant = compute_determinant(matrix)
            print("Определитель матрицы:", determinant)
        else:
            print("Ошибка: матрица не была введена.")
    elif choice == '5':
        if matrix is not None:
            inverted_matrix = invert_matrix(matrix)
            print("Обратная матрица:")
            print(inverted_matrix)
        else:
            print("Ошибка: матрица не была введена.")
    elif choice == '6':
        experiment1()
    elif choice == '7':
        experiment2()
    elif choice == '8':
        experiment3()
    elif choice == '9':
        experiment4()
    elif choice == '0':
        print("Программа завершена.")
        break
    else:
        print("Ошибка: некорректный выбор.")
