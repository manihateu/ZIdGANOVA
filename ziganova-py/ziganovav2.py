import random
import time
import numpy as np

def lu_factorization(A):
    n = len(A)
    LU = np.copy(A)
    P = np.eye(n)
    for j in range(n):
        row_max = j + np.argmax(np.abs(LU[j:, j]))
        if j != row_max:
            LU[[j, row_max], j:] = LU[[row_max, j], j:]
            P[[j, row_max], :] = P[[row_max, j], :]
        pivot = float(LU[j, j])
        LU = LU.astype(type(pivot))
        for i in range(j+1, n):
            LU[i,j:] -= LU[i,j]*LU[j,j:]/ pivot
            LU[i,j] /= pivot
    return P, LU

def determinant_LU(A):
    n = len(A)
    L = np.eye(n)
    U = np.copy(A)
    LU = lu_factorization(A)
    L = L.astype(np.float64)
    U = U.astype(np.float64)
    sign = 1

    for j in range(n):
        row_max = j + np.argmax(np.abs(U[j:, j]))
        if j != row_max:
            U[[j, row_max], j:] = U[[row_max, j], j:]
            L[[j, row_max], :j] = L[[row_max, j], :j]
            sign = -sign
        if U[j, j] == 0:
            return 0
        for i in range(j + 1, n):
            mult = U[i, j] / U[j, j]
            L[i, j] = mult
            U[i, j:] -= mult * U[j, j:]

    det = sign * np.prod(np.diag(U))

    return det

def solve_system_lu(A, b):
    n = len(A)
    A = A.astype(np.float64)
    b = b.astype(np.float64)
    for k in range(n):
        # Выбор главного элемента
        max_row = k
        for i in range(k + 1, n):
            if abs(A[i, k]) > abs(A[max_row, k]):
                max_row = i
        if max_row != k:
            A[[k, max_row], k:] = A[[max_row, k], k:]
            b[k], b[max_row] = b[max_row], b[k]

        # Нормирование первого столбца
        for i in range(k + 1, n):
            c = A[i, k] / A[k, k]
            A[i, k:] -= c * A[k, k:]
            b[i] -= c * b[k]

    # Обратный ход метода Гаусса
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]

    return x

def invert_matrix_elem(A):
    n = len(A)
    A_aug = np.hstack((A, np.eye(n)))  # Augmented matrix [A | I]

    # Perform row operations to transform A_aug to [I | B]
    for j in range(n):
        # Divide row j by A_jj to make A_jj equal to 1
        pivot = A_aug[j, j]
        if pivot != 0:
            A_aug[j, :] /= pivot

        # Subtract A_ij times row j from all other rows to make all A_ij below the diagonal equal to 0
        for i in range(n):
            if i != j:
                factor = A_aug[i, j]
                A_aug[i, :] -= factor * A_aug[j, :]

    # Extract the inverse matrix from the augmented matrix
    A_inv = A_aug[:, n:]
    return A_inv

def invert_matrix_AXE(A):
    n = len(A)
    E = np.eye(n)
    X = np.linalg.solve(A, E)
    return X

def menu():
    print("Лабораторный проект №1. Вариант 4. Разложение на основе гауссова исключения по столбцам")
    print("Авторы: Нуштаев Никита Петрович, Яин Максим Викторович, ПРИ-О-21/1")
    print("")
    print("2) Факторизация")
    print("3) Решение СЛАУ")
    print("4) Определитель")
    print("5) Обращение через элементарные преобразования")
    print("6) Обращение через AX=E")
    print("7) Эксперимент 1")
    print("8) Эксперимент 2")
    print("9) Эксперимент 3")
    print("10) Выход")
    print("Введите выбор - ")
    sch = int(input())
    return sch


def solve_system(A, b, n, auto_fill):
    if auto_fill:
        # Автоматическое заполнение
        if n is None:
            n = np.random.randint(2, 6)
        A = np.random.randint(-9, 10, size=(n, n))
        while np.linalg.det(A) == 0:
            A = np.random.randint(-9, 10, size=(n, n))
        b = np.random.randint(-9, 10, size=n)
    else:
        # Ручное заполнение
        if A is None:
            print("Введите матрицу A")
            A = np.array([list(map(float, input().split())) for _ in range(n)])
        if b is None:
            print("Введите вектор b")
            b = np.array(list(map(float, input().split())))
    return A, b

def print_matrix_with_b(A, b):
    n = len(A)
    for i in range(n):
        eqn = ""
        for j in range(n):
            eqn += f"{A[i, j]:7.2f} x{j + 1} "
            if j < n - 1:
                eqn += "+ "
        eqn += f"= {b[i]:7.2f}"
        print(eqn)


def pretty_print(arr):
    for row in arr:
        row_str = " ".join([str(elem) for elem in row])
        print(row_str)

#first exprtiment
def generate_system(n):
    """
    Generates a random linear system Ax = b with n equations.
    Returns A and b.
    """
    A = np.random.rand(n, n)
    b = np.random.rand(n)
    return A, b

def run_experiment():
    """
    Runs the experiment and prints the results table.
    """
    print("  n     Время    Погрешность      Теор. числ.    Реал. число")
    print("------------------------------------------------------------")
    for n in range(5, 101, 5):
        A, b = generate_system(n)
        
        start_time = time.time()
        x = solve_system_lu(A, b)
        end_time = time.time()
        time_elapsed = end_time - start_time
        
        theoretical_ops = n**3 / 3
        actual_ops = solve_system_lu.num_operations
        
        error = np.linalg.norm(A @ x - b)
        
        print(f"{n:4}  {time_elapsed:7.5f}  {error:9.2e}  {theoretical_ops:15.2e}  {actual_ops:8}")
    
        solve_system_lu.num_operations = 0
    print()
    print("решение СЛАУ")
    for i in range(len(x)):
            print(x[i], " x[", i, "]")    

print("введите уравнения - ")
print("Кол-во уравнений - ")
lenn = int(input())
A = np.zeros((lenn, lenn), dtype=np.float64)
b = np.zeros(lenn, dtype=np.float64)
print("1) Автоматически 2) Вручную")
auto = int(input())
A, b = solve_system(A, b, lenn, auto == 1)
print("получившаяся матрица - ")
print_matrix_with_b(A, b)
while 1:
    sch = menu()
    if sch == 2:
        print("матрица перестановок - ")
        P, LU = lu_factorization(A)
        print("p - ")
        pretty_print(P)
        print("LU-")
        pretty_print(LU)
    if sch == 3:
        print("корни уравнения - ")
        x = solve_system_lu(A, b)
        for i in range(len(x)):
            print(x[i], " x[", i, "]")
    if sch == 4:
        print("определитель матрицы - ")
        print(determinant_LU(A))
        print(np.linalg.det(A))
    if sch == 5:
        print(" Обращение через элементарные преобразования - ")
        A_inv = invert_matrix_elem(A)
        pretty_print(A_inv)
    if sch == 6:
        print(" Обращение через AX=E - ")
        pretty_print(invert_matrix_AXE(A))
    if sch == 7:
        run_experiment()

