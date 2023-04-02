import random

import numpy as np


def solve_gauss_with_pivoting(A, b):
    P, L, U = lu_factorization_with_pivoting(A)
    b = P @ b
    n = len(L)
    y = np.zeros(n)
    for i in range(n):
        y[i] = b[i] - np.dot(L[i, :i], y[:i])
    x = np.zeros(n)
    for i in reversed(range(n)):
        x[i] = (y[i] - np.dot(U[i, i + 1:], x[i + 1:])) / U[i, i]
    return x


def lu_factorization_with_pivoting(A):
    n = len(A)
    P = np.eye(n)
    L = np.zeros((n, n))
    U = np.copy(A)
    for j in range(n):
        row_max = j + np.argmax(np.abs(U[j:, j]))
        if j != row_max:
            U[[j, row_max], j:] = U[[row_max, j], j:]
            P[[j, row_max], :] = P[[row_max, j], :]
            if j >= 1:
                L[[j, row_max], :j] = L[[row_max, j], :j]
        pivot = U[j, j]
        L[j, j] = 1
        if pivot != 0:
            U[j, j:] /= pivot
        for i in range(j + 1, n):
            L[i, j] = U[i, j]
            U[i, j:] -= L[i, j] * U[j, j:]
    return P, L, U

def determinant(A):
    P, L, U = lu_factorization_with_pivoting(A)
    n = len(A)
    det = np.prod(np.diag(U)) * np.prod(np.diag(P))
    return det


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
    """
    Inverts a square matrix A by solving the linear system AX = E for X, where E is the identity matrix.
    Returns the inverse matrix A_inv.
    """
    n = len(A)
    E = np.eye(n)
    X = np.linalg.solve(A, E)
    return X

def menu():
    print("Лабораторный проект №1. Вариант 4. Разложение на основе гауссова исключения по столбцам")
    print("Авторы: Нуштаев Никита Петрович, Яин Максим Викторович, ПРИ-О-21/1")
    print("1) Ввод данных")
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
    sch = input()
    return sch


def InputData(A, b):
    print("введите кол-ва уравнений - ")
    lenn = int(input())
    print("как заполним? 1) автоматически 2) вручную")
    sch = int(input())
    if sch == 1:
        sch = int(input())
        for i in range(lenn):
            for j in range(lenn):
                A[i][j] = random.random()
            b[i] = random.random()
    if sch == 2:
        sch = int(input())
        for i in range(lenn):
            for j in range(lenn):
                A[i][j] = float(input())
            b[i] = float(input())

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

def main_func():
    while 1:
        sch = menu()
        A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
        b = np.array([1, 2, 3])
        if sch == 1:
            InputData(A, b)
            print("получившаяся матрица - ")
            print_matrix_with_b(A, b)
        if sch == 2:
            print("матрица перестановок - ")
            print(lu_factorization_with_pivoting(A))
        if sch == 3:
            print("корни уравнения - ")
            print(solve_gauss_with_pivoting(A, b))
        if sch == 4:
            print("определитель митрицы - ")
            print(determinant(A))
        if sch == 5:
            print(" Обращение через элементарные преобразования - ")
            print(invert_matrix_elem(A))
        if sch == 6:
            print(" Обращение через AX=E - ")
            print(invert_matrix_AXE(A))
main_func()
