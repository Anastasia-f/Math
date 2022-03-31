import numpy as np
import sympy as sp
from scipy import linalg
from scipy.integrate import quad
import matplotlib.pyplot as plt

left = -1
right = 1


def print_a(n, M, V):
    print("\nСЛАУ в матричном виде для нахождения коэфициентов a1, ..., an")
    for i in range(n):
        for j in range(n):
            print("{:.5f}".format(M[i][j]), end=' ')
        print("    |     {:.5f}".format(V[i]))


def error(F, LF, left, right):
    def integrand1(x_):
        x = sp.symbols('x')
        try:
            F_ = F.subs([(x, x_)])
        except:
            F_ = F
        return (F_) ** 2

    def integrand2(x_):
        x = sp.symbols('x')
        try:
            LF_ = LF.subs([(x, x_)])
        except:
            LF_ = LF
        return (LF_) ** 2
    E = (quad(integrand1, left, right)[0] / (right - left)) / (quad(integrand2, left, right)[0] / (right - left))
    return np.sqrt(E)


def collocation_method(P, Q, F, n, left, right, fi):
    X = np.zeros(n)
    for i in range(n):
        X[i] = left + (right - left) / (2 * n) + i * (right - left) / n + (
                    i * (right - left) / n / 1000000)
    print("Узлы коллокации:", X)

    M = np.zeros((n, n))
    V = np.zeros(n)

    for i in range(n):
        for j in range(n):
            x = sp.symbols('x')
            diff2_fi = sp.diff(fi(x, j + 1), x, x)
            diff1_fi = sp.diff(fi(x, j + 1), x)
            try:
                Q_ = Q.subs([(x, X[i])])
            except:
                Q_ = Q
            try:
                P_ = P.subs([(x, X[i])])
            except:
                P_ = P
            M[i][j] = diff2_fi.subs([(x, X[i])]) + diff1_fi.subs([(x, X[i])]) * P_ + \
                      fi(X[i], j + 1) * Q_
            try:
                V[i] = F.subs([(x, X[i])])
            except:
                V[i] = F

    print_a(n, M, V)
    A = list(linalg.solve(M, V))
    print("\nРешение системы: вектор A = {0}".format(A))
    return A


def galerkin_method(P, Q, F, n, left, right, fi):

    M = np.zeros((n, n))
    V = np.zeros(n)

    for i in range(n):
        for j in range(n):
            x = sp.symbols('x')
            diff2_fi = sp.diff(fi(x, j+1), x, x)
            diff1_fi = sp.diff(fi(x, j+1), x)
            def integrand_M(x_):
                try:
                    Q_ = Q.subs([(x, x_)])
                except:
                    Q_ = Q
                try:
                    P_ = P.subs([(x, x_)])
                except:
                    P_ = P
                return (diff2_fi.subs([(x, x_)]) + P_*diff1_fi.subs([(x, x_)]) + Q_*fi(x_, j+1))\
                   * fi(x_, i+1)
            M[i][j] = quad(integrand_M, left, right)[0]
        def integrand_V(x_):
            try:
                F_ = F.subs([(x, x_)])
            except:
                F_ = F
            return F_ * fi(x_, i+1)
        V[i] = quad(integrand_V, left, right)[0]
    print_a(n, M, V)

    A = list(linalg.solve(M, V))
    print("\nРешение системы: вектор A = {0}".format(A))
    return A


def integral_LSM(P, Q, F, n, left, right, fi):
    M = np.zeros((n, n))
    V = np.zeros(n)

    for i in range(n):

        x = sp.symbols('x')
        diff2_fi_i = sp.diff(fi(x, i + 1), x, x)
        diff1_fi_i = sp.diff(fi(x, i + 1), x)

        for j in range(n):
            diff2_fi = sp.diff(fi(x, j + 1), x, x)
            diff1_fi = sp.diff(fi(x, j + 1), x)

            def integrand_M(x_):
                try:
                    Q_ = Q.subs([(x, x_)])
                except:
                    Q_ = Q
                try:
                    P_ = P.subs([(x, x_)])
                except:
                    P_ = P
                return (diff2_fi.subs([(x, x_)]) + P_ * diff1_fi.subs([(x, x_)]) + Q_ * fi(x_, j + 1)) \
                       * (diff2_fi_i.subs([(x, x_)]) + P_ * diff1_fi_i.subs([(x, x_)]) + Q_ * fi(x_, i + 1))

            M[i][j] = quad(integrand_M, left, right)[0]

        def integrand_V(x_):
            try:
                F_ = F.subs([(x, x_)])
            except:
                F_ = F
            try:
                P_ = P.subs([(x, x_)])
            except:
                P_ = P
            try:
                Q_ = Q.subs([(x, x_)])
            except:
                Q_ = Q
            return F_ * (diff2_fi_i.subs([(x, x_)]) + P_ * diff1_fi_i.subs([(x, x_)]) + \
                         Q_ * fi(x_, i + 1))

        V[i] = quad(integrand_V, left, right)[0]

    print_a(n, M, V)
    A = list(linalg.solve(M, V))
    print("\nРешение системы: вектор A = {0}".format(A))
    return A


def discrete_LSM(P, Q, F, n, N, left, right, fi):
    X = np.zeros(N)
    for i in range(N):
        X[i] = left + ((i + 1) - 0.5) / N * (right + 1)
    M = np.zeros((n, n))
    V = np.zeros(n)

    for i in range(n):
        x = sp.symbols('x')
        diff2_fi_i = sp.diff(fi(x, i + 1), x, x)
        diff1_fi_i = sp.diff(fi(x, i + 1), x)
        for j in range(n):
            diff2_fi = sp.diff(fi(x, j + 1), x, x)
            diff1_fi = sp.diff(fi(x, j + 1), x)
            sum_M = 0
            for k in range(N):
                try:
                    Q_ = Q.subs([(x, X[k])])
                except:
                    Q_ = Q
                try:
                    P_ = P.subs([(x, X[k])])
                except:
                    P_ = P
                sum_M += (diff2_fi.subs([(x, X[k])]) + diff1_fi.subs([(x, X[k])]) * P_ + \
                          fi(X[k], j + 1) * Q_) * \
                         (diff2_fi_i.subs([(x, X[k])]) + diff1_fi_i.subs([(x, X[k])]) * P_ + \
                          fi(X[k], i + 1) * Q_)
            M[i][j] = sum_M
        sum_V = 0
        for k in range(N):
            try:
                F_ = F.subs([(x, X[k])])
            except:
                F_ = F
            try:
                Q_ = Q.subs([(x, X[k])])
            except:
                Q_ = Q
            try:
                P_ = P.subs([(x, X[k])])
            except:
                P_ = P

            sum_V += F_ * (diff2_fi_i.subs([(x, X[k])]) + diff1_fi_i.subs([(x, X[k])]) * P_ + \
                           fi(X[k], i + 1) * Q_)
        V[i] = sum_V
    print_a(n, M, V)
    A = list(linalg.solve(M, V))
    print("\nРешение системы: вектор A = {0}".format(A))
    return A


def fi(x, i):
    if i == 0:
        return 0
    if i > 0:
        return x**(i-1) * (1 - x**2)
    raise ValueError('Индекс функций базисной системы не может быть отрицательным')


def print_res(n, A, fi, x):
    print("\nПриближенное решение ДУ при n = {:}: \ny = ".format(n), end='')
    for i in range(len(A)):
        if round(A[i], 5) == 0:
            continue
        elif A[i] > 0 and i != 0:
            print("+{:.5f}".format(A[i]), "*", fi(x, i + 1), end=' ')
        else:
            print("{:.5f}".format(A[i]), "*", fi(x, i + 1), end=' ')
    print()


def solution(A, x):
    y = 0
    for i in range(len(A)):
        y += A[i] * fi(x, i + 1)
    return y

def y_fun1(x, A1):
    return solution(A1, x)

def y_fun2(x, A2):
    return solution(A2, x)

def y_fun3(x, A3):
    return solution(A3, x)

def y_fun4(x, A4):
    return solution(A4, x)

def print_solution(A1, A2, A3, A4):
    xlist = np.linspace(left, right, 20)
    ylist1 = [y_fun1(x, A1) for x in xlist]
    ylist2 = [y_fun2(x, A2) for x in xlist]
    ylist3 = [y_fun3(x, A3) for x in xlist]
    ylist4 = [y_fun4(x, A4) for x in xlist]

    fig = plt.figure()

    plt.plot(xlist, ylist1, 'r--', label='Coll')
    plt.plot(xlist, ylist2, 'g--', label='Gal')
    plt.plot(xlist, ylist3, 'b--', label='Int')
    plt.plot(xlist, ylist4, 'black', linestyle='--', label='Dis')

    plt.legend(loc='best')
    plt.show()


def L(F, P, Q, x):
    diff2 = sp.diff(F, x, x)
    diff1 = sp.diff(F, x)
    res = diff2 + diff1 * P + F * Q
    return res


def my_F(A, fi, P, Q):
    x = sp.symbols('x')
    Y = 0
    for i in range(len(A)):
        Y += A[i] * fi(x, i+1)
    return L(Y, P, Q, x)


def main():
    k = 7
    a = np.sin(k)
    b = np.cos(k)
    P = 0
    x = sp.symbols('x')
    Q = (1 + b * x ** 2) / a
    F = -1 / a
    left = -1
    right = 1
    n = 5
    N = n * 2

    print("Условие:")
    print("sin({0}) * y\" + (1 + cos({0}) * x**2)*y = -1".format(k))

    print("\nБазисная система:")
    for i in range(n + 1):
        print('fi[{:}] = {:}'.format(i, fi(x, i)))

    print("\n\nМетод коллокаций:\n")
    A1 = collocation_method(P, Q, F, n, left, right, fi)
    print_res(n, A1, fi, x)

    print()

    print("\n\nМетод Галеркина:\n")
    A2 = galerkin_method(P, Q, F, n, left, right, fi)
    print_res(n, A2, fi, x)

    print()

    print("\n\nИнтегральный метод наименьших квадратов:\n")
    A3 = integral_LSM(P, Q, F, n, left, right, fi)
    print_res(n, A3, fi, x)

    print()

    print("\n\nДискретный метод наименьших квадратов:\n")
    print("N = ", N)
    A4 = discrete_LSM(P, Q, F, n, N, left, right, fi)
    print_res(n, A4, fi, x)

    print()
    ###################################################################################################################
    print("###########################################################################################################")
    print("\nСравнение методов\n")
    print("Метод коллокаций A =", A1)
    print("Метод Галеркина A =", A2)
    print("Интегральный метод наименьших квадратов A =", A3)
    print("Дискретный метод наименьших квадратов A =", A4)

    A = [A1, A2, A3, A4]
    X = [1 / 4, 1 / 2, 3 / 4]
    length_X = len(X)

    res = np.zeros((length_X, 4))
    for j in range(4):
        for k in range(length_X):
            for i in range(len(A[j])):
                if A[j][i] == 0:
                    continue
                res[k][j] += A[j][i] * fi(X[k], i + 1)

    print("\nЗначения в точках: ", X, "\n  Кол.     Гал.     ИМНК     ДМНК")
    np.set_printoptions(suppress=True, precision=5, floatmode="fixed")
    print(res)

    razn_F = F - my_F(A1, fi, P, Q)
    E1 = error(razn_F, L(razn_F, P, Q, x), left, right)
    print("\nПогрешность метода коллокаций: {:}".format(E1))

    razn_F = F - my_F(A2, fi, P, Q)
    E2 = error(razn_F, L(razn_F, P, Q, x), left, right)
    print("Погрешность метода Галеркина:  {:}".format(E2))

    razn_F = F - my_F(A3, fi, P, Q)
    E3 = error(razn_F, L(razn_F, P, Q, x), left, right)
    print("Погрешность интегрального метод наименьших квадратов: {:}".format(E3))

    razn_F = F - my_F(A4, fi, P, Q)
    E4 = error(razn_F, L(razn_F, P, Q, x), left, right)
    print("Погрешность дискретного метод наименьших квадратов:   {:}".format(E4))

    print_solution(A1, A2, A3, A4)


if __name__=="__main__":
    main()

