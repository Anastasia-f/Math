import numpy as np 
import sympy as sp 
import matplotlib.pyplot as plt
from scipy import linalg


def plot_res(n, X, y, label): 
    '''печать n графков по известным массивам и названиям'''
    fig = plt.figure()
    for i in range(n):
        plt.plot(X, y[i], label=label[i])
    plt.legend(loc='best')
    plt.show()


def findSolution(SET, f, K, C, A, B, a, b, needPrint):
    '''Аналитическое решение 
       SET - номер набора, f, K, C, A, B - функции из условия, a, b - границы по х,
       bool needPrint: 0 - ничего не печатать, 1 - печать основного, 2 - подробная печать'''
    x = sp.symbols('x')
    c1 = sp.symbols('c1')
    c2 = sp.symbols('c2')

    u = -sp.integrate( (  (sp.integrate(f(x), x) - c1) / K(x, SET)  ), x ) + c2
    if (needPrint == 2):
        print('Проинтегрируем исходное уравнение дважды (С = {:}):\nu(x) = {:}'.format(SET, u))
    elif (needPrint == 1):
        print('Набор параметров №{:}'.format(SET))

    u1 = u.subs(x, a)
    u2 = u.subs(x, b)
    if (needPrint == 2):
        print('Подставим граничные условия в полученное уравнение:')
        print('u({:}) = {:} = {:}'.format(a, u1, A(SET)))
        print('u({:}) = {:} = {:}'.format(b, u2, B(SET)))

    f1 = -float(u1.args[0]) + A(SET)
    c21 = float(u1.args[1].subs(c2,1))
    try:
        c11 = float(u1.args[2].subs(c1,1))
    except:
        c11 = 0

    f2 = -float(u2.args[0]) + B(SET)
    c22 = float(u2.args[1].subs(c2,1))
    c12 = float(u2.args[2].subs(c1,1))

    M = np.array([[c11, c21], [c12, c22]]) 
    v = np.array([f1, f2]) 
    #print(M)
    #print(v)
    c = list(linalg.solve(M, v))
    if (needPrint in [1, 2]):
        print('Решение системы линейных уравнений: \nс1 = {:},    c2 = {:}'.format(c[0], c[1]))

    u = u.subs(c1, c[0])
    u = u.subs(c2, c[1])
    if (needPrint in [1, 2]):
        print('Аналитическое решение:\n u(x) = ', u, '\n')
    return u


def task1():
    # условие -d/dx(K(x)*du/dx) = f   K(x) = K(x,C)

    # наборы данных (7 шт.)
    def C(set):
        if set in [1, 4, 5, 6, 7]:
            return 1
        if set == 2:
            return 2
        if set == 3:
            return 0.1

    def K(x, set):
        if set in [1, 2, 3, 5, 6, 7]:
            return C(set) * k(x)
        if set == 4:
            return 1 / k(x)

    def A(set):
        if set in [1, 2, 3, 4, 6]:
            return A_
        if set in [5, 7]:
            return -A_

    def B(set):
        if set in [1, 2, 3, 4, 5]:
            return B_
        if set in [6, 7]:
            return -B_

    # # Вариант 3
    # a = 0.5
    # b = 1.5
    # A_ = 2
    # B_ = 6
    #
    # def k(x):
    #     return x ** (-2)
    #
    # def f(x):
    #     return -2 * x ** 2 - 2 * x

    # Вариант 7
    a = 1.0
    b = 2
    A_ = 3
    B_ = 3

    def k(x):
        return x

    def f(x):
       return 3 * x + x**2

    # Тестовый пример (Вариант 1)
    '''a = 1
    b = 2
    A_ = 3
    B_ = 0

    def k(x):
        return x**(3)
 
    def f(x):
       return 10 * x**(1/4)'''

   
    print("Задание 1.\nПромоделировать стационарные процессы теплопроводности стержня в зависимости" +
        " от входных данных задачи")
    # решение уравнения для всех наборов и значения у для графика для всех наборов
    u = []
    for i in range(7):
        u.append(findSolution(i+1, f, K, C, A, B, a, b, 1))
    X = np.linspace(a, b, 100)
    y = [[], [], [], [], [], [], []]
    for j in range(7):
        for i in range(len(X)):
            y[j].append(u[j].subs(sp.symbols('x'), X[i]))

    # Пункт 1, 2
    print('\n\n', '~' * 10, 'Пункт 1-2', '~' * 10)
    findSolution(1, f, K, C, A, B, a, b, 2)

    # Пункт 3
    print('\n\n', '~' * 10, 'Пункт 3', '~' * 10)
    for i in range(3):
        findSolution(i + 1, f, K, C, A, B, a, b, 1)

    # Пункт 4
    print('\n\n', '~' * 10, 'Пункт 4', '~' * 10, '\nПостроение графиков')
    plot_res(3, X, [y[0], y[1], y[2]], ['Набор №1', 'Набор №2', 'Набор №3'])

    # Пункт 5
    print('\n\n', '~' * 10, 'Пункт 5', '~' * 10)
    findSolution(4, f, K, C, A, B, a, b, 2)
    plot_res(2, X, [y[0], y[3]], ['Набор №1', 'Набор №4'])

    # Пункт 6
    print('\n\n', '~' * 10, 'Пункт 6', '~' * 10)
    for i in range(5, 8):
        findSolution(i, f, K, C, A, B, a, b, 1)
    plot_res(3, X, [y[4], y[5], y[6]], ['Набор №5', 'Набор №6', 'Набор №7'])
import numpy as np 
import sympy as sp 
import matplotlib.pyplot as plt, matplotlib.pylab


def print_array(y):
    '''печать массива с 6 знаками после запятой'''
    for i in range(len(y)):
        print("{:.6f}".format(y[i]), end="  ")
    print()

def plot_res(n, X, y, label, title): 
    '''печать n графков по известным массивам и названиям, y - двумерный массив, X, label - одномерные,
        n - размерность label, X и первая размерность y'''
    fig = plt.figure()
    for i in range(n):
        plt.plot(X, y[i], label=label[i])
    plt.legend(loc='best')
    plt.title(title)
    plt.show()


def plot_res_pro(n, X, y, label, title):
    '''печать n графков на 2 сист. координат(пополам) по известным массивам и названиям, 
        y - двумерный массив, X, label - одномерные,
        n - размерность label, X и первая размерность y'''
    plt.figure()
    ax1 = plt.subplot(1,2,1)
    ax2 = plt.subplot(1,2,2)
    plt.sca(ax1)
    for i in range(int(n/2)):
        plt.plot(X, y[i], label=label[i])
    plt.legend(loc='best')

    plt.sca(ax2)
    for i in range(int(n/2), n):
        plt.plot(X, y[i], label=label[i])
    plt.legend(loc='best')
    plt.title(title)
    plt.show()

def findSolution(k, f, h, n, X, A_, B_, a_, b_, c1, c2, x0, C):
    '''метод прогонки для решения стационарного уравнения теплопроводности методом баланса (2k -> 1 точка разб.)
        с гр. усл. 1 рода ( u(a) = Ua ) 
        k - массив из 3 const
        f, - функция от х в условии дифф уравнения
        h, n, X, - шаг, кол-во разбиений отрезка, массив точек разбиения
        A_, B_ - из граничных условий, 
        a_, b_, c1, c2 - границы отрезка и точки третьих частей отрезка
        x0 - массив, 1 или 2 точки расположения источников тепла
        С - массив, мощность источника тепла (соотв х0 1 или 2 значения)'''

    # массивы для коэффициентов трехдиагональной матрицы
    a = np.zeros(n-1)
    b = np.zeros(n-1)
    c = np.zeros(n-1)
    d = np.zeros(n-1)

    # коэффициенты трехдиагональной матрицы
    # a = coef_a, b = -(coef_a+coef_b), c = coef_b, d = -f(x)*h^2
    def coef_a(x, k1, k2, k3, c1, c2):
        if x < c1:
            return k1
        elif x-h < c1 < x:
            return h / ((c1 - x + h) / k1 + (x - c1) / k2)
        elif x <= c2:
            return k2
        elif x-h < c2 < x:
            return h / ((c2 - x + h) / k2 + (x - c2) / k3)
        else:
            return k3
    
    def coef_b(x, k1, k2, k3, c1, c2):
        if x < c1:
            return k1
        elif x < c1 < x+h:
            return h / ((c1 - x) / k1 + (x+h - c1) / k2)
        elif x+h <= c2:
            return k2
        elif x < c2 < x+h:
            return h / ((c2 - x) / k2 + (x+h - c2) / k3)
        else:
            return k3

    for i in range(n-1):
        a[i] = coef_a(X[i], k[0], k[1], k[2], c1, c2)
        b[i] = -(coef_a(X[i], k[0], k[1], k[2], c1, c2) + coef_b(X[i], k[0], k[1], k[2], c1, c2))
        c[i] = coef_b(X[i], k[0], k[1], k[2], c1, c2)
        sum = 0
        for j in range(len(x0)):
            sum += f(X[i], x0[j], C[j])
        d[i] = - h * sum

    # подстановка у0 и уn из граничных условий
    d[0] -= a[0] * A_
    d[-1] -= c[-1] * B_

    # сам метод прогонки
    n = n-1
    alpha = np.zeros(n)
    beta = np.zeros(n)

    y = np.zeros(n)
    y[0] = b[0]
    alpha[0] = -c[0] / y[0]
    beta[0] = d[0] / y[0]
    for i in range(1, n-1):
        y[i] = b[i] + a[i]*alpha[i-1]
        alpha[i] = -c[i] / y[i]
        beta[i] = (d[i] - a[i]*beta[i-1]) / y[i]
    y[n-1] = b[n-1] + a[n-1]*alpha[n-2]
    beta[n-1] = (d[n-1] - a[n-1]*beta[n-2]) / y[n-1]
    # обратный ход метода прогонки
    x = np.zeros(n)
    x[n-1] = beta[n-1]
    for i in range(n-2, -1, -1):
        x[i] = alpha[i]*x[i+1] + beta[i]

    # нахождение y0 и yn из граничных условий
    x = np.insert(x, 0, A_)
    x = list(x)
    x.append(B_)
    return x


def punkt4a(f, h, n, X, a, b, A_, B_, k1, k2, x0, C):
    medium = 0.5*(b + a)
    #заглушки, т.к. только 2 материала
    medium2 = b + h 
    k3 = list(np.zeros(len(k1)))
    u, label = [], []
    for i in range(len(k1)):
        u.append(findSolution([k1[i], k2[i], k3[i]], f, h, n, X, A_, B_, a, b, medium, medium2, x0, C))
        label.append('k1={:}, k2={:}'.format(k1[i], k2[i]))
    #print("\nВектор-решение задачи:")
    #print_array(u)
    plot_res(len(u), X, u, label, 'Мощн. ист. С={:}, коорд. x0={:}, k1{:}k2'.format(
        C, x0, '>>' if k1>k2 else '<<'))
    #plot_res_pro(len(u), X, u, label, 'Мощн. ист. С={:}, коорд. x0={:}, k1{:}k2'.format(
    #C, x0, '>>' if k1>k2 else '<<'))


def punkt4b(f, h, n, X, a, b, A_, B_, k1, k2, k3, x0, C):
    c1 = a + (b - a)/3
    c2 = a + 2*(b - a)/3
    u, label = [], []
    for i in range(len(k1)):
        u.append(findSolution([k1[i], k2[i], k3[i]], f, h, n, X, A_, B_, a, b, c1, c2, x0, C))
        label.append('k1={:}, k2={:}, k3={:}'.format(k1[i], k2[i], k3[i]))
    if k1>k2>k3: text = 'k1>k2>k3'
    elif k1<k2<k3: text = 'k1<k2<k3'
    elif k1==k3<k2: text = 'k1=k3=k k2=2k'
    elif k1==k3>k2: text = 'k1=k2=20k k2=k'
    plot_res(len(u), X, u, label, 'Мощн. ист. С={:}, коорд. x0={:}, {:}'.format(C, x0, text))
    #plot_res_pro(len(u), X, u, label, 'Мощн. ист. С={:}, коорд. x0={:}, {:}'.format(C, x0, text))


def task2():

    test = False

    # условие -d/dx(k(x)*du/dx) = f   

    if not test:
        # Вариант 7
        a = 1.0
        b = 2
        A_ = 3
        B_ = 3
    else:
        # Тестовый пример (Вариант 1)
        a = 1
        b = 2
        A_ = 3
        B_ = 0

    def f(x, x0, C):
            if abs(x - x0) - h/2 < 1e-5:
                return C/2
            elif x - h/2 < x0 < x + h/2:
                return C
            else:
                return 0

    print("Задание 2.\nПромоделировать стационарные процессы теплопроводности стержня в зависимости" +
          " от входных данных задачи – переменного коэффициента теплопроводности k(x) и плотности" + 
          " источников тепла f(x)")

    h = (b - a) / 150

    print("\nМетодом баланса (с использованием метода прогонки)")
  
    n = int((b - a) / h)
    X = np.zeros(n+1)
    for i in range(n+1):
        X[i] = a + i * h
    
    # пункт 4a
    if not test:
        k1 = [[0.7, 0.4, 0.1], 
              [50, 10, 2]]
        k2 = [[50, 10, 2],
              [0.7, 0.4, 0.1]]
    else:
        k1, k2 = [[0.1, 0.4], [2, 10]], [[2, 10], [0.1, 0.4]]

    # пункт 5а
    # источник тепла x0 в центре отрезка C - его мощность
    x0 = [(b + a)/2]
    if not test:
        C = [[10], [50], [200]]
    else:
        C = [[50]]
    for i in range(len(k1)):
        for j in range(len(C)):
            punkt4a(f, h, n, X, a, b, A_, B_, k1[i], k2[i], x0, C[j])
    # пункт 5б
    # 2 источникa тепла x0 (симметроично центра), C - мощности (одинаковые)
    x0 = [(b + a)/2-0.3, (b + a)/2+0.3]
    if not test:
        C = [[10, 10], [50, 50], [200, 200]]
    else:
        C = [[50, 50]]
    for i in range(len(k1)):
        for j in range(len(C)):
            punkt4a(f, h, n, X, a, b, A_, B_, k1[i], k2[i], x0, C[j])
    # пункт 5в
    # 2 источникa тепла x0(симметроично центра), C - мощности (разные)
    x0 = [(b + a)/2-0.3, (b + a)/2+0.3]
    if not test:
        C = [[50, 200], [200, 50]]
    else:
        C = [[50, 200]]
    for i in range(len(k1)):
        for j in range(len(C)):
            punkt4a(f, h, n, X, a, b, A_, B_, k1[i], k2[i], x0, C[j])
    # пункт 5г
    # 2 источникa тепла x0(свое расположение (1.3)), C - мощности (разные)
    x0 = [a + (b-a)/8, (b + a)/2 + (b-a)/8]
    if not test:
        C = [[50, 200], [200, 50], [50, 50], [200, 200]]
    else:
        C = [[50, 50], [50, 200]]
    for i in range(len(k1)):
        for j in range(len(C)):
            punkt4a(f, h, n, X, a, b, A_, B_, k1[i], k2[i], x0, C[j])
    
    # пункт 4б
    if not test:
        k1 = [[10, 0.75, 58, 0.6],
              [12, 2.6, 103, 0.9],
              [0.6, 3, 32, 0.9],
              [12, 58, 34, 18]]
        k2 = [[11, 1.2, 85, 0.7],
              [11, 1.2, 85, 0.7],
              [1.2, 6, 64, 1.8],
              [0.6, 2.9, 1.7, 0.9]]
        k3 = [[12, 2.6, 103, 0.9],
              [10, 0.75, 58, 0.6],
              [0.6, 3, 32, 0.9],
              [12, 58, 34, 18]]
    else:
        k1 = [[10, 0.75],
              [12, 2.6],
              [0.6, 3],
              [12, 58]]
        k2 = [[11, 1.2],
              [11, 1.2],
              [1.2, 6],
              [0.6, 2.9]]
        k3 = [[12, 2.6],
              [10, 0.75],
              [0.6, 3],
              [12, 58]]

    # пункт 5а
    # источник тепла x0 в центре отрезка C - его мощность
    x0 = [(b + a)/2]
    if not test:
        C = [[10], [50], [200]]
    else:
        C = [[50]]
    for i in range(len(k1)):
        for j in range(len(C)):
             punkt4b(f, h, n, X, a, b, A_, B_, k1[i], k2[i], k3[i], x0, C[j])

    # пункт 5б
    # 2 источникa тепла x0 (симметроично центра), C - мощности (одинаковые)
    x0 = [(b + a)/2-0.3, (b + a)/2+0.3]
    if not test:
        C = [[10, 10], [50, 50], [200, 200]]
    else:
        C = [[50, 50]] 
    for i in range(len(k1)):
        for j in range(len(C)):
            punkt4b(f, h, n, X, a, b, A_, B_, k1[i], k2[i], k3[i], x0, C[j])

    # пункт 5в
    # 2 источникa тепла x0(симметроично центра), C - мощности (разные)
    x0 = [(b + a)/2-0.3, (b + a)/2+0.3]
    if not test:
        C = [[50, 200], [200, 50]]
    else:
        C = [[50, 200]]
    for i in range(len(k1)):
        for j in range(len(C)):
            punkt4b(f, h, n, X, a, b, A_, B_, k1[i], k2[i], k3[i], x0, C[j])
    
    # пункт 5г
    # 2 источникa тепла x0(свое расположение (1.3)), C - мощности (разные)
    x0 = [a + (b-a)/8, (b + a)/2 + (b-a)/8]
    if not test:
        C = [[50, 200], [200, 50], [50, 50], [200, 200]]
    else:
        C = [[50, 200], [50, 50]]
    for i in range(len(k1)):
        for j in range(len(C)):
            punkt4b(f, h, n, X, a, b, A_, B_, k1[i], k2[i], k3[i], x0, C[j])
    
    

        import numpy as np
import sympy as sp 
import matplotlib.pyplot as plt
from celluloid import Camera


def findSolution(t, h, a, b, A_, B_, T, k, f, fi):
    N = int((b - a)/h) + 1
    M = int(T / t) + 1

    M_ = np.zeros((M, N))
    X = np.linspace(a, b, N)
    
    M_[0,:] = np.array([fi(i) for i in np.linspace(a, b, N)])
    M_[:,0] = np.array([A_ for i in np.linspace(fi(a), fi(b), M)])
    M_[:,-1] = np.array([B_ for xi in np.linspace(fi(a), fi(b), M)])
    for i in range(1, M):
        for j in range(1, N-1):
            M_[i][j] = k(X[j] + h/2)*M_[i-1][j+1]*t/(h*h) + \
                (1 - (k(X[j] + h/2) + k(X[j] - h/2))*t/(h*h)) * M_[i-1][j] + \
                k(X[j] - h/2)*M_[i-1][j-1]*t/(h*h) + t*f(X[j])*(1 - np.exp(-i*t))

    print('t = {:}, h = {:}, M = {:}, N = {:}'.format(t, h, M, N))
    print('fi = {:}, k = {:}, f = {:}\nX = {:}\nMatrix:\n{:}'.format(
        fi(sp.symbols('x')), k(sp.symbols('x')), f(sp.symbols('x')), X, M_))
    return M_, X


def task3():
    # Вариант 3
    a = 0.5
    b = 1.5
    A_ = 2.0
    B_ = 6.0
    l = b - a

    def k(x):
        return x**(-2)
    def f(x):
       return -2.0 * x**2 - 2.*x
    def fi(x):
       return (B_ - A_)*(x - 0.5)/l + A_

    T = abs(fi(b)-fi(a))

    
    # Тестовый пример (Вариант 1)
    '''a = 1
    b = 2
    A_ = 3.0
    B_ = 0
    l = b - a
    def k(x):
        return x**(3)
    def f(x):
       return 10 * x**(1/4)
    def fi(x):
       return (B_ - A_)*(x - a)/l + A_
    T = abs(fi(b)-fi(a))
    '''

    print("Задание 3.\nПромоделировать нестационарные процессы теплопроводности в зависимости от" +
         "входных данных задачи − коэффициента теплопроводности и начальной температуры")

    t, h = 0.0005, 0.1
    #t, h = 0.05, 0.01
    
    if fi(b) > fi(a):
        if fi(b)*t/h**2 > 1./2.:
            print('Wrong h or t')
            return
    else:
        if fi(a)*t/h**2 > 1./2.:
            print('Wrong h or t')
            return    

    Matrix, X = findSolution(t, h, a, b, A_, B_, T, k, f, fi)

    print(Matrix)
    plt.plot(X, Matrix[0], label='T = 0t')
    plt.plot(X, Matrix[5], label='T = 5t')
    plt.plot(X, Matrix[20], label='T = 20t')
    plt.plot(X, Matrix[200], label='T = 200t')
    plt.legend(loc='best')
    plt.show()

    # animation
    fig = plt.figure()
    camera = Camera(fig)
    for i in range(0, len(Matrix), 50):  #400 для рисунка в статике (50 для анимации)
        plt.plot(X, Matrix[i], label=i)
        camera.snap()
    animation = camera.animate()
    animation.save('task3.gif', writer = 'imagemagick')
    #animation.save('task3-test.gif', writer = 'imagemagick')
    
    # решение из задания 1
    import Task1
    C = lambda set: 1
    K = lambda x, set: C(set) * k(x)
    Y2 = (Task1.findSolution(1, f, K, C, lambda set: A_, lambda set: B_, a, b, 0))
    X2 = np.linspace(a, b, 100)
    plt.plot(X2, [Y2.subs(sp.symbols('x'), X2[i]) for i in range(len(X2))], 'r:')
    #plt.legend(loc='best')
    plt.show()

    
    Matrix_2, X_2 = findSolution(t, h, a, b, A_, B_, T, k, f, lambda x: 2*(B_ - A_)*(x - a)/l + A_)
    Matrix_3, X_3 = findSolution(t, h, a, b, A_, B_, T, k, f, lambda x: 0.5*(B_ - A_)*(x - a)/l + A_)

    t_ = [100, 300, 500, 1000]
    for i in range(len(t_)):
        plt.plot(X, Matrix[t_[i]], label='fi = (B_ - A_)*(x - a)/l + A_')
        plt.plot(X_2, Matrix_2[t_[i]], label='fi = 2*(B_ - A_)*(x - a)/l + A_')
        plt.plot(X_3, Matrix_3[t_[i]], label='fi = 0.5*(B_ - A_)*(x - a)/l + A_')
        plt.title('T = {:}t'.format(t_[i]))
        plt.legend(loc='best')
        plt.show()import numpy as np
import scipy as sp
from scipy.misc import derivative 
from scipy import sparse
from scipy import integrate
import matplotlib.pyplot as plt
'''
a = 0
b = 1
ua = lambda t: 5*t
ub = lambda t: 5*t
k = 0.2
f = lambda x: 0
phi = lambda x: 1 - x*x
h = (b - a) / 50
t = 0.5 * h*h /k
T = 0.2


a = 0
b = 1
ua = lambda t: 2
ub = lambda t: 6
k = 0.1
f = lambda x: 0
phi = lambda x: x*(1 - x)
h = (b - a) / 50
t = 0.5 * h*h /k
T = 0.5


nK = int((b - a)/h) + 1
nT = int(T/t)+1
def solve4():
    M = np.zeros(shape=(nT, nK))
    M[:,0] = np.array([ua(i) for i in np.linspace(0, T, nT)])
    M[:,-1] = np.array([ub(i) for i in np.linspace(0, T, nT)])
    M[0,:] = np.array([phi(xi) for xi in np.linspace(a, b, nK)])
    #print(M)
    for i in range(1, nT):
        for j in range(1, nK-1):
            xj = a + j*h
            M[i][j] = k * M[i-1][j+1]*t/(h*h) + \
                      (1 - 2*k*t/(h*h)) * M[i-1][j] + \
                      k*M[i-1][j-1]*t/(h*h) + \
                      t*f(xj)
    return M
matrix = solve4()
X = np.linspace(a, b, nK)

for i in range(0, nT, 50):
    plt.plot(X, matrix[i], label=i)
plt.legend()
plt.show()
'''

import numpy as np
import sympy as sp 
import matplotlib.pyplot as plt
from celluloid import Camera


def findSolution(t, h, a, b, A_, B_, T, k, f, fi):
    N = int((b - a)/h) + 1
    M = int(T / t) + 1

    M_ = np.zeros((M, N))
    X = np.linspace(a, b, N)
    
    M_[0,:] = np.array([fi(i) for i in np.linspace(a, b, N)])
    M_[:,0] = np.array([A_(i) for i in np.linspace(0, T, M)])
    M_[:,-1] = np.array([B_(i) for i in np.linspace(0, T, M)])
    for i in range(1, M):
        for j in range(1, N-1):
            M_[i][j] = t/(h*h)*k*M_[i-1][j+1] + (1 - 2*k*t/(h*h))*M_[i-1][j] + \
                t/(h*h)*k*M_[i-1][j-1] + t*f(X[j])
    
    print('a = {:}, b = {:}, Ua = {:}, Ub = {:}'.format(a, b, A_(sp.symbols('x')), B_(sp.symbols('x'))))
    print('t = {:}, h = {:}, M = {:}, N = {:}'.format(t, h, M, N))
    print('fi = {:}, k = {:}, f = {:}\nX = {:}\nMatrix:'.format(
        fi(sp.symbols('x')), k, f(sp.symbols('x')), X))
    # print matrix 
    for i in range(len(M_)):
        for j in range(len(M_[i])):
            print('{:.4f}'.format(M_[i][j]), end=' ')
        print()

    return M_, X


def task4():
    a = 0
    b = 1
    A_ = lambda t: 3
    B_ = lambda t: 3
    T = 0.02
    k = 2
    def f(x):
       return x*(1-x)
    def fi(x):
       return 0
    
    # Тестовый пример (Вариант 1)
    '''a = 0
    b = 1
    A_ = lambda t: 3
    B_ = lambda t: 0
    T = 0.5
    k = 1
    def f(x):
       return x
    def fi(x):
       return 0
    '''

    print("Задание 4.\nПромоделировать нестационарные процессы теплопроводности в зависимости от" +
         "входных данных задачи. Найти приближенное решение начально-краевой задачи для уравнения" +
         "теплопроводности")

    n = 10
    h = (b - a)/n
    t = 1/2*(h*h/k)
    Matrix, X = findSolution(t, h, a, b, A_, B_, T, k, f, fi)

    for i in range(0, len(Matrix), 2):
        plt.plot(X, Matrix[i], label='{:}t'.format(i))
    plt.legend(loc='best')
    plt.title('h = {:.4f}'.format(h))
    plt.show()

    # animation
    h2 = (b - a)/30
    t2 = 1/2*(h2*h2/k)
    Matrix2, X2 = findSolution(t2, h2, a, b, A_, B_, T, k, f, fi)
    fig = plt.figure()
    camera = Camera(fig)
    for i in range(0, len(Matrix2)):
        plt.plot(X2, Matrix2[i])
        camera.snap()
    animation = camera.animate()
    animation.save('task4.gif', writer = 'imagemagick')
    plt.title('h = {:.4f}'.format(h2))
    plt.show()

