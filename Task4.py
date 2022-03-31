import math
import numpy as np
import sympy as sp 
import PrintMatrix as PM


def explicit(t, h, a, b, g1, g2, T, k, f, fi, gy, print_res):
    N = int((b - a)/h) + 1
    M = math.ceil(T / t) + 1

    M_ = np.zeros((M, N))
    X = np.linspace(a, b, N)
    
    M_[0,:] = np.array([fi(i) for i in np.linspace(a, b, N)])
    M_[:,0] = np.array([g1(i) for i in np.linspace(0, T, M)])

    for i in range(1, M):
        for j in range(1, N-1):
            M_[i][j] = k*t/(h*h) * M_[i-1][j-1] + (1 - 2*k*t/(h*h)) * M_[i-1][j] + \
                k*t/(h*h) * M_[i-1][j+1] + t*f(X[j], i*t)
        if (gy == 1):
            M_[i][-1] = h*g2(i*t) + M_[i][-2]
        elif (gy == 2):
            M_[i][-1] = M_[i][-3] + 2 * h * g2(t * i)
        else:
            print("Wrong gy")
            exit()

    if (print_res):
        print('t = {:}, h = {:}, M = {:}, N = {:}'.format(t, h, M, N))
        print('fi = {:}, k = {:}, f = {:}\nX = {:}\nMatrix:'.format(
            fi(sp.symbols('x')), k, f(sp.symbols('x'), sp.symbols('t')), X))
        PM.print_matrix(M_)
    return M_, X

import math
import numpy as np
import sympy as sp 
import PrintMatrix as PM


def implicit(t, h, a, b, g1, g2, T, k, f, fi, gy, print_res):
    N = int((b - a)/h) + 1
    M = math.ceil(T / t) + 1

    M_ = np.zeros((M, N))
    X = np.linspace(a, b, N)

    M_[0,:] = np.array([fi(i) for i in np.linspace(a, b, N)])
    M_[:,0] = np.array([g1(i) for i in np.linspace(0, T, M)])
    
    for i in range(1, M):
        M2 = np.zeros(shape=(N-2, N-2))
        Y = np.zeros(N-2)
        M2[0][0] = -(1 + 2*t*k / (h*h))
        M2[0][1] = t*k / (h*h)

        Y[0] = -(M_[i-1][1] + t*f(1,i*t) + t*k/(h*h) * g1(i*t))

        for j in range(1, N-3):
            M2[j][j-1] = t*k / (h*h)
            M2[j][j] = -(1 + 2*t*k / (h*h))
            M2[j][j+1] = t*k / (h*h)
            Y[j] = -(M_[i-1][j+1] + t*f(X[j], i*t))

        if (gy == 1):
            M2[-1][-1] = -(1 + t*k/(h*h))
            M2[-1][-2] = t*k/(h*h)
            Y[-1] = -(M_[i-1][-2] + t*f(X[-2], i*t) + t*k/(h*h) * h*g2(i*t))
        elif (gy == 2):
            M2[-1][-1] = -(1 + 2*t*k/(h*h))
            M2[-1][-2] = 2*t*k/(h*h)
            Y[-1] = -(M_[i-1][-2] + t*f(X[-2], i*t) + 2*t*k/(h*h) * h*g2(i*t))
        else:
            print("Wrong gy")
            exit()
        

        res = list(np.linalg.solve(M2, Y))
        if (gy == 1):
            res.append(res[-1] + h*g2(i*t))
        else:
            res.append(res[-2] + 2*h*g2(i*t))

        #print('\n\nres \n', res)

        M_[i, 1:] = res
        
        #print('\n\nM2 \n', M2)
        #print('\n\nY \n', Y)
    
    if (print_res):
        print('t = {:}, h = {:}, M = {:}, N = {:}'.format(t, h, M, N))
        print('fi = {:}, k = {:}, f = {:}\nX = {:}\nMatrix:'.format(
            fi(sp.symbols('x')), k, f(sp.symbols('x'), sp.symbols('t')), X))
        PM.print_matrix(M_)
    return M_, X

import ExplicitSchema as ES
import ImplicitSchema as IS
import PrintTable as PT
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

var = 3
print("ЛР 14. Аппроксимации граничных условий второго рода в методе конечных разностей на примере" +
      " уравнения теплопроводности.\n")
print("\t\t\tВариант {0}".format(var))
print("-" * 50)


# Вариант 3
a = 0
b = 1
g1 = lambda t: 5*t
g2 = lambda t: 5*t
T = 0.5
k = 0.1
def f(x, t):
    return 0
def fi(x):
    return x*(1 - x)


n = 10
h = (b - a) / n
t_cond = 1
if (t_cond == 1):
    t = (h*h) / (2*k)
else:
    t = (h*h) / 6

# step of grafic
num = 2

e11 = PT.table(a, b, T, g1, g2, k, f, fi, gy=1, t_cond=1, schema_cond='explicit')
i1 = PT.table(a, b, T, g1, g2, k, f, fi, gy=1, t_cond=1, schema_cond='implicit')
e12 = PT.table(a, b, T, g1, g2, k, f, fi, gy=1, t_cond=2, schema_cond='explicit')
e21 = PT.table(a, b, T, g1, g2, k, f, fi, gy=2, t_cond=1, schema_cond='explicit')
i2 = PT.table(a, b, T, g1, g2, k, f, fi, gy=2, t_cond=1, schema_cond='implicit')
e22 = PT.table(a, b, T, g1, g2, k, f, fi, gy=2, t_cond=2, schema_cond='explicit')

# table for comparison
tab = np.zeros(shape=(4, 6))
tab[:,0] = [r[2] for r in e11]
tab[:,1] = [r[2] for r in i1]
tab[:,2] = [r[2] for r in e12]
tab[:,3] = [r[2] for r in e21]
tab[:,4] = [r[2] for r in i2]
tab[:,5] = [r[2] for r in e22]
df = pd.DataFrame(tab, columns=['e11', 'i1', 'e12', 'e21', 'i2', 'e22'])
df[:-1]
print(df)

Matrix, X = ES.explicit(t, h, a, b, g1, g2, T, k, f, fi, 1, print_res=True)
for i in range(0, len(Matrix), num):
    plt.plot(X, Matrix[i], label='{:}t'.format(i))
plt.legend(loc='best')
plt.title('Явная схема, 1 способ. h = {:.4f}'.format(h))
plt.show()

Matrix, X = IS.implicit(t, h, a, b, g1, g2, T, k, f, fi, 1, print_res=True)
for i in range(0, len(Matrix), num):
    plt.plot(X, Matrix[i], label='{:}t'.format(i))
plt.legend(loc='best')
plt.title('Неявная схема, 1 способ. h = {:.4f}'.format(h))
plt.show()

Matrix, X = ES.explicit(t, h, a, b, g1, g2, T, k, f, fi, 2, print_res=True)
for i in range(0, len(Matrix), num):
    plt.plot(X, Matrix[i], label='{:}t'.format(i))
plt.legend(loc='best')
plt.title('Явная схема, 2 способ. h = {:.4f}'.format(h))
plt.show()

Matrix, X = IS.implicit(t, h, a, b, g1, g2, T, k, f, fi, 2, print_res=True)
for i in range(0, len(Matrix), num):
    plt.plot(X, Matrix[i], label='{:}t'.format(i))
plt.legend(loc='best')
plt.title('Неявная схема, 2 способ. h = {:.4f}'.format(h))
plt.show()


'''
PT.table(a, b, T, g1, g2, k, f, fi, gy=1, t_cond=1, schema_cond='explicit')
PT.table(a, b, T, g1, g2, k, f, fi, gy=1, t_cond=1, schema_cond='implicit')
PT.table(a, b, T, g1, g2, k, f, fi, gy=1, t_cond=2, schema_cond='explicit')
PT.table(a, b, T, g1, g2, k, f, fi, gy=2, t_cond=1, schema_cond='explicit')
PT.table(a, b, T, g1, g2, k, f, fi, gy=2, t_cond=1, schema_cond='implicit')
PT.table(a, b, T, g1, g2, k, f, fi, gy=2, t_cond=2, schema_cond='explicit')
'''
  

def print_matrix(M):

    for i in range(len(M)):
        for j in range(len(M[i])):
            print('{:.4f}'.format(M[i][j]), end='   ')
        print()

import numpy as np
import ExplicitSchema as ES
import ImplicitSchema as IS
import pandas as pd
import pylab


def table(a, b, T, g1, g2, k, f, fi, gy, t_cond, schema_cond):

    ns = [10, 20, 40, 80]#, 160]

    table = np.zeros(shape=(len(ns), 8))
    for i in range(0, len(ns)):

        h = (b - a) / ns[i]
        N = ns[i]

        if (t_cond == 1):
            t = (h*h) / (2*k)
        else:
            t = (h*h) / 6

        M = int(T / t) + 1
        # 1/4 and center of tau coordinate 
        tn1 = int(M / 4)
        tn2 = int(M / 2)
    
        if (schema_cond == 'explicit'):
            M1 = ES.explicit(t, h, a, b, g1, g2, T, k, f, fi, gy, print_res=False)[0]
            M2 = ES.explicit(t, h*2, a, b, g1, g2, T, k, f, fi, gy, print_res=False)[0]
        elif (schema_cond == 'implicit'):
            M1 = IS.implicit(t, h, a, b, g1, g2, T, k, f, fi, gy, print_res=False)[0]
            M2 = IS.implicit(t, h*2, a, b, g1, g2, T, k, f, fi, gy, print_res=False)[0]
        else:
            printf('Wrong schema name')
            exit()
        
        sol1l = M1[tn1]
        sol1c = M2[tn1]
        sol2l = M1[tn2]
        sol2c = M2[tn2]  

        table[i][0] = ns[i]
        table[i][1] = t
        table[i][2] = np.sqrt(sum((sol1c - sol1l[::2]) ** 2)/len(sol1c))
        table[i][3] = np.sqrt(sum((sol2c - sol2l[::2]) ** 2)/len(sol1c))
        table[i][4] = max(abs(sol1c - sol1l[::2]))
        table[i][5] = max(abs(sol2c - sol2l[::2]))
        table[i][6] = tn1*t
        table[i][7] = tn2*t

    #print_table(t_cond, gy, schema_cond, table, ns, a, b, k)

    return table



def print_table(t_cond, gy, schema_cond, table, ns, a, b, k):
    SUB = str.maketrans("12", "₁₂")
    df = pd.DataFrame(table, columns=['N', 't', 's(t=t1)'.translate(SUB), 's(t=t2)'.translate(SUB), 
                                      'max|t1|'.translate(SUB), 'max|t2|'.translate(SUB), 't1', 't2'])
    df[:-1]

    if (t_cond == 1):
        formula = 't = (h*h) / (2*k)'
    else:
        formula = 't = (h*h) / 6' 
    print('\n\n{:} schema. GY variant {:}. formula {:}\n'.format(schema_cond, gy, formula))
    print(df)
    
    error = [r[2] for r in table]
    h = [(b - a) / i for i in ns]
    t = [(hi*hi) / (2*k) for hi in h]
    #print(error)
    #print(h)

    figure = pylab.figure()
    ax1 = figure.add_subplot (1, 2, 1)
    pylab.plot (h, error, 'r')
    #ax1.set_xscale ('log')
    ax1.set_yscale ('log')
    ax1.set_xlabel('h')
    ax1.set_ylabel('error')

    ax2 = figure.add_subplot (1, 2, 2)
    pylab.plot (t, error)
    #ax2.set_xscale ('log')
    ax2.set_yscale ('log')
    ax2.set_xlabel('t')
    ax2.set_ylabel('error')
    pylab.show()

