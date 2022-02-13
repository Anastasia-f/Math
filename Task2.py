import numpy as np
import logging as log
import matplotlib.pyplot as plt


def ndiag_forward(A, shift):
    h, w = A.shape
    if shift <= 0: # shift_v
        return
    B = np.concatenate((A, np.zeros((shift - 1, w), dtype=float)))
    indices = [[*range(shift - i, w - 1 - i), w - 1] for i in range(shift + 1)]
    cells = ([(n,) for n in range(shift)], indices[1:])
    for row in range(0, len(A) - 1):
        q = B[row + 1: row + shift + 1][cells][:,0] / B[row, shift]
        B[row + 1: row + shift + 1][cells] -= np.outer(q, B[row][indices[0]])
    np.copyto(A, B[:h])


def ndiag_backward(A, shift):
    h, w = A.shape
    offset = w - shift - 2
    if offset <= 0: # shift_h
        return
    B = np.concatenate((np.zeros((offset - 1, w), dtype=float), A))
    indices = [[w - 2 - i, w - 1] for i in range(offset + 1)]
    cells = ([(n,) for n in range(offset)], indices[:-1])
    for row in range(-1, -h, -1):
        q = B[row - offset: row][cells][:,0] / B[row, shift]
        B[row - offset: row][cells] -= np.outer(q, B[row][indices[-1]])
    np.copyto(A, B[-h:])


def ndiag_method(A, shift):
    ndiag_forward(A, shift)
    ndiag_backward(A, shift)
    A[:,-1] /= A[:,shift]
    A[:,shift] = 1
    return A[:,-1]


def dschema_common(aprxs, aprx_a, aprx_b, funcs, a, b, shift, eps, n=5):
    result = []
    while(len(result) < 2 or np.abs(result[-1][1][1:-1:2] - (result[-2][1][:-1] + result[-2][1][1:]) / 2).max() >= eps):
        h = (b-a)/n
        xs = np.linspace(a, b, n+1)
        aprxs_h = np.array([aprx(h) for aprx in aprxs], dtype=float)
        A = np.array(
            [aprx_a(h)] +
            [sum([ aprxs_h[i] * funcs[i](x) for i in range(len(aprxs_h)) ]) for x in xs[1:-1]] +
            [aprx_b(h)]
        )
        ys = ndiag_method(A, shift)
        result.append((xs, ys))
        n *= 2
    return result


def make_aprxs(shift_v=1, shift_h=1):
    aprxs = np.zeros((4, shift_v + shift_h + 2), dtype=float)
    base_aprxs = np.array([
        [0, 0, 0, 1],
        [0, 1, 0, 0],
        [-1, 0, 1, 0],
        [1, -2, 1, 0],
    ], dtype=float)
    aprxs[:, shift_v - 1: shift_v + 2] = base_aprxs[:, :-1]
    aprxs[:, -1:] = base_aprxs[:, -1:]
    return [lambda h: aprxs[0] * (h * h),
            lambda h: aprxs[1] * (h * h),
            lambda h: aprxs[2] * (h / 2),
            lambda h: aprxs[3]]


def make_edge_aprxs_kind1(qa, fa, qb, fb):
    return [lambda h: [0, qa, 0, fa],
            lambda h: [0, qb, 0, fb]]


def make_edge_aprxs_kind2(pa, qa, fa, pb, qb, fb):
    return [lambda h: [0, 0, -3 * pa + qa * h * 2, 4 * pa, -pa, fa * h * 2],
            lambda h: [pb, -4 * pb, 3 * pb + qb * h * 2, 0, 0, fb * h * 2]]


def make_edge_aprxs_balanced(pa, qa, fa, int_rka, int_qa, int_fa,
                             pb, qb, fb, int_rkb, int_qb, int_fb):
    return [lambda h: [0, int_qa(h) - 1 / int_rka(h) + qa / pa, 1 / int_rka(h), int_fa(h) + fa / pa],
            lambda h: [1 / int_rkb(h), int_qb(h) - 1 / int_rkb(h) - qb / pb, 0, int_fb(h) - fb / pb]]


def show_result(pairs, f_exact=None):
    m = len(pairs) - 1
    for i in range(m + 1):
        plt.plot(pairs[i][0],
                 pairs[i][1],
                 color=np.clip((4 * i / m - 2, -4 * i / m + 2, 2 - abs(4 * i / m - 2)), 0., 1.),
                 label="h={}".format( (pairs[i][0][-1] - pairs[i][0][0]) / (len(pairs[i][0]) - 1) ))
    if f_exact:
        xs = np.linspace(pairs[0][0][0], pairs[0][0][-1], 1000)
        plt.plot(xs, [f_exact(x) for x in xs], color=(1,0.5,0), label="exact")
    plt.legend()
    plt.grid()
    plt.show()


def show_eps(pairs):
    m = len(pairs) - 1
    for i in range(1, m + 1):
        ys = pairs[i][1][1:-1:2] - (pairs[i-1][1][:-1] + pairs[i-1][1][1:]) / 2
        plt.plot(np.hstack((pairs[i][0][0], pairs[i][0][1:-1:2], pairs[i][0][-1])),
                 np.abs(np.hstack((pairs[i][1][0] - pairs[i-1][1][0], ys, pairs[i][1][-1] - pairs[i-1][1][-1]))),
                 color=np.clip((4 * i / m - 2, -4 * i / m + 2, 2 - abs(4 * i / m - 2)), 0., 1.),
                 label="h={}".format( (pairs[i][0][-1] - pairs[i][0][0]) / (len(pairs[i][0]) - 1) ))
    plt.legend()
    plt.grid()
    plt.show()


def start_test(variant):
    np.set_printoptions(precision=3, suppress=True)
    log.basicConfig(level=log.INFO, format="%(message)s")

    log.info("TASK1")
    result = dschema_common(
        make_aprxs(1, 1),
        *make_edge_aprxs_kind1(1, 0, 1, 0),
        [lambda x: -1,
         lambda x: 1 + np.cos(variant) * x * x,
         lambda x: 0,
         lambda x: np.sin(variant)],
        -1., 1., 1, 0.001
    )
    show_result(result)
    show_eps(result)

    log.info("TASK2")
    result = dschema_common(
        make_aprxs(1, 1),
        *make_edge_aprxs_kind1(1, 4, 1, 0),
        [lambda x: 8 * x * x,
         lambda x: np.exp(-x) * (8 + x * x),
         lambda x: np.log(1 + x * x),
         lambda x: 1],
        0., 2., 1, 0.01
    )
    show_result(result)
    show_eps(result)

    log.info("TASK3")
    result = dschema_common(
        make_aprxs(2, 2),
        *make_edge_aprxs_kind2(1, 0, 0, -3, 1, 2),
        [lambda x: 2. * x,
         lambda x: 5.,
         lambda x: -4 * x,
         lambda x: 1],
        2., 4., 2, 0.07
    )
    show_result(result)
    show_eps(result)

    log.info("TASK4")
    result = dschema_common(
        make_aprxs(1, 1),
        *make_edge_aprxs_balanced(
          1, 0.5, 0,
          lambda h: h / -1.5, lambda h: 8.3 * h*0.5, lambda h: -14 * (np.exp(-0.5*0.5*h) - np.exp(-0.5*0.)),
          -1, 0.5, 0,
          lambda h: h / -0.6, lambda h:  12 * h*0.5, lambda h: -14 * (np.exp(-0.5*3.) - np.exp(-0.5*(3.-0.5*h))) ),
        [lambda x: 7*np.exp(-0.5 * x),
         lambda x: 8.3 if x < 1.875 else 12,
         lambda x: 0,
         lambda x: -1.5 if x < 1.875 else -0.6],
        0., 3., 1, 0.0001
    )
    show_result(result)
    show_eps(result)

    return 0


if __name__=="__main__":
    try:
        start_test(7)
    except ValueError as e:
        log.error("Exception:\n%s\nProgram terminated.", e)
    except KeyboardInterrupt:
        log.error("Program interrupted by user.")
    else:
        log.info("Test is done.")

