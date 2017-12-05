#!/usr/bin/env python
"""
    Portierung der Matlab-Skripts zu Python-Funktionen
    Matlab-Code Prof. Dr. Uwe Dettmar, uwe.dettmar@th-koeln.de
    Portierung zu Python 8/2017 Joschka Wirges, joschka.wirges@smail.th-koeln.de

    Teil des Praktikums zur Vorlesung Quellen- und Kanalcodierung

"""

import numpy as np, matplotlib.pyplot as plt
import scipy.special
from matplotlib import cm
from scipy import integrate

#converting the error function to the q-function
qfunc = lambda x: 0.5-0.5*scipy.special.erf(x/np.sqrt(2))

def my_range(start, end, step):
    """
    range for a loop with decimal incrementation
    stepsize step in range of [start, end]
    """
    while start <= end:
        yield start
        start += step

def ld(x):
    """
    "safe" workaround for np.log2(0)=0
    H=-sum(p*ld(p)) with p=0 equals 0
    """
    for i in range (0,2):
        if (x[i]>0.0):
            x[i]=np.log2(x[i])
    return x

def entropy(p):
    """
    H=ENTROPY(P) returns the entropy function of
    the probability vector P.
    """
    # if (len(np.nonzero(p < 0)) != 0):
    #     print('Error: Not a prob. vector, negative component(s)')
    #    return 0  # else log2 will be fed with negative numbers, entropy will be false, output 0
    # if (abs(sum(p) - 1) > 10e-10):
    #     print('Error: Not a prob. vector, components do not add up to 1')
    #     return 0  # entropy will be false, output 0
    if type(p) is dict:
        p=np.array(list(p.values()))
    h = sum(-p * ld(p))
    return h

def ip_08_07():
    """
    MATLAB script
    for Illustrative Problem 7, Chapter 8.
    """
    ep = 0.3
    p = np.zeros(61)
    for i in range(1, 61, 2):
        for j in my_range(int((i + 1) / 2), i, 1):
            p[i] = p[i] + np.math.factorial(i) / (np.math.factorial(j) * np.math.factorial(i - j)) * ep ** j * (
                                                                                                               1 - ep) ** (
                                                                                                               i - j)
    t = np.arange(1, 62)
    plt.stem(t, p)
    plt.xlabel('n')
    plt.ylabel('pe')
    plt.title('Error probability as a function of n in simple repetition code')
    plt.show()
    return

def ip_08_08():
    """
    MATLAB script for Illustrative Problem 8, Chapter 8.
    Generate U, denoting all information sequences.
    """
    k = 4
    u = np.zeros((16, 4))
    for i in range(0, 2 ** k):
        for j in range((k - 1), -1, -1):
            if np.remainder(i, 2 ** (-j + k)) >= 2 ** (-j + k - 1):
                u[i, j] = 1
    # Define G, the generator matrix.
    g = np.array([[1, 0, 0, 1, 1, 1, 0, 1, 1, 1],
                  [1, 1, 1, 0, 0, 0, 1, 1, 1, 0],
                  [0, 1, 1, 0, 1, 1, 0, 1, 0, 1],
                  [1, 1, 0, 1, 1, 1, 1, 0, 0, 1]])
    # Generate codewords
    c = np.remainder(u @ g, 2)
    # Find the minimum distance.
    w_min = min(sum((c[2:2 ** k, :]).T))
    print('w_min=', w_min)
    return

def ip_08_09():
    """
    MATLAB script for Illustrative Problem 2, Chapter 8.
    echo on
    """
    k = 11
    u = np.zeros((2048, 11))
    for i in range(0, 2 ** k):
        for j in range((k - 1), -1, -1):
            if np.remainder(i, 2 ** (-j + k)) >= 2 ** (-j + k - 1):
                u[i, j] = 1

    g = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
                  [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
                  [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
                  [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
                  [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1],
                  [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1],
                  [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0],
                  [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1],
                  [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]])

    c = np.remainder(u @ g, 2)
    print('Die ersten fünf Codeworte:\n', c[0:5, :])
    w_min = min(sum((c[2:2 ** k, :]).T))
    print('w_min=', w_min)
    return

def p_e_hd_a(gamma_db_l,gamma_db_h,k,n,d_min):
    """
    function for computing the error probability in
    hard decision decoding of a linear block code
    when antipodal signaling is used.
    [p_err,gamma_db]=p_e_hd_a(gamma_db_l,gamma_db_h,k,n,d_min)
    gamma_db_l=lower E_b/N_0
    gamma_db_h=higher E_b/N_0
    k=number of information bits in the code
    n=code block length
    d_min=minimum distance of the code
    """
    gamma_db=np.zeros(21)
    gamma_db=np.arange(gamma_db_l,gamma_db_h,(gamma_db_h-gamma_db_l)/20)
    gamma_b=10**(gamma_db/10)
    R_c=k/n
    p_b=qfunc(np.sqrt(2*R_c*gamma_b))
    p_err=(2**k-1)*(4*p_b*(1-p_b))**(d_min/2)
    return p_err,gamma_db

def p_e_hd_o(gamma_db_l,gamma_db_h,k,n,d_min):
    """
    function for computing the error probability in
    hard decision decoding of a linear block code
    when orthogonal signaling is used.
    [p_err,gamma_db]=p_e_hd_o(gamma_db_l,gamma_db_h,k,n,d_min)
    gamma_db_l=lower E_b/N_0
    gamma_db_h=higher E_b/N_0
    k=number of information bits in the code
    n=code block length
    d_min=minimum distance of the code
    """
    gamma_db=np.zeros(21)
    gamma_db=np.arange(gamma_db_l,gamma_db_h,(gamma_db_h-gamma_db_l)/20)
    gamma_b=10**(gamma_db/10)
    R_c=k/n
    p_b=qfunc(np.sqrt(R_c*gamma_b))
    p_err=(2**k-1)*(4*p_b*(1-p_b))**(d_min/2)
    return p_err,gamma_db

def p_e_sd_o(gamma_db_l,gamma_db_h,k,n,d_min):
    """
    function for computing the error probability in
    soft decision decoding of a linear block code
    when orthogonal signaling is used.
    [p_err,gamma_db]=p_e_sd_o(gamma_db_l,gamma_db_h,k,n,d_min)
    gamma_db_l=lower E_b/N_0
    gamma_db_h=higher E_b/N_0
    k=number of information bits in the code
    n=code block length
    d_min=minimum distance of the code
    """
    gamma_db=np.zeros(21)
    gamma_db=np.arange(gamma_db_l,gamma_db_h,(gamma_db_h-gamma_db_l)/20)
    gamma_b=10**(gamma_db/10)
    R_c=k/n
    p_err=(2**k-1)*qfunc(np.sqrt(d_min*R_c*gamma_b/2))
    return p_err,gamma_db

def p_e_sd_a(gamma_db_l,gamma_db_h,k,n,d_min):
    """
    function for computing the error probability in
    soft decision decoding of a linear block code
    when antipodal signaling is used.
    [p_err,gamma_db]=p_e_sd_a(gamma_db_l,gamma_db_h,k,n,d_min)
    gamma_db_l=lower E_b/N_0
    gamma_db_h=higher E_b/N_0
    k=number of information bits in the code
    n=code block length
    d_min=minimum distance of the code
    """
    gamma_db=np.zeros(21)
    gamma_db=np.arange(gamma_db_l,gamma_db_h,(gamma_db_h-gamma_db_l)/20)
    gamma_b=10**(gamma_db/10)
    R_c=k/n
    p_err=(2**k-1)*qfunc(np.sqrt(d_min*R_c*gamma_b))
    return p_err,gamma_db

def ip_08_12():
    """
    MATLAB script for Illustrative Problem 12, Chapter 8.
    """
    [p_err_ha,gamma_b]=p_e_hd_a(10,16,11,15,3)
    [p_err_ho,gamma_b]=p_e_hd_o(10,16,11,15,3)
    [p_err_so,gamma_b]=p_e_sd_o(10,16,11,15,3)
    [p_err_sa,gamma_b]=p_e_sd_a(10,16,11,15,3)

    line1, = plt.semilogy(gamma_b,p_err_ha,label="HD antipodal signalling")
    line2, = plt.semilogy(gamma_b,p_err_ho,label="HD orthogonal signalling")
    line3, = plt.semilogy(gamma_b,p_err_so,label="SD orthogonal signalling")
    line4, = plt.semilogy(gamma_b,p_err_sa,label="SD antipodal signalling")
    plt.legend(handles=[line1, line2, line3, line4,], loc=3)
    plt.show()
    return

def entropy2(p):
    """
    h=entropy2(p) Returns the binary entropy function of the components
    of the vevtor p.
    """
    n=len(p)
    h=np.zeros(n)
    for i in range(0,n):
        p1=np.array([p[i], 1-p[i]])
        h[i]=entropy(p1);
    return h

def ip_08_01():
    """
    MATLAB script for Illustrative Problem 1, Chapter 8.
    """
    gamma_db = np.arange(-20, 20, 0.1)
    gamma = 10 ** (gamma_db / 10)
    p_error = qfunc(np.sqrt(2 * gamma))
    capacity = 1 - entropy2(p_error)
    fig1 = plt.semilogx(gamma, p_error)
    plt.xlabel('SNR/bit')
    plt.title('Error probability versus SNR/bit')
    plt.ylabel('Error Prob.')
    plt.grid(True, which="both", ls="-")
    plt.show()
    fig2 = plt.semilogx(gamma, capacity)
    plt.xlabel('SNR/bit')
    plt.title('Channel capacity versus SNR/bit')
    plt.ylabel('Channel capacity')
    plt.grid(True, which="both", ls="-")
    plt.show()
    return

def ip_08_05():
    """
    MATLAB script for Illustrative Problem 5, Chapter 8.
    """
    w = np.array([np.arange(1, 20, 5), np.arange(25, 100, 20), np.arange(130, 300, 50)]).ravel()
    w = np.append(w, np.arange(400, 1100, 100))
    w = np.append(w, np.arange(1250, 5250, 250))
    w = np.append(w, np.arange(5500, 10500, 500))
    pn0_db = np.arange(-20, 31, 1)
    pn0 = 10 ** (pn0_db / 10)
    c = np.zeros((45, 51))
    for i in range(0, 45):
        for j in range(0, 51):
            c[i, j] = w[i] * np.log2(1 + pn0[j] / w[i])
    k = np.array([0.9, 0.8, 0.5, 0.6])
    s = np.array([-70, 35])

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    w, pn0_db = np.meshgrid(w, pn0_db)
    surf = ax.plot_surface(w, pn0_db, c.T, cmap=cm.coolwarm,
                           antialiased=True)  # https://matplotlib.org/examples/color/colormaps_reference.html
    plt.title('Capacity in Infobit per second vs. bandwidth W and P/N_0');
    ax.set_xlabel('single sided bandwidth W')
    ax.set_ylabel(r'$P/N_0$')
    ax.set_zlabel('Capacity in Bits per s')
    ax.view_init(45, 225)
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
    return

def ip_08_04():
    """
    MATLAB script for Illustrative Problem 4, Chapter 8.
    """
    a_db = np.arange(-13, 13.5, 0.5)
    a = 10 ** (a_db / 10)
    a_hard = a.copy()
    c_hard = 1 - entropy2(qfunc(a_hard))
    f = np.zeros(53)
    g = np.zeros(53)
    c_soft = np.zeros(53)

    il3_8fun = lambda x, p: 1 / np.sqrt(2 * np.pi) * np.exp((-(x - p) ** 2) / 2) * np.log2(2 / (1 + np.exp(-2 * x * p)))

    for i in range(0, 53):
        f[i] = integrate.quadrature(il3_8fun, a[i] - 5, a[i] + 5, args=(a[i],), tol=1e-3)[0]  # ,1e-3,[],a[i]
        g[i] = integrate.quadrature(il3_8fun, -a[i] - 5, -a[i] + 5, args=(-a[i],), tol=1e-3)[0]  # ,1e-3
        c_soft[i] = 0.5 * f[i] + 0.5 * g[i]
    plt.title('Capacity for BPSK transmisison on the AWGN channel for Hard and Soft Decision')
    plt.xlabel(r'A/$\sigma$')
    plt.ylabel('Capacity in bits Per channel Use')
    # plt.grid(True,which="both",ls="-")
    line1, = plt.semilogx(a, c_soft, label="Soft Decision")
    line2, = plt.semilogx(a_hard, c_hard, label="Hard Decision")
    plt.legend(handles=[line1, line2, ], loc=4)
    plt.show()
    return