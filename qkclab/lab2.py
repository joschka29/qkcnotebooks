#!/usr/bin/env python
"""
    Portierung der Matlab-Skripts zu Python-Funktionen
    Matlab-Code Prof. Dr. Uwe Dettmar, uwe.dettmar@th-koeln.de
    Portierung zu Python 8/2017 Joschka Wirges, joschka.wirges@smail.th-koeln.de

    Teil des Praktikums zur Vorlesung Quellen- und Kanalcodierung

"""

import numpy as np, matplotlib.pyplot as plt
#from ipywidgets import interact
#import ipywidgets as widgets

def my_range(start, end, step):
    """
    range for a loop with decimal incrementation
    stepsize step in range of [start, end]
    """
    while start <= end:
        yield start
        start += step

def ip_08_07():
    """
    MATLAB script
    for Illustrative Problem 7, Chapter 8.
    """
    ep = 0.3
    p=np.zeros(61)
    for i in range(1,61,2): #0-60 statt 1-61
        for j in my_range(int((i+1)/2),i,1):
            p[i] = p[i] + np.prod(np.arange(1,i)) / (np.prod(np.arange(1,j))*np.prod(np.arange(1, (i - j)))*ep**j * (1 - ep) **(i - j))
    t=np.arange(1,62)
    print('länge von t',len(t))
    print('länge von p', len(p))
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
    k=4
    for i in range(1,2**k):
      for j in range(k,1,-1):
        if rem(i-1,2^(-j+k+1))>=2^(-j+k)
          u(i,j)=1
        else
          u(i,j)=0
    # Define G, the generator matrix.
    g=[1 0 0 1 1 1 0 1 1 1,
       1 1 1 0 0 0 1 1 1 0,
       0 1 1 0 1 1 0 1 0 1,
       1 1 0 1 1 1 1 0 0 1]
    # Generate codewords
    c=rem(u*g,2)
    # Find the minimum distance.
    w_min=min(sum((c[2:2**k,:]).T))
    print('w_min=',w_min)
    return

def ip_08_09():
    """
    MATLAB script for Illustrative Problem 2, Chapter 8.
    echo on
    """
    k=11
    for i in range(1,2**k):
      for j in range(k,1,-1):
        if rem(i-1,2^(-j+k+1))>=2^(-j+k)
          u(i,j)=1
        else
          u(i,j)=0
    g=[1 0 0 0 0 0 0 0 0 0 0 1 1 0 0;
       0 1 0 0 0 0 0 0 0 0 0 0 1 1 0;
       0 0 1 0 0 0 0 0 0 0 0 0 0 1 1;
       0 0 0 1 0 0 0 0 0 0 0 1 0 1 0;
       0 0 0 0 1 0 0 0 0 0 0 1 0 0 1;
       0 0 0 0 0 1 0 0 0 0 0 0 1 0 1;
       0 0 0 0 0 0 1 0 0 0 0 1 1 1 0;
       0 0 0 0 0 0 0 1 0 0 0 0 1 1 1;
       0 0 0 0 0 0 0 0 1 0 0 1 0 1 1;
       0 0 0 0 0 0 0 0 0 1 0 1 1 0 1;
       0 0 0 0 0 0 0 0 0 0 1 1 1 1 1]
    c=rem(u*g,2)
    print('All Codewords:',c)
    print('only the first 5 Codewords',c[1:5,:])
    w_min=min(sum((c[2:2**k,:]).T))
    print('w_min=', w_min)
    return

def p_e_hd_a(gamma_db_l,gamma_db_h,k,n,d_min):
    """
    p_e_hd_a.m 	Matlab function for computing error probability in
    hard decision decoding of a linear block code
    when antipodal signaling is used.
    [p_err,gamma_db]=p_e_hd_a(gamma_db_l,gamma_db_h,k,n,d_min)
    gamma_db_l=lower E_b/N_0
    gamma_db_h=higher E_b/N_0
    k=number of information bits in the code
    n=code block length
    d_min=minimum distance of the code
    """
    gamma_db=[gamma_db_l:(gamma_db_h-gamma_db_l)/20:gamma_db_h]
    gamma_b=10**(gamma_db/10)
    R_c=k/n
    p_b=q(sqrt(2*R_c*gamma_b))
    p_err=(2**k-1)*(4*p_b*(1-p_b))**(d_min/2)
    return p_err,gamma_db

def p_e_hd_o(gamma_db_l,gamma_db_h,k,n,d_min):
    """
    p_e_hd_o.m 	Matlab function for computing error probability in
                	hard decision decoding of a linear block code
                	when orthogonal signaling is used.
      		[p_err,gamma_db]=p_e_hd_o(gamma_db_l,gamma_db_h,k,n,d_min)
      		gamma_db_l=lower E_b/N_0
      		gamma_db_h=higher E_b/N_0
    		k=number of information bits in the code
      		n=code block length
      		d_min=minimum distance of the code
    """
    gamma_db=[gamma_db_l:(gamma_db_h-gamma_db_l)/20:gamma_db_h]
    gamma_b=10**(gamma_db/10)
    R_c=k/n
    p_b=q(sqrt(R_c*gamma_b))
    p_err=(2**k-1)*(4*p_b*(1-p_b))**(d_min/2)
    return p_err,gamma_db

def p_e_sd_o(gamma_db_l,gamma_db_h,k,n,d_min):
    """
     p_e_sd_o.m Matlab function for computing error probability in
                soft decision decoding of a linear block code
               	when orthogonal signaling is used.
     		[p_err,gamma_db]=p_e_sd_o(gamma_db_l,gamma_db_h,k,n,d_min)
     		gamma_db_l=lower E_b/N_0
    		gamma_db_h=higher E_b/N_0
    		k=number of information bits in the code
    	    n=code block length
     		d_min=minimum distance of the code
    """
    gamma_db=[gamma_db_l:(gamma_db_h-gamma_db_l)/20:gamma_db_h]
    gamma_b=10**(gamma_db/10)
    R_c=k/n
    p_err=(2**k-1)*q(sqrt(d_min*R_c.*gamma_b/2))
    return p_err,gamma_db

def p_e_sd_a(gamma_db_l,gamma_db_h,k,n,d_min):
    """
     p_e_sd_a.m 	Matlab function for computing error probability in
                	soft decision decoding of a linear block code
                	when antipodal signaling is used.
      		[p_err,gamma_db]=p_e_sd_a(gamma_db_l,gamma_db_h,k,n,d_min)
      		gamma_db_l=lower E_b/N_0
      		gamma_db_h=higher E_b/N_0
      		k=number of information bits in the code
      		n=code block length
      		d_min=minimum distance of the code
    """
    gamma_db=[gamma_db_l:(gamma_db_h-gamma_db_l)/20:gamma_db_h]
    gamma_b=10**(gamma_db/10)
    R_c=k/n
    p_err=(2**k-1)*q(sqrt(d_min*R_c*gamma_b))
    return p_err,gamma_db

def ip_08_12():
    """
    MATLAB script for Illustrative Problem 12, Chapter 8.
    """
    [p_err_ha,gamma_b]=p_e_hd_a(10,16,11,15,3)
    [p_err_ho,gamma_b]=p_e_hd_o(10,16,11,15,3)
    [p_err_so,gamma_b]=p_e_sd_o(10,16,11,15,3)
    [p_err_sa,gamma_b]=p_e_sd_a(10,16,11,15,3)
    semilogy(gamma_b,p_err_sa,gamma_b,p_err_so,gamma_b,p_err_ha,gamma_b,p_err_ho)
    return