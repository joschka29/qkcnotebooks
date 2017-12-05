#!/usr/bin/env python
"""
    Portierung der Matlab-Skripts zu Python-Funktionen
    Matlab-Code Prof. Dr. Uwe Dettmar, uwe.dettmar@th-koeln.de
    Portierung zu Python 8/2017 Joschka Wirges, joschka.wirges@smail.th-koeln.de

    Teil des Praktikums zur Vorlesung Quellen- und Kanalcodierung

"""
import numpy as np, matplotlib.pyplot as plt
from ipywidgets import interact
import ipywidgets as widgets

def my_range(start, end, step):
    """
    range for a loop with decimal incrementation
    stepsize step in range of [start, end]
    """
    while start <= end:
        yield start
        start += step


def reconprk_old():
    """
    Version ohne ipython-Widget
    """
    T = 1  # T ist die Zeitkonstante der raised cosine function
    T_s = 0.5 * T  # für alpha=0.5 beträgt die Bandbreite zwar nur 1.5/T, es wird aber trotzdem mit T_s=0.5T (über-)abgetastet.
    alpha = 0.499  # Roll-Off Faktor
    amp_si = (1 + alpha) / T * T_s  # Amplitudenfaktor =2f_g*T_s
    t = np.linspace(-4, 4, 81)  # Betrachteter Ausschnitt des Signals
    x = np.sinc(t / T) * np.cos(np.pi * alpha * t / T) / (1 - 4. * alpha ** 2 * t ** 2 / T ** 2)
    plt.grid(True)
    plt.plot(t, x)
    plt.title('Raised cosine function')
    plt.ylabel(r'f(t)=sinc(\frac{t}{T}) cos($\pi \alpha$ \frac{t}{T})/(1-4$\alpha^2t^2/T^2$)')
    plt.xlabel('t/T')
    plt.text(2, 0.95, r'$\alpha$ = 0.5')
    plt.show()
    # Interpolationszwischenergebnis und Ausgabe
    # Beschriftung der Plots
    plt.subplot(2, 1, 1)
    plt.xlabel('t/T')
    plt.subplot(2, 1, 2)
    plt.xlabel('t/T')
    plt.tight_layout()  # Platz zwischen Subplots
    # Ausgabe
    inti = np.zeros(len(t))
    y = np.zeros((1, len(t)))
    for i in my_range(-4, 4, 0.5):
        y = x[(4 - np.int(i)) * 10] * amp_si * np.sinc((t - i) * (1 + alpha) / T)
        inti = inti + y
        plt.subplot(2, 1, 1)
        plt.plot(t, y, 'r--')
        plt.subplot(2, 1, 2)
        plt.plot(t, inti, 'b')
        plt.show()
    plt.subplot(2, 1, 2)
    plt.xlabel('t/T')
    plt.ylabel(r'f(t)=sinc(t/T) cos($\pi \alpha$ t/T)/(1-4 $\alpha$ ^2t^2/T^2)')
    plt.text(2, 1, r'$\alpha$ = 0.5')
    plt.text(2, 0.5, r'$T_s$ = 0.5T')
    plt.grid(True)
    return


def reconprk():
    """
    Demonstration der Rekonstruktion eines abgetasteten Signals
    als Zeitfunktion wird der RRC mit alpha=0.5 verwendet.
    Bandbreite des RRC Pulses: 2f_g=3/T
    Die Abtastzeit T_s  soll halb so groß wie T sein,
    d.h.BT_s=0.5$
    Der Rekonstruktionstiefpass ist ein idealer TP mit der Bandbreite
    B_id=2f_g=3/T
    Die nächste Interpolationsschritt wird durch Drücken der
    Returntaste durchgeführt. So baut sich das Signal erst langsam auf

    UwD 13.4.99
    """
    T = 1  # T ist die Zeitkonstante der raised cosine function
    T_s = 0.5 * T  # für alpha=0.5 beträgt die Bandbreite zwar nur 1.5/T, es wird aber trotzdem mit T_s=0.5T (über-)abgetastet.
    alpha = 0.499  # Roll-Off Faktor
    amp_si = (1 + alpha) / T * T_s  # Amplitudenfaktor =2f_g*T_s
    t = np.linspace(-4, 4, 81)  # Betrachteter Ausschnitt des Signals
    x = np.sinc(t / T) * np.cos(np.pi * alpha * t / T) / (1 - 4. * alpha ** 2 * t ** 2 / T ** 2)
    plt.grid(True)
    plt.plot(t, x)
    plt.title('Raised cosine function')
    plt.ylabel(r'f(t)=sinc($\frac{t}{T}$) cos($\pi \alpha \frac{t}{T}$)/(1-4$\alpha^2t^2/T^2$)')
    plt.xlabel(r'$\frac{t}{T}$')
    plt.text(2, 0.95, r'$\alpha$ = 0.5')
    plt.show()
    print('Steuern Sie über den Slider die Rekonstruktion des Raised Cosine von t=-4 bis t=4.')
    slider = widgets.FloatSlider(
        value=-4.0,
        min=-4 - 0,
        max=4.0,
        step=0.5,
        description='t',
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f',
    )

    def update(w):
        f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.autoscale(False, 'both', None)
        ax1.set_title('Rekonstruktion des RRC')
        ax1.set_xlim([-4, 4])
        ax1.set_ylim([-0.5, 1.25])
        inti = np.zeros(len(t))
        y = np.zeros(len(t))
        for i in my_range(-4, w, 0.5):
            y = 0
            y = x[int((4 - i) * 10)] * amp_si * np.sinc((t - (i)) * (1 + alpha) / T)
            inti = inti + y
            ax1.plot(t, y, '--r')
        ax2.plot(t, inti, 'b')
        plt.show()

    interact(update, w=slider);
    return

def u_pcm_old(a, n):
    """
    U_PCM  	Uniform PCM encoding of a sequence.
    [SQNR,A_QUAN,CODE]=U_PCM(A,N)
    a=input sequence.
    n=number of quantization levels (even).
    sqnr=output SQNR (in dB).
	a_quan=quantized output before encoding.
	code=the encoded output.
    """
    amax = max(abs(a))
    a_quan = a / amax
    b_quan = a_quan.copy()
    d = 2 / n
    q = d * np.linspace(0, (n - 1), n)
    q = q - ((n - 1) / 2) * d
    for i in range(1, n):
        a_quan[np.nonzero((q[i - 1] - d / 2 <= a_quan) & (a_quan <= q[i - 1] + d / 2))] = \
            q[i - 1] * np.ones(len(np.nonzero((q[i - 1] - d / 2 <= a_quan) & (a_quan <= q[i - 1] + d / 2))))
        b_quan[np.nonzero((a_quan == q[i - 1]))] = (i - 2) * np.ones((1, len(np.nonzero(a_quan == q[i - 1]))))
    a_quan = a_quan * amax
    nu = np.int(np.ceil(np.log2(n)))
    code = np.zeros((len(a), nu))
    for i in range(0, len(a)):
        for j in range(nu, 0, -1):
            if (np.fix(b_quan[i] / (2 ** j)) == 1):
                code[i, (nu - j)] = 1
                b_quan[i] = b_quan[i] - 2 ** j
    sqnr = 20 * np.log10(np.linalg.norm(a) / np.linalg.norm(a - a_quan))
    print(a_quan)
    return sqnr, a_quan, code


def u_pcm(a, n):
    """
    U_PCM  	Uniform PCM encoding of a sequence.
    [SQNR,A_QUAN,CODE]=U_PCM(A,N)
    a=input sequence.
    n=number of quantization levels (even).
    sqnr=output SQNR (in dB).
	a_quan=quantized output before encoding.
	code=the encoded output.
    """
    amax = max(a)
    a_quan = a / amax
    b_quan = a_quan.copy()
    d = 2 / n
    q = d * np.linspace(0, (n - 1), n)
    q = q - ((n - 1) / 2) * d
    for i in range(0, n):
        # print('i=',i,q[i] * np.ones(len(np.argwhere((q[i] - d / 2 <= a_quan) & (a_quan <= q[i] + d / 2))))) # PROB CORRECT
        # print('i=',i,a_quan[np.where((q[i] - d/2 <= a_quan) & (a_quan <= q[i] + d / 2))]) #CORRECT
        a_quan[np.where((q[i] - d / 2 <= a_quan) & (a_quan <= q[i] + d / 2))] = q[i] * np.ones(
            len(np.argwhere((q[i] - d / 2 <= a_quan) & (a_quan <= q[i] + d / 2))))
        b_quan[np.nonzero((a_quan == q[i]))] = (i - 1) * np.ones(len(np.nonzero(a_quan == q[i])))
    a_quan = a_quan * amax
    nu = np.int(np.ceil(np.log2(n)))
    code = np.zeros((len(a), nu))
    for i in range(0, len(a) - 1):
        for j in range(nu, 0, -1):
            if (np.fix(b_quan[i] / (2 ** j)) == 1):
                code[i, (nu - j)] = 1
                b_quan[i] = b_quan[i] - 2 ** j
    sqnr = 20 * np.log10(np.linalg.norm(a) / np.linalg.norm(a - a_quan))
    return sqnr, a_quan, code


def ip_04_09():
    """
    MATLAB script for Illustrative Problem 9, Chapter 4
    """
    t = np.linspace(0, 10, 1001)
    a = np.sin(t)
    [sqnr8, aquan8, code8] = u_pcm(a, 8)
    [sqnr16, aquan16, code16] = u_pcm(a, 16)
    [sqnr32, aquan32, code32] = u_pcm(a, 32)
    print('sqnr8=', sqnr8, 'dB')
    print('sqnr16=', sqnr16, 'dB')
    print('sqnr32=', sqnr32, 'dB')
    line1, = plt.plot(t, a, '-b', label="unquantized")
    line2, = plt.plot(t, aquan8, '--r', label="8 levels")
    line3, = plt.plot(t, aquan16, '-.g', label="16 levels")
    line4, = plt.plot(t, aquan32, 'm', label="32 levels")
    plt.plot(t, np.zeros(len(t)), 'k')
    plt.legend(handles=[line1, line2, line3, line4], loc=4)
    plt.show()
    return


def ip_04_10():
    """
    MATLAB script for Illustrative Problem 10, Chapter 4

    """
    a = np.random.randn(500)
    [sqnr, a_quan, code] = u_pcm(a, 64)
    print('sqnr=', sqnr)
    print('First five input values:', a[1:5])
    print('first five quantized values:', a_quan[1:5])
    print('first five codewords', code[1:5, :])
    return


def ip_04_10otherlevels():
    """
     finding the SQNR for a gaussian random variable N(0,1) for uniform
     quantiasation with 8, 16 and 32 stages
    """
    a = np.random.randn(500)
    [sqnr8, a_quan8, code8] = u_pcm(a, 8)
    [sqnr16, a_quan16, code16] = u_pcm(a, 16)
    [sqnr32, a_quan32, code32] = u_pcm(a, 32)
    print('sqnr8 =', sqnr8)
    print('sqnr16 =', sqnr16)
    print('sqnr32 =', sqnr32)
    return


def saege():
    """
    MATLAB script for Illustrative Problem 9, Chapter 4
    """
    t = np.linspace(0, 10, 1001)
    a = t - 5
    [sqnr8, aquan8, code8] = u_pcm(a, 8)
    [sqnr16, aquan16, code16] = u_pcm(a, 16)
    [sqnr32, aquan32, code32] = u_pcm(a, 32)
    print('sqnr8=', sqnr8)
    print('sqnr16=', sqnr16)
    print('sqnr32=', sqnr32)
    line1, = plt.plot(t, t - 5 - aquan8, '-', label="8 levels")
    line2, = plt.plot(t, t - 5 - aquan16, '-r', label="16 levels")
    line3, = plt.plot(t, t - 5 - aquan32, '-g', label="32 levels")
    plt.plot(t, np.zeros(len(t)), 'k')
    plt.legend(handles=[line1, line2, line3, ], loc=4)
    plt.show()
    return


def signum(x):
    """
    #SIGNUM	finds the signum of a vector.
    #Y=SIGNUM(X)
    #X=input vector
    """
    y = x.copy()
    y[np.nonzero(x > 0)] = np.ones(len(np.nonzero(x > 0)))
    y[np.nonzero(x < 0)] = -np.ones(len(np.nonzero(x < 0)))
    y[np.nonzero(x == 0)] = np.zeros(len(np.nonzero(x == 0)))
    return y


def mulaw(x, mu):
    """
    function [y,a]=mulaw(x,mu)
    %MULAW		mu-law nonlinearity for non-uniform PCM.
    %		Y=MULAW(X,MU)
    %		X=input vector.
    """
    a = abs(max(x))
    y = (np.log10(1 + mu * abs(x / a)) / np.log10(1 + mu)) * np.sign(x)  # *signum(x)
    return y, a


def invmulaw(y, mu):
    """
    function x=invmulaw(y,mu)
    INVMULAW		The inverse of mu-law nonlinearity
    X=INVMULAW(Y,MU)	Y=Normalized output of the mu-law nonlinearity
    """
    x = (((1 + mu) ** (abs(y)) - 1) / mu) * signum(y)
    return x


def mula_pcm(a, n, mu):
    """
    function [sqnr,a_quan,code]=mula_pcm(a,n,mu)
    %MULA_PCM 	mu-law PCM encoding of a sequence.
    %       	[SQNR,A_QUAN,CODE]=MULA_PCM(A,N,MU)
    %       	a=input sequence.
    %       	n=number of quantization levels (even).
    %       	sqnr=output SQNR (in dB).
    %		a_quan=quantized output before encoding.
    %		code=the encoded output.
    """
    [y, maximum] = mulaw(a, mu)
    [sqnr, y_q, code] = u_pcm(y, n)
    a_quan = invmulaw(y_q, mu)
    a_quan = maximum * a_quan
    sqnr = 20 * np.log10(np.linalg.norm(a) / np.linalg.norm(a - a_quan))
    return sqnr, a_quan, code


def ip_04_10mu_otherlevels():
    """
    % MATLAB script for Illustrative Problem 10, Chapter 4
    % finding the SQNR for a random gaussian variable N(0,1) with non-uniform
    % quantisation (u-law) for 8,16, and 32 stages
    """
    au = np.random.randn(500)
    [sqnru8, au_quan, codeu] = mula_pcm(au, 8, 255)
    [sqnru16, au_quan, codeu] = mula_pcm(au, 16, 255)
    [sqnru32, au_quan, codeu] = mula_pcm(au, 32, 255)
    print('sqnru8=', sqnru8)
    print('sqnru16=', sqnru16)
    print('sqnru32=', sqnru32)
    return


def ip_04_10modified():
    """
    % MATLAB
    script
    for Illustrative Problem 10, Chapter 4
    """
    ag = 0.1 * np.random.randn(50000)
    ag = np.append(ag, [-1, 1])
    count, bins, ignored = plt.hist(ag, 30)
    plt.plot(bins, np.ones_like(bins))
    plt.show()
    mu = 255
    n = 32
    [sqnrg, ag_quan, codeg] = u_pcm(ag, n)
    [sqnru, au_quan, codeu] = mula_pcm(ag, n, mu)
    print('gleichsqnr=', sqnrg)
    print('ungleichsqnr=', sqnru)
    return


def saegemu():
    """
    % MATLAB script for Illustrative Problem 9, Chapter 4

    """
    t = np.linspace(0, 10, 1001)
    a = np.linspace(0, 10, 1001) - 5
    [sqnr8, aquan8, code8] = mula_pcm(a, 8, 255)
    [sqnr16, aquan16, code16] = mula_pcm(a, 16, 255)
    [sqnr32, aquan32, code32] = mula_pcm(a, 32, 255)
    print('sqnr8=', sqnr8)
    print('sqnr16=', sqnr16)
    print('sqnr32=', sqnr32)
    plt.clf()
    line1, = plt.plot(t, t - 5 - aquan8, 'b-', label="8 levels")
    line2, = plt.plot(t, t - 5 - aquan16, 'g-', label="16 levels")
    line3, = plt.plot(t, t - 5 - aquan32, '-r', label="32 levels")
    plt.plot(t, np.zeros(len(t)), 'k')
    plt.legend(handles=[line1, line2, line3], loc=4)
    plt.show()
    return

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

"""
def huffman(p):
    #HUFFMAN 	Huffman code generator.
	#[h,l]=huffman(p), Huffman code generator
	#returns h the Huffman code matrix, and l the
    #average codeword length for a source with
	#probability vector p.
    
   # if (len(np.nonzero(p < 0)) != 0):
   #     print('Error:Not a prob. vector, negative component(s)')
   # if ((abs(sum(p) - 1)) > 10e-10):
   #     error('Error: Not a prob. vector, components do not add up to 1')
    n = len(p)
    q = p.copy()
    m = np.zeros((n - 1, n))
    for i in range(1, n):
        l = np.argsort(q, axis=-1)
        q = np.sort(q,axis=0)
        l+=1
        m[i-1, :] = np.append(l[0:n - i + 1], np.zeros((1, i - 1)))
        q = np.append(q[0] + q[1], q[2:n])
        q=np.append(q,1)
    #for i=1:n-1 #c(i,:)=blanks(n*n);
    c = np.empty((n-1,n*n), dtype='str')
    c[:] = ' '
    c[n - 2, n-1] = 0
    c[n - 2, 2 * (n-1)] = 1
    print(m)
    for i in range(2, n):
        #print((np.arange(np.ravel(n*(1+np.argwhere(m[(n - i ), :] ==1))-(n-2),),1+ n  * (1+np.argwhere(m[n - i, :] == 1)))-1).shape)
        a=np.ravel(n*(1+np.argwhere(m[(n - i ), :] ==1))-(n-2))
        print(a)
        b=np.ravel(1+ n  * (1+np.argwhere(m[n - i, :] == 1)))
        print(b)
        #c[n - i-1, 0:(n-2)] = c[n - i, a:b]
        #print(c[n - i, :(n-2)])
        #print(n - i)
        print(np.arange(np.ravel(n*(1+np.argwhere(m[(n - i ), :] ==1))-(n-2)),1+ n  * (1+np.argwhere(m[n - i, :] == 1)))-1)
        #print('i=',i,'n=',n)
        #print('vor minus',np.ravel(n*(1+np.argwhere(m[(n - i ), :] ==1)))) #CORRECT
        #print('nach minus',np.arange((n - 2),1+ n  * (1+np.argwhere(m[n - i, :] == 1)))) #CORRECT
        #print(c[n - i, n * (np.argwhere(m[n - i, : ] == 1)) - (n - 2): n * (np.argwhere(m[ n - i, : ] == 1))])
        #print('gesamt',np.ravel(n*(1+np.argwhere(m[(n - i ), :] ==1)))-np.arange((n - 2),1+ n  * (1+np.argwhere(m[n - i, :] == 1)))) # 7 6 5 4 3 2 1 0 statt 7 8 9 10
        print('gesamt',np.arange(np.ravel(n*(1+np.argwhere(m[(n - i ), :] ==1))-(n-2),),1+ n  * (1+np.argwhere(m[n - i, :] == 1)))) #CORRECT
        """
"""
        c[n - i-1, n-1] = 0
        c[n - i-1, n-1 + 1:2 * (n - 2)] = c[n - i-1, 1:n - 2]
        c[n - i-1, 2 * (n-1)] = 0
        for j in range(1, i - 1):
            for k in range(1, n):
                c_index2 = k * np.nonzero(m[n - i + 1, :] == j + 1)
            c[n - i, (j + 1) * n + 1: (j + 2) * n] = c[
                n - i + 1, n * (np.nonzero(m[n - i + 1, :] == j + 1) - 1) + c_index2]
    for i in range(1, n):
        h[i, 1:n] = c[1, n * (np.nonzero(m[1, :] == i)) - 1]
        for i in range(1, n):
            j_max = np.nonzero(m[1, :] == i)
            for j in range(1, j_max):
                h_index2 = j * n
            h[i, 1: n] = c[1, n * (np.nonzero(m[1, :] == i) - 1) + h_index2]
            l1[i] = len(np.nonzero(abs(h[i, :]) != 32))
    l = sum(p * l1)
    
    return h, l
    """


# https://gist.github.com/mreid/fdf6353ec39d050e972b
def huffman(p):
    '''Return a Huffman code for an ensemble with distribution p.'''
    # assert(sum(p.values()) == 1.0) # Ensure probabilities sum to 1

    # Base case of only two symbols, assign 0 or 1 arbitrarily
    if (len(p) == 2):
        return dict(zip(p.keys(), ['0', '1']))
    # if input p is a ndarray convert it to a dict
    if (type(p) == type(np.array([1]))):
        return huffman(dict(enumerate(p.flatten(), 1)))

    # Create a new distribution by merging lowest prob. pair
    codelength = 0
    p_prime = p.copy()
    a1, a2 = lowest_prob_pair(p)
    p1, p2 = p_prime.pop(a1), p_prime.pop(a2)
    p_prime[a1 + a2] = p1 + p2
    # Recurse and construct code on new distribution
    c = huffman(p_prime)
    ca1a2 = c.pop(a1 + a2)
    c[a1], c[a2] = ca1a2 + '0', ca1a2 + '1'
    return c


# https://gist.github.com/mreid/fdf6353ec39d050e972b
def lowest_prob_pair(p):
    '''Return pair of symbols from distribution p with lowest probabilities.'''
    # assert(len(p) >= 2)
    # Ensure there are at least 2 symbols in the dist.
    sorted_p = sorted(p.items(), key=lambda p: p[1])
    return sorted_p[0][0], sorted_p[1][0]

   # l = sum(len(i) for i in np.array(list(c.values()))) / len(np.array(list(c.values())))
def codelength(c,p):
    if (type(p) == type(np.array([1]))):
        p = dict(enumerate(p.flatten(), 1))
    l = 0.0
    sorted_p = sorted(p.values()) #sorted_p is a list with ascending p-values
    codewords = list(c.values())
    for i in range(0, len(p)):
        l+= sorted_p.pop(len(sorted_p)-1)*len(codewords.pop(0))
    return l

def p_ger_lan():
    """
    Aufrittswahrscheinlichkeiten für deutsche Sprache
    UwD  99 / 4 / 29
    returns a dictionary pgl with keys 'a'-'z' and '-'(space) and their corresponding p-values
    """
    pgl = {'a': 0.0549, 'b': 0.0138, 'c': 0.0255, 'd': 0.0546, 'e': 0.1440, 'f': 0.0078, 'g': 0.0236, 'h': 0.0361,
           'i': 0.0628, 'j': 0.0028, 'k': 0.0071, 'l': 0.0345, 'm': 0.0172,
           'n': 0.0865, 'o': 0.0211, 'p': 0.0067, 'q': 0.0005, 'r': 0.0622, 's': 0.0646, 't': 0.0536, 'u': 0.0422,
           'v': 0.0079, 'w': 0.0113, 'x': 0.0008, 'y': 0.0001, 'z': 0.0092,
           '-': 0.1486}
    return pgl