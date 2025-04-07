import numpy as np 
import matplotlib.pyplot as plt 


"""
Integration with Guass Quadrature 
"""
# when you (know the function but no known points)

def g(t):
    A = 0
    B = 1
    x = 0.5*(A * (1-t) + B*(t + 1)) # A is lower limit, B is upper
    return x * np.exp(x**2)

def GuassQuad(g,n): # n is degree of precision
    nodes = [] # = tj
    weights = [] # = wj
    I = 0 # just to init

    if n == 7:
        # define the nodes, tj = nodes
        nodes.append(((3/7) + (2/7)*(6/5)**0.5)**0.5)
        nodes.append(-1 * ((3/7) + (2/7)*(6/5)**0.5)**0.5)
        nodes.append(((3/7) - (2/7)*(6/5)**0.5)**0.5)
        nodes.append(-1 * ((3/7) - (2/7)*(6/5)**0.5)**0.5)

        # define the weights:
        weights.append((18 - 30**0.5)/36)
        weights.append((18 - 30**0.5)/36)
        weights.append((18 + 30**0.5)/36)
        weights.append((18 + 30**0.5)/36)

        # run the integration
        for j in range(0,4):
            I += weights[j]*g(nodes[j])
        I *= 0.5
        return I
    else:
        return

print("Guass solution:", GuassQuad(g,7))


"""
Integration of f(x) for a spacing of h between A and B with known points
"""
A = 0
B = 1
N = 10000 # total number of known points
h = (B-A)/N

x_known = np.linspace(A,B+h,N)

def f(x_known):
    return x_known * np.exp(x_known**2)

y_known = f(x_known)

"""
Simpsons 1/3 method
"""

def Simp13(y_known, N, h):
    I_even = 0 # init
    I_odd = 0 
    I = y_known[0] + y_known[-1]
    for n in range(0,N):
        if n%2 == 0: # n is even
            I_even += y_known[n]
        else:
            I_odd += y_known[n]
    I = h/3 * (I + 4*I_odd + 2*I_even)
    return I

print("simp", Simp13(y_known,N,h))

"""
Adaptive Simpsons 1/3 method
"""
# if you need to get I to be more accurate, reduce h:

def Simp13adaptive(y_known, N, h,tol):

    I_even = 0 # init
    I_odd = 0 
    I = y_known[0] + y_known[-1]
    for n in range(0,N):
        if n%2 == 0: # n is even
            I_even += y_known[n]
        else:
            I_odd += y_known[n]
    I = h/3 * (I + 4*I_odd + 2*I_even)

    act = GuassQuad(g,7) # as its the same function
    err = (act - I)
    print("err",err)
    if abs(err) > tol:
        N += 100 # increase number of points
        x_known = np.linspace(A,B,N) # recalc x
        y_known = f(x_known) # refind more y values
        h = x_known[1] - x_known[0]
        print("h: ",h)
        return Simp13adaptive(y_known,N,h,tol)
    return I

print("simpadapt", Simp13adaptive(y_known,N,h,0.00005), "vs guass quad: ", GuassQuad(g,7))


"""
Trapezium Rule
"""

def TrapInt(A,B,h,N,y_known):
    I = 0 # init
    for n in range(1,N):
        I += y_known[n] 
    I = h*(I + f(A)/2 + f(B)/2)
    return I

print("Trap result: ",TrapInt(A,B,h,N,y_known))