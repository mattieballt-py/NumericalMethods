import numpy as np
import matplotlib.pyplot as plt  

"""
Lagrangian Interpolation
"""
# want to find f(xp) at a known value xp from a set of x_knowns

def f(x):
    return x**2

x_known = np.arange(0,10,1)
y_known = f(x_known)

def Langrangian(x_known, y_known, xp):
    n = len(x_known) # the degree of P polynomial
    print("n",n)
    P = 0 # init this is the Langrangian Polynomial
    Li = 1
    for i in range(0,n):
        print("i",i)
        for j in range(0,n):
            if i != j:
                print("j",j)
                Li *= (xp - x_known[j])/(x_known[i]-x_known[j])
        P += y_known[i] * Li
        Li = 1
    return P

print(Langrangian(x_known,y_known,0.4))