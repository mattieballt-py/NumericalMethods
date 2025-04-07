# libraries:
import numpy as np
"""
Jacobian Numerical
"""
def numJacobian(F, x, h=1e-8): # F is a function f(x)
    n = len(x)
    J = np.zeros((n, n))
    Fx = F(x)
    for i in range(n):
        x_h = x.copy()
        x_h[i] += h
        J[:, i] = (F(x_h) - Fx) / h
    return J

"""
Golden Search Method to find max/min of a function 
"""
# want to find maximum of a function using only one func eval per iteration

# the function we want to search for a max in
def df(x): # will find min/max of df/dx which corresponds to where f''x = 0 => inflection point
    return -2*x*np.exp((-x**2))

#not recursively
def goldsrch(xl,xu,tol):
    R = (np.sqrt(5) - 1) / 2  # Golden ratio
    e = 100 # set high initial error to start off loop
    i = 0
    # just to start it off:
    x1 = xl + (xu-xl)*R # compute x1
    x2 = xu - (xu-xl)*R # compute x2 
    fx1 = df(x1)
    fx2 = df(x2)
    while abs(e) > tol:
        xn = (xl + xu)/2 # before iteration best guess
        print("xl:",xl,"xu:",xu)
        print("df(xl):",df(xl),"df(xu):",df(xu))
        if fx1 < fx2:
            i += 1
            xl = x2
            x2 = x1
            #x2' = x1
            x1 = xl + (xu-xl)*R # compute x1
            fx2 = fx1 # and no need to recalculate fx2  
            fx1 = df(x1)            
        else: # fx2 < fx1:
            i += 1
            xu = x1
            x1 = x2
            #x1' = x2
            fx1 = fx2 # and no need to recalculate fx1
            x2 = xu - (xu-xl)*R # compute x2 
            fx2 = df(x2)
        e = (((xl+xu)/2) - xn)/((xl+xu)/2) # compute width => error (xn+1 - xn)/xn
        #print("new best guess at (n + 1), i: ",i,"is: ", (xl+xu)/2)
    print("best guess under tol where derivative=0:",(xl+xu)/2)
    return (xl+xu)/2
        
# or

# recursively:
"""
Needs fixing
def goldsearch(xl,xu,i):
    xl = xl
    xu = xu
    fxl = df(xl) 
    fxu = df(xu)
    if fxl * fxu > 0:
        return print("choose a new sub interval")
    R = (np.sqrt(5) - 1) / 2  # Golden ratio
    x1 = xl + R * (xu - xl)
    x2 = xu - R * (xu - xl)

    if abs(xu - xl) < 0.01:  # If the interval is small enough, return result
        return (xl + xu) / 2
    else:
        x1 = xl + d
        x2 = xu - d
        fx1 = df(x1)
        fx2 = df(x2)
        if fx2 < fx1:
            i += 1
            print(i)
            return goldsearch(xl, x1,i)  # Shrink the interval towards x2
        else:
            i += 1
            print(i)
            return goldsearch(x2, xu,i)  # Shrink the interval towards x1
"""

#print(goldsrch(0,1,0.01))

"""
Newton Raphson Method (same as for RootFinding but f'=0 not f=0)
"""
def g(x):
    return x**2 + 5

# could numerically find g' or:
def gder(x):
    return 2*x

def gder2nd(x):
    return 2

# if linear system of eqns, could solve with x = A-1 linalgsolve with b
# for non linear syst eqns like x^2 + y^2 = 2 etc, use:
def NewtonRaph(x0,tol):
    xn = x0
    err = 10 # init error
    while err > tol:
        xn_1 = xn - gder(xn)/gder2nd(xn)
        err = abs((xn_1 - xn)/xn_1)
        xn = xn_1
    return xn

"""
Steepest Ascent Method (Optimal so choosing distance to move strategically)
"""
