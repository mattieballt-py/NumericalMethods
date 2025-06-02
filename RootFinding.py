# basically algorithms to find where f(x) = 0
import numpy as np
import matplotlib.pyplot as plt 

"""
Bisection Method (1 of 2 Bracketing Methods)
"""
R = 2 # radius of a cone or something, an oval
def f(x):
    return (3.14 * R * x **2) - ((3.14/3)*x**3) - 25


# bisection method recursively
def bisecrecurs(a,b,tol,i):
    if f(a)*f(b)>0: # check valid interval
        print("choose new interval")
        return None

    n = (a + b)/2 #new point, bang in the middle of a and b
    print('midpoint', n)

    if f(a)*f(n)<0: # function changes sign between a and n
        # best guess before was n, now n + 1 = (a + n)/2
        # err = (new - old) /new
        errnew = (((a + n)/2) - n)/((a + n)/2)
        if abs(errnew) < tol:
            return ((a + n)/2) # if interval less than min, the root is at centre of it ish
        else:
            i += 1
            print(i)
        return bisecrecurs(a,n,tol,i)
  
    else:  # function changes sign between n and b
        errnew = (((n + b)/2) - n)/((n + b)/2)
        if abs(errnew) < tol:
            return ((n + b)/2) # if interval less than min, the root is at centre of it ish
        i += 1
        print(i)
        return bisecrecurs(n,b,tol,i)


print(bisecrecurs(0,2*R,0.001,1))
# a, b is the interval to evaluate inside


"""
False Position Method (2 of 2 Bracketing Methods)
"""
# also diff recursive in playing 
def g(x):
    return x**2 + (x-2)**3 - 4

def falsepos(a,b,tol,i):
    if f(a)*f(b)>0: # check valid interval
        print("choose new interval")
        return None
    
    n = b - (f(b)*(a-b))/(f(a)-f(b)) # new point n
    print('midpoint', n, "fxn", f(n))
    if f(a)*f(n)<0:
        print("in range a to n")
        # then sign changes betweem a and n
        
        xnew = n - (f(n)*(a-n))/(f(a)-f(n)) # new point n
        xold = n
        err = (xnew-xold)/xnew
        print("xnew",xnew,"xold",xold)
        print("err",err)
        if abs(err) < tol:
            return (a + n)/2
        else:
            i += 1
            print(i)
            return falsepos(a,n,tol,i)
    else:
        print("in range n to b")
        # then sign changes betweem n and b
        xnew = b - (f(b)*(n-b))/(f(n)-f(b)) # new point n
        xold = n
        err = (xnew-xold)/xnew
        print("xnew",xnew,"xold",xold)
        print("err",err)
        if abs(err) < tol:
            return (b + n)/2
        i += 1
        print(i)
        return falsepos(n,b,tol,i)

#print("from here on, doing falsepos funct")
#print(falsepos(0,2*R,0.000001,2))


"""
Newton Raphson Method (1 of 2 open methods)
"""
# useful for finding e.g a/b efficiently
def g(x):
    return x**2 + 5

# could numerically find g' or:
def gder(x):
    return 2*x

# if linear system of eqns, could solve with x = A-1 linalgsolve with b
# for non linear syst eqns like x^2 + y^2 = 2 etc, use:
def NewtonRaph(x0,tol):
    xn = x0
    err = 100 # init error
    while err > tol:
        xn_1 = xn - g(xn)/gder(xn)
        err = abs((xn_1 - xn)/xn_1)
        xn = xn_1
    return xn

"""
Secant Method (2 of 2 open methods)
"""

# 2 of 2 open method for when can't f'(x)
def secant(x0,x1,tol): # two initial guesses
    xn = x0 # init
    xn_1 = x1
    err = 100 # just to init
    while err > tol:
        xn_2 = xn_1 - (f(xn_1)*(xn - xn_1))/(f(xn)-f(xn_1))
        err = abs((xn_2 - xn_1)/xn_2)
        xn,xn_1= xn_1,xn_2
    return xn_2

"""
Bisection method with langrangian interpolated points
"""

with open('/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/txt_Data/Temperatures.txt', 'r') as d:
    T = d.readlines()

# Initialize lists to store cleaned float values
Temps_known = []
for item in T:
    term = item.strip()
    try:
        clean_term = float(term)
        Temps_known.append(round(clean_term))
    except ValueError:
        continue

Temp_known = np.array(Temps_known)
print("Temps", Temp_known)

# Corresponding time array
t_known = np.arange(0, 6, 0.5)
print("time", t_known)

# Interpolation using Lagrangian method
def Langrangian(x_known, y_known, xp):
    n = len(x_known)
    P = 0
    for i in range(n):
        Li = 1
        for j in range(n):
            if i != j:
                Li *= (xp - x_known[j]) / (x_known[i] - x_known[j])
        P += y_known[i] * Li
    return P

def InterpPoints(x_known, y_known, InterpMethod, xp_points):
    yp_points = xp_points.copy()
    for i in range(len(xp_points)):
        yp_points[i] = InterpMethod(x_known, y_known, xp_points[i])
    return yp_points

# Bracketing interval between t = 3 and t = 3.5
t_interp = np.arange(3, 3.5, 0.005)
Temp_interp = InterpPoints(t_known, Temp_known, Langrangian, t_interp)

# Plotting interpolated section
plt.figure(figsize=(8, 4))
plt.scatter(t_interp, Temp_interp, color='black')
plt.title("Interpolated points for root finding")
plt.xlabel("t")
plt.ylabel("T(t)")
plt.grid(True)
plt.tight_layout()
plt.show()

# Bisection method for langrangian interpolated points
def bisecrecurs(a, b, tol, i):
    fa = Langrangian(t_known, Temp_known, a)
    fb = Langrangian(t_known, Temp_known, b)

    if fa * fb > 0:
        print("No sign change — choose a new interval.")
        return None

    mid = (a + b) / 2
    fmid = Langrangian(t_known, Temp_known, mid)
    print(f"Iteration {i}, midpoint t = {mid:.5f}, T = {fmid:.5f}")

    if abs(fmid) < tol:
        print(f"Root found at t = {mid:.5f}")
        return mid

    if fa * fmid < 0:
        return bisecrecurs(a, mid, tol, i + 1)
    else:
        return bisecrecurs(mid, b, tol, i + 1)

# Call fixed bisection function over time interval [3, 3.5]
root_time = bisecrecurs(3.0, 3.5, 0.01, 0)
print(f"Estimated root at time t = {root_time:.4f} s")