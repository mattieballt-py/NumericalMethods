# cid: 02413400

#importing needed libraries:
import numpy as np
import matplotlib.pyplot as plt 
import math as mt 


"""
Task A: 2D Simpson’s Rule Integration over triangular domain A
"""

# Define function: f(x, y) = sin(xy) * cos(x + y)
def f_2d(x, y):
    return np.sin(x * y) * np.cos(x + y)

# Grid spacings
dx = 0.05
dy = 0.04

# x from 0 to π
x = np.arange(0, mt.pi + dx, dx)

# Ensure odd number of x‐points
if len(x) % 2 == 0:
    x = x[:-1]

# Simpson’s Rule Integration over triangular domain
def simpsons_over_triangle(dx, dy):
    I = 0  # total integral

    for i in range(len(x)):
        xi = x[i]
        ymax = mt.pi - xi
        y = np.arange(0, ymax + dy, dy)

        # Ensure odd number of y‐points
        if len(y) % 2 == 0:
            y = y[:-1]

        #Simpson’s rule in y for fixed x = xi
        Iy = 0
        Ny = len(y)

        for j in range(Ny):
            yj = y[j]                   # current y‐value
            fij = f_2d(xi, yj)          # f(xi, yj)

            # determine Simpson coefficient in y‐direction
            if j == 0 or j == Ny - 1:
                coef_y = 1
            elif j % 2 == 1:
                coef_y = 4
            else:
                coef_y = 2

            Iy += coef_y * fij

        Iy *= dy / 3                     # finish Simpson’s 1/3 in y

        # Simpson’s coefficient in x‐direction
        if i == 0 or i == len(x) - 1:
            coef_x = 1
        elif i % 2 == 1:
            coef_x = 4
        else:
            coef_x = 2

        I += coef_x * Iy

    I *= dx / 3                           # finish Simpson’s 1/3 in x
    return I

# Compute and print result
I_result = simpsons_over_triangle(dx, dy)
print("Integral of sin(xy)*cos(x+y) over triangular region A:", I_result)


"""
Task B: Lagrnagian interp of two functions
"""
# open set 1 values
with open('/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/txt_Data/Sety1.txt', 'r') as d:
    T1 = d.readlines()
# set 2 values
with open('/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/txt_Data/Sety2.txt', 'r') as d:
    T2 = d.readlines()

# Initialize lists to store cleaned float values
T1_known = []
for item in T1:
    term = item.strip()
    try:
        clean_term = float(term)
        T1_known.append(round(clean_term))
    except ValueError:
        continue

T1n= np.array(T1_known)
print("T1", T1n)

# Initialize lists to store cleaned float values for set 2
T2_known = []
for item in T2:
    term = item.strip()
    try:
        clean_term = float(term)
        T2_known.append(round(clean_term))
    except ValueError:
        continue

T2n= np.array(T2_known)
print("T2", T2n)


# Corresponding xn's
xn1 = np.arange(0, 10+0.5, 0.5)
xn2 = np.arange(0,10+0.6,0.6)
print("xn1",xn1)
print("xn2",xn2)
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
x_interp = np.arange(2, 8, 0.05)

T1_interp = InterpPoints(xn1,T1n, Langrangian, x_interp)
T2_interp = InterpPoints(xn2,T2n, Langrangian, x_interp)

Z = T1_interp+ T2_interp # interp first then adding so x's align

# Plotting interpolated section
plt.figure(figsize=(8, 4))
plt.plot(x_interp, Z, color='black')
plt.scatter(xn1,T1n)
plt.scatter(xn2,T2n)
plt.title("Interpolated points Z")
plt.xlabel("t")
plt.ylabel("T(t)")
plt.grid(True)
plt.tight_layout()
plt.show()

"""
Task C: solve forward Euler first order ODE in dy/dx
"""

def func(x, Y):
    y0, y1 = Y
    dy0 = y1
    dy1 = (2/x)*y1 - 2*y0/(x**2)
    return np.array([dy0, dy1])

# Forward Euler Method for systems
def FwEuler(x0, h, x_final,a):
    # Initial conditions
    y0_init = 10
    y1_init = a

    Y0 = np.array([y0_init, y1_init])
    x_values = [x0]
    y_values = [Y0]
    x = x0
    Y = Y0.copy()

    while x < x_final:
        Y = Y + h * func(x, Y)
        x = x + h

        x_values.append(x)
        y_values.append(Y.copy())

    return np.array(x_values), np.array(y_values)


# x domain
x0 = 4
x_final = 20
h = 0.1 # dx

# Solve using FWD Euler method
a = [0,1,2,4,8]
plt.figure(figsize=(10, 4))
for i in range(0,len(a)):
    x_vals, y_vals = FwEuler(x0, h, x_final,a[i])
    # Plot y(x) for each a 
    plt.plot(x_vals, y_vals[:, 0], label=f'a = {a[i]}')  # plot y0
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Numerical Solution y(x) using Forward Euler')
    plt.grid()
    plt.legend()
    plt.tight_layout()
plt.show()


"""
Task D: bisection method to find intersection of two curves
"""
# two functions 
# for circle y1 = sqrt(r^2 - (x-4)^2)
# for line y2 = x - 4
# so want to solve Z = y1 - y2 and find where this =0 as this is when y1 = y2 bosh then sub in to find (x,y) of intersection


def f1(x):
    return (4**2 - (x-4)**2)**0.5
def f2(x):
    return (x - 4)

def f(x):
    return f1(x) - f2(x)

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


y = bisecrecurs(4,8,0.001,1) # choosing bracket from center of circle to the right, so dont get confused with root on lhs
x = y + 4 # as using the line one
# a, b is the interval to evaluate inside
