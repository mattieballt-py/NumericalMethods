# cid = 02413400
# the one on icecreams

#importing needed libraries:
import numpy as np
import matplotlib.pyplot as plt 
import math as mt 

"""
Task A: Finding root using bisection method
"""
# two functions 
# for circle y1 = sqrt(r^2 - (x-4)^2)
# for line y2 = x - 4
# so want to solve Z = y1 - y2 and find where this =0 as this is when y1 = y2 bosh then sub in to find (x,y) of intersection


# Geometry setup
R = 0.030  # Radius in meters
m = np.tan(np.pi / 3)  # Slope = tan(60°) = √3

# Function for the scoop: gives x in terms of z (your y)
def f1(y):
    return np.sqrt(2 * R * y - y**2)

# Function for the wafer: x = (y - R)/m, rearranged from z = -m*x + R
def f2(y):
    return (y - R) / m

# Difference of the two functions — root occurs at intersection
def f(y):
    return f1(y) - f2(y)

# Recursive Bisection Method
def bisecrecurs(a, b, tol, i):
    if f(a) * f(b) > 0:
        print("Choose new interval")
        return None

    n = (a + b) / 2
    print('Midpoint:', n)

    if f(a) * f(n) < 0:
        errnew = abs(((a + n) / 2 - n) / ((a + n) / 2))
        if errnew < tol:
            return (a + n) / 2
        return bisecrecurs(a, n, tol, i + 1)
    else:
        errnew = abs(((n + b) / 2 - n) / ((n + b) / 2))
        if errnew < tol:
            return (n + b) / 2
        return bisecrecurs(n, b, tol, i + 1)


# Call bisection to find z of intersection
zA = bisecrecurs(0, 2 * R, 0.001, 1)
xA = f1(zA)  # or f2(zA), both give x

print("intersection at xA,zA",xA,zA)