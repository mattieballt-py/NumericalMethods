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

