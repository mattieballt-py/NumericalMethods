import numpy as np
import matplotlib.pyplot as plt  

"""
Lagrangian Interpolation
"""
# want to find f(xp) at a known value xp from a set of x_knowns

def f(x):
    return x**2 + np.cos(x) + np.exp(-3*np.sin(x))

x_known = np.linspace(1,2,3)
y_known = f(x_known)

def Langrangian(x_known, y_known, xp):
    n = len(x_known) # the degree of P polynomial
    #print("n",n)
    P = 0 # init this is the Langrangian Polynomial
    Li = 1
    for i in range(0,n):
        #print("i",i)
        for j in range(0,n):
            if i != j:
                #print("j",j)
                Li *= (xp - x_known[j])/(x_known[i]-x_known[j])
        P += y_known[i] * Li
        Li = 1
    return P

def InterpPoints(x_known,y_known, InterpMethod, xp_points):
    yp_points = xp_points.copy()
    for i in range(0,len(xp_points)):
        yp_points[i] = InterpMethod(x_known,y_known,xp_points[i])
    return yp_points

xp = np.arange(0,3,0.05)
yp = InterpPoints(x_known,y_known,Langrangian,xp)

# for the plot
# same intervals, different degree of polynomial interp
xn_2 = np.linspace(1,2,2) # n = 2
yn_2 = f(xn_2)
xn_3 = np.linspace(1,2,3) # n = 3
yn_3 = f(xn_3)
xn_4 = np.linspace(1,2,4) # n = 4
yn_4 = f(xn_4)
xn_7 = np.linspace(1,2,7) # n = 7
yn_7 = f(xn_7)
xn_12 = np.linspace(1,2,12) # n = 7
yn_12 = f(xn_12)


yp_3 = InterpPoints(xn_3,yn_3,Langrangian,xp)
yp_4 = InterpPoints(xn_4,yn_4,Langrangian,xp)
yp_7 = InterpPoints(xn_7,yn_7,Langrangian,xp)
yp_12 = InterpPoints(xn_12,yn_12,Langrangian,xp)
plt.scatter(xn_2,yn_2, color="red", label="known")
plt.scatter(xn_3,yn_3, color="orange", label="known")
plt.scatter(xn_7,yn_7, color="orange", label="known")
plt.scatter(xp,yp_3, color="blue", label = "n = 3")
plt.scatter(xp,yp_4, color="turquoise", label = "n = 4")
plt.scatter(xp,yp_7, color="deepskyblue", label = "n = 7")
plt.scatter(xp,yp_12, color="indigo", label = "n = 12")
plt.legend()
plt.grid(True)
#plt.show()

"""
2D Lagrangian
"""
import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# 1D Lagrange Basis Polynomial
# ------------------------

def lagrange_basis(x, xi, i):
    """
    Computes the i-th Lagrange basis polynomial at x for given xi.

    Args:
        x (float): Point to evaluate the basis polynomial at
        xi (list): List of known x coordinates
        i (int): Index of the basis polynomial to evaluate

    Returns:
        float: Value of the Lagrange basis polynomial L_i(x)
    """
    L = 1.0
    for j in range(len(xi)):
        if j != i:
            L *= (x - xi[j]) / (xi[i] - xi[j])
    return L

# ------------------------
# 2D Bilinear Interpolation Using Lagrange
# ------------------------

def bilinear_lagrange_interpolation(x, y, x_vals, y_vals, f_vals):
    """
    Computes the bilinear interpolation using Lagrange polynomials in 2D.

    Args:
        x (float): x-coordinate of interpolation point
        y (float): y-coordinate of interpolation point
        x_vals (list): List of two x values [x0, x1]
        y_vals (list): List of two y values [y0, y1]
        f_vals (2x2 array): Function values at the 4 corners of the rectangle:
                            [[f(x0,y0), f(x1,y0)],
                             [f(x0,y1), f(x1,y1)]]

    Returns:
        float: Interpolated value at (x, y)
    """
    interp = 0.0
    for i in range(2):
        for j in range(2):
            Li = lagrange_basis(x, x_vals, i)
            Lj = lagrange_basis(y, y_vals, j)
            interp += f_vals[j][i] * Li * Lj
    return interp

# ------------------------
# Example
# ------------------------

# Define grid
x_vals = [0.0, 1.0]
y_vals = [0.0, 1.0]

# Define function z = f(x, y) at the grid points
# Example: f(x, y) = x + y at (x0,y0), (x1,y0), (x0,y1), (x1,y1)
f_vals = [
    [0.0, 1.0],  # f(x0,y0), f(x1,y0)
    [1.0, 2.0]   # f(x0,y1), f(x1,y1)
]

# Interpolation point
x_interp = 0.3
y_interp = 0.6

# Interpolate
z_interp = bilinear_lagrange_interpolation(x_interp, y_interp, x_vals, y_vals, f_vals)
print(f"Interpolated value at ({x_interp}, {y_interp}) = {z_interp:.4f}")

# ------------------------
# Visualization (Optional)
# ------------------------

def f_actual(x, y):
    return x + y  # Example function used to generate f_vals

# Grid for plotting
X = np.linspace(0, 1, 20)
Y = np.linspace(0, 1, 20)
Z = np.zeros((len(Y), len(X)))

for i, yi in enumerate(Y):
    for j, xi in enumerate(X):
        Z[i, j] = bilinear_lagrange_interpolation(xi, yi, x_vals, y_vals, f_vals)

# Plot
plt.figure(figsize=(6, 5))
plt.contourf(X, Y, Z, levels=20, cmap='viridis')
plt.colorbar(label='Interpolated z')
plt.scatter([x_interp], [y_interp], color='red', label='Interpolated point')
plt.title('2D Bilinear Interpolation using Lagrange Polynomials')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()

### 3D Lagrangian in Playing

"""
Barycentric Interpolation
"""






"""
2D Interpolation
"""
def f(coord):
    x,y = coord
    return x*y

r = np.array([[0, 0], [0, 4], [5, 0]])
f = np.array([f(pt) for pt in r]) # creates array cycling through tuples of points in r
p = np.array([1, 1])
print(r,f,p)

def inverse_dist_triangle(r,func,p):
    return




# this 2d Interpolation done on images in playing