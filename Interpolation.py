import numpy as np
import matplotlib.pyplot as plt  
# this 2d Interpolation done on images in playing

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

def lagrange_basis(x, points, i):
    #Computes the i-th Lagrange basis polynomial L_i(x) for given grid points.
    # x (float): Evaluation point
    # points (list): Grid coordinates (e.g., [x0, x1])
    # i (int): Index for basis function

    L = 1.0
    for j in range(len(points)):
        if j != i:
            L *= (x - points[j]) / (points[i] - points[j])
    return L

def bilinear_lagrange_interpolation(x, y, x_vals, y_vals, f_vals):
    """
    Bilinear interpolation using Lagrange basis functions in 2D.

    Args:
        x (float): x-coordinate to interpolate
        y (float): y-coordinate to interpolate
        x_vals (list): [x0, x1]
        y_vals (list): [y0, y1]
        f_vals (2x2 array): Function values at grid points:
            [[f(x0,y0), f(x1,y0)],
             [f(x0,y1), f(x1,y1)]]

    Returns:
        float: Interpolated value at (x, y)
    """
    interpolated = 0.0
    for i in range(2):
        Li = lagrange_basis(x, x_vals, i)
        for j in range(2):
            Lj = lagrange_basis(y, y_vals, j)
            interpolated += f_vals[j][i] * Li * Lj
    return interpolated

# Define grid
x_vals = [0.0, 1.0]
y_vals = [0.0, 1.0]

# Define the actual function
def f_actual(x, y):
    return np.exp(x * np.sinh(y))  # Example function

# Compute function values at corners:
# [[f(x0,y0), f(x1,y0)],
#  [f(x0,y1), f(x1,y1)]]
f_vals = np.array([
    [f_actual(x_vals[0], y_vals[0]), f_actual(x_vals[1], y_vals[0])],
    [f_actual(x_vals[0], y_vals[1]), f_actual(x_vals[1], y_vals[1])]
])

# Interpolation point
x_interp = 0.32
y_interp = 0.63

# Interpolate
z_interp = bilinear_lagrange_interpolation(x_interp, y_interp, x_vals, y_vals, f_vals)
print(f"Interpolated value at ({x_interp:.2f}, {y_interp:.2f}) = {z_interp:.4f}")

# Generate grid for plotting
X = np.linspace(0, 1, 20)
Y = np.linspace(0, 1, 20)
Z = np.zeros((len(Y), len(X)))

for i, y in enumerate(Y):
    for j, x in enumerate(X):
        Z[i, j] = bilinear_lagrange_interpolation(x, y, x_vals, y_vals, f_vals)

# Plot
plt.figure(figsize=(6, 5))
plt.contourf(X, Y, Z, levels=20, cmap='viridis')
plt.colorbar(label='Interpolated z')
plt.scatter([x_interp], [y_interp], color='red', label='Interpolated Point')
plt.title('Bilinear Interpolation using Lagrange Polynomials')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.tight_layout()
#plt.show()

### 3D Lagrangian in Playing

"""
Barycentric Interpolation
"""

# Define the actual function
def f_actual(pt):
    x,y = pt
    return np.exp(x * np.cos(y))  # Example function

def Barycentric(r,func,p):
    lam = np.zeros(3, dtype=float) # lambda 1  = area U/(U+V+W) # find areas with cross product of vectors between corners/2
    lam[0] = ((r[1][1] - r[2][1])*(p[0] - r[2][0]) + (r[2][0] - r[1][0])*(p[1] - r[2][1]))/((r[1][1] - r[2][1])*(r[0][0] - r[2][0]) + (r[2][0] - r[1][0])*(r[0][1] - r[2][1]))
    lam[1] = ((r[2][1] - r[0][1])*(p[0] - r[2][0]) + (r[0][0] - r[2][0])*(p[1] - r[2][1]))/((r[1][1] - r[2][1])*(r[0][0] - r[2][0]) + (r[2][0] - r[1][0])*(r[0][1] - r[2][1]))
    lam[2] = 1 - lam[0] - lam[1]
    for i in range(0,2):
        if lam[i] < 0:
            return print("p is not inside the three points in r")
    return lam[0]*func(r[0]) + lam[1]*func(r[1]) + lam[2]*func(r[2])

r = np.array([[1, 0], [0, 0], [0.5, 1]])
p = np.array([0.45, 0.42])

# Just calculate interpolated value at p
Bary = Barycentric(r, f_actual, p)
print(f"Barycentric interpolation at {p}: {Bary:.4f}")

# Grid for plotting
x_vals = np.linspace(0, 1, 300)
y_vals = np.linspace(0, 1.1, 300)

xs, ys, zs = [], [], []

for x in x_vals:
    for y in y_vals:
        pt = np.array([x, y])
        z = Barycentric(r, f_actual, pt)
        if z is not None:
            xs.append(x)
            ys.append(y)
            zs.append(z)

# Scatter plot the interpolated surface
plt.figure(figsize=(6, 6))
sc = plt.scatter(xs, ys, c=zs, cmap='coolwarm', s=2)
plt.plot(*r[[0, 1, 2, 0]].T, c='k', linewidth=2)  # triangle outline
plt.scatter(*p, color='black', label='Interpolated point')
plt.colorbar(sc, label='Interpolated Value')
plt.title('Barycentric Interpolation (Simple Visual Test)')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.legend()
plt.grid(True)
#plt.show()



"""
2D Interpolation
"""
def f(coord):
    x,y = coord
    return x*y

r = np.array([[1, 1], [1, 4], [5, 1]])
func = np.array([f(pt) for pt in r]) # creates array cycling through tuples of points in r
p = np.array([1.2, 1.3])
print(r,f,p)

def inverse_dist_triangle(r,func,p): # uses 1/distance between point and triangle corners to find weighting
    # weightings are 1/distance from p to corners r
    # f(p) = these weightings times the function values at the rs
    w = np.zeros(3, dtype=float)
    d = np.zeros(3, dtype=float)

    print("w",w, "r",r)
    for i in range(0,len(r)): # going through each point
        d[i] = ((r[i][0] - p[0])**2 + (r[i][1] - p[1])**2)**0.5
        print("d",d)
        w[i] = 1/d[i]
    return (w[0] * func[0] + w[1] * func[1] + w[2] * func[2])/(w[0] + w[1] + w[2])

#print(inverse_dist_triangle(r,func,p))

