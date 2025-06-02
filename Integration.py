import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D  # for 3D plotting
import math as mt  

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

# or for discrete y values
def TrapIntDisc(A,B,x_known,y_known):
    I = 0 # init
    h = x_known[1] - x_known[0]
    for n in range(1,len(y_known)):
        I += y_known[n] 
    I = h*(I + y_known[0] + y_known[-1]/2)
    return I

"""
3D Integration – Cone and Torus
"""


#1. Volume of a cone using Trapezium Rule


# define parameters for the cone
A = 0      # lower limit (base)
B = 3      # upper limit (height)
N = 1000   # number of intervals
h = (B - A) / N

# profile function for cone (linear slope)
def f_cone(x):
    return x  # straight line from base to tip

# known x values
x_known = np.linspace(A, B, N+1)
y_known = f_cone(x_known)

# V = π ∫[a to b] (f(x))² dx – volume around x-axis
def VolumeTrapX(A, B, h, N, y_known):
    integrand = np.pi * y_known**2
    I = 0
    for n in range(1, N):
        I += integrand[n]
    I = h * (I + integrand[0]/2 + integrand[-1]/2)
    return I

print("Volume of cone (around x-axis):", VolumeTrapX(A, B, h, N, y_known))

# 3D surface of revolution plot (cone)
theta = np.linspace(0, 2 * np.pi, 100)
X, T = np.meshgrid(x_known, theta)
R = f_cone(X)
Y = R * np.cos(T)
Z = R * np.sin(T)

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='plasma', edgecolor='none', alpha=0.9)
ax.set_title('3D Surface of Cone (Revolved around x-axis)')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.tight_layout()
plt.show()

"""
2. Volume of torus defined by:
(R - sqrt(x² + y²))² + z² = r²
"""

# define torus parameters
R = 10  # distance from center to tube
r = 4   # radius of the tube

# define grid in x and y (centered at 0)
x = np.arange(-R - r, R + r, 0.01)
y = np.arange(-R - r, R + r, 0.01)
X, Y = np.meshgrid(x, y)

# create mask for points inside the torus
inside = ((R - np.sqrt(X**2 + Y**2))**2) < r**2

# solve for z using: z² = r² - (R - sqrt(x² + y²))²
Z = np.zeros_like(X)
Z[inside] = np.sqrt(r**2 - (R - np.sqrt(X[inside]**2 + Y[inside]**2))**2) 

# V = ∭ dV over enclosed volume
# integrate z twice (positive and negative halves)
def VolumeTorusMasked(Z, dx, dy):
    return 2* np.sum(Z) * dx * dy  # *2 for full height (±z)

dx = x[1] - x[0]
dy = y[1] - y[0]
V_torus = VolumeTorusMasked(Z, dx, dy)
print("Volume of torus (numerical double integration):", V_torus)

#3. 3D Plot of Torus Surface of Revolution

# Define torus parameters
R = 10  # distance from center of hole to tube center
r = 4   # radius of the tube

# Create parameter grid
theta = np.linspace(0, 2 * np.pi, 100)  # around tube
phi = np.linspace(0, 2 * np.pi, 100)    # around center
theta, phi = np.meshgrid(theta, phi)

# Parametric equations for torus
X = (R + r * np.cos(theta)) * np.cos(phi)
Y = (R + r * np.cos(theta)) * np.sin(phi)
Z = r * np.sin(theta)


# Plot torus with same scale axis
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='plasma', edgecolor='none', alpha=0.9)

# Compute axis limits for equal aspect
x_limits = [X.min(), X.max()]
y_limits = [Y.min(), Y.max()]
z_limits = [Z.min(), Z.max()]

x_range = x_limits[1] - x_limits[0]
y_range = y_limits[1] - y_limits[0]
z_range = z_limits[1] - z_limits[0]
max_range = max(x_range, y_range, z_range)

x_middle = np.mean(x_limits)
y_middle = np.mean(y_limits)
z_middle = np.mean(z_limits)

ax.set_xlim(x_middle - max_range/2, x_middle + max_range/2)
ax.set_ylim(y_middle - max_range/2, y_middle + max_range/2)
ax.set_zlim(z_middle - max_range/2, z_middle + max_range/2)

ax.set_title('3D Parametric Torus')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.tight_layout()
plt.show()


"""
Generic 3D Integration using the Trapezium Rule
Use: I = ∫∫ f(x, y) dx dy over grid
"""

# Example function: f(x, y) = x * y
def f_3d(x, y):
    return x * y

# Define 2D grid
x = np.linspace(0, 1, 100)
y = np.linspace(0, 2, 100)
X, Y = np.meshgrid(x, y)
Z = f_3d(X, Y)

def Trap3D(X, Y, Z):
    dx = X[0,1] - X[0,0]
    dy = Y[1,0] - Y[0,0]
    I = 0
    for i in range(1, Z.shape[0]-1):
        for j in range(1, Z.shape[1]-1):
            I += Z[i,j]
    I *= dx * dy
    return I

print("Generic 3D integral (x*y over 0≤x≤1, 0≤y≤2):", Trap3D(X, Y, Z))


"""
2D Simpson’s Rule Integration over triangular domain Double Integral of f(x,y) dxdy
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

        # ==== Simpson’s rule in y for fixed x = xi ====
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

        # ==== Simpson’s coefficient in x‐direction ====
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
