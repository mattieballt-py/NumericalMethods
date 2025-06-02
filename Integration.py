import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D  # for 3D plotting

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


"""
3D Integration – Cone and Torus
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
1. Volume of a cone using Trapezium Rule
"""

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
    return 2 * np.sum(Z) * dx * dy  # *2 for full height (±z)

dx = x[1] - x[0]
dy = y[1] - y[0]
V_torus = VolumeTorusMasked(Z, dx, dy)
print("Volume of torus (numerical double integration):", V_torus)

"""
3. 3D Plot of Torus Surface of Revolution
"""

# define upper semi-circle of torus profile
x_profile = np.linspace(R - r, R + r, 200)
y_profile = np.sqrt(r**2 - (x_profile - R)**2)

# define rotation mesh
theta = np.linspace(0, 2 * np.pi, 200)
X, T = np.meshgrid(x_profile, theta)
Y = np.tile(y_profile, (len(theta), 1))
Z = Y * np.sin(T)
Y = Y * np.cos(T)

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='plasma', edgecolor='none', alpha=0.9)
ax.set_title('3D Surface of a Torus')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.tight_layout()
plt.show()
