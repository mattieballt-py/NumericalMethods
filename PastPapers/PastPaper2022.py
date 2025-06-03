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


"""
Task B: Heat Diffusion PDE
"""

# this specific solution is from the ice cream Q, adjust paramaters, bc's etc as needed:
# spatial domain (in mm)
R = 30
dx = 0.5 # spatial increment
# greate the spatial grid points
x = np.arange(-R,R+dx,dx)
Nx = len(x)

# temporal domain
tend = 3600
dt = 1 # temporal increment
# create the temporal grid points
t = np.arange(0,tend+dt,dt)
Nt = len(t)

# create the solution matrix
T = np.ndarray((Nt,Nx))

# set the physics
# boundary conditions (fixed temperature at boundaries)
Ta = 35
Tb = 35
# set the initial T for watermelon
T[0,:int(Nx/2)] = 2
# set the initial T for pistachio
T[0,int(Nx/2):] = 6

# set the thermal diffusivity
alpha = np.ndarray(len(x))
# set the diffusivity T for watermelon
alpha[:int(Nx/2)] = 4.5e-2
# set the diffusivity for pistachio
alpha[int(Nx/2):] = 8.6e-2

averagep = np.average(T[0,int(Nx/2):])
averagew = np.average(T[0,:int(Nx/2)])
avp = [averagep]
avw = [averagew]
Tmelt = 15

c = alpha * dt / dx**2

# compute the solution incrementally at subsequent time steps
k = 1
while averagep < Tmelt:
    # compute at time step p, i.e. t = k * dt
    # do it for every node in the spatial grid
    # start with the boundaries
    T[k,0] = Ta
    T[k,-1] = Tb
    # do the interior nodes
    for i in range(1,Nx-1):
        # apply the discretised equation
        T[k,i] = c[i] * ( T[k-1,i+1] + T[k-1,i-1] ) + (1 - 2*c[i]) * T[k-1,i]
    averagep = np.average( T[k,int(Nx/2):])
    averagew = np.average( T[k,:int(Nx/2)])
    avp += [averagep]
    avw += [averagew]
    k += 1


plt.plot(x,T[k-1,:])
plt.grid()
plt.show()


plt.plot(t[:k],avp,c='Green')
plt.plot(t[:k],avw,c='Red')
plt.grid()
plt.show()

print('Ice cream melted after minutes:')
print((k-1)*dt/60)

print('Courant:')
print(c[-1],c[0])

"""
Task C: Lagrangian Interpolation of randomly spaced co ords for a logo
"""

# opening and cleaning x vals of the logo
with open('/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/txt_Data/LogoXn.txt', 'r') as d:
    xvals = d.readlines()

# Initialize lists to store cleaned float values
xnknown = []
for item in xvals:
    term = item.strip()
    try:
        clean_term = float(term)
        xnknown.append((clean_term))
    except ValueError:
        continue

xn = np.array(xnknown)
print("xn", xn)

with open('/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/txt_Data/LogoYn.txt', 'r') as d:
    yvals = d.readlines()

# Initialize lists to store cleaned float values
ynknown = []
for item in yvals:
    term = item.strip()
    try:
        clean_term = float(term)
        ynknown.append(clean_term)
    except ValueError:
        continue
yn = np.array(ynknown)
print("yn", yn)


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
x_interp = np.arange(-19, 19+0.5, 0.5)
Temp_interp = InterpPoints(xn, yn, Langrangian, x_interp)

# Plotting interpolated section
plt.figure(figsize=(8, 4))
plt.plot(x_interp, Temp_interp, color='black')
plt.scatter(xn, yn, color='red',s = 50) # make smaller so can see interp points
plt.title("Interpolated points for root finding")
plt.xlabel("x")
plt.ylabel("y(x)")
plt.grid(True)
plt.tight_layout()
plt.show()

"""
3D Cone Integration with trapezium rule
"""

#1. Volume of a cone using Trapezium Rule in mm


# define parameters for the cone
A = 0      # lower limit (base)
B = 150    # upper limit (height)
h = 5 # dz
N = int((B-A)/h)  # number of intervals


# profile function for cone (linear slope)
def f_cone(x):
    return 0.2*x  # straight line from base to tip

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


# or

# in cartesian: (discretization in z too, much more accurate. for when dx,dy and dz given)
# function trapzeqd: compute numerical integration with trapezium rule with equidistant nodes
def trapzeqd(x,y):
    # get the interval h: distance between any two consecutives nodes
    h = x[1] - x[0]
    # get the number of intervals
    N = len(x) - 1  # obviously x and y must have same length
    
    # compute the integral
    
    # compute the sum for the intermediate points
    S = 0.0
    for n in range(1,N):
        S += y[n]  # add the current calue of y
    # add first and last points: see the formula for trapezoidal method
    I = h * (y[0]/2 + S + y[-1]/2 )
    
    # an alternative approach, with slicing and the function np.sum(), the integral can be computed within one line
    I = h * (y[0]/2 + np.sum(y[1:-1]) + y[-1]/2 )
    
    return I
# set dimensions (in mm)
H = 150
R = 30

# set discretising steps
dx = 0.05
dy = 0.05
dz = 5

# integrate for every z
# create the z range
zr = np.arange(dz,H+dz,dz)
# determine all the circular sections surfaces
Iz = []
for z in zr:
    # integrate along x
    Rz = 0.2*z
    # define the x range
    xr = np.arange(-Rz+dx,Rz,dx)
    Ix = []
    for x in xr:
        # integrate along y
        # define the y range
        m = np.sqrt(Rz**2-x**2)
        yr = np.arange(-m+dy,m,dy)
        # define the function to be integrated (all ones)
        F = np.ones(len(yr))
        # integrate with trapezium rule along y
        Ix += [trapzeqd(yr,F)]
    Ix = np.array(Ix)
    # integrate with trapezium rule along x
    Iz += [trapzeqd(xr,Ix)]        
# integrate with trapezium rule along x
Iz = np.array(Iz)
I = trapzeqd(zr,Iz)      

print(I)

# double check with analytical formula
I = np.pi*R**2*H/3
print(I)