# using algorithms from the other py files in this folder
# putting them here as I often mess them up with different parameters and alterations

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg
import csv # for writing
import math as mt 

"""
Dealing with csv's
"""

with open('txt_Data/Rocket.txt', 'r') as d:
    Rckt = d.readlines()

# Initialize lists to store cleaned float values
RHeights = []

# Process each line in the Rocket file
for item in Rckt:
    # Strip leading/trailing whitespace
    term = item.strip()
    # Try to convert the term to a float
    try:
        clean_term = float(term)
        RHeights.append(round(clean_term)) # Round to the nearest integer if needed
    except ValueError:
        # If conversion fails, skip this item
        continue

# Example 2D array
array = np.array([[1.2, 2.5, 3.8],
                  [4.0, 5.1, 6.3],
                  [7.7, 8.8, 9.9]])

# Save to CSV
with open('txt_data/output.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(array)


"""
Wave Eqn Explicit Non-Dimensionalised, Newman LHS and RHS BC
"""
# main algorithm from PDE_Explicit.py
# explicit solution of d^2T/dx^2 - d^2T/dt^2 = 0  non dimensionalised

dx = 0.2 # space step
dt = 0.2 # time step
r = 1 # r = 1 for this
#print("convergence factor r",r)

x = np.arange(0,1+dx,dx) # mesh in x here gives 0, 0.2, 0.4, 0.6, 0.8, 1
t = np.arange(0,2+dt,dt) # mesh in time up to 2 seconds
U = np.zeros((len(t),len(x))) # initialise u array, time, space

# BC's:
U[:,0] = 0
U[0,:] = np.sin((np.pi)*x) # at t=0, u = sinwx
U[:,-1] = 0
U[1,:] = U[0,:]# because newman bc du/dt = 0 at t=0

def explicitWave(x,t,U,r): # with BC's set up as tutorial question
    
    for i in range(2,len(t)): # so at each time step   
        for j in range(1,len(x)-1): # from LHS to RHS stepping along x
            U[i,j] = U[i-1,j-1] + U[i-1,j+1] - U[i-2,j]

        # reapply BC's to next step:
        U[:,0] = 0
        U[0,:] = np.sin((np.pi)*x) # at t=0, u = sinwx
        U[:,-1] = 0
        U[1,:] = U[0,:]# because newman bc du/dt = 0 at t=0

    return U

diff = explicitWave(x,t,U,r)
print(diff)

"""
# Plot results
plt.figure(figsize=(10, 6))
for n in [0, int(len(t)/4), int(len(t)/2), len(t)-1]:
    plt.plot(x, diff[n, :], label=f't = {t[n]:.2f}')
plt.xlabel('x')
plt.ylabel('wave position')
plt.title('Wave Equation Solution - Explicit Method')
plt.legend()
plt.grid(True)
#plt.show()
"""

"""
2D Interpolation of Images
"""
# Provide the correct relative or absolute path to your image
image_path = '/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/img_Data/Flower.jpg'
img_matrix = mpimg.imread(image_path) # Load the image as a NumPy array

# shrink the flower image
n = 1
img_shrunk = img_matrix[::n, ::n] # skips every n as start,stop,step then , for next direction ::n
img_sml = plt.imsave('shrunk_flower.jpg',img_shrunk)
#print(img_matrix)

# Nearest Neighbour interpolation of image
def Neighb2d_Interp(A,n):
     # receive image matrix A, and scale to interpolate to, n
    orig_height, orig_width, channels = A.shape
    new_height, new_width = orig_height * n, orig_width * n
    Larger_Img = np.zeros((new_height, new_width, channels), dtype=A.dtype)
    Larger_Img[::n,::n,:] = A # fill in known values
    for i in range(0,new_height):
        for j in range(0,new_width):
            if (i % n != 0) or (j % n != 0):  # Skip known pixels as % is remainder from dividing
                nearest_i = round(i / n)
                nearest_j = round(j / n)

                # Clamp to image bounds to avoid overflow at the edges
                nearest_i = min(nearest_i, orig_height - 1)
                nearest_j = min(nearest_j, orig_width - 1)

                # Copy the pixel value
                Larger_Img[i, j] = A[nearest_i, nearest_j]
    return Larger_Img

# Bilinear Interpolation of image
def Interp2d_Bilinear(A,n):
    # receive image matrix A, and scale to interpolate to, n
    orig_height, orig_width, channels = A.shape
    new_height, new_width = orig_height * n, orig_width * n

    Larger_Img = np.zeros((new_height, new_width, channels), dtype=A.dtype)
    Larger_Img[::n,::n,:] = A # fill in known values
    for i in range(new_height):
        for j in range(new_width):
            if (i % n != 0) or (j % n != 0):  # Skip known pixels   
                i_l = (i//n)*n # top left hand corner (i,j)
                j_t = (j//n)*n
                i_r = min(i_l + n, new_height - 1) # min of +n or if at edge
                j_b = min(j_t + n, new_width - 1)
                # now we have the square points, (i_tlh,j_tlh) (i_tlh,j_brh) etc
                dx = (i-i_l)/n
                dy = (j-j_t)/n
                for c in range(channels): # retrieve corner values
                    f00 = Larger_Img[i_l, j_t, c]
                    f10 = Larger_Img[i_r, j_t, c]
                    f01 = Larger_Img[i_l, j_b, c]
                    f11 = Larger_Img[i_r, j_b, c]

                    p = ((1-dx)*(1-dy)*f00 + dx*(1-dy)*f10 + (1-dx)*dy*f01 + dx*dy*f11)

                    Larger_Img[i,j,c] = p
    return Larger_Img


#img_Interp= Interp2d_Bilinear(img_shrunk,2)
#img_nearestneigh = Neighb2d_Interp(img_shrunk,2)
# Display the image
#plt.imshow(img_Interp)
#plt.axis('off')  # Hide axes
#plt.show()

# img_matrix is now a NumPy array representing the image
#print("Image shape:", img_matrix.shape)


"""
3D Lagrangian
"""

def trilinear_interp(x, y, z, x_vals, y_vals, z_vals, f_vals):
    """
    Perform trilinear interpolation at point (x, y, z).

    Args:
        x, y, z: Coordinates where the function is interpolated
        x_vals: [x0, x1] - known x coordinates (must be sorted)
        y_vals: [y0, y1] - known y coordinates (must be sorted)
        z_vals: [z0, z1] - known z coordinates (must be sorted)
        f_vals: 2x2x2 numpy array containing function values at corners

    Returns:
        Interpolated scalar value at (x, y, z)
    """
    # Normalize x, y, z into [0, 1] local coordinates within the cube
    xd = (x - x_vals[0]) / (x_vals[1] - x_vals[0])
    yd = (y - y_vals[0]) / (y_vals[1] - y_vals[0])
    zd = (z - z_vals[0]) / (z_vals[1] - z_vals[0])

    # Interpolate along x
    c00 = f_vals[0, 0, 0] * (1 - xd) + f_vals[1, 0, 0] * xd
    c01 = f_vals[0, 0, 1] * (1 - xd) + f_vals[1, 0, 1] * xd
    c10 = f_vals[0, 1, 0] * (1 - xd) + f_vals[1, 1, 0] * xd
    c11 = f_vals[0, 1, 1] * (1 - xd) + f_vals[1, 1, 1] * xd

    # Interpolate along y
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    # Interpolate along z
    c = c0 * (1 - zd) + c1 * zd

    return c

"""
Teapot Barycentric Interpolation of a triangulated mesh hob and teapot
"""
# have two meshes, want to find if function values overlap
# if  they do, interpolate the teapot temps from the hob temps
data = {}

def read_lines(filename):
    with open(filename, 'r') as f:
        return f.readlines()

# Helper: clean & convert strings to float (rounded), skip invalids
def clean_floats(lines):
    cleaned = []
    for line in lines:
        parts = line.strip().replace(',', ' ').split()  # convert commas to spaces first
        row = []
        for part in parts:
            try:
                row.append((float(part)))
            except ValueError:
                continue
        if row:
            cleaned.append(row)
    return cleaned


# File mapping: variable name -> path
file_map = {
    'Hb_Elements': 'txt_Data/Hob.Elements.txt',
    'Hb_Nodes': 'txt_Data/Hob.Nodes.txt',
    'Hb_Temps': 'txt_Data/Hob.Temperatures.txt',
    'Tp_Elements': 'txt_Data/TeaPot.Elements.txt',
    'Tp_Nodes': 'txt_Data/TeaPot.Nodes.txt'
}


# Load and clean files into variables
for var, path in file_map.items():
    raw_lines = read_lines(path)
    #print(f"{var} raw:", raw_lines[:3])  # Preview first 3 lines

    cleaned_vals = clean_floats(raw_lines)
    #print(f"{var} cleaned:", cleaned_vals[:3])  # Preview cleaned data

    data[var] = cleaned_vals

Hb_Elements = data['Hb_Elements']
Hb_Nodes = data['Hb_Nodes']
Hb_Temps = data['Hb_Temps']
Tp_Elements = data['Tp_Elements']
Tp_Nodes = data['Tp_Nodes']



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


"""
Root Finding to find Water Tank min Volume
"""
# given volume equation in regards to parameters
# using 2 of 2 bracketing methods: False Position Method

R = 2 # radius of water tank
def f(x):
    return (3.14 * R * x **2) - ((3.14/3)*x**3) - 25


# false position method recursively
def falseposition(a,b,tol,i):
    if f(a)*f(b)>0: # check valid interval
        print("choose new interval")
        return None

    n = b - (f(b)*(a-b))/(f(a) - f(b)) #new x point, at intersection with 0 instead of bang in middle like bisection method
    #print('new n', n)

    if f(a)*f(n)<0: # function changes sign between a and n
        # best guess before was n, now n + 1 = (a + n)/2
        # err = (new - old) /new
        n_plus1 = n - (f(n)*(a-n))/(f(a) - f(n)) #new x point, at intersection with 0 instead of bang in middle like bisection method
        errnew = ((n_plus1 - n)/n_plus1)
        #print("ner error: ",errnew)
        if abs(errnew) < tol:
            return ((a + n)/2) # if interval less than min, the root is at centre of it ish
        else: # replace b with n
            i += 1
            #print(i)
        return falseposition(a,n,tol,i)
  
    else:  # function changes sign between n and b
        n_plus1 = b - (f(b)*(n-b))/(f(n) - f(b)) #new x point, at intersection with 0 instead of bang in middle like bisection method
        errnew = ((n_plus1 - n)/n_plus1)
        #print("ner error: ",errnew)
        if abs(errnew) < tol:
            return ((n + b)/2) # if interval less than min, the root is at centre of it ish
        i += 1
        #print(i)
        return falseposition(n,b,tol,i)


#print(falseposition(0,2*R,0.001,1))
# a, b is the interval to evaluate inside

### with Newton Raphson:

# need to differentiate f(x) aka V(h) (central as more accurate)
# Numerical derivative using central difference
def numerical_derivative(h, delta):
    return (f(h + delta) - f(h - delta)) / (2 * delta)

# if linear system of eqns, could solve with x = A-1 linalgsolve with b
# for non linear syst eqns like x^2 + y^2 = 2 etc, use:
def NewtonRaph(x0,tol):
    xn = x0 # guess
    err = 10 # init error
    while err > tol:
        delta=1e-5
        xn_1 = xn - f(xn)/numerical_derivative(xn,delta)
        err = abs((xn_1 - xn)/xn_1)
        xn = xn_1
        #print("tryin again")
    return xn

#print(NewtonRaph(1.99*R,0.01))


# plotting the water tank function
x = np.arange(-2,5,0.01)
y = f(x)
# Scatter plot the interpolated surface
plt.figure(figsize=(6, 6))
sc = plt.scatter(x, y)
plt.title('Barycentric Interpolation (Simple Visual Test)')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.show()

## with secant method:

# 2 of 2 open method for when can't f'(x)
def secant(x0,x1,tol): # two initial guesses
    xn = x0 # init
    xn_1 = x1
    err = 100 # just to init
    while err > tol:
        xn_2 = xn_1 - (f(xn_1)*(xn - xn_1))/(f(xn)-f(xn_1))
        err = abs((xn_2 - xn_1)/xn_2)
        xn,xn_1= xn_1,xn_2
        #print("tried again secant")
    return xn_2

print(secant(R,1.99*R,0.01))


"""
Newton Raphson for Multiple Roots
"""
# Numerical derivative analytically
def f(x):
    return x**3 - 2*x**2 - 4*x + 8

def df(x):
    return 3*x**2 - 4*x - 4

# if linear system of eqns, could solve with x = A-1 linalgsolve with b
# for non linear syst eqns like x^2 + y^2 = 2 etc, use:
def NewtonRaph(x0,tol):
    xn = x0 # guess
    err = 10 # init error
    while err > tol:
        delta=1e-5
        xn_1 = xn - 2*f(xn)/df(xn)
        err = abs((xn_1 - xn)/xn_1)
        xn = xn_1
        #print("tryin again")
    return xn

#print(NewtonRaph(1.2,0.01))


"""
2D Integration of Royal Albert Hall
"""

# function trapz: compute numerical integration with trapezium rule, for nodes at any distance
def trapz(x,y):
    # get the number of subintervals
    N = len(x) - 1
    # compute the integral
    # set range for the trapezia: there are as many trapezia as the number of intervals
    R = range(0,N)
    S = 0
    for i in R:
        # compute the area of this single trapezium (remind yourself the area of a trapezium)
        S += 0.5 * (y[i+1] + y[i]) * (x[i+1] - x[i])
    return S

a = 67 # major axis
b = 56 # minor axis
h = 25 # minor axis

# set the step intervals in x and y
dx = 0.5
dy = 0.5

# set the x range, not including the boundaries
x = np.arange(-a+dx,a,dx)
N = len(x)
# the y range depends of the various values of x, and cannot be fixed here

# integrate in dy, for all the value of x, i.e. find G(x)

G = np.zeros(N)
# for every x
for i in range(0,N):
    # determine the boundaries m and p for this x
    mx = np.sqrt(b**2*(1-x[i]**2/a**2))
    px = mx
    # set the y points for this x, not including the boundaries
    y = np.arange(-mx+dy,px,dy)
    z = np.zeros(len(y))
    # determine the values of the function z(x,y)
    for j in range(0,len(y)):
        z[j] = np.sqrt(h**2*(1-x[i]**2/a**2-y[j]**2/b**2)) 
    
    # integrate in dy from cx to dx (for this specific x)
    G[i] = trapz(y,z) # G(x)

# integrate G(x) in dx
I = trapz(x,G)

# for an emisphere the volume is:
print((4/3*np.pi*a*b*h)/2)

"""
3D Integration of an aerofoil
"""

import numpy as np
import matplotlib.pyplot as plt

data = {}

def read_lines(filename):
    with open(filename, 'r') as f:
        return f.readlines()

# Helper: clean & convert strings to float (rounded), skip invalids
def clean_floats(lines):
    cleaned = []
    for line in lines:
        parts = line.strip().replace(',', ' ').split()  # convert commas to spaces first
        row = []
        for part in parts:
            try:
                row.append(float(part))
            except ValueError:
                continue
        if row:
            cleaned.append(row)
    return cleaned

# File mapping: variable name -> path
file_map = {
    'Aero': 'txt_Data/Aerofoil.txt',
}

# Load and clean files into variables
for var, path in file_map.items():
    raw_lines = read_lines(path)
    cleaned_vals = clean_floats(raw_lines)
    data[var] = cleaned_vals

AeroAll = data['Aero']
print("Aeroall:", AeroAll)
print("Just tn=12:", AeroAll[12])

x_terms = [item[0] for item in AeroAll]
y_terms = [item[1] for item in AeroAll]
Zt_terms = [item[2] for item in AeroAll]  # top surface
Zb_terms = [item[3] for item in AeroAll]  # bottom surface

# Create a 3D surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# You need to reshape or simulate a meshgrid, but for now, just scatter top surface
ax.plot_trisurf(x_terms, y_terms, Zt_terms, cmap='Oranges')

# Set labels and title
ax.set_title('Top of Aerofoil')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.view_init(15, 60)

# function trapz: compute numerical integration with trapezium rule, for nodes at any distance
def trapz(x, y):
    N = len(x) - 1
    R = range(0, N)
    S = 0
    for i in R:
        S += 0.5 * (y[i + 1] + y[i]) * (x[i + 1] - x[i])
    return S

a = 67  # major axis
b = 56  # minor axis
h = 25  # height

dx = 0.5
dy = 0.5

x = np.arange(-a + dx, a, dx)
N = len(x)

G = np.zeros(N)

for i in range(N):
    mx = np.sqrt(b**2 * (1 - x[i]**2 / a**2))
    y = np.arange(-mx + dy, mx, dy)
    z = np.zeros(len(y))
    for j in range(len(y)):
        term = 1 - x[i]**2 / a**2 - y[j]**2 / b**2
        z[j] = np.sqrt(h**2 * term) if term >= 0 else 0
    G[i] = trapz(y, z)

I = trapz(x, G)
print("Integrated volume I =", I)

plt.show()


"""
Langrangian for CID no. Points
"""

xn = np.arange(1,9) # xn points 1,2,3,4,5,6,7,8
yn = np.array([0,2,4,1,3,4,0,0])

dx = 0.1
xp = np.arange(1,8+dx,dx)# from 1 to 8 with dx=0.1 steps

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

yp = InterpPoints(xn,yn,Langrangian,xp)
print("interpolated lagrangian points: ",yp)

# Plot results
plt.figure(figsize=(10, 6))
plt.scatter(xn,yn,c="red", s=20) # known (s = controls size of points)
plt.scatter(xp,yp,c="blue", s=5) # interpolated
plt.xlabel('x')
plt.ylabel('y')
plt.title('Lagrnagian interpolation 0 to 8 with my CID')
plt.legend()
plt.grid(True)
plt.show()


"""
Solving a 2nd Order ODE
"""
# Set parameters and initial conditions
x0 = 0
x_end = 15
dx = 0.02
x = np.arange(x0, x_end + dx, dx)

# Initial conditions (use your CID digits as required)
y0 = 4      # y(0) = 3rd digit of CID
dy0 = 3     # y'(0) = 5th digit of CID

# Initialize arrays
y = np.zeros_like(x)
dy = np.zeros_like(x)

# Apply initial values
y[0] = y0
dy[0] = dy0

# --------------------------------------
# Define your second-order ODE here:
# dy2 = f(x, y, dy)
# --------------------------------------
def f(x, y, dy):
    return -5*x*dy - (x + 7)*np.sin(x)

# --------------------------------------
# Forward Euler Solver
# --------------------------------------
for i in range(1, len(x)):
    y[i] = y[i-1] + dx * dy[i-1]
    dy[i] = dy[i-1] + dx * f(x[i-1], y[i-1], dy[i-1])

# --------------------------------------
# Plotting
# --------------------------------------
plt.figure(figsize=(10, 5))
plt.plot(x, y, label='y(x)', linewidth=2)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Forward Euler Solution to Second-Order ODE")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


"""
Integrate Area between circle and a function
"""
dx = 0.01
xn = np.arange(0,13+dx,dx)

def circ(xn): # (x - 3)^2 + y^2 = 10^2
    return (100 - (xn-3)**2)**0.5


def TrapInt(h,N,y_known):
    I = 0 # init
    for n in range(1,N):
        I += y_known[n] 
    I = h*(I + y_known[0]/2 + y_known[-1]/2)
    return I

N = len(xn)
yn = circ(xn)
I1 = TrapInt(dx,N,yn) # this gives area under circle from x = 0 to 13

# Red curve function, clipped to only keep positive y-values
def yfunc(xn, n):
    y = 5 * np.sin((2 * mt.pi * n * xn) / 13) * np.exp(-xn / 10)
    y[y < 0] = 0
    return y
ns = np.array([0,2,4,1,3,4,0,0])

I=[] # list to store integral values of area for each n
for i in range(0,len(ns)):
    I +=[TrapInt(h,N,yfunc(xn,ns[i]))]
         
Itot = [] # list of final intervals for each n
for j in range(0,8):
    Itot += [I1 - I[j]]

print("I",I,"Itot",Itot,"I1",I1)
plt.figure(figsize=(10, 5))
plt.scatter(ns, Itot, label='Int', linewidth=2)
plt.xlabel("n")
plt.ylabel("I")
plt.title("Integral area with trap rule for various n's")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# clearer: 
# Define the function y = 5sin((2π/13) * n * x) * exp(-x/10)
def red_function(x, n):
    return 5 * np.sin((2 * np.pi * n * x) / 13) * np.exp(-x / 10)

# Define the upper boundary of the circle (since the circle is symmetric about x-axis)
def circle_upper(x):
    return np.sqrt(100 - (x-3)**2)  # from x^2 + y^2 = 100 => y = sqrt(100 - x^2)

# Trapezoidal rule implementation
def trapezoidal_rule(x, y):
    return np.trapz(y, x)

# Set up the x range and step size
x = np.arange(0, 13, 0.01)

# Values of n to evaluate (1st to 8th digit of student CID)
n_values = ns
areas = []

# Calculate area between circle and red curve for each n
for n in n_values:
    y_curve = red_function(x, n)
    y_circle = circle_upper(x)

    # Ensure we're only integrating where the curve is below the circle
    diff = y_circle - y_curve
    diff[diff < 0] = 0  # Ignore negative area (i.e., where curve is above circle)

    area = trapezoidal_rule(x, diff)
    areas.append(area)

    print(f"n = {n}, Area = {area:.4f}")

# Plotting the results
plt.figure(figsize=(10, 5))
plt.scatter(n_values, areas, marker='o', linestyle='-', color='green')
plt.title("Area Between Circle and Red Curve vs n")
plt.xlabel("n (digit of CID)")
plt.ylabel("Area (using Trapezoidal Rule)")
plt.grid(True)
plt.show()


#prof answer int between two lines:
#Generate domain of integration and eqn for circle:
R = 10
dx = 0.01
x = np.arange(0,13+dx,dx)
yup = np.sqrt(R**2-(x-3)**2)
#Area of circle:
Sup = TrapInt(dx,len(x),yup)
#Loop for every digit of the CID:
S = []
Rn = [0,1,2,3,4,5,6,7,8,9]
for n in Rn:
    yt = 5*np.sin(2*np.pi/13*n*x)*np.exp(-x/10)
    ydown = np.zeros(len(x))
    ydown[yt>=0] = yt[yt>=0]
    Sdown = TrapInt(dx,len(x),ydown) 

    S += [Sup-Sdown]
print("profs s",S)