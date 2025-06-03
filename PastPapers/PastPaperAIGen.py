# cid = 02413400
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg


"""
Modifying an image scaling and pixel playing
"""
# Provide the correct relative or absolute path to your image
image_path = '/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/img_Data/FilterImage.png'
img_matrix = mpimg.imread(image_path) # Load the image as a NumPy array
img_mod = img_matrix.copy()  # creates a writable version

# 1) Display the image:

plt.figure()
plt.imshow(img_matrix)
plt.title("Original Image")



Nscale = 2 #Scale down size

Nx = img_mod.shape[1] # number of columns
Ny = img_mod.shape[0] # number of rows


# shrink the image, by taking a pixel every Nscale pixels
ImShrunk = np.zeros((int(Ny/Nscale),int(Nx/Nscale))) #initialise Shrink image array 

ImShrunk = img_mod[0:Ny:Nscale,0:Nx:Nscale] #including pixels at every Nscale'th p
plt.imshow(ImShrunk) 
plt.show()

"""
Task B ODE Diffusion
"""

# Newton's Law of Cooling: dT/dt = -k(T - T_env)
def func(t, T):
    k = 0.03  # per minute
    Tenv = 22  # ambient temp
    return -k * (T - Tenv)

# Forward Euler method
def FwEuler(T0, t0, dt, t_final):
    t_values = [t0]
    T_values = [T0]

    T = T0
    t = t0

    while t < t_final:
        T += dt * func(t, T)
        t += dt
        t_values.append(t)
        T_values.append(T)

    return np.array(t_values), np.array(T_values)

# Initial condition
T0 = 90
t0 = 0
t_final = 60  # 60 minutes
dt = 1        # step size in minutes

# Solve ODE
t_vals, T_vals = FwEuler(T0, t0, dt, t_final)

# Plot
plt.figure(figsize=(10, 4))
plt.plot(t_vals, T_vals, label='Temperature (Euler Method)', color='brown')
plt.xlabel('Time (minutes)')
plt.ylabel('Temperature (°C)')
plt.title("Cooling of Coffee Over Time")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


"""
Task C interpolating of roast quality and optimisation
"""

# opening and cleaning x vals of the logo
with open('/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/txt_Data/RoastTemps.txt', 'r') as d:
    rvals = d.readlines()

# Initialize lists to store cleaned float values
rnknown = []
for item in rvals:
    term = item.strip()
    try:
        clean_term = float(term)
        rnknown.append((clean_term))
    except ValueError:
        continue

rn = np.array(rnknown)
print("xn", rn)

# Step 1: Compute second derivatives (M values) for natural cubic spline
def compute_spline_coefficients(x, y):
    n = len(x)
    h = np.diff(x)

    # Set up matrix A and RHS vector b
    A = np.zeros((n, n))
    b = np.zeros(n)

    # Natural boundary conditions
    A[0, 0] = 1
    A[-1, -1] = 1

    for i in range(1, n-1):
        A[i, i-1] = h[i-1]
        A[i, i]   = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
        b[i] = 6 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])

    # Solve the linear system A * M = b
    M = np.linalg.solve(A, b)
    return M

# Step 2: Evaluate spline at a single point xp
def CubicSplineInterp(x, y, M, xp):
    n = len(x)
    for i in range(n - 1):
        if x[i] <= xp <= x[i+1]:
            h = x[i+1] - x[i]
            A = (x[i+1] - xp) / h
            B = (xp - x[i]) / h
            S = (A * y[i] + B * y[i+1] +
                 ((A**3 - A) * M[i] + (B**3 - B) * M[i+1]) * h**2 / 6)
            return S
    return None  # if xp is out of bounds

# Step 3: Interpolation driver
def InterpPoints(x_known, y_known, InterpMethod, xp_points, M=None):
    yp_points = np.zeros_like(xp_points)
    for i in range(len(xp_points)):
        yp_points[i] = InterpMethod(x_known, y_known, M, xp_points[i])
    return yp_points

# Example known data
x_known = np.array([0, 5, 10, 15, 20])
y_known = np.array([22, 24, 31, 45, 60])

# Compute spline coefficients
M_coeffs = compute_spline_coefficients(x_known, y_known)

# Interpolate
x_interp = np.arange(3, 17, 0.1)
y_interp = InterpPoints(x_known, y_known, CubicSplineInterp, x_interp, M_coeffs)

# Plotting
plt.figure(figsize=(8, 4))
plt.plot(x_interp, y_interp, label='Cubic Spline', color='black')
plt.scatter(x_known, y_known, color='red', label='Known Points')
plt.title("Cubic Spline Interpolation (Numerical)")
plt.xlabel("x")
plt.ylabel("y(x)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
