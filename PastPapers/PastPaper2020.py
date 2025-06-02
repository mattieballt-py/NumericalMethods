# CID: 02413400 

# import required libraries
import numpy as np
import matplotlib.pyplot as plt  
import matplotlib.image as mpimg
import csv # for writing

"""
Task A: Modify Banksy image colours
"""
# Provide the correct relative or absolute path to your image
image_path = '/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/img_Data/BanksyEx.jpg'
img_matrix = mpimg.imread(image_path) # Load the image as a NumPy array
img_mod = img_matrix.copy()  # creates a writable version

# 1) Display the image:
plt.imshow(img_matrix)
plt.axis('off')  # Hide axes
#plt.show()

# 2) Render blue and red
# img_matrix is 333x500 matrix and rgb, top triangle only b, bottom only r
slope = (np.shape(img_matrix)[0])/(np.shape(img_matrix)[1])
# looping through each pixel
for i in range(0,(np.shape(img_matrix)[0])): # going down height
    for j in range(0,(np.shape(img_matrix)[1])): # going along width half way
        # Top triangle → BLUE (remove red and green)
        if i < slope * j and i < (np.shape(img_matrix)[0]) - slope * j:
            img_mod[i, j, 0] = 0  # red off
            img_mod[i, j, 1] = 0  # green off

        # Bottom triangle → RED (remove green and blue)
        elif i > slope * j and i > (np.shape(img_matrix)[0]) - slope * j:
            img_mod[i, j, 1] = 0  # green off
            img_mod[i, j, 2] = 0  # blue off

plt.imsave('BanksyMod.jpg',img_mod)


"""
Task B: interpolate langrangian CID
"""
xn = np.arange(1,9)
print(xn)
yn = np.array([0,2,4,1,3,4,0,0])

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
x_interp = np.arange(1, 8, 0.1)
f_interp = InterpPoints(xn,yn, Langrangian, x_interp)

# Plotting interpolated section
plt.figure(figsize=(8, 4))
plt.scatter(x_interp, f_interp, color='black')
plt.scatter(xn,yn, color="red", s = 10) # s for marker size, smaller so can see as well as interp points
plt.title("Interpolated points for root finding")
plt.xlabel("x")
plt.ylabel("cid")
plt.grid(True)
plt.tight_layout()
plt.show()



"""
Task C: Solve an ODE Explicit Fwd Euler
"""

def func(x, Y):
    y0, y1 = Y
    dy0 = y1
    dy1 = -5*y1 - (x + 7)*np.sin(x)
    return np.array([dy0, dy1])

# Forward Euler Method for systems
def FwEuler(Y0, x0, h, x_final):
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

# Initial conditions
y0_init = 4
y1_init = 3

Y0 = np.array([y0_init, y1_init])

# x domain
x0 = 0
x_final = 15
h = 0.02 # dx

# Solve using FWD Euler method
x_vals, y_vals = FwEuler(Y0, x0, h, x_final)

# Extract relevant components
y_x = y_vals[:, 0]         # y(t)
y2_x = y_vals[:, 1]        # d²y/dt²

# Plot y(x) 
plt.figure(figsize=(10, 4))
plt.plot(x_vals, y_x, label='y(x)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Numerical Solution y(x) using Forward Euler')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

# Plot d²y/dx²
plt.figure(figsize=(10, 4))
plt.plot(x_vals, y2_x, label='d²y/dx²', color='orange')
plt.xlabel('x')
plt.ylabel('d²y/dx²')
plt.title('Second Derivative of y(x)')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

