# using algorithms from the other py files in this folder
# putting them here as I often mess them up with different parameters and alterations

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg
import csv # for writing

"""
Dealing with csv's
"""

with open('txt_Data/Hob.Elements.txt', 'r') as d:
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
#print(diff)

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
2D Interpolation of Images
"""
# Provide the correct relative or absolute path to your image
image_path = '/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/img_Data/Flower.jpg'

# Load the image as a NumPy array
img_matrix = mpimg.imread(image_path)

M = np.ndarray((10,10 )) # empty 10 by 10 matrix
#Image_Matrix = plt.imread() # import an image into a matrix
#Matrix_to_Image = plt.imshow() # convert a matrix to an image
#Save_matrix2image = plt.imsave() # display matrix as an image? 
#size_as_tuple = M.shape() # should give size of matrix M

# Display the image
plt.imshow(img_matrix)
plt.axis('off')  # Hide axes
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
                row.append(round(float(part)))
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
    print(f"{var} cleaned:", cleaned_vals[:3])  # Preview cleaned data

    data[var] = cleaned_vals



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


