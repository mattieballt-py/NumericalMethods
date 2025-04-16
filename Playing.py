# using algorithms from the other py files in this folder
# putting them here as I often mess them up with different parameters and alterations

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg

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

# Plot results
plt.figure(figsize=(10, 6))
for n in [0, int(len(t)/4), int(len(t)/2), len(t)-1]:
    plt.plot(x, diff[n, :], label=f't = {t[n]:.2f}')
plt.xlabel('x')
plt.ylabel('wave position')
plt.title('Wave Equation Solution - Explicit Method')
plt.legend()
plt.grid(True)
plt.show()


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
plt.show()

# img_matrix is now a NumPy array representing the image
print("Image shape:", img_matrix.shape)

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
