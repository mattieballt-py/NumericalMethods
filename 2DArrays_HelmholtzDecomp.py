import numpy as np
import matplotlib.pyplot as plt  
from mpl_toolkits import mplot3d

"""
Generate Multidimensional Arrays and Grids
"""
# define two multivariable functions:
def f(x,y):
    return np.sin(x) * np.cos(y)

def g(x,y):
    return np.sin(x) * np.cos(y)

# set up arrays for domain to evaluate the functions
dx = 0.1
x = np.arange(-2, 2, dx)
y = np.arange(-1, 2, dx)

Nx = len(x)
Ny = len(y)

Xg, Yg = np.meshgrid(x,y) # creates two meshgrids of x and y, each max NxN

# compute two functions of f and g
s_xy = f(Xg,Yg) + g(Xg,Yg)
p_xy = f(Xg,Yg)*g(Xg,Yg)
print("sxy to compare",s_xy)
# surface plot functions on x,y domain 
ax = plt.axes(projection='3d') 
ax.plot_surface(Xg,Yg,s_xy) 
ax.plot_wireframe(Xg,Yg,s_xy)
ax.plot_surface(Xg,Yg,p_xy) 
ax.plot_wireframe(Xg,Yg,p_xy)
#plt.show() # this is plotting wireframe and surface of the two functions on the domain x by y


"""
plotting a 3 variable function (with time)
"""

# defining r(x,y,t):
def r(x,y,t):
    return f(x,y)*np.exp(-0.5*t) # not sure if it needs to be Xg or x

#adding array for time, t:
dt = 0.05
t = np.arange(0,10,dt)
Xg,Yg,Tg = np.meshgrid(x,y,t)

# calc the r(x,y,t):
r_xyt = r(Xg,Yg,t)
r_t0 = r_xyt[:,:,0]

# surface plot functions on x,y domain 
ax = plt.axes(projection='3d') 
ax.plot_surface(Xg,Yg,r_t0) 
ax.plot_wireframe(Xg,Yg,s_xy)
ax.plot_surface(Xg,Yg,p_xy) 
ax.plot_wireframe(Xg,Yg,p_xy)
#plt.show() # this is plotting wireframe and surface of the two functions on the domain x by y
