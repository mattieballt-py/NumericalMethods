# CID: 02413400

# importing libraries:
import numpy as np
import matplotlib.pyplot as plt  
import math as mt

"""
Task A: Signals and fourier transform
"""

t = np.arange(0, 2*mt.pi,0.01)

yts = [5*np.sin(t*mt.pi/2), 2*np.sin(t*mt.pi), 0*np.sin(t*mt.pi/2), 3*np.sin(t), 1*np.sin(3*t*mt.pi/2), 8*np.sin(t*mt.pi), 2*np.sin(t*0), 3*np.sin(3*t*mt.pi/2), 7*np.sin(t*mt.pi), 6*np.sin(t*0)]
for i in range(0,len(yts)):
    pl1 = plt.scatter(t,yts[i],s=2)
plt.title('Plot 1, components of y(t)')
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()

# task b) plot sum of all components
yt = 5*np.sin(t*mt.pi/2) + 2*np.sin(t*mt.pi) + 0*np.sin(t*mt.pi/2) + 3*np.sin(t) + 1*np.sin(3*t*mt.pi/2) + 8*np.sin(t*mt.pi) + 2*np.sin(t*0) + 3*np.sin(3*t*mt.pi/2) + 7*np.sin(t*mt.pi) + 6*np.sin(t*0)
pl2 = plt.scatter(t,yt,s=3)
plt.title('Plot 2, yt sum of components')
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()

# c) discrete fourier transform: