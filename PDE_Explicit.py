import numpy as np
import matplotlib.pyplot as plt 

"""
Explicit solution Diffusion Eqn Non-Dimensionalised
"""
# explicit solution of d^2T/dx^2 - dT/dt = 0  non dimensionalised

h = 0.25
r = 0.5
k = r * h**2

x = np.arange(0,1,h) # mesh in x
t = np.arange(0,1,k) # mesh in time
T = np.zeros((len(t),len(x)))

def explicit(x,t,T):
    for k in range(0,len(t)): # so at each time step
        # BC's:
        T[0,:] = 1
        T [:,-1] = 0
        for h in range(0,len(x)): # from LHS to RHS stepping along x
            # Nueman BC:
            T[k,h+1] = T[k,h-1]
            T[k,h] = 0.5 * T[k-1,h+1] + 0.5*T[k-1,h-1]
    return T



