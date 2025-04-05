import numpy as np
import matplotlib.pyplot as plt 

"""
Explicit solution Diffusion Eqn Non-Dimensionalised
"""
# explicit solution of d^2T/dx^2 - dT/dt = 0  non dimensionalised

h = 0.25 # space step
k = 0.05 # time step
r = k/h**2 # max it can be while converging
#print("convergence factor r",r)

x = np.arange(0,1+h,h) # mesh in x
t = np.arange(0,1+k,k) # mesh in time
T = np.zeros((len(t),len(x))) # initialise Temp arra
# BC's:
T[0,:] = 1


def explicit(x,t,T,r): # with BC's set up as tutorial question
    
    for i in range(0,len(t)-1): # so at each time step
        # BC's:
        T[i,0] = T[i,1]
        T [i,-1] = 0    
        for j in range(1,len(x)-1): # from LHS to RHS stepping along x
            # Nueman BC:
            T[i+1,j] = r * T[i,j+1] + r*T[i,j-1]

        # reapply BC's to next step:
        T[i+1,0] = T[i+1,1] 
        T[i+1,-1] = 0
    return T

diff = explicit(x,t,T,r)


# Plot results
plt.figure(figsize=(10, 6))
for n in [0, int(len(t)/4), int(len(t)/2), len(t)-1]:
    plt.plot(x, diff[n, :], label=f't = {t[n]:.2f}')
plt.xlabel('x')
plt.ylabel('Temperature')
plt.title('Heat Equation Solution - Explicit Method')
plt.legend()
plt.grid(True)
#plt.show()