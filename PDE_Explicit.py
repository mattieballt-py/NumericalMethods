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


def explicitDiff(x,t,T,r): # with BC's set up as tutorial question
    
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

diff = explicitDiff(x,t,T,r)


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


"""
Explicit solution Wave Eqn Non-Dimensionalised
"""
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

