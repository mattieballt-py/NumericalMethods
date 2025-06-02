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


"""
Two Var FWD Euler
"""
# if u have initial value problem d^2y/dx^2 + dy/dxy + x = 0
# make w = dy/dx to have two vars, w and dw/dx
# then:

def func1(x,y):
    f = y[1] # w = dy/dt
    return f

def func2(x,y):
    f = -5*x*y[1]-(x+7)*np.sin(x)
    return f

# set initial conditions:
y0 = np.ndarray(2)
y0[0] = 4 # initial y
y0[1]=3 # initial dy/dx

def FwEulerTwo(Y0,t0,tend,h):
 # compose nodal times
 t = np.arange(t0,tend+h,h)
 # determine the number of time steps
 N = len(t)
 # allocate output array
 Y = np.ndarray((2,N))
 # initialise the solution
 t[0] = t0
 Y[0,0] = Y0[0]
 Y[1,0] = Y0[1]
 # compute the solution incrementally at subsequent time steps
 for n in range(1,N):
    Y[0,n] = Y[0,n-1] + func1(t[n-1],Y[:,n-1]) * h
    Y[1,n] = Y[1,n-1] + func2(t[n-1],Y[:,n-1]) * h
 return (t,Y)

(x,y) = FwEulerTwo(y0,0,15,0.02)


"""
Heat Diffusion PDE
"""

# this specific solution is from the ice cream Q, adjust paramaters, bc's etc as needed:
# spatial domain (in mm)
R = 30
dx = 0.5 # spatial increment
# greate the spatial grid points
x = np.arange(-R,R+dx,dx)
Nx = len(x)

# temporal domain
tend = 3600
dt = 1 # temporal increment
# create the temporal grid points
t = np.arange(0,tend+dt,dt)
Nt = len(t)

# create the solution matrix
T = np.ndarray((Nt,Nx))

# set the physics
# boundary conditions (fixed temperature at boundaries)
Ta = 35
Tb = 35
# set the initial T for watermelon
T[0,:int(Nx/2)] = 2
# set the initial T for pistachio
T[0,int(Nx/2):] = 6

# set the thermal diffusivity
alpha = np.ndarray(len(x))
# set the diffusivity T for watermelon
alpha[:int(Nx/2)] = 4.5e-2
# set the diffusivity for pistachio
alpha[int(Nx/2):] = 8.6e-2
# ================================================

averagep = np.average(T[0,int(Nx/2):])
averagew = np.average(T[0,:int(Nx/2)])
avp = [averagep]
avw = [averagew]
Tmelt = 15

c = alpha * dt / dx**2

# compute the solution incrementally at subsequent time steps
k = 1
while averagep < Tmelt:
    # compute at time step p, i.e. t = k * dt
    # do it for every node in the spatial grid
    # start with the boundaries
    T[k,0] = Ta
    T[k,-1] = Tb
    # do the interior nodes
    for i in range(1,Nx-1):
        # apply the discretised equation
        T[k,i] = c[i] * ( T[k-1,i+1] + T[k-1,i-1] ) + (1 - 2*c[i]) * T[k-1,i]
    averagep = np.average( T[k,int(Nx/2):])
    averagew = np.average( T[k,:int(Nx/2)])
    avp += [averagep]
    avw += [averagew]
    k += 1


plt.plot(x,T[k-1,:])
plt.grid()
plt.show()


plt.plot(t[:k],avp,c='Green')
plt.plot(t[:k],avw,c='Red')
plt.grid()
plt.show()

print('Ice cream melted after minutes:')
print((k-1)*dt/60)

print('Courant:')
print(c[-1],c[0])