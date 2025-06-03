import numpy as np
import matplotlib.pyplot as plt 

"""
Implicit crank nicholson, solution of Diffusion Eqn non dimensionalised
"""
# implicit solution of d^2T/dx^2 - dT/dt = 0  non dimensionalised

h = 0.25 # space step
k = 0.5 # time step so r = 1
r = k/h**2 # max it can be while converging = 0.5 for explicit
#print("convergence factor r",r)
print(r)

x = np.arange(0,1+h,h) # mesh in x
t = np.arange(0,1+k,k) # mesh in time
T = np.zeros((len(t),len(x))) # initialise Temp array


# BC's:
T[0,:] = 1

# Matrix A (left-hand side) for interior points only (x1, x2, x3)
N = len(x) - 2  # number of interior points
A = np.zeros((N, N))

# Fill A matrix: tridiagonal with 4 on the diagonal and -1 on off-diagonals
for i in range(N):
    A[i, i] = 4
    if i > 0:
        A[i, i - 1] = -1
    if i < N - 1:
        A[i, i + 1] = -1

def crankynick(x,t,T,r,A): # with BC's set up as tutorial question
    N = len(x) - 2  # interior points
    for n in range(len(t) - 1):  # time loop
        b = np.zeros((N,))

        # Apply Neumann BC: T[n,0] = T[n,1]
        T[n, 0] = T[n, 1]
        # Dirichlet BC: T[n,-1] = 0
        T[n, -1] = 0

        # Construct RHS vector b
        for j in range(1, N + 1):  # j = 1 to N (interior nodes)
            b[j - 1] = T[n, j - 1] + 2 * T[n, j] + T[n, j + 1]

        # Solve the system A u_{n+1} = b
        u_new = np.linalg.solve(A, b)

        # Update T at next time step
        T[n + 1, 1:-1] = u_new

        # Re-apply BCs to next time step
        T[n + 1, 0] = T[n + 1, 1]  # Neumann
        T[n + 1, -1] = 0  # Dirichlet

    return T

# Run the method
T = crankynick(x,t,T,r,A)

# Plot results
plt.figure(figsize=(10, 6))
for n in [0, int(len(t) / 4), int(len(t) / 2), len(t) - 1]:
    plt.plot(x, T[n, :], label=f't = {t[n]:.2f}')
plt.xlabel('x')
plt.ylabel('Temperature')
plt.title('Crank-Nicolson Solution for Heat Equation')
plt.legend()
plt.grid(True)
plt.show()



"""
BWD Euler Implicit
"""
#Example question: 
    # d^2y/dt^2 + c/m (dy/dt) + (g/L) sin(y) = 0

"""To alter: the equation solve for yn"""
def func(yn,ynp,ynp2,tn,h):
    """Put all terms on one side = f"""
    f = (yn-2*ynp+ynp2)/h**2+c/m*(yn-ynp)/h+g/L*np.sin(yn)
    #f = yn*(1+h**2*c/m) + h**2*g/L*np.sin(yn) - 2*ynp + ynp2 
    return f

def mybisection(a,b,ynp,ynp2,tn,h,eps):
    """To change: any input parameters to mybisection"""
    # repeat the split of teh interval until the bracketing intervla becomes smaller than the accuracy
    while abs(a-b)>eps:
        # calculate the mid point
        xm = (a + b) / 2
        # establish in which subinterval the solution lies
        # compute f(a) * f(xm)
        ff = func(a,ynp,ynp2,tn,h) * func(xm,ynp,ynp2,tn,h)
        if ff < 0: 
            # the solution lies in the left interval
            # set the upper bracket as xm
            b = xm
        else:
            # the solution lies in the right interval
            # set the lower bracket as xm
            a = xm
            
    # the true solution is bracketed within the latest interval [a,b]
    # we can approximate it with the midpoint
    sol = (a + b) / 2
    
    return sol

"""Define constants and domain"""
g = 9.8 #Gravity ms-2
L = 2 #Length m
m = 10 #Mass kg
c = 0.001 #damping Ns/m

dt = 0.1 #time step
tend = 20 #end time
y0 = 45.0 #initial angle in degrees
vinit = 0 # initial velocity
t = np.arange(0,tend,dt)
y = np.ndarray(len(t))
y[0] = y0*np.pi/180 # initial condition on position
y[1] = vinit*dt+y[0] # intial condition on velocity
ynp = y[1] #y_n-1
ynp2 = y[0] #y_n-2

# iterate from the third time step
for i in range(2,len(t)):
    # determine the solution at this time step
    """To change, input parameters to "mybisection"""
    y[i] = mybisection(-np.pi/2,np.pi/2,ynp,ynp2,t[i],dt,0.001) # the angle of solution will not exceed -pi/2 +pi/2
    ynp = y[i] # make the current solution as the solution at previous time step
    ynp2 = y[i-1] # update the solution at two timesteps ago
    
plt.plot(t,y*180/np.pi,'b')
plt.grid()
plt.show()

# determine the velocity
v = np.ndarray(len(t)-1)
v = (y[1:] - y[:-1])/dt
plt.plot(t[:-1],v,'b')
plt.grid()
plt.show()