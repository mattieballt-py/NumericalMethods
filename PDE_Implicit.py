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
