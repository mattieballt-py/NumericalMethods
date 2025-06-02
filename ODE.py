import numpy as np
import matplotlib.pyplot as plt

"""
First Order ODE FWD Euler
"""

# Define the function f(t, y) = dy/dt
def func(t, y):
    return -2 * y - 2 * np.exp(-t)

# Forward Euler method
def FwEuler(y0, t0, h, t_final):
    t_values = [t0]
    y_values = [y0]

    t = t0
    y = y0

    while t < t_final:
        y = y + h * func(t, y)
        t = t + h

        t_values.append(t)
        y_values.append(y)

    return np.array(t_values), np.array(y_values)

"""
Runge Kutter ODE 1st Order
"""
# Runge-Kutta 4th Order method
def ODERK4(y0, t0, h, t_final):
    t_values = [t0]
    y_values = [y0]

    t = t0
    y = y0

    while t < t_final:
        k1 = func(t, y)
        k2 = func(t + h/2, y + h * k1 / 2)
        k3 = func(t + h/2, y + h * k2 / 2)
        k4 = func(t + h, y + h * k3)

        y = y + (h / 6) * (k1 + 2*k2 + 2*k3 + k4)
        t = t + h

        t_values.append(t)
        y_values.append(y)

    return np.array(t_values), np.array(y_values)

"""
ODE 2nd Order Fwd Euler
"""

x0 = 0
x_end = 15
dx = 0.02
x = np.arange(x0, x_end + dx, dx)

# Initial conditions (use your CID digits as required)
y0 = 3      # y(0) = 3rd digit of CID
dy0 = 5     # y'(0) = 5th digit of CID

# Initialize arrays
y = np.zeros_like(x)
dy = np.zeros_like(x)

# Apply initial values
y[0] = y0
dy[0] = dy0

# Define your second-order ODE here:
# dy2 = f(x, y, dy)

def f(x, y, dy):
    return -5*x*dy - (x + 7)*np.sin(x)

for i in range(1, len(x)):
    y[i] = y[i-1] + dx * dy[i-1]
    dy[i] = dy[i-1] + dx * f(x[i-1], y[i-1], dy[i-1])


# Plotting
plt.figure(figsize=(10, 5))
plt.plot(x, y, label='y(x)', linewidth=2)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Forward Euler Solution to Second-Order ODE")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()



"""
BWD Euler ODE Solver
"""

# Define the function f(t, y) = dy/dt
def func(t, y):
    return -2 * y - 2 * np.exp(-t)

# Derivative of f w.r.t. y
def dfdy(t, y):
    return -2

# Newton-Raphson solver to find y_{n+1}
def newton_solver(f, dfdy, t_next, y_guess, y_prev, h, tol=1e-8, max_iter=50):
    y = y_guess
    for _ in range(max_iter):
        G = y - y_prev - h * f(t_next, y)
        dG = 1 - h * dfdy(t_next, y)
        y_new = y - G / dG
        if abs(y_new - y) < tol:
            return y_new
        y = y_new
    raise RuntimeError("Newton-Raphson did not converge")


# bisection method recursively
def bisecrecurs(a,b,tol,i):
    if f(a)*f(b)>0: # check valid interval
        print("choose new interval")
        return None

    n = (a + b)/2 #new point, bang in the middle of a and b
    print('midpoint', n)

    if f(a)*f(n)<0: # function changes sign between a and n
        # best guess before was n, now n + 1 = (a + n)/2
        # err = (new - old) /new
        errnew = (((a + n)/2) - n)/((a + n)/2)
        if abs(errnew) < tol:
            return ((a + n)/2) # if interval less than min, the root is at centre of it ish
        else:
            i += 1
            print(i)
        return bisecrecurs(a,n,tol,i)
  
    else:  # function changes sign between n and b
        errnew = (((n + b)/2) - n)/((n + b)/2)
        if abs(errnew) < tol:
            return ((n + b)/2) # if interval less than min, the root is at centre of it ish
        i += 1
        print(i)
        return bisecrecurs(n,b,tol,i)


print(bisecrecurs(0,2*0.3,0.001,1))

# Backward Euler method using Newton-Raphson
def BWDEuler(y0, t0, h, t_final):
    t_values = [t0]
    y_values = [y0]

    t = t0
    y = y0

    while t < t_final:
        t_next = t + h
        y_guess = y  # initial guess for Newton-Raphson
        y_next = newton_solver(func, dfdy, t_next, y_guess, y, h)
        t = t_next
        y = y_next

        t_values.append(t)
        y_values.append(y)

    return np.array(t_values), np.array(y_values)

# Example usage
t_values, y_values = BWDEuler(0, 0, 0.1, 10)

# Plotting (optional)
C = 0 - 1  # since y0 = 0, exact solution C = y0 - 1
y_exact = 1 - np.exp(-t_values) + C * np.exp(-2 * t_values)

plt.plot(t_values, y_values, 'bo-', label='Backward Euler')
plt.plot(t_values, y_exact, 'k--', label='Analytical Solution')
plt.xlabel("t")
plt.ylabel("y(t)")
plt.legend()
plt.grid(True)
plt.title("Backward Euler with Newton-Raphson")
plt.show()

"""
FWD Euler for system of equations
"""

import numpy as np
import matplotlib.pyplot as plt

# Define system of ODEs: dy/dt = f(t, y)
def func(t, Y):
    y0, y1, y2, y3 = Y  # unpack vector

    # dy0/dt = y1, dy1/dt = y2, dy2/dt = y3
    dy0 = y1
    dy1 = y2
    dy2 = y3
    dy3 = -2*y2 - y0 + 3*np.sin(t) - 5*np.cos(t)

    return np.array([dy0, dy1, dy2, dy3])

# Forward Euler Method for systems
def FwEuler(Y0, t0, h, t_final):
    t_values = [t0]
    y_values = [Y0]

    t = t0
    Y = Y0.copy()

    while t < t_final:
        Y = Y + h * func(t, Y)
        t = t + h

        t_values.append(t)
        y_values.append(Y.copy())

    return np.array(t_values), np.array(y_values)

# --- Initial conditions ---
y0_init = 5
y1_init = 8
y2_init = 3
y3_init = 10

Y0 = np.array([y0_init, y1_init, y2_init, y3_init])

# --- Time domain ---
t0 = 0
t_final = 15
h = 0.001

# --- Solve using Euler method ---
t_vals, y_vals = FwEuler(Y0, t0, h, t_final)

# --- Extract relevant components ---
y_t = y_vals[:, 0]         # y(t)
y2_t = y_vals[:, 2]        # d²y/dt²

# --- Plot y(t) ---
plt.figure(figsize=(10, 4))
plt.plot(t_vals, y_t, label='y(t)')
plt.xlabel('t')
plt.ylabel('y')
plt.title('Numerical Solution y(t) using Forward Euler')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

# --- Plot d²y/dt² ---
plt.figure(figsize=(10, 4))
plt.plot(t_vals, y2_t, label='d²y/dt²', color='orange')
plt.xlabel('t')
plt.ylabel('d²y/dt²')
plt.title('Second Derivative of y(t)')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
