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

