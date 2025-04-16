import numpy as np
import matplotlib.pyplot as plt  


def choose(a, b):
    """
    Computes the binomial coefficient "a choose b" = a! / (b! * (a - b)!)
    """
    # Edge case: choose(a, 0) = choose(a, a) = 1
    if b == 0 or b == a:
        return 1

    # Compute factorials safely
    def factorial(n):
        result = 1
        for i in range(1, n+1):
            result *= i
        return result

    return factorial(a) / (factorial(b) * factorial(a - b))

"""
K-th Order Derivative
"""

# numerical result of d/dx eg
# function that returns central diff kth derivative


def central_diff_derivative_general(x_known, y_known, k, n, step_multiple=1):
    """
    Compute the k-th derivative at point index n using central difference with
    spacing of h * step_multiple (e.g., 2h for step_multiple=2)
    """
    h = x_known[1] - x_known[0]
    h_step = h * step_multiple

    # Number of stencil points needed (minimum: k+1, better: odd number around center)
    stencil_size = k + 1 if k % 2 == 0 else k + 2
    if stencil_size % 2 == 0:
        stencil_size += 1  # Make it odd for central symmetry

    half = stencil_size // 2

    # Create stencil offsets (e.g., [-2, -1, 0, 1, 2] for step=1, or [-4, -2, 0, 2, 4] for step=2)
    offsets = [i * step_multiple for i in range(-half, half + 1)]
    
    # Check bounds
    idxs = [n + i for i in range(-half, half + 1)]
    if min(idxs) < 0 or max(idxs) >= len(y_known):
        raise IndexError("Not enough points for this stencil.")

    # Build Vandermonde system to solve for coefficients
    A = np.array([[offset**i for offset in offsets] for i in range(stencil_size)])
    b = np.zeros(stencil_size)
    b[k] = np.math.factorial(k)

    coeffs = np.linalg.solve(A, b)

    # Apply to data
    D = sum(c * y_known[n + i] for c, i in zip(coeffs, range(-half, half + 1)))
    return D / (h_step**k)




def forward_diff_derivative(x_known, y_known, k, n):
    h = x_known[1] - x_known[0]  # Step size
    D = 0  # Initialize result
    
    # Check if we have enough forward points to compute kth derivative
    if n + k >= len(y_known):
        raise IndexError("Not enough forward points to compute this derivative.")
    
    # Forward difference formula:
    # D ≈ (1/h^k) * sum_{i=0}^{k} (-1)^i * C(k, i) * f(x + i*h)
    for i in range(0, k + 1):
        sign = (-1)**i
        coeff = choose(k, i)
        D += sign * coeff * y_known[n + i]
    
    return D / (h**k)

# Example use:
x_known = [0, 1, 2]
y_known = [0, 1, 4]  # f(x) = x^2, so f'(x) = 2x, f''(x) = 2

# Compute first derivative at x = 1
print("Central 1st derivative at x=1:", central_diff_derivative(x_known, y_known, 1, 1))  # Should be close to 2
print("Forward 1st derivative at x=0:", forward_diff_derivative(x_known, y_known, 1, 0))  # Should be close to 1

