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

def central_diff_derivative(x_known, y_known, k, n):
    
    h = x_known[1] - x_known[0]  # Uniform step size
    D = 0  # Final result

    # For odd k, use spacing of 2h
    if k % 2 != 0:
        # Central difference over 2h spacing
        half_stencil = (k + 1) // 2  # How far out to go
        if n - half_stencil < 0 or n + half_stencil >= len(y_known):
            raise IndexError("Not enough points for central difference at this index.")

        for i in range(-half_stencil, half_stencil + 1):
            coeff = central_diff_coeff(k, i * 2)  # Spacing = 2h
            D += coeff * y_known[n + i]

        D /= (2 * h)**k

    else:
        # Even-order derivative with 1h spacing
        half_stencil = k // 2
        if n - half_stencil < 0 or n + half_stencil >= len(y_known):
            raise IndexError("Not enough points for central difference at this index.")

        for i in range(-half_stencil, half_stencil + 1):
            coeff = central_diff_coeff(k, i)
            D += coeff * y_known[n + i]

        D /= h**k

    return D



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

