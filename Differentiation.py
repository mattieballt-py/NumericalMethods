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
    
    h = x_known[1] - x_known[0]  # Assume uniform spacing
    D = 0  # To accumulate result
    
    # Check: do we have enough data points around index `n`?
    half_k = k // 2
    if n - half_k < 0 or n + half_k >= len(y_known):
        raise IndexError("Not enough points for central difference at this index.")
    
    # Compute central difference formula (based on binomial coefficients)
    for i in range(0, k + 1):
        sign = (-1)**i
        coeff = choose(k, i)
        index = int(n + half_k - i)  # Shift around center point
        D += sign * coeff * y_known[index]
    
    return D / (h**k)



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

