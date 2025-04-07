import numpy as np 
import matplotlib.pyplot as plt 


"""
Integration with Guass Quadrature
"""

def f(t):
    A = 0
    B = 1
    x = 0.5*(A * (1-t) + B*(t + 1))
    return x * np.exp(x**2)

def GuassQuad(f,n): # n is degree of precision
    nodes = [] # = tj
    weights = [] # = wj
    I = 0 # just to init

    if n == 7:
        # define the nodes, tj = nodes
        nodes.append(((3/7) + (2/7)*(6/5)**0.5)**0.5)
        nodes.append(-1 * ((3/7) + (2/7)*(6/5)**0.5)**0.5)
        nodes.append(((3/7) - (2/7)*(6/5)**0.5)**0.5)
        nodes.append(-1 * ((3/7) - (2/7)*(6/5)**0.5)**0.5)

        # define the weights:
        weights.append((18 - 30**0.5)/36)
        weights.append((18 - 30**0.5)/36)
        weights.append((18 + 30**0.5)/36)
        weights.append((18 + 30**0.5)/36)
        print(nodes, weights)
        # run the integration
        for j in range(0,4):
            I += weights[j]*f(nodes[j])
        I *= 0.5
        return I
    else:
        return

print("Guass solution:", GuassQuad(f,7))
"""
Integration with Simpson's 1/3 rule
"""


