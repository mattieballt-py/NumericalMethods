# libraries:
import numpy as np

"""
Golden Search Method
"""
# want to find maximum of a function using only one func eval per iteration

# the function we want to search for a max in
def df(x):
    return -2*x*np.exp((-x**2))

# search interval
xl = 0 # lower boundary
xu = 10 # upper boundary

def goldsearch(xl,xu,i):
    xl = xl
    xu = xu
    fxl = df(xl) 
    fxu = df(xu)
    if fxl * fxu > 0:
        return print("choose a new sub interval")
    R = (np.sqrt(5) - 1) / 2  # Golden ratio
    x1 = xl + R * (xu - xl)
    x2 = xu - R * (xu - xl)

    if abs(xu - xl) < 0.01:  # If the interval is small enough, return result
        return (xl + xu) / 2
    else:
        x1 = xl + d
        x2 = xu - d
        fx1 = df(x1)
        fx2 = df(x2)
        if fx2 < fx1:
            i += 1
            print(i)
            return goldsearch(xl, x1,i)  # Shrink the interval towards x2
        else:
            i += 1
            print(i)
            return goldsearch(x2, xu,i)  # Shrink the interval towards x1



print(goldsearch(0,1,0))

