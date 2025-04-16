import numpy as np
import matplotlib.pyplot as plt  

"""
Lagrangian Interpolation
"""
# want to find f(xp) at a known value xp from a set of x_knowns

def f(x):
    return x**2 + np.cos(x) + np.exp(-3*np.sin(x))

x_known = np.linspace(1,2,3)
y_known = f(x_known)

def Langrangian(x_known, y_known, xp):
    n = len(x_known) # the degree of P polynomial
    print("n",n)
    P = 0 # init this is the Langrangian Polynomial
    Li = 1
    for i in range(0,n):
        print("i",i)
        for j in range(0,n):
            if i != j:
                print("j",j)
                Li *= (xp - x_known[j])/(x_known[i]-x_known[j])
        P += y_known[i] * Li
        Li = 1
    return P

def InterpPoints(x_known,y_known, InterpMethod, xp_points):
    yp_points = xp_points.copy()
    for i in range(0,len(xp_points)):
        yp_points[i] = InterpMethod(x_known,y_known,xp_points[i])
    return yp_points

xp = np.arange(0,3,0.05)
yp = InterpPoints(x_known,y_known,Langrangian,xp)

# for the plot
# same intervals, different degree of polynomial interp
xn_2 = np.linspace(1,2,2) # n = 2
yn_2 = f(xn_2)
xn_3 = np.linspace(1,2,3) # n = 3
yn_3 = f(xn_3)
xn_4 = np.linspace(1,2,4) # n = 4
yn_4 = f(xn_4)
xn_7 = np.linspace(1,2,7) # n = 7
yn_7 = f(xn_7)
xn_12 = np.linspace(1,2,12) # n = 7
yn_12 = f(xn_12)


yp_3 = InterpPoints(xn_3,yn_3,Langrangian,xp)
yp_4 = InterpPoints(xn_4,yn_4,Langrangian,xp)
yp_7 = InterpPoints(xn_7,yn_7,Langrangian,xp)
yp_12 = InterpPoints(xn_12,yn_12,Langrangian,xp)
plt.scatter(xn_2,yn_2, color="red", label="known")
plt.scatter(xn_3,yn_3, color="orange", label="known")
plt.scatter(xn_7,yn_7, color="orange", label="known")
plt.scatter(xp,yp_3, color="blue", label = "n = 3")
plt.scatter(xp,yp_4, color="turquoise", label = "n = 4")
plt.scatter(xp,yp_7, color="deepskyblue", label = "n = 7")
plt.scatter(xp,yp_12, color="indigo", label = "n = 12")
plt.legend()
plt.grid(True)
plt.show()


"""
Barycentric Interpolation
"""






"""
2D Interpolation
"""

# this 2d Interpolation done on images in playing