# CID: 02413400
# Fourier Transform Task A

import numpy as np
import matplotlib.pyplot as plt
"""
# (a) Plot each component y_k(t)
t = np.arange(0, 2*np.pi, 0.01)
N = len(t)

A_k = [5, 2, 0, 3, 1, 8, 2, 3, 7, 6]  # amplitudes
phi_k = [np.pi/2, np.pi, 0, 0, 3*np.pi/2, np.pi, 0, 3*np.pi/2, np.pi, 0]  # phases

components = []

plt.figure(figsize=(10, 6))
for k in range(1, 11):
    y_k = A_k[k-1] * np.sin(k * t + phi_k[k-1])
    components.append(y_k)
    plt.plot(t, y_k, label=f"$y_{k}(t)$")

plt.title("(a) Components $y_k(t)$")
plt.xlabel("t")
plt.ylabel("Amplitude")
plt.grid(True)
plt.legend(fontsize='small')
plt.tight_layout()
plt.show()

# (b) Plot total signal y(t)
y_total = sum(components)

plt.figure(figsize=(8, 4))
plt.plot(t, y_total, color='black')
plt.title("(b) Total Signal $y(t)$")
plt.xlabel("t")
plt.ylabel("y(t)")
plt.grid(True)
plt.tight_layout()
plt.show()

# (c) Compute DFT and frequency data

def DFT(yn):
    # y: values of the function, in time domain
    N = len(yn)
    w = 2*np.pi/N
    FTk = np.zeros(N,dtype=complex)
    for k in range(0,N):
        for n in range(0,N):
            FTk[k] += np.exp(-1j*k*w*n)*yn[n]
    return FTk

Y = DFT(y_total)
mag = np.abs(Y)
phase = np.angle(Y)
freq = np.fft.fftfreq(N, d=0.01)
half = N // 2
Δf = freq[1] - freq[0]

# (c) Plot DFT Magnitude – bar chart
plt.figure(figsize=(8, 4))
plt.bar(freq[:half], mag[:half], width=Δf, color='blue', edgecolor='black')
plt.title("(c) DFT Magnitude Spectrum")
plt.xlabel("Frequency [Hz]")
plt.ylabel("|Y(f)|")
plt.grid(True)
plt.tight_layout()
plt.show()

# (c) Plot DFT Phase – bar chart
plt.figure(figsize=(8, 4))
plt.bar(freq[:half], phase[:half], width=Δf, color='orange', edgecolor='black')
plt.title("(c) DFT Phase Spectrum")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Phase [rad]")
plt.grid(True)
plt.tight_layout()
plt.show()

# (d) Print frequency resolution info
f_min = np.min(freq[:half])
f_max = np.max(freq[:half])
print(f"(d) Frequency step Δf = {Δf:.4f} Hz")
print(f"    Minimum frequency = {f_min:.2f} Hz")
print(f"    Maximum frequency = {f_max:.2f} Hz")
"""
# I hate Task A

# Task C

with open('/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/txt_Data/Temperatures.txt', 'r') as d:
    T = d.readlines()

# Initialize lists to store cleaned float values
Temps_known = []
for item in T:
    term = item.strip()
    try:
        clean_term = float(term)
        Temps_known.append(round(clean_term))
    except ValueError:
        continue

Temp_known = np.array(Temps_known)
print("Temps", Temp_known)

# Corresponding time array
t_known = np.arange(0, 6, 0.5)
print("time", t_known)

# Interpolation using Lagrangian method
def Langrangian(x_known, y_known, xp):
    n = len(x_known)
    P = 0
    for i in range(n):
        Li = 1
        for j in range(n):
            if i != j:
                Li *= (xp - x_known[j]) / (x_known[i] - x_known[j])
        P += y_known[i] * Li
    return P

def InterpPoints(x_known, y_known, InterpMethod, xp_points):
    yp_points = xp_points.copy()
    for i in range(len(xp_points)):
        yp_points[i] = InterpMethod(x_known, y_known, xp_points[i])
    return yp_points

# Bracketing interval between t = 3 and t = 3.5
t_interp = np.arange(3, 3.5, 0.005)
Temp_interp = InterpPoints(t_known, Temp_known, Langrangian, t_interp)

# Plotting interpolated section
plt.figure(figsize=(8, 4))
plt.scatter(t_interp, Temp_interp, color='black')
plt.title("Interpolated points for root finding")
plt.xlabel("t")
plt.ylabel("T(t)")
plt.grid(True)
plt.tight_layout()
plt.show()

# Bisection method for langrangian interpolated points
def bisecrecurs(a, b, tol, i):
    fa = Langrangian(t_known, Temp_known, a)
    fb = Langrangian(t_known, Temp_known, b)

    if fa * fb > 0:
        print("No sign change — choose a new interval.")
        return None

    mid = (a + b) / 2
    fmid = Langrangian(t_known, Temp_known, mid)
    print(f"Iteration {i}, midpoint t = {mid:.5f}, T = {fmid:.5f}")

    if abs(fmid) < tol:
        print(f"Root found at t = {mid:.5f}")
        return mid

    if fa * fmid < 0:
        return bisecrecurs(a, mid, tol, i + 1)
    else:
        return bisecrecurs(mid, b, tol, i + 1)

# Call fixed bisection function over time interval [3, 3.5]
root_time = bisecrecurs(3.0, 3.5, 0.005, 0)
print(f"Estimated root at time t = {root_time:.4f} s")
