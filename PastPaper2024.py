# CID: 02413400
# Fourier Transform Task A

import numpy as np
import matplotlib.pyplot as plt

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

