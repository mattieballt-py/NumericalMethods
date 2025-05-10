import numpy as np
import matplotlib.pyplot as plt 

"""
Discrete Fourier Series
"""

# Question 4: Discrete Fourier Transform at f = 0.3684
N = 1024  # Sample points
t = np.linspace(0, 1, N, endpoint=False)
y = np.sin(2 * np.pi * t) + np.sin(6 * np.pi * t)
Y = fft(y)
frequencies = fftfreq(N, d=1/N)
idx = np.argmin(np.abs(frequencies - 0.3684))
value_at_f = Y[idx].real  # Taking real part for amplitude
#print("Question 4 Answer:", round(value_at_f, 2))

# Question 5: Filtering a Noisy Signal
def low_pass_filter(signal, cutoff_freq, sampling_rate):
    freqs = np.fft.fftfreq(len(signal), 1/sampling_rate)
    spectrum = np.fft.fft(signal)
    spectrum[np.abs(freqs) > cutoff_freq] = 0
    return np.fft.ifft(spectrum).real

sampling_rate = 100
cutoff_freq = 0.5
time = np.linspace(0, 10, sampling_rate * 10)
noisy_signal = np.sin(2 * np.pi * 1.0 * time) + 0.5 * np.random.randn(len(time))  # Simulated noisy signal
clean_signal = low_pass_filter(noisy_signal, cutoff_freq, sampling_rate)
clean_value_t5 = np.interp(5, time, clean_signal)
#print("Question 5 Answer:", round(clean_value_t5, 2))



# for running scripts you're unsure will work:
# try:
    # print something?
# except:
    # something else it runs

"""
Task B: Anologue Filters and Bode Plots
"""
inc = np.arange(1,5) # 1,2,3,4
x = np.log10(10**inc)
#print(x)

#basic H and phase
w = np.arange(0,100000,1)
H = 1/(1 + 1j* 0.1*w)
Hamp = np.abs(H)
HampdB = 20*np.log10(Hamp)
Hphase = np.arctan2(H.imag,H.real)*180/np.pi
    
#plottign the basic H and phase

"""
plt.figure()
plt.semilogx(w, H)  # Logarithmic x-axis, amplitide vs w
plt.xlabel('w')
plt.ylabel('|H(w)|')
plt.title('Magnitude of H(w)')
plt.grid()
plt.show()

plt.figure()
plt.semilogx(w, HampdB)  # Logarithmic x-axis, amplitide dB vs w
plt.xlabel('w')
plt.ylabel('|H(w)| dB')
plt.title('Magnitude of H(w) DdB')
plt.grid()
plt.show()
"""
# plt.figure()
# plt.semilogx(w, Hphase)  # Logarithmic x-axis, phase vs w
# plt.xlabel('w')
# plt.ylabel('phase')
# plt.title('phase')
# plt.grid()
#plt.show()

#ii)
# impedance of a resitor in parallel with a capaior- Zeq = ZrZc/(Zr+Zc)
#Zr = R
#Zc = 1/jwc (w = 2pif)

R1 = 1000
R2 = 2000
C1 = 1E-3
C2 = 2E-3

w = np.arange(0.001,100,0.01)
ZR1 = R1
ZR2 = R2
ZC1 = 1/(1j*w*C1)
ZC2 = 1/(1j*w*C2)

Zh =(ZR1*ZC1)/(ZR1 + ZC1)    #Impedance 1st resistor and capcitor
Zv = (ZR2*ZC2)/(ZR2 + ZC2)   #Impedance 2nd resistor and capcitor
H2= 1/(1 + Zh/Zv)            #H = v0/v1 = zv/(zv+zh)
Hamp2 = np.abs(H2)
HampdB2 = 20*np.log10(Hamp2)
Hphase2 = np.arctan2(H2.imag,H2.real)*180/np.pi

# for w = 0.5
w05 = 0.5
ZC105 = 1/(1j*w05*C1)
ZC205 = 1/(1j*w05*C2)
Zh05 =(ZR1*ZC105)/(ZR1 + ZC105)    #Impedance 1st resistor and capcitor
Zv05 = (ZR2*ZC205)/(ZR2 + ZC205) 
H205= 1/(1 + Zh05/Zv05)  
Hamp205 = np.abs(H205)
HampdB205 = 20*np.log10(Hamp205)
#print('this answer q2', H205, Hamp205, HampdB205)

"""
plt.figure()
plt.semilogx(w, HampdB2)  # Logarithmic x-axis, amplitide dB vs w
plt.xlabel('w')
plt.ylabel('|H(w)| dB')
plt.title('Magnitude of H(w) DdB')
plt.grid()
plt.show()
"""
# plt.figure()
# plt.semilogx(w, Hphase2)  # Logarithmic x-axis, phase vs w
# plt.xlabel('w')
# plt.ylabel('phase')
# plt.title('phase')
# plt.grid()
# plt.show()


#iii)
#iv)
w = np.arange(0.0001,10000,0.01)
Hlp = 1/(1  + 1j*0.01*w) #low pass filter
Hhp = (1j*40*w)/(1 + 1j*40*w) #high pass filter
#these numbers were from c1=20mF, c2 = 10uF, R2 = 1000, R1 = 2000
#cascade means its yimes together

Hcasc = Hlp*Hhp
HampCasc = np.abs(Hcasc)
HampdBCasc = 20*np.log10(HampCasc)
HphaseCasc = np.arctan2(Hcasc.imag,Hcasc.real)*180/np.pi

"""
plt.figure()
plt.semilogx(w, HampdBCasc)  # Logarithmic x-axis, amplitide dB vs w
plt.xlabel('w')
plt.ylabel('|H(w)| dB')
plt.title('Magnitude of H(w) casc DdB')
plt.grid()
plt.show()

plt.figure()
plt.semilogx(w, HphaseCasc)  # Logarithmic x-axis, phase vs w
plt.xlabel('w')
plt.ylabel('phase casc')
plt.title('phase cas')
plt.grid()
plt.show()
"""

"""
Task C: Fourier Series
"""

def FourierSqr(t, N, T):
    sum_terms = np.zeros_like(t)  # Initialize sum array
    
    for n in range(1, N+1, 2):  # Only odd values of n
        sum_terms += (np.sin((2 * n * np.pi * t) / T)) / n

    y = (4 / np.pi) * sum_terms  # Apply Fourier scaling factor
    return y

# Define time range
T = 5
t = np.linspace(0, 2*T, 1000)  # Smooth curve with 1000 points

# Evaluate for N = [2,6,50]
N_values = [2, 6, 50]
plt.figure(figsize=(10, 5))

for N in N_values:
    y = FourierSqr(t, N, T)
    plt.plot(t, y, label=f'N={N}')

"""
plt.xlabel('Time (t)')
plt.ylabel('y(t)')
plt.title('Fourier Series Approximation of Square Wave')
plt.legend()
plt.grid()
plt.show()
"""

"""
Task D: Discrete Fourier Transform
"""


def DFourTrans(yn):
    N = len(yn)
    DFT = np.zeros(N, dtype=complex)  # Store complex results

    for k in range(N):
        for n in range(N):
            DFT[k] += yn[n] * np.exp(-2j * np.pi * k * n / N)
    return DFT


# Example Input
yn = np.arange(10)

# Compute DFT
DFT_result = DFourTrans(yn)
#print(DFT_result)


#ii 

def DFTinv(DFT_result):
    N = len(DFT_result)
    y = np.zeros(N, dtype=complex)  # Store complex results

    for n in range(N):
        for k in range(N):
            y[n] += (DFT_result[k]) * np.exp((2j * np.pi * k * n )/ N)
    y = y/N
    return y.real

# testing
# iii
tn = np.arange(0,6*mt.pi,1)
yn = np.exp(-((tn-5)**2)/4)
yaxis = DFourTrans(yn)
yaxis = np.abs(yaxis)

# Frequency axis
freqs = np.fft.fftfreq(len(tn), d=1)  # d is the time step

"""
# Plot results
plt.figure(figsize=(8,4))
plt.plot(freqs, yaxis)
plt.title("Discrete Fourier Transform of Gaussian Function")
plt.xlabel("Frequency")
plt.ylabel("Magnitude")
plt.grid(True)
plt.show()
"""

"""
Task E: Signal Processing
"""
from scipy.signal import butter, filtfilt

with open('Vibration.txt', 'r') as d:
    Vib = d.readlines()
with open('Noisy.txt', 'r') as d:
    Nos = d.readlines()

# Initialize lists to store cleaned float values
Vibration = []
Noisy = []

# Process each line in the file
for item in Vib:
    # Strip leading/trailing whitespace
    term = item.strip()
    Vibration += [float(term)]

# Process each line in the file
for item in Nos:
    # Strip leading/trailing whitespace
    term = item.strip()
    # Try to convert the term to a float
    Noisy += [float(term)]

dt = 0.01  # Sampling interval (seconds)
n = len(Vibration)  # Number of samples

# Time and frequency vectors
tn = np.arange(0, n * dt, dt)
freqs = np.fft.fftfreq(n, d=dt)  # Correct frequency calculation

# Compute the Fourier Transform
DFourTrans = np.fft.fft(Vibration)
yaxis = np.abs(DFourTrans)  # Magnitude of the FFT

# Only take the positive frequencies (first half of data)
half_n = n // 2
plt.figure(figsize=(8, 4))
plt.plot(freqs[:half_n], yaxis[:half_n])  # Ensure same dimensions

plt.title("Discrete Fourier Transform of Vibrating Beam Response")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.grid(True)
plt.show()


# Convert data to float array
Noisy = np.array([float(item.strip()) for item in Nos])

# Given parameters
sampling_rate = 20  # samples per second
dt = 1 / sampling_rate  # time step
n = len(Noisy)  # Number of samples
t = np.arange(0, n * dt, dt)  # Time vector

# Define the low-pass filter (Butterworth)
def lowpass_filter(data, cutoff, fs, order=4):
    nyquist = 0.5 * fs  # Nyquist frequency
    normal_cutoff = cutoff / nyquist  # Normalize cutoff frequency
    b, a = butter(order, normal_cutoff, btype='low', analog=False)  # Butterworth filter
    filtered_data = filtfilt(b, a, data)  # Apply filter
    return filtered_data

# Apply filter with cutoff frequency of 0.5 Hz
cutoff_frequency = 0.5  # Hz
Cleaned_Signal = lowpass_filter(Noisy, cutoff_frequency, sampling_rate)

# Find the value of the clean signal at t = 5 sec
t_target = 5  # seconds
idx = np.argmin(np.abs(t - t_target))  # Find closest index to t = 5 sec
clean_value_at_5s = Cleaned_Signal[idx]

# Display the result
print(f"Clean signal value at t = 5 sec: {clean_value_at_5s:.4f}")

# Plot the original and filtered signals
plt.figure(figsize=(10, 5))
plt.plot(t, Noisy, label="Noisy Signal", alpha=0.5)
plt.plot(t, Cleaned_Signal, label="Filtered Signal (0.5 Hz cutoff)", linewidth=2)
plt.axvline(x=5, color='r', linestyle='--', label="t = 5 sec")
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Signal Amplitude")
plt.title("Noisy vs. Filtered Signal")
plt.grid(True)
plt.show()


"""
Own, active low pass filter
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

# Given circuit parameters (modify these as needed)
R1 = 1e3  # 1kΩ
R2 = 2e3  # 2kΩ (Given R2/R1 = 2)
C = 10e-6 # 10µF

# Define cutoff frequency
fc = 1 / (2 * np.pi * R2 * C)  # Cutoff frequency (Hz)
wc = 1 / (R2 * C)
print("cut off frequency: ", fc, "or in w" ,wc)

# Define frequency range (log scale)
frequencies = np.logspace(0, 5, 500)  # From 1Hz to 100kHz
w = 2 * np.pi * frequencies  # Convert to rad/s

# Compute magnitude response
H_mag = (R2 / R1) / np.sqrt(1 + (w * R2 * C) ** 2)
H_mag_db = 20 * np.log10(H_mag)  # Convert to dB

# Compute phase response
H_phase = 180 - np.arctan(w * R2 * C) * (180 / np.pi)  # Convert to degrees

# Plot Bode magnitude response
plt.figure(figsize=(10, 5))
plt.subplot(2, 1, 1)
plt.semilogx(frequencies, H_mag_db, label="Magnitude (dB)")
plt.axvline(fc, color='r', linestyle="--", label=f"Cutoff = {fc:.2f} Hz")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Gain (dB)")
plt.title("Bode Plot: Magnitude Response")
plt.legend()
plt.grid()

# Plot Bode phase response
plt.subplot(2, 1, 2)
plt.semilogx(frequencies, H_phase, label="Phase (degrees)", color='g')
plt.axvline(fc, color='r', linestyle="--")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Phase (degrees)")
plt.title("Bode Plot: Phase Response")
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()
