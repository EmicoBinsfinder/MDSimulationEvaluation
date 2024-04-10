import sys, argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot
from scipy import integrate
from scipy.constants import Boltzmann
from scipy.optimize import curve_fit
from tqdm import trange
from os.path import join as join

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define ACF
def acf_(data):
    steps = data.shape[0]
    size = steps // 2

    autocorrelation = np.zeros(size, dtype=float)
    for shift in trange(size, ncols=100, desc='Progress'):
            autocorrelation[shift] = np.mean( (data[:steps-shift]) * (data[shift:]) )

    return autocorrelation

# Define ACF using FFT
def acf(data):
    steps = data.shape[0]
    lag = steps // 2

    # Nearest size with power of 2 (for efficiency) to zero-pad the input data
    size = 2 ** np.ceil(np.log2(2 * steps - 1)).astype('int')

    # Compute the FFT
    FFT = np.fft.fft(data, size)

    # Get the power spectrum
    PWR = FFT.conjugate() * FFT

    # Calculate the auto-correlation from inverse FFT of the power spectrum
    COR = np.fft.ifft(PWR)[:steps].real

    autocorrelation = COR / np.arange(steps, 0, -1)

    return autocorrelation[:lag]

df = pd.read_csv('Stress_AVGOnespd_squalane_T313KP1atm.out')
# print(df.columns)
# df = df.drop(columns=['TimeStep'])
# print(df)

unit = 'atm' #Pressure, bar or temperature
datafile = 'Stress_AVGOnespd_squalane_T313KP1atm.out'
diag = False # Use the diagonal components for viscosity prediction
steps = len(df) # Num steps to read from the pressure tensor file
timestep = 1 # What timestep are you using in the pressure tensor file
temperature = 313 #System temp
volume = 176163.49
each = 100 # Sample frequency
plot = True # Create plot for autocorrelation function

# Conversion ratio from atm/bar to Pa
if unit == 'Pa':
    conv_ratio = 1
elif unit == 'atm':
    conv_ratio = 101325
elif unit == 'bar':
    conv_ratio = 100000

# Calculate the kBT value
kBT = Boltzmann * temperature

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initiate the pressure tensor component lists
Pxx, Pyy, Pzz, Pxy, Pxz, Pyz = [], [], [], [], [], []

# Read the pressure tensor elements from data file
print('\nReading the pressure tensor data file')
with open(datafile, "r") as file:
    next(file)
    for _ in trange(steps, ncols=100, desc='Progress'):
        line = file.readline()
        step = list(map(float, line.split()))
        Pxx.append(step[1]*conv_ratio)
        Pyy.append(step[2]*conv_ratio)
        Pzz.append(step[3]*conv_ratio)
        Pxy.append(step[4]*conv_ratio)
        Pxz.append(step[5]*conv_ratio)
        Pyz.append(step[6]*conv_ratio)

# Convert lists to numpy arrays
Pxx = np.array(Pxx)
Pyy = np.array(Pyy)
Pzz = np.array(Pzz)
Pxy = np.array(Pxy)
Pxz = np.array(Pxz)
Pyz = np.array(Pyz)

# Generate the time array
end_step = steps * timestep
Time = np.linspace(0, end_step, num=steps, endpoint=False)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Viscosity from Einstein relation
def einstein(timestep):

    Pxxyy = (Pxx - Pyy) / 2
    Pyyzz = (Pyy - Pzz) / 2

    '''
    Calculate the viscosity from the Einstein relation 
    by integrating the components of the pressure tensor
    '''
    timestep = timestep * 10**(-12)

    Pxy_int = integrate.cumtrapz(y=Pxy, dx=timestep, initial=0)
    Pxz_int = integrate.cumtrapz(y=Pxz, dx=timestep, initial=0)
    Pyz_int = integrate.cumtrapz(y=Pyz, dx=timestep, initial=0)

    Pxxyy_int = integrate.cumtrapz(y=Pxxyy, dx=timestep, initial=0)
    Pyyzz_int = integrate.cumtrapz(y=Pyyzz, dx=timestep, initial=0)

    integral = (Pxy_int**2 + Pxz_int**2 + Pyz_int**2 + Pxxyy_int**2 + Pyyzz_int**2) / 5

    viscosity = integral[1:] * ( volume * 10**(-30) / (2 * kBT * Time[1:] * 10**(-12)) )

    return viscosity

viscosity = einstein(timestep=timestep)

print(f"\nViscosity (Einstein): {round((viscosity[-1] * 1000), 2)} [mPa.s]")

# Plot the running integral of viscosity
if plot:
    pyplot.figure()
    pyplot.plot(Time[:viscosity.shape[0]], viscosity[:]*1000, label='Viscosity')
    pyplot.xlabel('Time (ps)')
    pyplot.ylabel('Viscosity (mPa.s)')
    pyplot.legend()
    pyplot.show()

# Save the running integral of viscosity as a csv file
df = pd.DataFrame({"time(ps)" : Time[:viscosity.shape[0]:each], "viscosity(Pa.s)" : viscosity[::each]})
df.to_csv("viscosity_Einstein.csv", index=False)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Viscosity from Green-Kubo relation
def green_kubo(timestep):
    # Calculate the ACFs
    Pxy_acf = acf(Pxy)
    Pxz_acf = acf(Pxz)
    Pyz_acf = acf(Pyz)

    # Calculate the shear components of the pressure tensor and their ACF
    if diag:
        Pxxyy = (Pxx - Pyy) / 2
        Pyyzz = (Pyy - Pzz) / 2
        Pxxzz = (Pxx - Pzz) / 2

        Pxxyy_acf = acf(Pxxyy)
        Pyyzz_acf = acf(Pyyzz)
        Pxxzz_acf = acf(Pxxzz)

    if diag:
        avg_acf = (Pxy_acf + Pxz_acf + Pyz_acf + Pxxyy_acf + Pyyzz_acf + Pxxzz_acf) / 6
    else:
        avg_acf = (Pxy_acf + Pxz_acf + Pyz_acf) / 3

    # Integrate the average ACF to get the viscosity
    timestep = timestep * 10**(-12)
    integral = integrate.cumtrapz(y=avg_acf, dx=timestep, initial=0)
    viscosity = integral * (volume * 10**(-30) / kBT)

    return avg_acf, viscosity

avg_acf, viscosity = green_kubo(timestep)

# Plot the normalized average ACF
if plot:
    norm_avg_acf = avg_acf / avg_acf[0]
    pyplot.figure()
    pyplot.plot(Time[:avg_acf.shape[0]], norm_avg_acf[:], label='Average')
    pyplot.xlabel('Time (ps)')
    pyplot.ylabel('ACF')
    pyplot.legend()
    pyplot.show()

# Save the normalized average ACF as a csv file
norm_avg_acf = avg_acf / avg_acf[0]
df = pd.DataFrame({"time (ps)" : Time[:avg_acf.shape[0]], "ACF" : norm_avg_acf[:]})
df.to_csv("avg_acf.csv", index=False)

print(f"Viscosity (Green-Kubo): {round((viscosity[-1] * 1000), 2)} [mPa.s]")

# Plot the time evolution of the viscosity estimate
if plot:
    pyplot.figure()
    pyplot.plot(Time[:viscosity.shape[0]], viscosity[:]*1000, label='Viscosity')
    pyplot.xlabel('Time (ps)')
    pyplot.ylabel('Viscosity (mPa.s)')
    pyplot.legend()
    pyplot.show()

# Save running integral of the viscosity as a csv file
df = pd.DataFrame({"time(ps)" : Time[:viscosity.shape[0]:each], "viscosity(Pa.s)" : viscosity[::each]})
df.to_csv("viscosity_GK.csv", index=False)
