import sys, argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import integrate
from scipy.constants import Boltzmann
from scipy.optimize import curve_fit
from tqdm import trange
from os.path import join as join
from os import chdir
from os import getcwd
from os import listdir
import traceback
import os

# Define ACF using FFT
def acf(data):
    steps = data.shape[0]
    lag = steps // 2

    size = 2 ** np.ceil(np.log2(2 * steps - 1)).astype('int')
    FFT = np.fft.fft(data, size)
    PWR = FFT.conjugate() * FFT
    COR = np.fft.ifft(PWR)[:steps].real
    autocorrelation = COR / np.arange(steps, 0, -1)

    return autocorrelation[:lag]

# Viscosity from Einstein relation
def einstein(timestep):
    Pxxyy = (Pxx - Pyy) / 2
    Pyyzz = (Pyy - Pzz) / 2

    timestep = timestep * 10**(-12)
    Pxy_int = integrate.cumtrapz(y=Pxy, dx=timestep, initial=0)
    Pxz_int = integrate.cumtrapz(y=Pxz, dx=timestep, initial=0)
    Pyz_int = integrate.cumtrapz(y=Pyz, dx=timestep, initial=0)

    integral = (Pxy_int**2 + Pxz_int**2 + Pyz_int**2) / 3
    viscosity = integral[1:] * (volume * 10**(-30) / (2 * kBT * Time[1:] * 10**(-12)))

    return viscosity

# Define a function to pad the viscosity array
def pad_viscosity(viscosity, min_length, target_length):
    if len(viscosity) >= min_length and len(viscosity) <= target_length:
        # Calculate how many elements to add
        elements_to_add = target_length - len(viscosity)
        # Pad the viscosity array with NaNs
        viscosity = np.pad(viscosity, (0, elements_to_add), mode='constant', constant_values=np.nan)
    return viscosity

chdir('/rds/general/ephemeral/user/eeo21/ephemeral/Benchmarking/LOPLS')
STARTDIR = getcwd()
Names = [x for x in listdir(getcwd()) if os.path.isdir(x)]
Temps = [313, 373]

# Initialize a list to store the final viscosity values for each material and temperature
final_viscosities = []

for Name in Names:
    print(Name)
    for Temp in Temps:
        print(Temp)
        chdir(join(STARTDIR, Name))
        Runs = [x for x in listdir(getcwd()) if os.path.isdir(x)]

        DataframeEinstein = pd.DataFrame()
        DataframeGK = pd.DataFrame()

        for Run in Runs:
            try:
                chdir(join(STARTDIR, Name, Run))

                df = pd.read_csv(f'Stress_AVGOnespd_{Name}_T{Temp}KP1atm.out')

                unit = 'atm'
                datafile = f'Stress_AVGOnespd_{Name}_T{Temp}KP1atm.out'
                diag = False
                steps = len(df) - 1
                timestep = 1
                temperature = Temp

                with open(f'logGKvisc_{Name}_T{Temp}KP1atm.out', "r") as file:
                    content = file.readlines()
                    for line in content:
                        linecontent = line.split(' ')
                        linecontent = [x for x in linecontent if x != '']
                        if len(linecontent) == 18:
                            try:
                                vol = linecontent[9]
                                volume = float(vol)
                            except:
                                pass

                each = 2
                plot = True

                if unit == 'Pa':
                    conv_ratio = 1
                elif unit == 'atm':
                    conv_ratio = 101325
                elif unit == 'bar':
                    conv_ratio = 100000

                kBT = Boltzmann * temperature

                Pxx, Pyy, Pzz, Pxy, Pxz, Pyz = [], [], [], [], [], []

                with open(datafile, "r") as file:
                    next(file)
                    next(file)

                    for _ in range(steps):
                        line = file.readline()
                        step = list(map(float, line.split()))
                        Pxx.append(step[1]*conv_ratio)
                        Pyy.append(step[2]*conv_ratio)
                        Pzz.append(step[3]*conv_ratio)
                        Pxy.append(step[4]*conv_ratio)
                        Pxz.append(step[5]*conv_ratio)
                        Pyz.append(step[6]*conv_ratio)

                Pxx = np.array(Pxx)
                Pyy = np.array(Pyy)
                Pzz = np.array(Pzz)
                Pxy = np.array(Pxy)
                Pxz = np.array(Pxz)
                Pyz = np.array(Pyz)

                end_step = steps * timestep
                Time = np.linspace(0, end_step, num=steps, endpoint=False)

                viscosity = einstein(timestep=timestep)
                viscE = round((viscosity[-1] * 1000), 2)

                if plot:
                    plt.figure()
                    plt.plot(Time[:viscosity.shape[0]], viscosity[:]*1000, label='Viscosity')
                    plt.xlabel('Time (ps)')
                    plt.ylabel('Einstein Viscosity (mPa.s)')
                    plt.legend([f'Viscosity Estimate: {viscE} [mPa.s]'])
                    plt.title(f'{Name}_{Run}_{Temp}_Einstein')
                    plt.xlim(0, 1)
                    plt.ylim(0, 35)
                    plt.savefig(join(STARTDIR, Name, f'{Name}_{Run}_{Temp}_Einstein.png'))
                    plt.close()

                # Add to DataframeEinstein only if the viscosity array length is at least 1000

                if len(viscosity) > 1000:
                    viscosity_padded = pad_viscosity(viscosity, min_length=1000, target_length=1500)
                    DataframeEinstein[f'Viscosity_{Run}'] = viscosity_padded * 1000

                def green_kubo(timestep):
                    Pxy_acf = acf(Pxy)
                    Pxz_acf = acf(Pxz)
                    Pyz_acf = acf(Pyz)
                    avg_acf = (Pxy_acf + Pxz_acf + Pyz_acf) / 3

                    timestep = timestep * 10**(-12)
                    integral = integrate.cumtrapz(y=avg_acf, dx=timestep, initial=0)
                    viscosity = integral * (volume * 10**(-30) / kBT)

                    return avg_acf, viscosity

                avg_acf, viscosity = green_kubo(timestep)

                if plot:
                    norm_avg_acf = avg_acf / avg_acf[0]
                    plt.figure()
                    plt.plot(Time[:avg_acf.shape[0]], norm_avg_acf[:], label='Average')
                    plt.xlabel('Time (ps)')
                    plt.ylabel('Green Kubo ACF')
                    plt.legend()
                    plt.title(f'{Name} {Temp}K Green Kubo ACF')
                    plt.grid(color='grey', linestyle='--', linewidth=0.5)
                    plt.grid(which="minor", linestyle='--', linewidth=0.2)
                    plt.minorticks_on()
                    plt.savefig(join(STARTDIR, Name, f'{Name}_{Run}_{Temp}_2nsGKACF.png'))
                    plt.close()

                # Add to DataframeGK only if the viscosity array length is at least 1000
                if len(viscosity) >= 500:
                    viscosity_padded = pad_viscosity(viscosity, min_length=500, target_length=750)
                    DataframeGK[f'Viscosity_{Run}'] = viscosity_padded * 1000

            except Exception as E:
                print(E)
                traceback.print_exc()
                pass

        try:
            DataframeEinstein = DataframeEinstein.dropna()
            DataframeGK = DataframeGK.dropna()

            DataframeGK['Average'] = DataframeGK.mean(axis=1)
            DataframeGK['STD'] = DataframeGK.std(axis=1)
            DataframeEinstein['Average'] = DataframeEinstein.mean(axis=1)
            DataframeEinstein['STD'] = DataframeEinstein.std(axis=1)

            # Save the viscosity data to CSV files
            DataframeEinstein.to_csv(join(STARTDIR, Name, f'{Name}_{Temp}K_Einstein_Viscosity.csv'), index=False)
            DataframeGK.to_csv(join(STARTDIR, Name, f'{Name}_{Temp}K_GreenKubo_Viscosity.csv'), index=False)

            # Extract the final viscosity values
            ViscosityAvList = DataframeGK['Average'].tolist()
            ViscosityAv = ViscosityAvList[-1]
            ViscosityAvListEinstein = DataframeEinstein['Average'].tolist()
            ViscosityAvEinstein = ViscosityAvListEinstein[-1]

            # print(ViscosityAv)
            # ViscosityAvEinstein = round((DataframeEinstein['Average'].iloc[-1]), 2)

            # Add the final viscosities to the master list
            final_viscosities.append({
                "Material": Name,
                "Temperature (K)": Temp,
                "Einstein Viscosity (mPa.s)": ViscosityAvEinstein,
                "Green-Kubo Viscosity (mPa.s)": ViscosityAv
            })

            # Continue with plotting and additional analysis
            DataframeGKViscList_Average = DataframeGK['Average'].to_list()
            DataframeGKViscList_AverageSTD = DataframeGK['STD'].to_list()
            ViscGK_UpperSTD = [a + b for a, b in zip(DataframeGKViscList_Average, DataframeGKViscList_AverageSTD)]
            ViscGK_LowerSTD = [a - b for a, b in zip(DataframeGKViscList_Average, DataframeGKViscList_AverageSTD)]
            DataframeEinsteinList_Average = DataframeEinstein['Average'].to_list()
            DataframeEinsteinList_AverageSTD = DataframeEinstein['STD'].to_list()
            Einstein_LowerSTD = [a + b for a, b in zip(DataframeEinsteinList_Average, DataframeEinsteinList_AverageSTD)]
            Einstein_UpperSTD = [a - b for a, b in zip(DataframeEinsteinList_Average, DataframeEinsteinList_AverageSTD)]

            step = list(range(0, len(DataframeGKViscList_Average)))
            step = [x/1000 for x in step]

            ViscPlt, Vplot = plt.subplots()
            Vplot.set_title(f'GK Viscosity Average - {Name} {Temp}K')
            Vplot.set_ylabel('Viscosity [mPa.S]')
            Vplot.set_xlabel('Time (ns)')
            Vplot.plot(step, DataframeGKViscList_Average)
            Vplot.legend([f'Viscosity Estimate: {round((DataframeGKViscList_Average[-1]), 2)} [mPa.s]']) 
            Vplot.fill_between(step, ViscGK_LowerSTD, ViscGK_UpperSTD, alpha=0.4)
            Vplot.grid(color='grey', linestyle='--', linewidth=0.5)
            Vplot.grid(which="minor", linestyle='--', linewidth=0.2)
            plt.minorticks_on()
            Vplot.set_xlim(0, 1)
            Vplot.set_ylim(0, 35)
            plt.savefig(join(STARTDIR, Name,  f'AvViscPlotGK_{Name}{Temp}K_2ns.png'))
            plt.close()

            step = list(range(0, len(DataframeEinsteinList_Average)))
            step = [x/1000 for x in step]

            ViscPlt, Vplot = plt.subplots()
            Vplot.set_title(f'Einstein Viscosity Average - {Name} {Temp}K')
            Vplot.set_ylabel('Viscosity [mPa.S]')
            Vplot.set_xlabel('Time (ns)')
            Vplot.plot(step, DataframeEinsteinList_Average)
            Vplot.legend([f'Viscosity Estimate: {round((DataframeEinsteinList_Average[-1]), 2)} [mPa.s]']) 
            Vplot.fill_between(step, Einstein_LowerSTD, Einstein_UpperSTD, alpha=0.4)
            Vplot.grid(color='grey', linestyle='--', linewidth=0.5)
            Vplot.grid(which="minor", linestyle='--', linewidth=0.2)
            plt.minorticks_on()
            Vplot.set_xlim(0, 1)
            Vplot.set_ylim(0, 35)
            plt.savefig(join(STARTDIR, Name,  f'EinsteinAvViscPlot_{Name}{Temp}K_2ns.png'))
            plt.close()

            print(f'Green Kubo Visc: {ViscosityAv}')
            print(f'Einstein Visc: {ViscosityAvEinstein}')

        except Exception as E:
            print(E)
            traceback.print_exc()
            pass

# After processing all materials and temperatures, save the final viscosities to a master CSV
master_df = pd.DataFrame(final_viscosities)
master_df.to_csv(join(STARTDIR, 'Final_Viscosities.csv'), index=False)
