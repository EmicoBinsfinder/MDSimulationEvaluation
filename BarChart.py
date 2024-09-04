import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Saving the updated provided string data as a CSV file
csv_data = """Material,ForceField,Temperature (K),Einstein Viscosity (mPa.s),Green-Kubo Viscosity (mPa.s),Experimental
1-cyclohexyl-2-cyclohexylmethyl-heptadecane,LOPLS,313,20.87,28.77,24.30
1-cyclohexyl-2-cyclohexylmethyl-heptadecane,LOPLS,373,4.05,4.95,4.23
1-methyl_Napthalene,LOPLS,313,2.53,2.73,2.43
1-methyl_Napthalene,LOPLS,373,0.74,0.89,0.8
13-Phenylpentacosane,LOPLS,313,10.83,14.51,17.63
13-Phenylpentacosane,LOPLS,373,2.68,3.12,3.63
13-Cyclohexylpentacosane,LOPLS,313,10.46,13.91,18.99
13-Cyclohexylpentacosane,LOPLS,373,2.58,3.23,3.74
DCMP,LOPLS,313,29.83,34.79,30.8
DCMP,LOPLS,373,3.86,4.62,5
1-cyclohexyl-2-cyclohexylmethyl-heptadecane,COMPASS,313,17.29,19.1,24.30
1-cyclohexyl-2-cyclohexylmethyl-heptadecane,COMPASS,373,3.45,3.64,4.23
1-methyl_Napthalene,COMPASS,313,1.34,1.26,2.43
1-methyl_Napthalene,COMPASS,373,0.8,0.82,0.8
13-Phenylpentacosane,COMPASS,313,10.12,14.16,17.63
13-Phenylpentacosane,COMPASS,373,1.92,1.3,3.63
13-Cyclohexylpentacosane,COMPASS,313,12.12,10.96,18.99
13-Cyclohexylpentacosane,COMPASS,373,2.49,2.79,3.74
DCMP,COMPASS,313,4.54,3.46,30.8
DCMP,COMPASS,373,1.67,1.87,5
"""

# Writing the string to a CSV file
file_path = 'viscosity_data_updated.csv'
with open(file_path, 'w') as file:
    file.write(csv_data)

# Loading the CSV into a DataFrame
df = pd.read_csv(file_path)

# Splitting the data by temperature
df_313K = df[df['Temperature (K)'] == 313]
df_373K = df[df['Temperature (K)'] == 373]

# Function to plot data for a specific temperature
def plot_viscosities(df, temperature):
    fig, ax = plt.subplots(figsize=(14, 8))

    # Defining bar width
    bar_width = 0.6
    positions = np.arange(len(df['Material'].unique())) * 5  # positions for each material group

    # Plotting for each material
    unique_materials = df['Material'].unique()
    
    for i, material in enumerate(unique_materials):
        df_material = df[df['Material'] == material]

        # Plot Einstein viscosities for LOPLS and COMPASS
        ax.bar(positions[i] - bar_width*2, df_material[df_material['ForceField'] == 'LOPLS']['Einstein Viscosity (mPa.s)'].values, bar_width, label='Einstein (LOPLS)', color='blue')
        ax.bar(positions[i] - bar_width, df_material[df_material['ForceField'] == 'COMPASS']['Einstein Viscosity (mPa.s)'].values, bar_width, label='Einstein (COMPASS)', color='purple')

        # Plot Green-Kubo viscosities for LOPLS and COMPASS
        ax.bar(positions[i], df_material[df_material['ForceField'] == 'LOPLS']['Green-Kubo Viscosity (mPa.s)'].values, bar_width, label='Green-Kubo (LOPLS)', color='green')
        ax.bar(positions[i] + bar_width, df_material[df_material['ForceField'] == 'COMPASS']['Green-Kubo Viscosity (mPa.s)'].values, bar_width, label='Green-Kubo (COMPASS)', color='orange')

        # Plot Experimental viscosities
        ax.bar(positions[i] + bar_width*2, df_material['Experimental'].values, bar_width, label='Experimental', color='Red')

    # Adding labels and title
    ax.set_xlabel('Base Oil')
    ax.set_ylabel('Viscosity (mPa·s)')
    ax.set_title(f'Comparison of Viscosities for Different Cyclic Base Oils and Force Fields at {temperature}K')
    ax.set_xticks(positions)
    plt.grid(color='grey', linestyle='--', linewidth=0.5)
    plt.grid(which="minor", linestyle='--', linewidth=0.2)
    plt.minorticks_on()
    ax.set_xticklabels(unique_materials, rotation=30)

    # Simplify the legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())

    plt.tight_layout()
    plt.savefig(f'Forcefield Comparison {temperature}')
    plt.show()

# Plotting for 313K
plot_viscosities(df_313K, 313)

# Plotting for 373K
plot_viscosities(df_373K, 373)
