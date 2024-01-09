import os
from os.path import join
from os import getcwd
from os import chdir
from rdkit import Chem
import subprocess
from copy import deepcopy
from os import listdir
import random as rnd

def runcmd(cmd, verbose = False, *args, **kwargs):
    #bascially allows python to run a bash command, and the code makes sure 
    #the error of the subproceess is communicated if it fails
    process = subprocess.run(
        cmd,
        text=True,
        shell=True)
    
    return process

def MakeLAMMPSFile(Name, CWD, GKRuntime, Temp, ID):
    
    if os.path.exists(f"{os.path.join(CWD, f'{Name}_system_{Temp}K_{ID}.lammps')}"):
        print('Specified LAMMPS file already exists in this location, overwriting.')
        os.remove(f"{os.path.join(CWD, f'{Name}_system_{Temp}K_{ID}.lammps')}")

    # Write LAMMPS file for 40C run
    with open(os.path.join(CWD, f'{Name}_system_{Temp}K_{ID}.lammps'), 'x') as file:
        file.write(f"""
# Setup parameters
variable       		T equal {Temp} # Equilibrium temperature [K]
log             	logEQM_{Name}_T${{T}}KP1atm.out

#include         "{Name}_system.in.init" Need to determine correct Kspace params

# Potential information
units           	real
dimension       	3
boundary        	p p p
atom_style      	full

pair_style      	lj/cut/coul/cut 12.0 12.0 
bond_style      	harmonic
angle_style     	harmonic
dihedral_style 		opls
improper_style     	harmonic
pair_modify 		mix geometric tail yes
special_bonds   	lj/coul 0.0 0.0 0.5

# Read lammps data file consist of molecular topology and forcefield info
read_data       	{Name}_system.data
neighbor        	2.0 bin
neigh_modify 		every 1 delay 0 check yes

include         "{Name}_system.in.charges"
include         "{Name}_system.in.settings"

# Define variables
variable        	eqmT equal $T			 			# Equilibrium temperature [K]
variable        	eqmP equal 1.0						# Equilibrium pressure [atm]
variable    		p equal 100						    # Nrepeat, correlation length
variable    		s equal 10       					# Nevery, sample interval
variable    		d equal $s*$p  						# Nfreq, dump interval
variable 			rho equal density

# Minimize system at target temperature using the default conjugate gradient method
velocity        	all create ${{eqmT}} 482648
fix             	min all nve
thermo          	10
thermo_style 		custom step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol
# dump            	1 all custom 10 min_w_{Name}_T${{T}}KP1atm.lammpstrj id mol type x y z mass q
# dump            	2 all custom 10 min_u_{Name}_T${{T}}KP1atm.lammpstrj id mol type xu yu zu mass q
# dump_modify     	1 sort id
# dump_modify     	2 sort id
minimize        	1.0e-16 1.06e-6 100000 500000
# undump          	1
# undump          	2
write_restart   	Min_{Name}_T${{T}}KP1atm.restart

unfix           	min
reset_timestep  	0
neigh_modify 		every 1 delay 0 check yes

# NVT at high temperature
fix             	nvt{Temp}K all nvt temp {Temp}.0 {Temp}.0 100.0
thermo			    $d
thermo_style 		custom step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol
fix 				thermo_print all print $d "$(step) $(temp) $(press) $(density) $(pxx) $(pyy) $(pzz) $(pxy) $(pxz) $(pyz) $(pe) $(ke) $(etotal) $(evdwl) $(ecoul) $(epair) $(ebond) $(eangle) $(edihed) $(eimp) $(emol) $(etail) $(enthalpy) $(vol)" &
					append thermoNVT{Temp}K_{Name}_T${{T}}KP1atm.out screen no title "# step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol"
# dump            	1 all custom $d NVT{Temp}K_u_{Name}_T${{T}}KP1atm.lammpstrj id mol type xu yu zu mass q
# dump_modify     	1 sort id
run            		250000
# undump          	1
unfix			    nvt{Temp}K
unfix               thermo_print
write_restart   	NVT_{Name}_T${{T}}KP1atm.restart

# NPT: Isothermal-isobaric ensemble to set the desired pressure; compute average density at that pressure
fix 				NPT all npt temp ${{eqmT}} ${{eqmT}} 100.0 iso ${{eqmP}} ${{eqmP}} 25.0
fix             	dave all ave/time $s $p $d v_rho ave running file eqmDensity_{Name}_T${{T}}KP1atm.out
thermo				$d
thermo_style 		custom step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol
fix 				thermo_print all print $d "$(step) $(temp) $(press) $(density) $(pxx) $(pyy) $(pzz) $(pxy) $(pxz) $(pyz) $(pe) $(ke) $(etotal) $(evdwl) $(ecoul) $(epair) $(ebond) $(eangle) $(edihed) $(eimp) $(emol) $(etail) $(enthalpy) $(vol)" &
					append thermoNPT_{Name}_T${{T}}KP1atm.out screen no title "# step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol"
# dump            	1 all custom $d NPT_u_{Name}_T${{T}}KP1atm.lammpstrj id mol type xu yu zu mass q
# dump_modify     	1 sort id
run					1000000
# undump          	1
unfix				NPT
unfix               thermo_print
write_restart  		NPT_{Name}_T${{T}}KP1atm.restart

# NVT: Canonical ensemble to deform the box to match increase in P in previous step
variable        	averho equal f_dave
variable        	adjustrho equal (${{rho}}/${{averho}})^(1.0/3.0) # Adjustment factor needed to bring rho to averge rho
unfix				dave
fix             	NVT all nvt temp ${{eqmT}} ${{eqmT}} 100.0	
fix             	adjust all deform 1 x scale ${{adjustrho}} y scale ${{adjustrho}} z scale ${{adjustrho}}
thermo         		$d
thermo_style 		custom step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol
fix 				thermo_print all print $d "$(step) $(temp) $(press) $(density) $(pxx) $(pyy) $(pzz) $(pxy) $(pxz) $(pyz) $(pe) $(ke) $(etotal) $(evdwl) $(ecoul) $(epair) $(ebond) $(eangle) $(edihed) $(eimp) $(emol) $(etail) $(enthalpy) $(vol)" &
					append thermoNVT_{Name}_T${{T}}KP1atm.out screen no title "# step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol"
# dump            	1 all custom $d NVT_u_{Name}_T${{T}}KP1atm.lammpstrj id mol type xu yu zu mass q
# dump_modify     	1 sort id
run					500000
# undump          	1
unfix				NVT
unfix           	adjust
unfix               thermo_print
write_restart  		NVT_{Name}_T${{T}}KP1atm.restart

# NVE: Microcanonical ensemble to explore the configuration space at constant T and V; relax
fix	       			NVE all nve
fix 				thermostat all langevin ${{eqmT}} ${{eqmT}} 100.0 39847 
thermo          	$d
thermo_style 		custom step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol

fix 			thermo_print all print $d "$(step) $(temp) $(press) $(density) $(pxx) $(pyy) $(pzz) $(pxy) $(pxz) $(pyz) $(pe) $(ke) $(etotal) $(evdwl) $(ecoul) $(epair) $(ebond) $(eangle) $(edihed) $(eimp) $(emol) $(etail) $(enthalpy) $(vol)" &
			append thermoNVE_{Name}_T${{T}}KP1atm.out screen no title "# step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol"

run             	250000

unfix           	NVE

unfix 			thermostat

unfix               	thermo_print

# Output the state generated that is needed to shear the molecules

write_restart  		state_{Name}_T${{T}}KP1atm.restart
write_data 		equilibrated.data

# Green-Kubo method via fix ave/correlate

log                 logGKvisc_{Name}_T${{T}}KP1atm.out

# Define variables
variable        	eqmT equal $T 				# Equilibrium temperature [K]
variable        	tpdn equal 3*1E6 			# Time for production run [fs]

variable    		dt equal 1.0				# time step [fs]
variable 		    V equal vol

# convert from LAMMPS real units to SI
variable    		kB equal 1.3806504e-23    	
variable    		atm2Pa equal 101325.0		
variable    		A2m equal 1.0e-10 			
variable    		fs2s equal 1.0e-15 			
variable			Pas2cP equal 1.0e+3			
variable    		convert equal ${{atm2Pa}}*${{atm2Pa}}*${{fs2s}}*${{A2m}}*${{A2m}}*${{A2m}}*${{Pas2cP}}

##################################### Viscosity Calculation #####################################################
timestep     		${{dt}}						# define time step [fs]

compute         	TT all temp
compute         	myP all pressure TT
fix             	1 all nve
fix             	2 all langevin ${{eqmT}} ${{eqmT}} 100.0 482648

variable        	myPxx equal c_myP[1]
variable        	myPyy equal c_myP[2]
variable       		myPzz equal c_myP[3]
variable     		myPxy equal c_myP[4]
variable     		myPxz equal c_myP[5]
variable     		myPyz equal c_myP[6]

fix             	3 all ave/time 1 1 1 v_myPxx v_myPyy v_myPzz v_myPxy v_myPxz v_myPyz ave one file Stress_AVGOne111_{Name}_T${{T}}KP1atm.out
fix             	4 all ave/time $s $p $d v_myPxx v_myPyy v_myPzz v_myPxy v_myPxz v_myPyz ave one file Stress_AVGOnespd_{Name}_T${{T}}KP1atm.out

variable    kB equal 1.3806504e-23    # [J/K] Boltzmann
variable    atm2Pa equal 101325.0
variable    A2m equal 1.0e-10
variable    fs2s equal 1.0e-15
variable    convert equal ${{atm2Pa}}*${{atm2Pa}}*${{fs2s}}*${{A2m}}*${{A2m}}*${{A2m}}

fix          SS all ave/correlate $s $p $d &
             v_myPxy v_myPxz v_myPyz type auto file S0St.dat ave running

variable     scale equal ${{convert}}/(${{kB}}*$T)*$V*$s*${{dt}}
variable     v11 equal trap(f_SS[3])*${{scale}}
variable     v22 equal trap(f_SS[4])*${{scale}}
variable     v33 equal trap(f_SS[5])*${{scale}}

thermo       		$d
thermo_style custom step temp press v_myPxy v_myPxz v_myPyz v_v11 v_v22 v_v33 vol
fix thermo_print all print $d "$(temp) $(press) $(v_myPxy) $(v_myPxz) $(v_myPyz) $(v_v11) $(v_v22) $(v_v33) $(vol)" &
    append thermoNVE_{Name}_T${{T}}KP1atm.out screen no title "# temp press v_myPxy v_myPxz v_myPyz v_v11 v_v22 v_v33 vol"

# Dump all molecule coordinates

dump         1 all custom $d All_u_{Name}_T${{T}}KP1atm.lammpstrj id mol type xu yu zu mass q
run          {GKRuntime}
variable     v equal (v_v11+v_v22+v_v33)/3.0
variable     ndens equal count(all)/vol
print        "Average viscosity: $v [Pa.s] @ $T K, ${{ndens}} atoms/A^3"

write_restart   	GKvisc_{Name}_T${{T}}KP1atm.restart
write_data          GKvisc_{Name}_T${{T}}KP1atm.data

""")
        
def MakeMoltemplateFile(Name, CWD, NumMols, BoxL):
    if os.path.exists(f"{os.path.join(CWD, f'{Name}_system.lt')}"):
        print('Specified Moltemplate file already exists in this location, overwriting.')
        os.remove(f"{os.path.join(CWD, f'{Name}_system.lt')}")

    with open(os.path.join(CWD, f'{Name}_system.lt'), 'x') as file:
                file.write(f"""
import "{Name}.lt"  # <- defines the molecule type.

# Periodic boundary conditions:
write_once("Data Boundary") {{
   0.0  {BoxL}.00  xlo xhi
   0.0  {BoxL}.00  ylo yhi
   0.0  {BoxL}.00  zlo zhi
}}

ethylenes = new {Name} [{NumMols}]
""")
                
def MakePackmolFile(Name, CWD, NumMols, Seed, BoxL):
    if os.path.exists(f"{os.path.join(CWD, f'{Name}.inp')}"):
        print('Packmol file already exists in this location, overwriting')
        os.remove(f"{os.path.join(CWD, f'{Name}.inp')}")

    with open(os.path.join(CWD, f'{Name}.inp'), 'x') as file:
        file.write(f"""
tolerance 2.0
output {Name}_PackmolFile.pdb

filetype pdb

seed {Seed}

structure {Name}.pdb
number {NumMols} 
inside cube 0. 0. 0. {BoxL}.
end structure""")

def CreateArrayJob(CWD, SimName, SimType, TopValue, BotValue, WALLTIME):
    #Create an array job for each separate simulation

    if os.path.exists(f"{join(CWD, f'{SimType}.pbs')}"):
        print(f'Specified file already exists in this location, overwriting')
        os.remove(f"{os.path.join(CWD, f'{SimType}.pbs')}")       

    with open(join(CWD, f'{SimType}.pbs'), 'w') as file:
        file.write(f"""#!/bin/bash
#PBS -l select=1:ncpus=32:mem=62gb
#PBS -l walltime={WALLTIME}
#PBS -J {BotValue}-{TopValue}

module load intel-suite/2020.2
module load mpi/intel-2019.6.166

cd /rds/general/user/eeo21/home/HIGHTHROUGHPUTSTUDIES/MDSIMULATIONEVLUATION//Run_${{PBS_ARRAY_INDEX}}
mpiexec ~/tmp/bin/lmp -in {SimName}_${{PBS_ARRAY_INDEX}}.lammps
""")
    # os.rename(join(CWD, f'{SimType}.pbs'),  f'{SimType}.pbs')}")

STARTINGDIR = deepcopy(getcwd())
PYTHONPATH = 'python3'
PYTHONPATH = 'C:/Users/eeo21/AppData/Local/Programs/Python/Python310/python.exe'
Molecules = [x for x in listdir(STARTINGDIR) if '.pdb' in x]
NumMols = 100
NumRuns = 100
RunList = list(range(1, NumRuns+1))
# Values for the array job
TopValue = RunList[-1] 
BotValue = RunList[0]
BoxL = 80
LOPLS = False
WALLTIME = '07:59:59'

# runcmd('mkdir Trajectory_Studies')

# for Molecule in Molecules:
#     FolderName = Molecule.split('.')[0]
#     Path = join(STARTINGDIR, Molecule)
#     MolObject = Chem.MolFromPDBFile(Path)

#     if Molecule == 'squalane.pdb':
#         SMILESString = 'CC(=CCCC(=CCCC(=CCCC=C(C)CCC=C(C)CCC=C(C)C)C)C)C'
#     else:
#         SMILESString = Chem.MolToSmiles(MolObject)
    
#     if LOPLS:
#         LTCOMMAND = f"{join(STARTINGDIR, 'rdlt.py')} --smi '{SMILESString}' -n {FolderName} -l -c"
#     else:
#         LTCOMMAND = f"{join(STARTINGDIR, 'rdlt.py')} --smi '{SMILESString}' -n {FolderName} -c"
    
#     #Create initial Moltemplate file
#     runcmd(f'{PYTHONPATH} {LTCOMMAND} > {STARTINGDIR}/{FolderName}.lt')
    
#     #Enter Trajectory studies directory
#     chdir(join(getcwd(), 'Trajectory_Studies'))

#     #Make Molecule Folder in Trajectory studies directory
#     runcmd(f'mkdir {FolderName}')

#     #Copy molecule pdb to molecule directory
#     runcmd(f'cp "{join(STARTINGDIR, Molecule)}" {join(getcwd(), FolderName)}')

#     #Enter molecule directory
#     chdir(join(getcwd(), FolderName))

#     ### Make files for short vs long run convergence
#     """
#     - So we could run 100 trajectories for each molecule up to 10ns with 100 molecules
#     - We could generate trajectories for 500ps, 1ns, 2ns, and 4ns by
#     truncating each of the runs to the respective times
#     - We could then take between 5 and 10 random samples from the 100 trajectories
#     at each time length for varying numbers of trajectories
#     - Take an average of the samples as the predicted value, and use the extreme values 
#     or calculate the standard deviation as the uncertainty (2std for 95% confidence)
#     """

#     for x in RunList:
#         runcmd(f'mkdir Run_{x}')
#         Trajectory = f'Run_{x}'
#         # Copy pdb
#         runcmd(f'cp"{join(STARTINGDIR, Molecule)}" {join(getcwd(), Trajectory, Molecule)}')
#         # Copy lt file
#         ltfile = f'{FolderName}.lt'
#         runcmd(f'cp "{join(STARTINGDIR, ltfile)}" {join(getcwd(), Trajectory, ltfile)}')
#         # Set random seed
#         Seed = rnd.randint(0, 1E6)

#         # Make Packmol Input file
#         MakePackmolFile(Name=FolderName, CWD=join(getcwd(), Trajectory), NumMols=NumMols, Seed=Seed, BoxL=BoxL)
#         # Make Moltemplate file
#         MakeMoltemplateFile(Name=FolderName, CWD=join(getcwd(), Trajectory), NumMols=NumMols, BoxL=BoxL)
#         # Make Lammps Files
#         MakeLAMMPSFile(Name=FolderName, ID=x, CWD=join(getcwd(), Trajectory), GKRuntime=10000000, Temp=313)
#         MakeLAMMPSFile(Name=FolderName, ID=x, CWD=join(getcwd(), Trajectory), GKRuntime=10000000, Temp=363)
#         # Make packmol coordinate file and LAMMPS data file
#         chdir(join(getcwd(), Trajectory))
#         if PYTHONPATH == 'python3':
#             runcmd(f'packmol < {FolderName}.inp')
#             runcmd(f'moltemplate.sh -pdb {FolderName}_PackmolFile.pdb {FolderName}_system.lt')

#         # Return to starting directory
#         chdir(join(STARTINGDIR, 'Trajectory_Studies', FolderName))
#     chdir(STARTINGDIR)

### Files showing long term variations in viscosity and thermal conductivity prediction
"""
- Let's run 10 simulations for each molecule up to 100ns, will be pretty robust for 
uncertainty calculations.
- Run at 40C and 100C 
"""
STARTINGDIR = deepcopy(getcwd())
PYTHONPATH = 'python3'
PYTHONPATH = 'C:/Users/eeo21/AppData/Local/Programs/Python/Python310/python.exe'
Molecules = [x for x in listdir(STARTINGDIR) if '.pdb' in x]
NumRuns = 10
RunList = list(range(1, NumRuns+1))
# Values for the array job
TopValue = RunList[-1]
BotValue = RunList[0]
BoxL = 80
LOPLS = False
WALLTIME = '07:59:59'

runcmd('mkdir LongTermVariations')

for Molecule in Molecules:
    FolderName = Molecule.split('.')[0]
    Path = join(STARTINGDIR, Molecule)
    MolObject = Chem.MolFromPDBFile(Path)

    if Molecule == 'squalane.pdb':
        SMILESString = 'CC(=CCCC(=CCCC(=CCCC=C(C)CCC=C(C)CCC=C(C)C)C)C)C'
    else:
        SMILESString = Chem.MolToSmiles(MolObject)
    
    if LOPLS:
        LTCOMMAND = f"{join(STARTINGDIR, 'rdlt.py')} --smi {SMILESString} -n {FolderName} -l -c"
    else:
        LTCOMMAND = f"{join(STARTINGDIR, 'rdlt.py')} --smi {SMILESString} -n {FolderName} -c"
    
    #Create initial Moltemplate file
    runcmd(f'{PYTHONPATH} {LTCOMMAND} > {STARTINGDIR}/{FolderName}.lt')
    
    #Enter Trajectory studies directory
    chdir(join(getcwd(), 'LongTermVariations'))

    #Make Molecule Folder in Trajectory studies directory
    runcmd(f'mkdir {FolderName}')

    #Copy molecule pdb to molecule directory
    runcmd(f'copy "{join(STARTINGDIR, Molecule)}" {join(getcwd(), FolderName)}')

    #Enter molecule directory
    chdir(join(getcwd(), FolderName))

    NumMols = [25, 50, 100, 250, 500]
    for NumMol in NumMols:
        runcmd(f'mkdir NumMols_{NumMol}')
        chdir(join(getcwd(), f'NumMols_{NumMol}'))

        for x in RunList:
            runcmd(f'mkdir Run_{x}')
            Run = f'Run_{x}'
            # Copy pdb
            runcmd(f'copy "{join(STARTINGDIR, Molecule)}" {join(getcwd(), Run, Molecule)}')
            # Copy lt file
            ltfile = f'{FolderName}.lt'
            runcmd(f'copy "{join(STARTINGDIR, ltfile)}" {join(getcwd(), Run, ltfile)}')
            # Set random seed
            Seed = rnd.randint(0, 1E6)

            # Make Packmol Input file
            MakePackmolFile(Name=FolderName, CWD=join(getcwd(), Run), NumMols=NumMol, Seed=Seed, BoxL=BoxL)
            # Make Moltemplate file
            MakeMoltemplateFile(Name=FolderName, CWD=join(getcwd(), Run), NumMols=NumMol, BoxL=BoxL)
            # Make Lammps Files
            MakeLAMMPSFile(Name=FolderName, ID=x, CWD=join(getcwd(), Run), GKRuntime=10000000, Temp=313)
            MakeLAMMPSFile(Name=FolderName, ID=x, CWD=join(getcwd(), Run), GKRuntime=10000000, Temp=363)
            # Make packmol coordinate file and LAMMPS data file
            chdir(join(getcwd(), Run))
            if PYTHONPATH == 'python3':
                runcmd(f'packmol < {FolderName}.inp')
                runcmd(f'moltemplate.sh -pdb {FolderName}_PackmolFile.pdb {FolderName}_system.lt')

            # Return to starting directory
            chdir(join(STARTINGDIR, 'LongTermVariations', FolderName, f'NumMols_{NumMol}'))
    
        chdir(join(STARTINGDIR, 'LongTermVariations', FolderName))
    
    chdir(STARTINGDIR)

# ### Files showing finite size effects
# """
# Run experiments at 25, 50, 100, 250, 500 molecules
# Run 10 simulations for each number of molecules for 10ns
# """