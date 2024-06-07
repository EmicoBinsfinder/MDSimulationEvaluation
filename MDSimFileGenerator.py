import os
from os.path import join
from os import getcwd
from os import chdir
from rdkit import Chem
import subprocess
from copy import deepcopy
from os import listdir
import random as rnd
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import re
from rdkit.Chem import AllChem

def runcmd(cmd, verbose = False, *args, **kwargs):
    #bascially allows python to run a bash command, and the code makes sure 
    #the error of the subproceess is communicated if it fails
    process = subprocess.run(
        cmd,
        text=True,
        shell=True)
    
    return process

def MakeLAMMPSFile(Name, CWD, GKRuntime, Temp, ID):
    
    VelNumber = rnd.randint(0, 1000000)

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

units real
atom_style full
bond_style harmonic
angle_style harmonic
dihedral_style hybrid opls multi/harmonic
improper_style harmonic
pair_style lj/cut/coul/long 12.0 #change to cut/coul/long for kspace calculations
pair_modify mix geometric tail yes
special_bonds lj/coul 0.0 0.0 0.5
kspace_style pppm 0.0001

# Read lammps data file consist of molecular topology and forcefield info
read_data       	{Name}_system.data
neighbor        	2.0 bin
neigh_modify 		every 1 delay 0 check yes

include         "{Name}_system.in.charges"
include         "settings.txt"

########################################################################################
#CHARGES

set type 80 charge -0.22    #CARBON CH3 (Sui)
set type 81 charge -0.148   #CARBON CH2 (Sui)
set type 88 charge -0.160   #CARBON CH (Sui)
set type 85 charge  0.074   #HYDROGEN HC_CH2&CH3 single bond (Sui)
set type 89 charge  0.160   #HYDROGEN HC_CH (hydrogen with double bonded carbon) (Sui)
set type 406 charge 0.75    #Ester -COOR >>>>> CARBON >>> C_2 (Pluhackova)
set type 407 charge -0.55   #Ester C=O >>>>> OXYGEN >>> O_2 (Pluhackova)
set type 408 charge -0.45   #Ester CO-O-R >>>>> OXYGEN >>> OS (Pluhackova)

########################################################################################
#NON BONDED INTERACTIONS

pair_coeff   80  80    0.066000000000000   3.500000     
pair_coeff   81  81    0.066000000000000   3.500000
pair_coeff   82  82    0.075710400000000   3.550000    
pair_coeff   85  85    0.026290631000000   2.500000    
pair_coeff   89  89    0.030000000000000   2.420000   
pair_coeff   406 406   0.105000000000000   3.187500
pair_coeff   407 407   0.167360000000000   3.108000
pair_coeff   408 408   0.169360000000000   2.550000

#########################################################################################
#DIHEDRALS

dihedral_coeff 40  multi/harmonic   0.728438     -0.864800   -1.23214     1.336379    0.000000    #CT CT C2 O2
dihedral_coeff 51  multi/harmonic   -0.448883    -1.178005   0.624517     0.982929    0.000000    #CT CT C2 OS
dihedral_coeff 52  multi/harmonic   0.217797     0.636519    -0.023530    -0.856219   0.000000    #HC CT CT OS
dihedral_coeff 68  multi/harmonic   5.435142     0.00000     -5.715862    0.00000     0.00000     #CT OS C2 O2 
dihedral_coeff 70  multi/harmonic   7.761123     -1.553674   -5.715855    0.00000     0.00000     #CT C2 OS CT
dihedral_coeff 183 multi/harmonic   -1.432205    2.078667	 0.624936     -1.239164   0.000000    #C2 CT CT CT 
dihedral_coeff 187 multi/harmonic   -0.174038    -0.440071   -0.017288    0.588167    0.000000    #C2 CT CT HC
dihedral_coeff 196 multi/harmonic   0.124000     -0.055020   0.214340     -0.356440   0.000000    #CT CT CT CT (Sui)
dihedral_coeff 265 multi/harmonic   -1.68071     1.796674    1.243688     -1.243040   0.000000    #C2 OS CT CT 
dihedral_coeff 799 multi/harmonic   0.636100     0.379845    1.020183     -2.006533   0.000000    #CT CT CT OS

##########################################################################################
#Dictionary 
#HC = Aliphatic hydrogen 
#OH = hydroxyl oxygen 
#OS = alkoxy oxygen
#O_2 = ester carbonyl oxygen 
#C_2 = ester carbonyl carbon
#HC = alpha methoxy carbon
#CT = alkoxy carbon 


# Define variables
variable        	eqmT equal $T			 			# Equilibrium temperature [K]
variable        	eqmP equal 1.0						# Equilibrium pressure [atm]
variable    		p equal 100						    # Nrepeat, correlation length (Check different sample lengths)
variable    		s equal 10       					# Nevery, sample interval
variable    		d equal $s*$p  						# Nfreq, dump interval
variable 			rho equal density

# Minimize system at target temperature using the default conjugate gradient method
velocity        	all create ${{eqmT}} {VelNumber}
fix             	min all nve
thermo          	10
thermo_style 		custom step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol
minimize        	1.0e-16 1.06e-6 100000 500000
write_restart   	Min_{Name}_T${{T}}KP1atm.restart

unfix           	min
reset_timestep  	0
neigh_modify 		every 1 delay 0 check yes

# NPT: Isothermal-isobaric ensemble to set the desired pressure; compute average density at that pressure
fix 				NPT all npt temp ${{eqmT}} ${{eqmT}} 100.0 iso ${{eqmP}} ${{eqmP}} 25.0
fix             	dave all ave/time $s $p $d v_rho ave running file eqmDensity_{Name}_T${{T}}KP1atm.out
thermo				$d
thermo_style 		custom step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol
fix 				thermo_print all print $d "$(step) $(temp) $(press) $(density) $(pxx) $(pyy) $(pzz) $(pxy) $(pxz) $(pyz) $(pe) $(ke) $(etotal) $(evdwl) $(ecoul) $(epair) $(ebond) $(eangle) $(edihed) $(eimp) $(emol) $(etail) $(enthalpy) $(vol)" &
					append thermoNPT_{Name}_T${{T}}KP1atm.out screen no title "# step temp press density pxx pyy pzz pxy pxz pyz pe ke etotal evdwl ecoul epair ebond eangle edihed eimp emol etail enthalpy vol"
run					1000000
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
run					500000
unfix				NVT
unfix           	adjust
unfix               thermo_print
write_restart  		NVT_{Name}_T${{T}}KP1atm.restart

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
variable            kCal2J equal 4186.0/6.02214e23
variable    		atm2Pa equal 101325.0		
variable    		A2m equal 1.0e-10 			
variable    		fs2s equal 1.0e-15 			
variable			PaS2mPaS equal 1.0e+3			
variable    		convert equal ${{atm2Pa}}*${{atm2Pa}}*${{fs2s}}*${{A2m}}*${{A2m}}*${{A2m}}*${{PaS2mPaS}}
variable            convertWk equal ${{kCal2J}}*${{kCal2J}}/${{fs2s}}/${{A2m}}

##################################### Viscosity Calculation #####################################################
timestep     		${{dt}}						# define time step [fs]

compute         	TT all temp
compute         	myP all pressure TT

###### Thermal Conductivity Calculations 

compute         myKE all ke/atom
compute         myPE all pe/atom
compute         myStress all stress/atom NULL virial

# compute heat flux vectors
compute         flux all heat/flux myKE myPE myStress
variable        Jx equal c_flux[1]/vol
variable        Jy equal c_flux[2]/vol
variable        Jz equal c_flux[3]/vol

fix             	1 all nve

variable        	myPxx equal c_myP[1]
variable        	myPyy equal c_myP[2]
variable       		myPzz equal c_myP[3]
variable     		myPxy equal c_myP[4]
variable     		myPxz equal c_myP[5]
variable     		myPyz equal c_myP[6]

fix             	3 all ave/time 1 1 1 v_myPxx v_myPyy v_myPzz v_myPxy v_myPxz v_myPyz ave one file Stress_AVGOne111_{Name}_T${{T}}KP1atm.out
fix             	4 all ave/time $s $p $d v_myPxx v_myPyy v_myPzz v_myPxy v_myPxz v_myPyz ave one file Stress_AVGOnespd_{Name}_T${{T}}KP1atm.out
fix             	FluxVector all ave/time $s $p $d v_Jx v_Jy v_Jz ave one file HeatFlux_AVGOnespd_{Name}_T${{T}}KP1atm.out

fix          SS all ave/correlate $s $p $d &
             v_myPxy v_myPxz v_myPyz type auto file S0St.dat ave running

variable     scale equal ${{convert}}/(${{kB}}*$T)*$V*$s*${{dt}}
variable     v11 equal trap(f_SS[3])*${{scale}}
variable     v22 equal trap(f_SS[4])*${{scale}}
variable     v33 equal trap(f_SS[5])*${{scale}}

fix          JJ all ave/correlate $s $p $d &
             c_flux[1] c_flux[2] c_flux[3] type auto &
             file profile.heatflux ave running

variable        scaleWk equal ${{convertWk}}/${{kB}}/$T/$T/$V*$s*${{dt}}
variable        k11 equal trap(f_JJ[3])*${{scaleWk}}
variable        k22 equal trap(f_JJ[4])*${{scaleWk}}
variable        k33 equal trap(f_JJ[5])*${{scaleWk}}

##### Diffusion Coefficient Calculations 

compute         vacf all vacf   #Calculate velocity autocorrelation function
fix             5 all vector 1 c_vacf[4]
variable        vacf equal 0.33*${{dt}}*trap(f_5)

thermo       		$d
thermo_style custom step temp press v_myPxy v_myPxz v_myPyz v_v11 v_v22 v_v33 vol v_Jx v_Jy v_Jz v_k11 v_k22 v_k33 v_vacf

fix thermo_print all print $d "$(temp) $(press) $(v_myPxy) $(v_myPxz) $(v_myPyz) $(v_v11) $(v_v22) $(v_v33) $(vol) $(v_Jx) $(v_Jy) $(v_Jz) $(v_k11) $(v_k22) $(v_k33) $(v_vacf)" &
    append thermoNVE_{Name}_T${{T}}KP1atm.out screen no title "# temp press v_myPxy v_myPxz v_myPyz v_v11 v_v22 v_v33 vol v_Jx v_Jy v_Jz v_k11 v_k22 v_k33 v_vacf"

# save thermal conductivity to file
variable     kav equal (v_k11+v_k22+v_k33)/3.0
fix          fxave1 all ave/time $d 1 $d v_kav file lamda.txt

# save viscosity to a file
variable     visc equal (v_v11+v_v22+v_v33)/3.0
fix          fxave2 all ave/time $d 1 $d v_visc file visc.txt

# save diffusion coefficient to a file
fix          fxave3 all ave/time $d 1 $d v_vacf file diff_coeff.txt

dump        LAMMPS all custom $d NVE_Prod_{Name}_${{T}}KP1atm.lammpstrj id mol type xu yu zu mass q

run          {GKRuntime}

variable     ndens equal count(all)/vol
print        "Average viscosity: ${{visc}} [Pa.s] @ $T K, ${{ndens}} atoms/A^3"

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

cd {join(CWD, f'Run_${{PBS_ARRAY_INDEX}}')}
mpiexec ~/tmp/bin/lmp -in {SimName}_${{PBS_ARRAY_INDEX}}.lammps
""")
        
    # os.rename(join(CWD, f'{SimType}.pbs'),  f'{SimType}.pbs')}")

def GetMolMass(mol):
    formula = CalcMolFormula(mol)

    parts = re.findall("[A-Z][a-z]?|[0-9]+", formula)
    mass = 0

    for index in range(len(parts)):
        if parts[index].isnumeric():
            continue

        atom = Chem.Atom(parts[index])
        multiplier = int(parts[index + 1]) if len(parts) > index + 1 and parts[index + 1].isnumeric() else 1
        mass += atom.GetMass() * multiplier
    return mass

def CalcBoxLen(MolMass, TargetDens, NumMols):
    # Very conservative implementation of Packmol volume guesser
    BoxL = (((MolMass * NumMols * 2)/ TargetDens) * 1.75) ** (1./3.)
    BoxLRounded = int(BoxL)
    return BoxLRounded

def GeneratePDB(SMILES, PATH, CONFORMATTEMPTS=10):
    """
    Function that generates PDB file from RDKit Mol object, for use with Packmol
    Inputs:
        - SMILES: SMILES string to be converted to PDB
        - PATH: Location that the PDB will stored
        - CONFORMATTEMPTS: Max number of tries (x5000) to find converged conformer for molecule
    """
    SMILESMol = Chem.MolFromSmiles(SMILES) # Create mol object
    SMILESMol = Chem.AddHs(SMILESMol) # Need to make Hydrogens explicit

    AllChem.EmbedMolecule(SMILESMol, AllChem.ETKDG()) #Create conformer using ETKDG method

    # Initial parameters for conformer optimisation
    MMFFSMILES = 1 
    ConformAttempts = 1
    MaxConformAttempts = CONFORMATTEMPTS

    # Ensure that script continues to iterate until acceptable conformation found
    while MMFFSMILES !=0 and ConformAttempts <= MaxConformAttempts: # Checking if molecule converged
        MMFFSMILES = Chem.rdForceFieldHelpers.MMFFOptimizeMolecule(SMILESMol, maxIters=5000)
        ConformAttempts += 1
        
    # Record parameterised conformer as pdb to be used with packmol later 
    Chem.MolToPDBFile(SMILESMol, f'{PATH}')

CopyCommand = 'cp'
STARTINGDIR = deepcopy(getcwd())
PYTHONPATH = 'python3'
#PYTHONPATH = 'C:/Users/eeo21/AppData/Local/Programs/Python/Python310/python.exe'
Molecules = [x for x in listdir(STARTINGDIR) if '.pdb' in x]
NumMols = 200
NumRuns = 20
RunList = list(range(1, NumRuns+1))
# Values for the array job
TopValue = RunList[-1] 
BotValue = RunList[0]
LOPLS = False
WALLTIME = '23:59:59'

runcmd('mkdir ValidationStudies12ACutoff_200mols_LOPLS_NoKSPACE')

for Molecule in Molecules:
    FolderName = Molecule.split('.')[0]
    Path = join(STARTINGDIR, Molecule)
    MolObject = Chem.MolFromPDBFile(Path)

    if Molecule == '1-methylnapthalene.pdb':
        SMILESString = 'CC1=CC=CC2=CC=CC=C12'
    elif Molecule == 'squalane.pdb':
         SMILESString = 'CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C'
    elif Molecule == 'bis_2_ethylhexyl_sebacate.pdb':
        SMILESString = 'CCCCC(CC)COC(=O)CCCCCCCCC(=O)OCC(CC)CCCC'
    else:
        SMILESString = Chem.MolToSmiles(MolObject)
    
    if LOPLS:
        LTCOMMAND = f"{join(STARTINGDIR, 'rdlt.py')} --smi '{SMILESString}' -n {FolderName} -l -c"
    else:
        LTCOMMAND = f"{join(STARTINGDIR, 'rdlt.py')} --smi '{SMILESString}' -n {FolderName} -c"
  
    runcmd(f'{PYTHONPATH} {LTCOMMAND} > {STARTINGDIR}/{FolderName}.lt')

    GeneratePDB(SMILESString, PATH=join(STARTINGDIR, f'{FolderName}.pdb'))

    #Enter Trajectory studies directory
    chdir(join(getcwd(), 'ValidationStudies12ACutoff_200mols_LOPLS_NoKSPACE'))

    #Make Molecule Folder in Trajectory studies directory
    runcmd(f'mkdir {FolderName}')

    #Copy molecule pdb to molecule directory
    runcmd(f'{CopyCommand} "{join(STARTINGDIR, Molecule)}" {join(getcwd(), FolderName)}')

    #Enter molecule directory
    chdir(join(getcwd(), FolderName))

    CreateArrayJob(CWD=getcwd(), SimName=f'{FolderName}_system_313K', SimType=f'{FolderName}_313K',
                   TopValue=TopValue, BotValue=BotValue, WALLTIME=WALLTIME)
    
    CreateArrayJob(CWD=getcwd(), SimName=f'{FolderName}_system_373K', SimType=f'{FolderName}_373K',
                   TopValue=TopValue, BotValue=BotValue, WALLTIME=WALLTIME)

    ### Make files for short vs long run convergence
    """
    - So we could run 100 trajectories for each molecule up to 10ns with 100 molecules
    - We could generate trajectories for 500ps, 1ns, 2ns, and 4ns by
    truncating each of the runs to the respective times
    - We could then take between 5 and 10 random samples from the 100 trajectories
    at each time length for varying numbers of trajectories
    - Take an average of the samples as the predicted value, and use the extreme values 
    or calculate the standard deviation as the uncertainty (2std for 95% confidence)
    """

    MolMass = GetMolMass(MolObject)
    BoxL = CalcBoxLen(MolMass=MolMass, TargetDens=0.8, NumMols=NumMols)

    for x in RunList:
        runcmd(f'mkdir Run_{x}')
        Trajectory = f'Run_{x}'
        # Copy pdb
        runcmd(f'{CopyCommand} "{join(STARTINGDIR, Molecule)}" {join(getcwd(), Trajectory, Molecule)}')
        # Copy lt file
        ltfile = f'{FolderName}.lt'
        runcmd(f'{CopyCommand} "{join(STARTINGDIR, ltfile)}" {join(getcwd(), Trajectory, ltfile)}')
        # Copy settings file
        settingsfile = 'settings.txt'
        runcmd(f'{CopyCommand} "{join(STARTINGDIR, settingsfile)}" {join(getcwd(), Trajectory, settingsfile)}')

        # Set random seed
        Seed = rnd.randint(0, 1E6)

        # Make Packmol Input file
        MakePackmolFile(Name=FolderName, CWD=join(getcwd(), Trajectory), NumMols=NumMols, Seed=Seed, BoxL=BoxL)
        # Make Moltemplate file
        MakeMoltemplateFile(Name=FolderName, CWD=join(getcwd(), Trajectory), NumMols=NumMols, BoxL=BoxL)
        # Make Lammps Files
        MakeLAMMPSFile(Name=FolderName, ID=x, CWD=join(getcwd(), Trajectory), GKRuntime=3000000, Temp=313)
        MakeLAMMPSFile(Name=FolderName, ID=x, CWD=join(getcwd(), Trajectory), GKRuntime=3000000, Temp=373)
        # Make packmol coordinate file and LAMMPS data file
        chdir(join(getcwd(), Trajectory))
        if PYTHONPATH == 'python3':
            runcmd(f'packmol < {FolderName}.inp')
            runcmd(f'moltemplate.sh -pdb {FolderName}_PackmolFile.pdb {FolderName}_system.lt')

        # Return to starting directory
        chdir(join(STARTINGDIR, 'ValidationStudies12ACutoff_200mols_LOPLS_NoKSPACE', FolderName))
    chdir(STARTINGDIR)

    
# ## Files showing finite size effects
# """
# Run experiments at 25, 50, 100, 250, 500 molecules
# Run 10 simulations for each number of molecules for 10ns
# """

# Molecules = [x for x in listdir(STARTINGDIR) if '.pdb' in x]
# NumRuns = 10
# RunList = list(range(1, NumRuns+1))
# # Values for the array job
# TopValue = RunList[-1]
# BotValue = RunList[0]
# LOPLS = True
# WALLTIME = '72:00:00'

# runcmd('mkdir FiniteSizeEffects')

# for Molecule in Molecules:
#     FolderName = Molecule.split('.')[0]
#     Path = join(STARTINGDIR, Molecule)
#     MolObject = Chem.MolFromPDBFile(Path)

#     if Molecule == '1-methylnapthalene.pdb':
#         SMILESString = 'CC1=CC=CC2=CC=CC=C12'
#     elif Molecule == 'squalene.pdb':
#          SMILESString = 'CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C'
#     else:
#         SMILESString = Chem.MolToSmiles(MolObject)

#     if LOPLS:
#         LTCOMMAND = f"{join(STARTINGDIR, 'rdlt.py')} --smi '{SMILESString}' -n {FolderName} -l -c"
#     else:
#         LTCOMMAND = f"{join(STARTINGDIR, 'rdlt.py')} --smi '{SMILESString}' -n {FolderName} -c"
   
#     #Enter Trajectory studies directory
#     chdir(join(getcwd(), 'FiniteSizeEffects'))

#     #Make Molecule Folder in Trajectory studies directory
#     runcmd(f'mkdir {FolderName}')

#     #Copy molecule pdb to molecule directory
#     runcmd(f'copy "{join(STARTINGDIR, Molecule)}" {join(getcwd(), FolderName)}')

#     #Enter molecule directory
#     chdir(join(getcwd(), FolderName))

#     NumMols = [25, 50, 100, 250, 500]
#     for NumMol in NumMols:
#         MolMass = GetMolMass(MolObject)
#         BoxL = CalcBoxLen(MolMass=MolMass, TargetDens=0.8, NumMols=NumMol)
#         runcmd(f'mkdir NumMols_{NumMol}')
#         chdir(join(getcwd(), f'NumMols_{NumMol}'))

#         CreateArrayJob(CWD=getcwd(), SimName=f'{FolderName}_system_313K', SimType=f'{FolderName}_313K',
#                 TopValue=TopValue, BotValue=BotValue, WALLTIME=WALLTIME)
        
#         CreateArrayJob(CWD=getcwd(), SimName=f'{FolderName}_system_373K', SimType=f'{FolderName}_373K',
#                 TopValue=TopValue, BotValue=BotValue, WALLTIME=WALLTIME)
        
#         for x in RunList:
#             runcmd(f'mkdir Run_{x}')
#             Run = f'Run_{x}'
#             # Copy pdb
#             runcmd(f'{CopyCommand} "{join(STARTINGDIR, Molecule)}" {join(getcwd(), Run, Molecule)}')
#             # Copy lt file
#             ltfile = f'{FolderName}.lt'
#             runcmd(f'{CopyCommand} "{join(STARTINGDIR, ltfile)}" {join(getcwd(), Run, ltfile)}')
#             # Set random seed
#             Seed = rnd.randint(0, 1E6)

#             # Make Packmol Input file
#             MakePackmolFile(Name=FolderName, CWD=join(getcwd(), Run), NumMols=NumMol, Seed=Seed, BoxL=BoxL)
#             # Make Moltemplate file
#             MakeMoltemplateFile(Name=FolderName, CWD=join(getcwd(), Run), NumMols=NumMol, BoxL=BoxL)
#             # Make Lammps Files
#             MakeLAMMPSFile(Name=FolderName, ID=x, CWD=join(getcwd(), Run), GKRuntime=5000000, Temp=313)
#             MakeLAMMPSFile(Name=FolderName, ID=x, CWD=join(getcwd(), Run), GKRuntime=5000000, Temp=373)
#             # Make packmol coordinate file and LAMMPS data file
#             chdir(join(getcwd(), Run))
#             if PYTHONPATH == 'python3':
#                 runcmd(f'packmol < {FolderName}.inp')
#                 runcmd(f'moltemplate.sh -pdb {FolderName}_PackmolFile.pdb {FolderName}_system.lt')

#             # Return to starting directory
#             chdir(join(STARTINGDIR, 'FiniteSizeEffects', FolderName, f'NumMols_{NumMol}'))
    
#         chdir(join(STARTINGDIR, 'FiniteSizeEffects', FolderName))
    
#     chdir(STARTINGDIR)
    
    
# ### Files showing long term variations in viscosity and thermal conductivity prediction
# """
# - Let's run 10 simulations for each molecule up to 100ns, will be pretty robust for 
# uncertainty calculations.
# - Run at 40C and 100C 
# """

# Molecules = [x for x in listdir(STARTINGDIR) if '.pdb' in x]
# NumMols = 100
# NumRuns = 5
# RunList = list(range(1, NumRuns+1))
# # Values for the array job
# TopValue = RunList[-1]
# BotValue = RunList[0]
# LOPLS = True
# WALLTIME = '72:00:00'

# runcmd('mkdir LongSimEffects')

# for Molecule in Molecules:
#     FolderName = Molecule.split('.')[0]
#     Path = join(STARTINGDIR, Molecule)
#     MolObject = Chem.MolFromPDBFile(Path)

#     if Molecule == '1-methylnapthalene.pdb':
#         SMILESString = 'CC1=CC=CC2=CC=CC=C12'
#     elif Molecule == 'squalane.pdb':
#         SMILESString = 'CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C'
#     else:
#         SMILESString = Chem.MolToSmiles(MolObject)

#     print(SMILESString)

#     if LOPLS:
#         LTCOMMAND = f"{join(STARTINGDIR, 'rdlt.py')} --smi '{SMILESString}' -n {FolderName} -l -c"
#     else:
#         LTCOMMAND = f"{join(STARTINGDIR, 'rdlt.py')} --smi '{SMILESString}' -n {FolderName} -c"

#     runcmd(f'{PYTHONPATH} {LTCOMMAND} > {STARTINGDIR}/{FolderName}.lt')

#     MolMass = GetMolMass(MolObject)
#     BoxL = CalcBoxLen(MolMass=MolMass, TargetDens=0.8, NumMols=NumMols)
  
#     #Enter Long Sim Effects studies directory
#     chdir(join(getcwd(), 'LongSimEffects'))

#     #Make Molecule Folder in Trajectory studies directory
#     runcmd(f'mkdir {FolderName}')

#     #Copy molecule pdb to molecule directory
#     runcmd(f'{CopyCommand} "{join(STARTINGDIR, Molecule)}" {join(getcwd(), FolderName)}')

#     #Enter molecule directory
#     chdir(join(getcwd(), FolderName))

#     CreateArrayJob(CWD=getcwd(), SimName=f'{FolderName}_system_313K', SimType=f'{FolderName}_313K',
#                    TopValue=TopValue, BotValue=BotValue, WALLTIME=WALLTIME)
    
#     CreateArrayJob(CWD=getcwd(), SimName=f'{FolderName}_system_373K', SimType=f'{FolderName}_373K',
#                    TopValue=TopValue, BotValue=BotValue, WALLTIME=WALLTIME)

#     for x in RunList:
#         runcmd(f'mkdir Run_{x}')
#         Trajectory = f'Run_{x}'
#         # Copy pdb
#         runcmd(f'{CopyCommand} "{join(STARTINGDIR, Molecule)}" {join(getcwd(), Trajectory, Molecule)}')
#         # Copy lt file
#         ltfile = f'{FolderName}.lt'
#         runcmd(f'{CopyCommand} "{join(STARTINGDIR, ltfile)}" {join(getcwd(), Trajectory, ltfile)}')
#         # Set random seed
#         Seed = rnd.randint(0, 1E6)

#         # Make Packmol Input file
#         MakePackmolFile(Name=FolderName, CWD=join(getcwd(), Trajectory), NumMols=NumMols, Seed=Seed, BoxL=BoxL)
#         # Make Moltemplate file
#         MakeMoltemplateFile(Name=FolderName, CWD=join(getcwd(), Trajectory), NumMols=NumMols, BoxL=BoxL)
#         # Make Lammps Files
#         MakeLAMMPSFile(Name=FolderName, ID=x, CWD=join(getcwd(), Trajectory), GKRuntime=5000000, Temp=313)
#         MakeLAMMPSFile(Name=FolderName, ID=x, CWD=join(getcwd(), Trajectory), GKRuntime=5000000, Temp=373)
#         # Make packmol coordinate file and LAMMPS data file
#         chdir(join(getcwd(), Trajectory))
#         if PYTHONPATH == 'python3':
#             runcmd(f'packmol < {FolderName}.inp')
#             runcmd(f'moltemplate.sh -pdb {FolderName}_PackmolFile.pdb {FolderName}_system.lt')

#         # Return to starting directory
#         chdir(join(STARTINGDIR, 'LongSimEffects', FolderName))
#     chdir(STARTINGDIR)