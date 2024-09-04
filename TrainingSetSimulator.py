"""
Author: Egheosa Ogbomo
Date: 9th September 2024

Script to Generate Test Dataset using random mutations and LOPLS and COMPSS

Plan:

1. Create pool of known molecules of different sizes
2. Perform a molecule crossover and or one of the other mutations
3. Sequentially check a single  master data file to prevent duplication
4. Submit one large array job with different agents 25 lots of 2000 simulations(?)
5. Simulate with COMPASS for linear and branched molecules and esters, LOPLS for cyclic molecules and nitrogen containing molecules(?)
6. Can we write a script to replace a molecule in an array job on the fly if it doesn't work, or replace it with molecules from the Kajita
dataset or should we just simulated a selection of the Kajita dataset and remove that from the large dataset, then generate like 500 N containing
molecules?

Constraints

- Allow a maximum of 1 N atom per molecule
- Max number of heavy atoms
- No more than 2 rings
- Limit the number of double bonds to 3 (after the production run)
- Preference for aromatics with low viscosity, high thermal conductivity allowing for less competition to 
surface 


"""