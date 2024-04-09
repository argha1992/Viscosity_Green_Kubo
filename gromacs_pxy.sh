#!/bin/bash

# This script automates the extraction of pressure data from molecular dynamics (MD) simulations
# using GROMACS, a widely used software package for simulating proteins, lipids, and nucleic acids.
# It is designed to work across multiple simulation runs organized in subdirectories named r1, r2, ..., r100.
# For each subdirectory, the script:
# 1. Changes into the subdirectory.
# 2. Executes the 'gmx_mpi energy' command from GROMACS to extract pressure data.
#    The command uses 'production.edr' as the input energy file (.edr) and 'production.tpr' for the run input file,
#    and outputs the pressure data to 'pressure.xvg'.
#    It assumes the option for extracting pressure data (pressure-xy or offdiagonal components) corresponds to '22', which might need adjustment
#    depending on the specific setup of the energy groups in the GROMACS energy file.
# 3. Returns to the parent directory before moving on to the next subdirectory.
# This script simplifies the process of extracting relevant simulation data for further analysis or post-processing.

# Define the list of subdirectories
subdirs=(r{1..100})

# Loop over each subdirectory
for subdir in "${subdirs[@]}"; do
    echo "Processing $subdir..."
    
    # Change to the subdirectory
    cd "$subdir"
        # Run the GROMACS energy extraction command
    # Echoing 'Pressure' option (here assumed to be 8, replace with correct number if different)
    echo 22 | gmx_mpi energy -f production.edr -s production.tpr -o pressure.xvg
    
    # Return to the parent directory
    cd ..
done
