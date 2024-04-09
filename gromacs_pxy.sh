#!/bin/bash

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
