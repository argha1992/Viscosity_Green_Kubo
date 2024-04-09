import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda

# Path for the new directory where the output files will be saved
output_directory = 'pxy_cor_cnt12_12'
# Create the output directory if it does not exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Open a log file in write mode
with open(os.path.join('pxy_cor_cnt12_12', 'process_log.txt'), 'w') as log_file:
    # Open a separate file for recording viscosities
    with open(os.path.join('pxy_cor_cnt12_12', 'viscosities.txt'), 'w') as viscosities_file:
    # List of the subdirectories
        subdirectories = [f"r{i}" for i in range(1, 101)]

        for subdir in subdirectories:
            log_message = f"Processing {subdir}..."
            print(log_message)
            log_file.write(log_message + '\n')  # Write the status to the log file

            # pressure.xvg, production.tpr, and production.xtc based on current subdirectory
            pressure_file = os.path.join(subdir, 'pressure.xvg')
            tpr_file = os.path.join(subdir, 'production.tpr')
            xtc_file = os.path.join(subdir, 'production.xtc')
            
            # To process pressures
            pressures = np.loadtxt(pressure_file, skiprows=24) 
            pressures = pressures[:6001]
            pressures = np.array(pressures[:, 1])
            times = len(pressures)
            pp_corr = []
            
            # To process pressure pressure correlation with time average
            for i in range(times//2):
                sum = 0                    # Initialize sum for each autocorrelation step
                count = 0
                for j in range(times):
                    if (j+i) >= times:
                        break
                    sum += pressures[j]*pressures[j+i]
                    count += 1
                sum /= count               # Compute average
                pp_corr.append(sum)
            integral_sum = np.sum(pp_corr) # multiplicaton of dt was done during viscosity calculation

            # MDAnalysis code to compute volumes
            u = mda.Universe(tpr_file, xtc_file)
            volumes = []
            for ts in u.trajectory[:6001]:
                box = ts.dimensions[:3]   # Get the dimensions of the current timestep
                volume = np.prod(box)     # Compute the volume
                volume /= 1000            # Convert A^3 to nm^3
                volumes.append(volume)
            average_volume = np.mean(volumes)


            # Calculation for viscosity at 310K
            k = 1.381*10**-23                                # Boltzmann constant in J/K
            T = 310                                          # Temperature in K
            delta_t = 0.02*(10**-12)                         # Time 1 ps = 10**-12 s 
            average_volume = average_volume*((10**-9)**3)    # nm to m conversion 1 nm = 10**-9 m 
            integral_sum = integral_sum*np.square(10**5)     # 1 bar = 10**5 Pa

            viscosity=(average_volume/ (k * T)) * delta_t * integral_sum
            print("Viscosity:", viscosity, "Pa·s")

            # Write the calculated viscosity into the viscosities file
            viscosities_file.write(f"{subdir}\t{viscosity}\n")

            # Write average_volume and pp_corr to a txt file
            output_file = os.path.join(output_directory, f'pp_corr_cnt12_12_{subdir}.txt')
            with open(output_file, 'w') as f:
                f.write(f"{average_volume}\n")
                for item in pp_corr:
                    f.write(f"{item}\n")
            
            completed_message = f"Completed {subdir}"
            print(completed_message)
            log_file.write(completed_message + '\n')  # Log completion of this subdirectory

        log_file.write("Processing completed for all directories.\n")  # Log completion of entire process



