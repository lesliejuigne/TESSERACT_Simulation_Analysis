# SimRunner Workflow Overview

The SimRunner allow us to run jobs using the TESSERACT container for G4 simulations.
The SimTester.py allow us to check the Stability of the G4 Simulation Code. 

The SimRunner workflow consists of four main scripts and a configuration file:

1. **`job_configs.py`**
   - Stores all job parameters: internals, rock, concrete, isotopes, bias options.
   - Central place to define what simulations need to run.

2. **`launch_jobs.py`**
   - Reads `job_configs.py`.
   - Loops over materials, isotopes, and bias options.
   - Submits SLURM jobs using `python launch_jobs.sh`.

3. **`submit_simulation.sh`**
   - SLURM array job script.
   - Receives arguments from `launch_jobs.py`: `JOB_TYPE`, `MATERIAL`, `ISOTOPE`, `BIAS`.
   - Determines log/output directories based on job type.
   - Calls `run_in_container.sh` to run the actual simulation.

4. **`run_in_container.sh`**
   - Executes inside the container.
   - Loads Geant4 environment and required data paths.
   - Runs the `TesseractSim` executable with the specified macro, bias, and output file.

## Flow Summary

- `launch_jobs.py` → submits jobs → `submit_simulation.sh` → calls → `run_in_container.sh` → runs simulation.
- Output files and logs are separated by material/isotope and job type, avoiding overwrites.