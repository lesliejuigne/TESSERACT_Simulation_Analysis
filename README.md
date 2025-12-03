# TESSSA  
**TESSERACT Simulation Analysis Framework**

TESSSA (TESSERACT Simulations Analysis) is a lightweight framework designed to run, process and analyze simulation data produced with **GEANT4** for the TESSERACT experiment.  

It provides a simple, reproducible pipeline to go from raw GEANT4 simulation outputs to **normalized datasets and plots** in physical units.

## Prerequisites

Before running TESSSA, make sure the following software and modules are installed:

### 1. Python dependencies
- Python 3.9+  
- Core Python packages (install with `pip install -r requirements.txt`)

### 2. ROOT
ROOT version 6.26 or higher recommended
ROOT must be installed and accessible in your shell $PATH

Check version with:

```bash
root --version
```

## Installation

Clone the repository:

```bash
git clone https://github.com/lesliejuigne/TESSSA.git
cd TESSSA
```
There are two main installation modes:

Development mode (recommended for contributors)
Use this if you plan to modify the code:

```bash
pip install -e .
```

User mode (recommended for end users)
Use this if you only want to run TESSSA without changing the source code:

```bash
pip install .
```

## Dependencies

Install all dependencies with:

```bash
pip install -r requirements.txt
```
---

## Workflow Overview

1. **Simulation (GEANT4)**  
   - Run detector and shielding simulations with GEANT4 on the container.  
   - Produce raw ROOT files, typically in formats such as:  
     - **Internals:** `layer_isotope_nº.root` (example: `Cu_K40_2.root`)  
     - **Rock gammas:** `layer_particle_isotope_nº.root` (example: `Rock_Gammas_K40_3.root`) 
     Please note that the current nomenclature is a personal choice, a new nomenclature can easily be implemented to match other simulations needs. 

2. **Filtering**  
   - Use TESSSA’s filtering codes to select only events in the **zone of interest** (mainly the *virtual detector*).  
   - Output is a reduced file with `_filtered.root` appended.  
   - Example:  
     ```
     Cu_K40_2.root  →  Cu_K40_2_filtered.root
     Rock_Gammas_K40_3.root  →  Rock_Gammas_K40_3_filtered.root
     ```

3. **Processing & Normalization**  
   - The filtered files are passed to `sim_processing.py`.  
   - Events are normalized into **Counts / (keV·kg·day)** using isotope activities and detector parameters.  
   - Output can be:  
     - A **pandas DataFrame** (used in the example notebook).

4. **Analysis**  
   - Use the example Jupyter notebooks or custom scripts to:  
     - Plot energy spectra, background rates, and shielding comparisons.  
     - Export datasets for further use.  

---

## Repository Structure
```text
TESSSA/
│
├── filtering/                           # Codes to filter raw ROOT files
│   ├── cpp/
│   │   ├── VD_single_filtering.cc       # Keep only events in detector zone (single file)
│   │   ├── VD_group_filtering.cc        # Filter + merge all files
│   │   └── ...
│   └── bash_script_templates/
│       ├── single_filtering_slurm.sh
│       └── group_filtering_slurm.sh
│
├── processing/                          # Normalization & post-filtering analysis
│   ├── tesssapy/
│   │   ├── sim_processing.py            # Build normalized DataFrames from filtered data
│   │   ├── get_norm_parameter.py        # Get mass & initial beamOn counts from simulations
│   │   └── ...
│   ├── cache/                           # Cached inputs & normalization parameters
│   │   ├── material_data.csv            # Internal shielding isotope activities & uncertainties
│   │   ├── rock_data.csv                # External shielding isotope activities & uncertainties
│   │   ├── SetStyle_mplstyle.txt        # Matplotlib plot configuration
│   │   └── ...
│   └── tesssa_demo.ipynb                # Example notebook (processing & plotting)
│
├── SimRunner/                           # Simulation launcher (GEANT4 pipeline)
│   ├── lil_readme.md                    # Explanation of the mini-pipeline
│   ├── job_config.py                    # Choose what to simulate
│   ├── launch_jobs.py                   # Orchestrates & runs simulations
│   ├── run_in_container.sh              # Run simulation inside container
│   ├── submit_jobs.sh                   # SLURM submission script
│   └── SimTester.py                     # Automatic test for GEANT4 installation stability
│
├── docs/                                # Documentation
│   ├── METHODOLOGY.md                   # Normalization formulas & methodology
│   └── TESSSA_v1_presentation           # v1 presentation (workflow & initial results)
│
└── README.md                            # Project overview (this file)

```
---

## Normalization Methodology
TESSSA converts GEANT4 counts into physical units using scaling formulas that depend on whether the background is **internal** or **external**:  

- **Internals (materials & components):**  
  Uses isotope activity ($A_{iso}$), detector mass, and simulated decays.  

- **Externals (rock gammas & neutrons at LSM):**  
  Uses measured fluxes, isotopic concentrations, and Geant4-derived emission rates.  

Full formulas, tables, and references are in [`docs/METHODOLOGY.md`](docs/METHODOLOGY.md).  

---

## Command Example
```bash
# 1. Run Simulations
python3 launch_jobs.py

# 2. Filter raw GEANT4 data
root -l -b -q "VD_single_filtering.cc("input_file.root","output_folder/")"

# 3. Process and normalize filtered data
python -c "import TESSSA.sim_processing as ssp; s = ssp.g4_sim_proc('layer', 'path/to/data', plots=True)"
