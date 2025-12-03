# TESSSA  
**TESSERACT Simulation Analysis Framework**

TESSSA (TESSERACT Simulations Analysis) is a lightweight framework designed to process and analyze simulation data produced with **GEANT4** for the TESSERACT experiment.  

It provides a simple, reproducible pipeline to go from raw GEANT4 simulation outputs to **normalized datasets and plots** in physical units.

## Prerequisites

Before running TESSSA, make sure the following software and modules are installed:

### 1. Python dependencies
- Python 3.9+  
- Core Python packages (install with `pip install -r requirements.txt`):

```txt
numpy>=1.21
pandas>=1.3
matplotlib>=3.4
uproot>=4.1
tqdm>=4.62

### 2. ROOT
ROOT version 6.26 or higher recommended
ROOT must be installed and accessible in your shell $PATH

Check version with:

```bash
root --version

## Installation

Clone the repository:

```bash
git clone https://github.com/lesliejuigne/TESSSA.git
cd TESSSA

There are two main installation modes:

Development mode (recommended for contributors)
Use this if you plan to modify the code:

```bash
pip install -e .

User mode (recommended for end users)
Use this if you only want to run TESSSA without changing the source code:

```bash
pip install .

## Dependencies

Install all dependencies with:

```bash
pip install -r requirements.txt

---

## 📌 Workflow Overview

1. **Simulation (GEANT4)**  
   - Run detector and shielding simulations with GEANT4.  
   - Produce raw ROOT files, typically in formats such as:  
     - **Internals:** `layer_isotope_nº.root` (example: `Cu_K40_2.root`)  
     - **Rock gammas:** `layer_particle_isotope_nº.root` (example: `Rock_Gammas_K40_3.root`) 
     Please note that the current nomenclature is a personal choice, a new nomenclature easily be implemented to match other simulations needs. 

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
     - (In development) stored as `.h5` files for larger-scale workflows.  

4. **Analysis**  
   - Use the example Jupyter notebooks or custom scripts to:  
     - Plot energy spectra, background rates, and shielding comparisons.  
     - Export datasets for further use.  

---

## 📂 Repository Structure
TESSSA/
│
├── filtering/ # Codes to filter raw ROOT files
│	├── cpp/
│		├── VD_single_filtering.cc # Keeps only events in the detector zone for 1 file
│		├── VD_group_filtering.cc # Keeps only events in the detector zone and merge all the files
│     └── ...
└── bash_script_templates/
│		├── single_filtering_slurm.sh
│		└── group_filtering_slurm.sh
│
├── processing/ # Normalization & post-filtering analysis
│ ├── tesssapy/ 
│		├──sim_processing.py # Converts filtered data into normalized DataFrames
│		└── ...
│	└── example_notebook.ipynb # Demo notebook for processing & plotting
│
├── cache/ # Cached input data and normalization parameters
│ ├── materials.json # Material properties
│ ├── material_data.csv # Isotope activities & uncertainties for the internal shielding
│ └── rock_data.csv # Isotope activities & uncertainties for the external shielding
│
├── docs/ # Documentation
│ ├── METHODOLOGY.md # Normalization formulas and methodology
│	└── TESSSA_v1_presentation # Presentation of the v1 | How it works and initial results 
│
└── README.md # Project overview (this file)

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
# 1. Filter raw GEANT4 data
root -l -b -q "VD_single_filtering.cc("input_file.root","output_folder/")"

# 2. Process and normalize filtered data
python -c "import TESSSA.sim_processing as ssp; s = ssp.g4_sim_proc('layer', 'path/to/data', plots=True)"🧩 Dependencies
