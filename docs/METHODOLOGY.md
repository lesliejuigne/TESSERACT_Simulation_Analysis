# TESSSA Methodology

This document describes the methodology used in **TESSSA** (TESSERACT Simulation Analysis) to convert raw **GEANT4 simulation outputs** into physically meaningful datasets and plots.

---

## 1. Introduction

TESSSA is designed to analyze GEANT4 simulations for the TESSERACT experiment. Its main goals are:

- Filter simulated events to retain only those in the **zone of interest** (the *virtual detector*).  
- Normalize the filtered events into physical units (**Counts / (keV·kg·day)**).  
- Generate spectra and other analysis plots for background and shielding studies.  

Simulations include:

- **Internal components** (e.g., detector materials, cryostat, shielding layers)  
- **External sources** (rock gammas and neutrons)  

## 2. Normalization Methodology

Normalization converts raw event counts into **physical units**. TESSSA treats **internal components** and **external sources** differently.

### 2.1 Internals Components

The scaling factor is:

$$
f = \frac{n_{bins} \cdot A_{iso} \cdot m_{vol}}{\Delta E \cdot N_{decays} \cdot M_{detector}} \cdot \frac{3600.24}{1 \text{ day}}
$$

Where:

| Symbol | Meaning |
|--------|---------|
| $n_{bins}/\Delta E$ | Binning factor |
| $A_{iso}$ | Activity of isotope (Bq/kg) |
| $m_{vol}$ | Mass of simulated volume |
| $N_{decays}$ | Number of simulated decays (`beamOn * Nfiles`) |
| $M_{detector}$ | Detector's mass |

- Activity values are from `cache/material_data.csv` (radiopurity.org and literature).  
- Conversion factor $3600.24/1\text{ day}$ converts seconds to days.  
- Output units: **Counts / (keV·kg·day)**

---

### 2.2 External Shield – Rock (Gammas & Neutrons)

The scaling factor for rock emissions is:

$$
f = \frac{n_{bins}}{\Delta E} \cdot \frac{\phi \cdot S}{N_{decays} \cdot M_{detector}} \cdot \frac{3600.24}{1 \text{ day}}
$$

Where:

- $S$ = Irradiating surface 
- $\phi$ = Total measured flux at LSM  
- $n_{decays}$ = Number of simulated decays  
- $M_{detector}$ = Detector's mass  

#### Isotope Normalization

Since the total flux $\phi$ is measured, TESSSA normalizes per isotope:

$$
n_{norm/iso} = \frac{C_{iso} \cdot r_{iso}}{\sum_i C_i \cdot r_i}, \quad r_{iso} = \frac{N_{\gamma/iso}}{N_{decays/iso}}
$$

- $C_{iso}$ = Relative concentration of isotope in rock  
- $r_{iso}$ = Emission probability per decay (from Geant4 simulations)  

Normalized flux per isotope:

$$
\phi_{iso} = \phi_{LSM} \cdot n_{norm/iso}
$$

Example table:

| Isotope | $C_{iso}$ (Bq/kg) | $r_{iso}$ (γ/decay) | $n_{norm/iso}$ | $\phi_{iso}$ (γ·cm⁻²·s⁻¹) |
|---------|-----------------|------------------|----------------|------------------|
| K-40    | 182             | 2.48e-2          | 0.06           | 0.035            |
| Th-232  | 10.2            | 3.779            | 0.49           | 0.306            |
| U-238   | 11.8            | 2.99             | 0.45           | 0.280            |

---

## 3. Parameter Management

- Cached data are stored in `cache/`:

  - `materials.json` – material densities, composition, and notes  
  - `material_data.csv` – isotope activities and uncertainties for the internal shielding
  - `rock_data.csv` – isotope activities and uncertainties for the external shielding

- These parameters are used automatically during normalization.

---

## 4. Uncertainty Estimation

At the moment the uncertainty estimation is made using the available data used for the isotope activities. 

---

## 5. Output & Analysis

- Normalized datasets (DataFrame or `.h5`)  
- Plots: energy spectra, background comparisons

---

## 6. References

1. C. de Dominicis & M. Settimo, *Simulations and background estimate for the DAMIC-M experiment*, Il Nuovo Cimento C, **45** (2022), 6. [https://doi.org/10.1393/ncc/i2022-22006-y](https://doi.org/10.1393/ncc/i2022-22006-y)  
2. Radiopurity database: [https://radiopurity.org](https://radiopurity.org)  
3. GEANT4 Collaboration, *GEANT4 – A Simulation Toolkit*, Nucl. Instrum. Meth. A 506 (2003) 250-303
