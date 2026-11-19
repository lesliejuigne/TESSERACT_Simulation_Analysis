import uproot
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import os
import glob
import pandas as pd
from tqdm import tqdm
import warnings
from .utils import get_cached_data, load_style_file
from .get_norm_param import GetNormParam

warnings.simplefilter("ignore")

materials = get_cached_data("materials_data.csv")
rock      = get_cached_data("rock_data.csv")

plt.style.use(load_style_file('SetStyle_mplstyle.txt'))

class g4_sim_proc:

    def __init__(self, compoment, folder_path, bias="boff",  plots= True):
        
        # ======== parameters ========
        self.compoment = compoment      # internals, rock, concrete      
        self.folder_path = folder_path  # path to the folder containing the root files
        print(f"Processing files in {self.folder_path}")
        
        # ======== constants ========
        self.detMass = 0.797336 
        self.SecPerDay = float(3600 * 24)
        self.xLow = 0
        self.xHigh = 5000
        self.Bins = 20 
        self.BinSize = (self.xHigh - self.xLow) / self.Bins
        self.Bin = float(1.0 / self.BinSize)
        
        # ======== variables ========
        self.data = {}
        self.data_counts = {}
        self.counts = {}
        self.counts_err = {}
        self.energy = {}
        self.bias = bias
        # ======== functions calling ========
        self.load_raw_data()
         
        self.normalize_data()
        self.get_totals()
        if plots == True : 
            print("Generating plots...")
            self.get_spectrum()
            print("Generating total spectrum plot...")
            self.get_spectrum_totals()
        self.print_simulation_summary()
        
    def get_root_tree(self, file_path):
        try:
            root_file = uproot.open(file_path)
            tree = root_file["events"]
            Params = tree.keys()
            
            if not tree:
                return []
            edep = np.array(tree.arrays(Params[16])[Params[16]])
            return edep
        
        except Exception as e:
            return []
    
    def hist_it(self, X):
        counts, edges = np.histogram(X, self.Bins,[self.xLow,self.xHigh]) #, bins = logbins)
        centers = (edges[:-1] + edges[1:]) / 2.0
        return centers, counts, edges
    
    def load_raw_data(self):
        
        # ======== functions ========
        def get_layers_and_isotopes(component):
            if component == "internals":
                layers = materials["Material"].unique().tolist()
                get_isotopes = lambda l: materials[materials["Material"] == l]["Isotope"].unique().tolist()
            else:
                material_type = "Rock" if component == "rock" else "Concrete"
                layers = rock[rock["Material"] == material_type]["Particule"].unique().tolist()
                get_isotopes = lambda l: rock[rock["Particule"] == l]["Isotope"].unique().tolist()
            return layers, get_isotopes

        def build_filepath(component, layer, iso, i,bias="boff"):
            if component == "internals":
                return f"{self.folder_path}/{layer}_{iso}_{i}_{bias}_filtered.root"

            prefix = "Rock" if component == "rock" else "Concrete"
            
            if layer == "Neutrons":
                return f"{self.folder_path}/{prefix}_{layer}_{i}_{bias}_filtered.root"
            
            return f"{self.folder_path}/{prefix}_{layer}_{iso}_{i}_{bias}_filtered.root"
        
        # ======== main code ========
        layers, get_isotopes = get_layers_and_isotopes(self.compoment)
        if self.compoment != "internals":
            self.particule = layers

        total_files = sum(len(glob.glob(f"{self.folder_path}/{layer}*.root")) for layer in layers)

        with tqdm(total=total_files, desc="Processing Files", unit="file") as pbar:
            for layer in layers:
                isotopes = get_isotopes(layer)
                self.data[layer] = {iso: [] for iso in isotopes}
                self.data_counts[layer] = {iso: 0 for iso in isotopes}
                for iso in isotopes:
                    for i in range(300):
                        file_path = build_filepath(self.compoment, layer, iso, i,self.bias)

                        if not os.path.exists(file_path):
                            continue
                        data = self.get_root_tree(file_path)
                        if data is not None and len(data) > 0:
                            self.data[layer][iso].extend(data)
                            self.data_counts[layer][iso] += 1
                        pbar.update(1)

        print("Data loading complete.", self.data_counts)

    def normalize_data(self):
        # Load shielding and particle names from macros instead of CSVs
        self.shielding = materials["Material"].unique().tolist()
        self.particle = rock["Particule"].unique().tolist()

        def compute_normalized_counts(X, Y, norm, norm_err):
            return np.multiply(Y, norm), np.multiply(Y, norm_err), X

        def get_parameters(component, layer, iso, file_path):
            """
            Return normalization parameters based on macros inside the ROOT file.

            Robust to the PureCu <-> PCu naming mismatch and to pandas selection issues.
            """

            # --- name mapping: what to use when querying geometry vs materials ---
            # geometry (GetNormParam) stores "PureCu"
            # your materials CSV uses "PCu"
            geo_name_map = {"PCu": "PureCu", "PureCu": "PureCu"}    # what geometry expects
            mat_name_map = {"PCu": "PCu", "PureCu": "PCu"}          # what materials table expects

            # Normalize input (strip whitespace)
            layer_in = layer.strip() if isinstance(layer, str) else layer
            iso_in = iso.strip() if isinstance(iso, str) else iso

            # Decide names to use for lookup
            geo_layer = geo_name_map.get(layer_in, layer_in)
            mat_layer = mat_name_map.get(layer_in, layer_in)

            if component == "internals":
                # --- Mass from geometryTable ---
                geo = GetNormParam(file_path, "geometryTable")
                masses = getattr(geo, "total_mass", {}) or {}
                mass_val = masses.get(geo_layer, 0.0)

                # --- BeamOn from runMacro ---
                run = GetNormParam(file_path, "runMacro")
                beamon = getattr(run, "beamon_number", 1) or 1

                # --- Lookup activity and sigma in materials dataframe ---
                
                sel = materials[
                    (materials["Material"] == mat_layer) & (materials["Isotope"] == iso_in)
                ]

                if sel.empty:
                    # no matching row
                    print(f"No data for {mat_layer} {iso_in}: materials selection returned empty")
                    return mass_val, beamon, None, None

                # If multiple rows, pick the first (or you can decide to aggregate)
                if len(sel) > 1:
                    print(f"Warning: multiple rows found for {mat_layer} {iso_in}, using the first one.")

                try:
                    # use iat to extract the scalar value from the column directly
                    activity = float(sel["Activity"].iat[0])
                    sigma = float(sel["Sigma"].iat[0])
                except Exception as e:
                    print(f"Error converting Activity/Sigma to float for {mat_layer} {iso_in}: {e}")
                    return mass_val, beamon, None, None

                return mass_val, beamon, activity, sigma

            else:
                # Rock or Concrete: beamOn also from macro
                run = GetNormParam(file_path, "runMacro")
                beamon = getattr(run, "beamon_number", 1) or 1

                mat = "Rock" if component == "rock" else "Concrete"
                
                sel = rock[
                    (rock["Material"] == mat) &
                    (rock["Particule"] == layer_in) &
                    (rock["Isotope"] == iso_in)
                ]

                if sel.empty:
                    print(f"No data for {mat} {layer_in} {iso_in}: rock selection returned empty")
                    return None, beamon, None, None

                if len(sel) > 1:
                    print(f"Warning: multiple rows found for {mat} {layer_in} {iso_in}, using the first one.")

                try:
                    surface = float(sel["Surface"].iat[0])
                    flux = float(sel["Flux"].iat[0])
                    sigma = float(sel["Sigma"].iat[0])
                except Exception as e:
                    print(f"Error converting Surface/Flux/Sigma to float for {mat} {layer_in} {iso_in}: {e}")
                    return None, beamon, None, None

                return surface, beamon, flux, sigma

        # --- loop as before ---
        if self.compoment not in ["internals", "rock", "concrete"]:
            print("Unknown component:", self.compoment)
            return

        layers = self.shielding if self.compoment == "internals" else self.particle

        for layer in layers:
            self.counts[layer], self.counts_err[layer], self.energy[layer] = {}, {}, {}
            isotopes = (
                materials[materials["Material"] == layer]["Isotope"].tolist()
                if self.compoment == "internals"
                else rock[rock["Particule"] == layer]["Isotope"].tolist()
            )
            for iso in isotopes:
                try:
                    X, Y, self.edges = self.hist_it(self.data[layer][iso])

                    # Grab the first file path we used for this layer/iso
                    file_path = f"{self.folder_path}/{layer}_{iso}_1_boff_filtered.root" \
                        if self.compoment == "internals" \
                        else f"{self.folder_path}/{('Rock' if self.compoment=='rock' else 'Concrete')}_{layer}_{iso}_1_boff_filtered.root"

                    p1, p2, p3, p4 = get_parameters(self.compoment, layer, iso, file_path)

                    normalization_factor = (
                        p1 * self.SecPerDay *
                        (1.0 / (p2 * self.data_counts[layer][iso])) *
                        self.Bin * (1.0 / self.detMass)
                    )
                    norm, norm_err = normalization_factor * p3, normalization_factor * p4

                    self.counts[layer][iso], self.counts_err[layer][iso], self.energy[layer][iso] = \
                        compute_normalized_counts(X, Y, norm, norm_err)

                except Exception as e:
                    print(f"No data for {layer} {iso}:", e)
                    continue

    def get_totals(self):
        self.counts['total'] = 0
        self.counts_err['total'] = 0
        for layer in self.counts:
            if layer == 'total':
                continue
            self.counts[layer]['total'] = np.sum([v for k, v in self.counts[layer].items() if k != 'total'], axis=0)
            self.counts_err[layer]['total'] = np.sqrt(np.sum([v**2 for k, v in self.counts_err[layer].items() if k != 'total'], axis=0))
            self.counts['total'] += self.counts[layer]['total']
            self.counts_err['total'] = np.sqrt(self.counts_err['total']**2 + self.counts_err[layer]['total']**2)
    
    def get_spectrum(self):
        def plot_individual_spectrum(energy, counts, counts_err, label, color, ax):
            ax.step(energy, counts, where='mid', color=color, label=label)
            ax.step(energy, counts_err, where='mid', color=color)
            if counts[0] != 0.:
                ax.step(energy, -counts_err, where='mid', color=color)
                ax.fill_between(energy, counts - counts_err, counts + counts_err, color=color, step='mid', alpha=0.4)
                ax.scatter(energy, counts, color=color, marker='.')
                ax.errorbar(energy, counts, yerr=np.sqrt(counts), fmt='.', color=color, capsize=2)
            else:
                ax.fill_between(energy, counts + counts_err, color=color, step='mid', alpha=0.2)
                ax.scatter(energy, counts_err, color=color, marker='.')
                ax.errorbar(energy, counts_err, yerr=np.sqrt(counts_err), fmt='.', color=color, capsize=2)

        try:
            for k in self.counts:
                energy_key = 'Cu' if self.compoment == 'internals' else 'Gammas'
                energy = self.energy[energy_key]['K40']
                if k == 'total':
                    continue
                fig, ax = plt.subplots(figsize=(12, 8))
                cmap = cm.get_cmap('viridis', len(self.counts[k]))

                for i, j in enumerate(self.counts[k]):
                    color = cmap(i)
                    plot_individual_spectrum(energy, self.counts[k][j], self.counts_err[k][j], j, color, ax)

                ax.set_yscale('log')
                ax.set_xlabel("Energy [keV]", fontsize=20)
                ax.set_ylabel("Counts / (keV.kg.day)", fontsize=20)
                ax.set_title(f"{self.compoment} {k}", fontsize=20)
                ax.grid()
                #ax.set_xlim(0,2000)
                ax.legend(loc='upper right', fontsize=20)
                plt.show()
        except Exception as e:
            print(f"Error in plotting spectrum: {e}")
            return

    def get_spectrum_totals(self):
        fig, ax = plt.subplots(figsize=(12, 8))
        def plot_total(energy, counts, counts_err, label, color, ax):
            ax.step(energy, counts, where='mid', color=color, label=label)
            ax.fill_between(energy, counts + counts_err, counts - counts_err,
                            step='mid', alpha=0.2, color=color)
            ax.errorbar(energy, counts, yerr=counts_err,
                        fmt='.', capsize=2, color=color)

        try:
            # --------- Energy selection ---------
            energy_key = 'Cu' if self.compoment == 'internals' else 'Gammas'
            energy = self.energy[energy_key]['K40']

            plot_total(energy, self.counts['total'],self.counts_err['total'],label='Total',color='black', ax=ax)
            cmap = cm.get_cmap('viridis', len(self.counts)-1)
            for layer in self.counts:
                if layer != 'total':
                    color = cmap(list(self.counts.keys()).index(layer)-1)
                    plot_total(energy,self.counts[layer]['total'],self.counts_err[layer]['total'],label=f"{layer} ",color=color ,ax=ax)
                
            # --------- Plot formatting ---------
            ax.set_yscale('log')
            ax.set_xlabel("Energy [keV]", fontsize=20)
            ax.set_ylabel("Counts / (keV·kg·day)", fontsize=20)
            ax.set_title(f"{self.compoment}", fontsize=20)
            ax.grid()
            ax.legend(loc='upper right', fontsize=16)

            plt.show()

        except Exception as e:
            print(f"Error in plotting spectrum totals: {e}")
    
    def print_simulation_summary(self):
        print("\nSimulation Summary")
        print("===========================================")
        print(f"Shielding Type: {self.compoment.capitalize()} \n")
        print("Total counts per Layer of Shielding:")
        print("____________________________________\n")
        for layer in self.counts:
            if layer == 'total':
                continue
            print(f"{layer}: R = {np.sum(self.counts[layer]['total']):.2e} ± {np.sum(self.counts_err[layer]['total']):.2e} counts/keV.day.kg")
        print("\nTotal Counts:")
        print("_______________\n")
        print(f"R = {np.sum(self.counts['total']):.2e} ± {np.sum(self.counts_err['total']):.2e} counts/keV.day.kg")  
