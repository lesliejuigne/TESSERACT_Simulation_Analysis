import uproot
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import os
import glob
import pandas as pd
from tqdm import tqdm
import warnings
import json
from .utils import get_cached_data, load_style_file
from .get_norm_param import GetNormParam

warnings.simplefilter("ignore")

materials = get_cached_data("materials_data_LUX+cosmo.csv")
rock      = get_cached_data("rock_data.csv")

plt.style.use(load_style_file('SetStyle_mplstyle.txt'))

class G4SimProcessor:
    def __init__(self, compoment, folder_path, bias="boff", target='LHe'):
        # ======== parameters ========
        self.compoment = compoment.lower()      
        self.folder_path = folder_path  
        self.bias = bias
        self.target = target
        
        # ======== constants ========
        self.detMass = 0.199866 
        self.SecPerDay = float(3600 * 24)
        self.xLow = 0
        self.xHigh = 4000
        self.Bins = 4000 
        self.BinSize = (self.xHigh - self.xLow) / self.Bins
        self.Bin = float(1.0 / self.BinSize)
        
        # ======== variables ========
        self.data = {}
        self.data_counts = {} # This will now store TOTAL events simulated
        self.counts = {}
        self.counts_err = {}
        self.energy = {}
        self.raw_norm = {}
        self.summary_data = []
        
        # ======== Execution Flow ========
        print(f"Processing {self.compoment} in {self.folder_path}...")
        self.load_raw_data()
        
        if self.compoment == "internals":
            self.process_internals_data() 
        if self.compoment in ["rock", "concrete"]:
            self.process_externals_data()
            
        self.print_summary()
         
    #======================================================#
    
    def _get_root_tree(self, file_path):
        try:
            with uproot.open(file_path) as root_file:
                tree = root_file["events"]
                keys = tree.keys()
                
                # Identify the energy deposition column
                edep_key = [k for k in keys if "edep" in k.lower() or "energy" in k.lower()][0]
                # Fetch edep, eventID, AND trackWeight
                data = tree.arrays([edep_key, "eventID", "trackWeight"], library="pd")
                
                # --- GROUP BY EVENT ---
                # 1. Sum the energy deposits for the event
                # 2. Take the mean of the trackWeight (assuming weight is constant per event)
                
                grouped = data.groupby("eventID").agg({
                edep_key: "sum",
                "trackWeight": "mean"
                })
                event_data = grouped[edep_key].values
                
                # --- WEIGHT HANDLING ---
                if self.bias == "boff":
                    # If bias is off, we ignore trackWeight and use 1.0
                    weights = np.ones_like(event_data)
                else:
                    # If bias is on, use the actual track weights from the simulation
                    weights = grouped["trackWeight"].values
                return event_data, weights
            
        except Exception as e:
            print(f"Erreur lors de la lecture du fichier {file_path}: {e}")
        return [], [] 
    
    def _hist_it(self, X):
        """
        Computes the binned count of energy depositions (Spectrum).

        Input: Array (X) of energy values.
        Output: Tuple of (Bin Centers, Bin Counts, Bin Edges).
        
        Note: Binning resolution and range are defined by class attributes (Bins, xLow, xHigh).
        """
        counts, edges = np.histogram(X, bins=self.Bins, range=[self.xLow, self.xHigh])
        centers = (edges[:-1] + edges[1:]) / 2.0
        return centers, counts, edges
    
    def load_raw_data(self):
        # 1. Determine layers and isotopes (Logic remains solid)
        if self.compoment == "internals":
            layers = materials["Material"].unique().tolist()
            get_iso_list = lambda l: materials[materials["Material"] == l]["Isotope"].unique().tolist()
        else:
            mat_type = "Rock" if self.compoment == "rock" else "Concrete"
            layers = rock[rock["Material"] == mat_type]["Particule"].unique().tolist()
            get_iso_list = lambda l: rock[rock["Particule"] == l]["Isotope"].unique().tolist()

        prefix = "" if self.compoment == "internals" else (f"{mat_type}_")

        for layer in layers:
            isotopes = get_iso_list(layer)
            self.data[layer] = {}
            self.data_counts[layer] = {}
            
            for iso in isotopes:
                # Initialize based on bias
                if self.bias == "bon":
                    all_energies = []
                    all_weights = []
                else:
                    all_energies = []
                
                expected_start = f"{prefix}{layer}_{iso}_"
                pattern = os.path.join(self.folder_path, f"{expected_start}*_{self.bias}_filtered.root")
                
                raw_file_list = glob.glob(pattern)
                file_list = [f for f in raw_file_list if os.path.basename(f).startswith(expected_start)]
                successful_loads = 0
                try: 
                    for f in file_list:
                        d, w = self._get_root_tree(f)
                        if len(d) == 0: continue 

                        d = np.array(d)
                        mask = d > 0
                        d_filtered = d[mask]
                        
                        if len(d_filtered) > 0:
                            all_energies.extend(d_filtered)
                            
                            # Handle weights only if we are in "bon" mode and weights exist
                            if self.bias == "bon" and w is not None:
                                w_filtered = np.array(w)[mask]
                                all_weights.extend(w_filtered)
                            
                            successful_loads += 1
                    
                    # 2. Store data only if we actually found events
                    if all_energies:
                        if self.bias == "bon":
                            self.data[layer][iso] = {
                                'energy': np.array(all_energies), 
                                'weight': np.array(all_weights)
                            }
                        else:
                            self.data[layer][iso] = np.array(all_energies)
                        
                        self.data_counts[layer][iso] = successful_loads
                    else:
                        self.data_counts[layer][iso] = 0
                
                except Exception as e:
                    print(f"Error processing {layer}-{iso}: {e}")
                    self.data_counts[layer][iso] = 0
                    
                # 3. Fixed Event Counting for the print statement
                def get_len(item):
                    if isinstance(item, dict):
                        return len(item['energy'])
                    return len(item)

            total_events = sum(get_len(self.data[layer][i]) for i in self.data[layer] if i in self.data[layer])
            print(f"Loaded {total_events} events across {sum(self.data_counts[layer].values())} files for {layer}")

        print("===========================================")
        print("Data Loading Complete.")
    
    def process_internals_data(self):
        
        """
        Normalizes internal background data using Activity [Bq/kg] and Layer Mass [kg].
        
        Inputs: self.data, materials DataFrame, and ROOT geometry metadata.
        Outputs: 
            - self.raw_norm[layer][iso]: Dictionary with weighted energies.
            - self.summary_data: List of integrated DRU values per layer for the 0-1 keV window.
        
        Note: Includes a check to skip daughters in secular equilibrium with parents.
        """
        
        decay_chains = {
            'U238': ['Ra226', 'Pb214', 'Bi214', 'Pb210', 'Th234', 'Pa234m'],
            'Th232': ['Th228', 'Ac228', 'Pb212', 'Tl208', 'Ra228', 'Ra224', 'Rn220']
        } 
        all_daughters = [item for sublist in decay_chains.values() for item in sublist]
        geo_name_mapping = {'PCu': 'PureCu'}
        win1_min, win1_max = 0, 1 # Narrow window 0-1 keV

        for layer in self.data.keys():
            
            self.raw_norm[layer] = {}
            geo_key = geo_name_mapping.get(layer, layer)
            layer_sum = 0
            
            for iso in self.data[layer].keys():
                # 1. Smart Double Counting Check
                if iso in all_daughters:
                    parent = [p for p, daughts in decay_chains.items() if iso in daughts][0]
                    parent_prefix = f"{layer}_{parent}"
                    if any(f.startswith(parent_prefix) for f in os.listdir(self.folder_path)):
                        iso_row = materials[(materials['Material'] == layer) & (materials['Isotope'] == iso)]
                        parent_row = materials[(materials['Material'] == layer) & (materials['Isotope'] == parent)]
                        
                        if not iso_row.empty and not parent_row.empty:
                            if iso_row.iloc[0]['Activity'] == parent_row.iloc[0]['Activity']:
                                continue

                # 2. Metadata & Scaling
                try:
                    prefix = f"{layer}_{iso}"
                    matching_files = [f for f in os.listdir(self.folder_path) if f.startswith(prefix) and f.endswith(".root")]
                    if not matching_files: continue
                    
                    full_path = os.path.join(self.folder_path, matching_files[0])
                    macro_info = GetNormParam(full_path, "runMacro")
                    geo_info = GetNormParam(full_path, "geometryTable")
                    
                    n_sim = macro_info.beamon_number
                    n_files = self.data_counts[layer][iso]
                    mass_layer = geo_info.total_mass[geo_key] 
                    mass_det = geo_info.total_mass[self.target]
                    
                    row = materials[(materials['Material'] == layer) & (materials['Isotope'] == iso)]
                    if row.empty: continue
                    
                    raw_act = row.iloc[0]['Activity']
                    raw_sigma = row.iloc[0]['Sigma']
                    
                    if raw_act > 0:
                        activity = raw_act
                        act_err = raw_sigma  
                    else:
                        activity = raw_sigma
                        act_err = 0.0
                    # Formula: (Flux * Surface * SecPerDay) / (N_sim * Mass_Detector * Energy_Width)
                    dru_factor = (activity * mass_layer * self.SecPerDay * self.Bin) / (n_sim * n_files * mass_det)
                    
                    energies = self.data[layer][iso]
                    stat_err = np.sqrt(len(energies)) * dru_factor
                    sys_err = dru_factor * (act_err / activity if activity > 0 else 0)
                    
                    dru_win1 = np.sum((energies >= win1_min) & (energies < win1_max)) * dru_factor
                    n_win1 = np.sum((energies >= win1_min) & (energies < win1_max))
                    
                    dru_win1_stat_err = np.sqrt(n_win1) * dru_factor if n_win1 > 0 else 0
                    
                    layer_sum += dru_win1
                    
                    self.raw_norm[layer][iso] = {'energy': energies, 
                                                'stat_err': stat_err,
                                                'sys_err': sys_err,
                                                'weight': dru_factor, 
                                                'dru_0_1': dru_win1,
                                                'dru_0_1_stat_err': dru_win1_stat_err
                                                 }

                except Exception:
                    continue

            self.summary_data.append({'Layer': layer, 'DRU_0-1_keV': layer_sum, 'Stat_err': dru_win1_stat_err, 'Sys_err': sys_err})

    def process_externals_data(self):
        """
        Normalizes external background data.
        Consistent Scaling: (Flux * Area * Time * Bin) / (N_sim * N_files * Mass)
        """
        win1_min, win1_max = 0, 1
        mat_type = "Rock" if self.compoment == "rock" else "Concrete"

        for layer in self.data.keys():
            self.raw_norm[layer] = {}
            layer_sum = 0
            layer_stat_err_sq = 0
            layer_sys_err_sq = 0
            
            for iso in self.data[layer].keys():
                try:
                    # 1. File Metadata & Selection
                    prefix = f"{mat_type}_{layer}_{iso}"
                    suffix = f"_{self.bias}_filtered.root"
                    matching_files = [f for f in os.listdir(self.folder_path) if f.startswith(prefix) and f.endswith(suffix)]
                    if not matching_files: continue
                    
                    full_path = os.path.join(self.folder_path, matching_files[0])
                    macro_info = GetNormParam(full_path, "runMacro")
                    geo_info = GetNormParam(full_path, "geometryTable")
                    
                    n_sim = macro_info.beamon_number
                    n_files = self.data_counts[layer][iso]
                    mass_det = geo_info.total_mass[self.target]

                    # 2. Parameter Lookup
                    sel = rock[(rock["Particule"] == layer) & (rock["Isotope"] == iso)]
                    if sel.empty: continue
                    
                    surface_area = float(sel["Surface"].iat[0])
                    raw_flux, raw_sigma = float(sel["Flux"].iat[0]), float(sel["Sigma"].iat[0])
                    flux, flux_err = (raw_flux, raw_sigma) if raw_flux > 0 else (raw_sigma, 0.0)
                    
                    
                    # 3. Unified Scaling Factor (Applying n_sim * n_files to all)
                    dru_factor = (flux * surface_area * self.SecPerDay * self.Bin) / (n_sim * n_files * mass_det)
                    
                    # 4. Mode-Specific Statistical Logic
                    if self.bias == "bon":
                        energies = self.data[layer][iso]['energy']
                        weights = self.data[layer][iso]['weight']
                        
                        mask = (energies >= win1_min) & (energies < win1_max)

                        # Weighted normalization
                        dru_win1 = np.sum(weights[mask]) * dru_factor
                        # Stat Error: factor * sqrt(sum of weights squared)
                        dru_win1_stat_err = np.sqrt(np.sum(weights[mask]**2)) * dru_factor if any(mask) else 0
                        
                        self.raw_norm[layer][iso] = {'energy': energies, 'weights': weights, 'weight': dru_factor}
                                                
                    if self.bias == "boff":
                        
                        energies = self.data[layer][iso]
                        n_win1 = np.sum((energies >= win1_min) & (energies < win1_max))
                        # Standard normalization
                        dru_win1 = n_win1 * dru_factor
                        # Stat Error: factor * sqrt(N)
                        dru_win1_stat_err = np.sqrt(n_win1) * dru_factor if n_win1 > 0 else 0
                        self.raw_norm[layer][iso] = {'energy': energies, 'weight': dru_factor}

                    # 5. Systematic Error & Dictionary Update
                    dru_win1_sys_err = dru_win1 * (flux_err / flux if flux > 0 else 0)
                    
                    self.raw_norm[layer][iso].update({
                        'dru_0_1': dru_win1,
                        'dru_0_1_stat_err': dru_win1_stat_err,
                        'sys_err': dru_win1_sys_err
                    })

                    # 6. Summing into Layer Totals
                    layer_sum += dru_win1
                    layer_stat_err_sq += dru_win1_stat_err**2
                    layer_sys_err_sq += dru_win1_sys_err**2

                except Exception as e:
                    print(f"Error processing external {layer}-{iso}: {e}")
                    continue
                
            # Append final layer summary to summary_data list
            self.summary_data.append({
                'Layer': layer, 
                'DRU_0-1_keV': layer_sum, 
                'Stat_err': np.sqrt(layer_stat_err_sq), 
                'Sys_err': np.sqrt(layer_sys_err_sq)
            })
            
    def print_summary(self):
        """
        Prints a formatted DataFrame summary of the DRU for the 0-1 keV window.
        
        Inputs: Uses self.summary_data generated during processing.
        Outputs: 
            - console: A formatted table showing DRU per layer.
            - console: The overall global background total in DRU.
        """
        
        if not self.summary_data:
            print("No data processed. Call process_data() first.")
            return

        df_summary = pd.DataFrame(self.summary_data)
        stat_err_sum = np.sqrt(np.sum(df_summary['Stat_err']**2))
        sys_err_sum = np.sqrt(np.sum(df_summary['Sys_err']**2))
        print("\n" + "="*40)
        print(f"{'LAYER SUMMARY (DRU)':^40}")
        print("="*40)
        print(df_summary.to_string(index=False, formatters={'DRU_0-1_keV': '{:,.4e}'.format}))
        print("="*40)
        print(f"OVERALL 0-1 keV TOTAL: {df_summary['DRU_0-1_keV'].sum():.4e} ± {stat_err_sum:.4e}(stat.) ± {sys_err_sum:.4e}(sys.) DRU\n")

    def plot_results(self, show_plot=True):
        """
        Generates normalized spectra with variable binning to smooth high-energy noise.
        Automatically normalizes by bin width to maintain correct DRU units.
        """
        if not self.raw_norm:
            print("No normalized data available to plot.")
            return

        # 1. Create Variable Bins: Fine (1keV) -> Medium (5keV) -> Coarse (20keV)
        # This reduces statistical spikes at high energy where counts are low.
        bins_low = np.arange(self.xLow, 100, 1)
        bins_med = np.arange(100, 500, 5)
        bins_high = np.arange(500, self.xHigh + 21, 20)
        bins = np.concatenate([bins_low, bins_med, bins_high])
        
        # Calculate widths for normalization (DRU = counts / keV)
        bin_widths = np.diff(bins)

        fig1, axes = plt.subplots(len(self.raw_norm), 1, figsize=(12, 4 * len(self.raw_norm)), sharex=True)
        if len(self.raw_norm) == 1: axes = [axes]
        
        fig2, ax_sum = plt.subplots(figsize=(12, 7))
        global_total_hist = np.zeros(len(bins)-1)

        for i, (layer, isotopes) in enumerate(self.raw_norm.items()):
            layer_total_hist = np.zeros(len(bins)-1)
            
            for iso, data in isotopes.items():
                energies = data['energy']
                dru_factor = data['weight']
                
                # --- Bias Handling ---
                if self.bias == "bon":
                    weights = data['weights'] * dru_factor
                else:
                    weights = np.full_like(energies, dru_factor)

                # 2. Histogram with weights
                hist_counts, _ = np.histogram(energies, bins=bins, weights=weights)
                
                # 3. NORMALIZE BY BIN WIDTH
                # Crucial: Since bins have different sizes, we must divide by width 
                # to keep the Y-axis as "Rate per keV".
                weighted_hist = hist_counts / bin_widths
                
                # Plot Isotope Contribution
                axes[i].stairs(weighted_hist, bins, label=iso, baseline=None)
                layer_total_hist += weighted_hist
            
            # Layer Plot Formatting
            axes[i].set_yscale('log')
            axes[i].set_xscale('log')
            axes[i].set_xlim(self.xLow if self.xLow > 0 else 1, self.xHigh)
            axes[i].set_title(f"Isotope Contributions in {layer} ({self.bias.upper()})")
            axes[i].set_ylabel("Rate DRU")
            axes[i].legend(loc='upper right', ncol=3, fontsize='x-small')
            axes[i].grid(True, which='both', ls='--', lw=0.5)

            # Add to Summary Plot
            ax_sum.stairs(layer_total_hist, bins, label=f"{layer}", linewidth=1.5, alpha=0.8)
            global_total_hist += layer_total_hist

        # Global Summary Plot Formatting
        ax_sum.stairs(global_total_hist, bins, label="Total Background", color='black', linewidth=2.5, zorder=10)
        ax_sum.set_yscale('log')
        ax_sum.set_xscale('log')
        ax_sum.grid(True, which='both', ls='--', lw=0.5)
        ax_sum.set_xlim(self.xLow if self.xLow > 0 else 1, self.xHigh)
        ax_sum.set_xlabel("Energy [keV]")
        ax_sum.set_ylabel("Rate DRU")
        ax_sum.legend(loc='lower left', ncol=2)

        if show_plot:
            plt.show()
            
    def save_results(self, filename_prefix="background_results"):
        """
        Goal: Exports normalized spectra and summary tables to persistent files.
        
        Inputs: 
            - filename_prefix: Base name for the output files.
        Outputs: 
            - .csv: A flat table of the layer-wise DRU summaries.
            - .json: The full raw_norm dictionary (energies and weights) for external plotting.
        """
        # 1. Save Summary Table to CSV
        if self.summary_data:
            df = pd.DataFrame(self.summary_data)
            csv_name = f"{filename_prefix}_summary.csv"
            df.to_csv(csv_name, index=False)
            print(f"Summary saved to: {csv_name}")

        # 2. Save Full Normalized Data to JSON
        # Note: We convert NumPy arrays to lists so they are JSON-serializable
        export_dict = {}
        for layer, isotopes in self.raw_norm.items():
            export_dict[layer] = {}
            for iso, data in isotopes.items():
                export_dict[layer][iso] = {
                    'energy': data['energy'].tolist(),
                    'weight': float(data['weight']),
                    'dru_0_1': float(data['dru_0_1'])
                }
        
        json_name = f"{filename_prefix}_data.json"
        with open(json_name, 'w') as f:
            json.dump(export_dict, f)
        print(f"Full data saved to: {json_name}")
        
    def load_results(self, filename_prefix="background_results"):
        """
        Goal: Reconstructs the class state (raw_norm and summary_data) from saved files.
        
        Inputs: 
            - filename_prefix: Base name used when the data was saved.
        Outputs: 
            - self.raw_norm: Re-populated with energy lists and weights.
            - self.summary_data: Re-populated with the layer-wise DRU table.
            
        Note: Converts lists back into NumPy arrays to maintain compatibility with plotting.
        """
        # 1. Load Summary Table from CSV
        csv_name = f"{filename_prefix}_summary.csv"
        if os.path.exists(csv_name):
            df = pd.read_csv(csv_name)
            self.summary_data = df.to_dict(orient='records')
            print(f"Summary loaded from: {csv_name}")
        else:
            print(f"Warning: {csv_name} not found.")

        # 2. Load Full Normalized Data from JSON
        json_name = f"{filename_prefix}_data.json"
        if os.path.exists(json_name):
            with open(json_name, 'r') as f:
                temp_dict = json.load(f)
            
            # Convert lists back to NumPy arrays for plotting
            self.raw_norm = {}
            for layer, isotopes in temp_dict.items():
                self.raw_norm[layer] = {}
                for iso, data in isotopes.items():
                    self.raw_norm[layer][iso] = {
                        'energy': np.array(data['energy']),
                        'weight': data['weight'],
                        'dru_0_1': data['dru_0_1']
                    }
            print(f"Full data loaded from: {json_name}")
        else:
            print(f"Warning: {json_name} not found.")


        
    
    
    