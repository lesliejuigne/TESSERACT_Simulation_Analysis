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

#materials = get_cached_data("internals_radiogenic_data_v0.csv")
rock      = get_cached_data("rock_data.csv")

plt.style.use(load_style_file('SetStyle_mplstyle.txt'))

class G4SimProcessor:
    def __init__(self, compoment, folder_path, bias="boff", target='LHe', materials=None):
        # ======== parameters ========
        self.compoment = compoment.lower()      
        self.folder_path = folder_path  
        self.bias = bias
        self.target = target
        self.materials = materials
        
        # ======== constants ========
        self.detMass = 0.199866 
        self.SecPerDay = float(3600 * 24)
        self.xLow = 0
        self.xHigh = 4000
        self.Bins = 4000 
        self.BinSize = (self.xHigh - self.xLow) / self.Bins
        self.Bin = float(1.0 / self.BinSize)
        
        # ======== variables ========
        # self.data[layer][iso] is ALWAYS {'energy': array, 'weight': array}.
        # With bias off the weights are all 1.0, so every downstream method can
        # use the same weighted formulas for "bon" and "boff".
        self.data = {}
        self.data_counts = {} # number of files that contributed events
        self.counts = {}
        self.counts_err = {}
        self.energy = {}
        self.raw_norm = {}
        self.summary_data = []
        self.mixed_weight_events = 0  # events whose hits carry different trackWeights
        
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
        """
        Reads one ROOT file and returns one entry per event.

        Output: Tuple of (summed energy deposit per event, weight per event).
        With bias off the weights are 1.0. With bias on they are the mean
        trackWeight of the event's hits; events whose hits carry different
        weights are counted in self.mixed_weight_events.
        """
        try:
            with uproot.open(file_path) as root_file:
                tree = root_file["events"]
                keys = tree.keys()
                
                # Identify the energy deposition column
                edep_key = [k for k in keys if "edep" in k.lower() or "energy" in k.lower()][0]

                if self.bias == "boff":
                    data = tree.arrays([edep_key, "eventID"], library="pd")
                    event_data = data.groupby("eventID")[edep_key].sum().values
                    weights = np.ones(len(event_data), dtype=float)
                    return event_data, weights

                # Bias on: fetch edep, eventID AND trackWeight
                data = tree.arrays([edep_key, "eventID", "trackWeight"], library="pd")
                grouped = data.groupby("eventID")
                event_data = grouped[edep_key].sum().values

                track_w = grouped["trackWeight"]
                weights = track_w.mean().values.astype(float)

                # The mean is only exact if the weight is constant within an event
                spread = (track_w.max() - track_w.min()).values
                self.mixed_weight_events += int(np.sum(spread > 1e-6 * np.abs(weights)))

                return event_data, weights
            
        except Exception as e:
            print(f"Erreur lors de la lecture du fichier {file_path}: {e}")
        return np.array([]), np.array([])
    
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

    @staticmethod
    def _window_rate(energies, weights, dru_factor, e_min, e_max):
        """
        Weighted rate and statistical error in the window [e_min, e_max).

        Rate     = dru_factor * sum(w)
        Stat err = dru_factor * sqrt(sum(w^2))   (reduces to sqrt(N) when w = 1)
        """
        in_win = (energies >= e_min) & (energies < e_max)
        w = weights[in_win]
        rate = np.sum(w) * dru_factor
        stat_err = np.sqrt(np.sum(w**2)) * dru_factor
        return rate, stat_err
    
    def load_raw_data(self):
        # 1. Determine layers and isotopes
        if self.compoment == "internals":
            layers = self.materials["Material"].unique().tolist()
            get_iso_list = lambda l: self.materials[self.materials["Material"] == l]["Isotope"].unique().tolist()
            
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
                all_energies = []
                all_weights = []
                
                expected_start = f"{prefix}{layer}_{iso}_"
                pattern = os.path.join(self.folder_path, f"{expected_start}*_{self.bias}_filtered.root")
                raw_file_list = glob.glob(pattern)
                file_list = [f for f in raw_file_list if os.path.basename(f).startswith(expected_start)]
                successful_loads = 0
                try: 
                    for f in file_list:
                        d, w = self._get_root_tree(f)
                        if len(d) == 0: continue 

                        d = np.asarray(d)
                        w = np.asarray(w, dtype=float)
                        mask = d > 0.01
                        
                        if np.any(mask):
                            # Energies and weights are filtered with the same mask
                            # so they always stay aligned event by event
                            all_energies.append(d[mask])
                            all_weights.append(w[mask])
                            successful_loads += 1
                    
                    # 2. Store data only if we actually found events
                    if all_energies:
                        self.data[layer][iso] = {
                            'energy': np.concatenate(all_energies), 
                            'weight': np.concatenate(all_weights)
                        }
                        self.data_counts[layer][iso] = successful_loads
                    else:
                        self.data_counts[layer][iso] = 0
                
                except Exception as e:
                    print(f"Error processing {layer}-{iso}: {e}")
                    self.data[layer].pop(iso, None)
                    self.data_counts[layer][iso] = 0

            total_events = sum(len(entry['energy']) for entry in self.data[layer].values())
            print(f"Loaded {total_events} events across {sum(self.data_counts[layer].values())} files for {layer}")

        if self.mixed_weight_events > 0:
            print(f"WARNING: {self.mixed_weight_events} events have hits with different trackWeights; "
                  f"the event weight used is the mean over the hits.")

        print("===========================================")
        print("Data Loading Complete.")
    
    def process_internals_data(self):
        """
        Normalizes internal background data using Activity [Bq/kg] and Layer Mass [kg].
        
        Inputs: self.data, materials DataFrame, and ROOT geometry metadata.
        Outputs: 
            - self.raw_norm[layer][iso]: Dictionary with energies, per-event bias
              weights and the DRU scaling factor.
            - self.summary_data: List of integrated DRU values per layer for the 0-1 keV window.
        
        Notes: 
            - Bias weights are applied to the rate and to the statistical error
              (they are all 1.0 when bias is off).
            - Includes a check to skip daughters in secular equilibrium with parents.
        """
        print("Processing internal background simulations...")
        
        decay_chains = {
            'U238': ['Ra226', 'Pb214', 'Bi214', 'Pb210', 'Th234', 'Pa234m'],
            'Th232': ['Th228', 'Ac228', 'Pb212', 'Tl208', 'Ra228', 'Ra224', 'Rn220']
        } 
        
        nsims_u238early = {
            'Cu': 10000000,
            'PCu': 10000000,
            'MuMetal':10000000,
            'Pb': 10000000,
            'PE': 10000000,
            'PEext': 100000000,
            'SSi': 100000000,
            'Ti': 50000000
        }
        
        all_daughters = [item for sublist in decay_chains.values() for item in sublist]
        geo_name_mapping = {'PCu': 'PureCu'}
        win1_min, win1_max = 0, 1 # Narrow window 0-1 keV
        suffix = f"_{self.bias}_filtered.root"
        folder_files = os.listdir(self.folder_path)

        for layer in self.data.keys():
            self.raw_norm[layer] = {}
            geo_key = geo_name_mapping.get(layer, layer)
            layer_sum = 0
            layer_stat_err_sq = 0
            layer_sys_err_sq = 0

            for iso in self.data[layer].keys():
                try:
                    # 1. Smart Double Counting Check
                    if iso in all_daughters:
                        parent = [p for p, daughts in decay_chains.items() if iso in daughts][0]
                        parent_prefix = f"{layer}_{parent}"
                        if any(f.startswith(parent_prefix) for f in folder_files):
                            iso_row = self.materials[(self.materials['Material'] == layer) & (self.materials['Isotope'] == iso)]
                            parent_row = self.materials[(self.materials['Material'] == layer) & (self.materials['Isotope'] == parent)]
                            
                            if not iso_row.empty and not parent_row.empty:
                                if iso_row.iloc[0]['Activity'] == parent_row.iloc[0]['Activity']:
                                    print(f"Skipping {layer}-{iso} due to secular equilibrium with {parent}.")
                                    continue

                    # 2. Metadata & Scaling
                    # Only files of the current bias mode are used, so n_sim always
                    # comes from the same kind of run as the loaded events
                    prefix = f"{layer}_{iso}_"
                    matching_files = [f for f in folder_files if f.startswith(prefix) and f.endswith(suffix)]
                    if not matching_files: 
                        continue
                    
                    full_path = os.path.join(self.folder_path, matching_files[0])
                    macro_info = GetNormParam(full_path, "runMacro")
                    geo_info = GetNormParam(full_path, "geometryTable")
                    
                    n_sim = macro_info.beamon_number
                    if iso in ("U238early", "U238late"):
                        n_sim = nsims_u238early.get(layer, n_sim)  # Override with predefined values if available
                    
                    n_files = self.data_counts[layer][iso]
                    mass_layer = geo_info.total_mass[geo_key] 
                    mass_det = geo_info.total_mass[self.target]
                    
                    row = self.materials[(self.materials['Material'] == layer) & (self.materials['Isotope'] == iso)]
                    if row.empty: continue
                    
                    raw_act = row.iloc[0]['Activity']
                    raw_sigma = row.iloc[0]['Sigma']
                    
                    if raw_act > 0:
                        activity = raw_act
                        act_err = raw_sigma  
                    else:
                        # Upper limit: use sigma as the activity, no systematic error
                        activity = raw_sigma
                        act_err = 0.0

                    # Formula: (Activity * Mass_Layer * SecPerDay) / (N_sim * N_files * Mass_Detector * Energy_Width)
                    dru_factor = (activity * mass_layer * self.SecPerDay * self.Bin) / (n_sim * n_files * mass_det)
                    
                    # 3. Weighted rate (bias weights; all 1.0 when bias is off)
                    energies = self.data[layer][iso]['energy']
                    weights = self.data[layer][iso]['weight']

                    dru_win1, dru_win1_stat_err = self._window_rate(energies, weights, dru_factor, win1_min, win1_max)

                    # Full-spectrum statistical error
                    stat_err = np.sqrt(np.sum(weights**2)) * dru_factor

                    # 4. Systematic error from the activity uncertainty
                    rel_sys = act_err / activity if activity > 0 else 0.0
                    dru_win1_sys_err = dru_win1 * rel_sys
                    
                    self.raw_norm[layer][iso] = {'energy': energies, 
                                                 'weights': weights,
                                                 'weight': dru_factor, 
                                                 'stat_err': stat_err,
                                                 'sys_err': dru_win1_sys_err,
                                                 'dru_0_1': dru_win1,
                                                 'dru_0_1_stat_err': dru_win1_stat_err
                                                 }

                    # 5. Summing into Layer Totals
                    layer_sum += dru_win1
                    layer_stat_err_sq += dru_win1_stat_err**2
                    layer_sys_err_sq += dru_win1_sys_err**2

                except Exception as e:
                    print(f"Error processing {layer}-{iso}: {e}")
                    continue

            self.summary_data.append({
                'Layer': layer, 
                'DRU_0-1_keV': layer_sum, 
                'Stat_err': np.sqrt(layer_stat_err_sq), 
                'Sys_err': np.sqrt(layer_sys_err_sq)
            })

    def process_externals_data(self):
        """
        Normalizes external background data.
        Consistent Scaling: (Flux * Area * Time * Bin) / (N_sim * N_files * Mass)

        Bias weights are applied to the rate and to the statistical error
        (they are all 1.0 when bias is off).
        """
        win1_min, win1_max = 0, 1
        mat_type = "Rock" if self.compoment == "rock" else "Concrete"
        suffix = f"_{self.bias}_filtered.root"
        folder_files = os.listdir(self.folder_path)

        for layer in self.data.keys():
            self.raw_norm[layer] = {}
            layer_sum = 0
            layer_stat_err_sq = 0
            layer_sys_err_sq = 0
            
            for iso in self.data[layer].keys():
                try:
                    # 1. File Metadata & Selection
                    prefix = f"{mat_type}_{layer}_{iso}"
                    matching_files = [f for f in folder_files if f.startswith(prefix) and f.endswith(suffix)]
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
                    
                    print("===== DEBUG INFO =====")
                    print(f"Layer: {layer}, Isotope: {iso}")
                    print(f"Flux: {flux}, Surface Area: {surface_area}, SecPerDay: {self.SecPerDay}, Bin: {self.Bin}")
                    print(f"N_Sim: {n_sim}, N_Files: {n_files}, Mass_Det: {mass_det}")
                    print(f"Computed DRU Factor: {dru_factor}")

                    # 4. Weighted rate (bias weights; all 1.0 when bias is off)
                    energies = self.data[layer][iso]['energy']
                    weights = self.data[layer][iso]['weight']

                    dru_win1, dru_win1_stat_err = self._window_rate(energies, weights, dru_factor, win1_min, win1_max)

                    # 5. Systematic Error & Dictionary Update
                    dru_win1_sys_err = dru_win1 * (flux_err / flux if flux > 0 else 0)
                    
                    self.raw_norm[layer][iso] = {
                        'energy': energies, 
                        'weights': weights, 
                        'weight': dru_factor,
                        'stat_err': np.sqrt(np.sum(weights**2)) * dru_factor,
                        'sys_err': dru_win1_sys_err,
                        'dru_0_1': dru_win1,
                        'dru_0_1_stat_err': dru_win1_stat_err
                    }

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

    def plot_results(self, show_plot=True, save_plot=False, file_name = None):
        """
        Generates normalized spectra with variable binning to smooth high-energy noise.
        Automatically normalizes by bin width to maintain correct DRU units.
        Each event enters with (bias weight * DRU factor).
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
                # Per-event bias weights (all 1.0 when bias is off, or when the
                # data comes from an old results file saved without them)
                event_weights = data.get('weights')
                if event_weights is None:
                    event_weights = np.ones(len(energies), dtype=float)
                weights = event_weights * dru_factor

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
        #ax_sum.grid(True, which='both', ls='--', lw=0.5)
        ax_sum.set_xlim(self.xLow if self.xLow > 0 else 1, self.xHigh)
        ax_sum.set_xlabel("Energy [keV]")
        ax_sum.set_ylabel("Rate DRU")
        ax_sum.legend(loc='upper left', ncol=2)

        if show_plot:
            plt.show()
        if save_plot:
            if file_name is None:
                file_name = f"{self.compoment}_spectra_{self.bias}"
            fig1.savefig(f"{file_name}_spectra.png", dpi=300)
            fig2.savefig(f"{file_name}_total_spectrum.png", dpi=300)
            print(f"Plots saved as '{file_name}_spectra.png' and '{file_name}_total_spectrum.png'.")
            
    def save_results(self, filename_prefix="background_results"):
        """
        Goal: Exports normalized spectra and summary tables to persistent files.
        
        Inputs: 
            - filename_prefix: Base name for the output files.
        Outputs: 
            - .csv: A flat table of the layer-wise DRU summaries.
            - .json: The full raw_norm dictionary (energies, per-event bias weights
              and DRU factor) for external plotting.
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
                event_weights = data.get('weights')
                if event_weights is None:
                    event_weights = np.ones(len(data['energy']), dtype=float)
                export_dict[layer][iso] = {
                    'energy': np.asarray(data['energy']).tolist(),
                    'weights': np.asarray(event_weights).tolist(),
                    'weight': float(data['weight']),
                    'dru_0_1': float(data['dru_0_1']),
                    'dru_0_1_stat_err': float(data.get('dru_0_1_stat_err', 0.0)),
                    'stat_err': float(data.get('stat_err', 0.0)),
                    'sys_err': float(data.get('sys_err', 0.0))
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
            - self.raw_norm: Re-populated with energies, bias weights and DRU factors.
            - self.summary_data: Re-populated with the layer-wise DRU table.
            
        Note: Converts lists back into NumPy arrays to maintain compatibility with plotting.
        Files saved before the bias weights were exported are still readable; their
        events get a weight of 1.0.
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
                    energies = np.array(data['energy'], dtype=float)
                    if 'weights' in data:
                        event_weights = np.array(data['weights'], dtype=float)
                    else:
                        event_weights = np.ones(len(energies), dtype=float)
                    self.raw_norm[layer][iso] = {
                        'energy': energies,
                        'weights': event_weights,
                        'weight': data['weight'],
                        'dru_0_1': data['dru_0_1'],
                        'dru_0_1_stat_err': data.get('dru_0_1_stat_err', 0.0),
                        'stat_err': data.get('stat_err', 0.0),
                        'sys_err': data.get('sys_err', 0.0)
                    }
            print(f"Full data loaded from: {json_name}")
        else:
            print(f"Warning: {json_name} not found.")