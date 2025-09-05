import uproot
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import os
import glob
import pandas as pd
from tqdm import tqdm
import warnings
from tesssapy.utils import get_cached_data
from tesssapy import get_h5_files as ghd

warnings.simplefilter("ignore")

materials = get_cached_data("materials_data_2.csv")
rock      = get_cached_data("rock_data.csv")

class g4_sim_proc:

    def __init__(self, compoment, folder_path, plots= True):
        self.compoment = compoment
        #self.t = thickness
        #self.geometry = geo         
        self.folder_path = folder_path
        print(f"Processing files in {self.folder_path}")
        #To choose between the rock or the shielding 
         # either "rock" or "internals"
        
        # ======== constants ========
        self.detMass = 0.797336 
        self.SecPerDay = float(3600 * 24)
        #self.xLow = 0.001
        self.xLow = 0
        self.xHigh = 5000
        #self.xHigh = 5000.001
        self.Bins = 20 
        self.BinSize = (self.xHigh - self.xLow) / self.Bins
        self.Bin = float(1.0 / self.BinSize)
        
        # ======== variables ========
        self.data = {}
        self.data_counts = {}
        self.counts = {}
        self.counts_err = {}
        self.energy = {}
        
        # ======== functions calling ========
        self.load_raw_data()
        #self.load_h5py_data()
         
        self.normalize_data()
        self.get_totals()
        if plots == True : 
            self.get_spectrum()
            self.get_spectrum_totals()
        self.print_simulation_summary()
    
    def get_log_bins(self):
        self.bins = np.geomspace(self.xLow, self.xHigh, self.Bins)  
        
    def get_root_tree(self, file_path):
        try:
            root_file = uproot.open(file_path)
            tree = root_file["events"]
            Params = tree.keys()
            
            if not tree:
                #print(f"Warning: No tree found in {file_path}")
                return []
            edep = np.array(tree.arrays(Params[5])[Params[5]])
            return edep
        
        except Exception as e:
            #print(f"Error reading {file_path}: {e}")
            return []
    
    def hist_it(self, X):
        #logbins=np.geomspace(self.xLow, self.xHigh, self.Bins)
        counts, edges = np.histogram(X, self.Bins,[self.xLow,self.xHigh]) #, bins = logbins)
        centers = (edges[:-1] + edges[1:]) / 2.0
        return centers, counts, edges
    
    def load_raw_data(self):

        def get_layers_and_isotopes(component):
            if component == "internals":
                layers = materials["Material"].unique().tolist()
                get_isotopes = lambda l: materials[materials["Material"] == l]["Isotope"].unique().tolist()
            else:
                material_type = "Rock" if component == "rock" else "Concrete"
                layers = rock[rock["Material"] == material_type]["Particule"].unique().tolist()
                get_isotopes = lambda l: rock[rock["Particule"] == l]["Isotope"].unique().tolist()
            return layers, get_isotopes

        def build_filepath(component, layer, iso, i):
            if component == "internals":
                #return f"{self.folder_path}/{layer}_{iso}_{i}_{t}_proc.root"
                return f"{self.folder_path}/{layer}_{iso}_{i}_proc.root"

            prefix = "Rock" if component == "rock" else "Concrete"
            #return f"{self.folder_path}/{prefix}_{layer}_{iso}_{i}_{t}_proc.root"
            return f"{self.folder_path}/{prefix}_{layer}_{iso}_{i}_proc.root"

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
                        #file_path = build_filepath(self.compoment, layer, iso, i, self.t)
                        file_path = build_filepath(self.compoment, layer, iso, i)

                        if not os.path.exists(file_path):
                            continue
                        data = self.get_root_tree(file_path)
                        if data is not None and len(data) > 0:
                            self.data[layer][iso].extend(data)
                            self.data_counts[layer][iso] += 1
                        pbar.update(1)

        print("Data loading complete.", self.data_counts)

    def load_h5py_data(self):
        self.h5_file = ghd.load_h5_file(self.folder_path)
        for layer in self.h5_file.keys():
            self.data[layer] = {}
            self.data_counts[layer] = {}
            for iso in self.h5_file[layer].keys():
                self.data[layer][iso] = [self.h5_file[layer][iso]['edep']]
                self.data_counts[layer][iso] = len(self.h5_file[layer][iso]['edep'])
                
        
    def normalize_data(self):
        
        self.shielding = materials["Material"].unique().tolist()
        self.particle = rock["Particule"].unique().tolist()
        
        def compute_normalized_counts(X, Y, norm, norm_err):
            return np.multiply(Y, norm), np.multiply(Y, norm_err), X

        def get_parameters(component, layer, iso):
            if component == "internals":
                m = "Mass_octa"
                df = materials[(materials["Material"] == layer) & (materials["Isotope"] == iso)]
                return float(df[m]), float(df["Count"]), float(df["Activity"]), float(df["Sigma"])
            else:
                mat = "Rock" if component == "rock" else "Concrete"
                df = rock[(rock["Material"] == mat) & (rock["Particule"] == layer) & (rock["Isotope"] == iso)]
                return float(df["Surface"]), float(df["Count"]), float(df["Flux"]), float(df["Sigma"])

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
                    p1, p2, p3, p4 = get_parameters(self.compoment, layer, iso)

                    normalization_factor = (
                        p1 * self.SecPerDay *
                        (1.0 / (p2 * self.data_counts[layer][iso])) *
                        self.Bin * (1.0 / self.detMass)
                    )
                    norm, norm_err = normalization_factor * p3, normalization_factor * p4
                    #norm, norm_err = 1.0,1.0
                    self.counts[layer][iso], self.counts_err[layer][iso], self.energy[layer][iso] = compute_normalized_counts(X, Y, norm, norm_err)
                
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
            energy_key = 'Cu' if self.compoment == 'internals' else 'Gammas'

            energy = self.energy[energy_key]['K40']

            for k in self.counts:
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
                ax.set_xlim(0,2000)
                ax.legend(loc='upper right', fontsize=20)
                plt.show()
        except Exception as e:
            print(f"Error in plotting spectrum: {e}")
            return


    def get_spectrum_totals(self):
        def plot_total(energy, counts, counts_err, label, color, ax):
            ax.step(energy, counts, where='mid', color=color, label=label)
            ax.fill_between(energy, counts + counts_err, counts - counts_err, step='mid', alpha=0.2, color=color)
            ax.scatter(energy, counts, color=color, marker='.')
            ax.errorbar(energy, counts, yerr=np.sqrt(counts), fmt='.', capsize=2, color=color)
        try:
            energy_key = 'Cu' if self.compoment == 'internals' else 'Gammas'
            energy = self.energy[energy_key]['K40']
            fig, ax = plt.subplots(figsize=(12, 8))

            plot_total(energy, self.counts['total'], self.counts_err['total'], 'total', 'k', ax)

            colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
            for h, k in enumerate(self.counts):
                if k == 'total':
                    continue
                for j in self.counts[k]:
                    if j == 'total':
                        color = colors[h % len(colors)]
                        plot_total(energy, self.counts[k][j], self.counts_err[k][j], k, color, ax)

            ax.set_yscale('log')
            ax.set_xlabel("Energy [keV]", fontsize=20)
            ax.set_ylabel("Counts / (keV.kg.day)", fontsize=20)
            ax.set_title(f" {self.compoment}", fontsize=20)
            ax.grid()
            ax.legend(loc='upper right', fontsize=20)
            plt.show()
        except Exception as e:
            print(f"Error in plotting spectrum totals: {e}")
            return
    
    def get_spectrum_totals2(self):
        def plot_total(energy, counts, counts_err, label, color, ax):
            # Replace zero values in counts with their corresponding error
            counts_clean = np.where(counts == 0, counts_err, counts)
            ax.step(energy, counts_clean, where='mid', color=color, label=label)
            ax.scatter(energy, counts_clean, color=color, marker='.')

        try:
            energy_key = 'Cu' if self.compoment == 'internals' else 'Gammas'
            energy = self.energy[energy_key]['K40']
            fig, ax = plt.subplots(figsize=(12, 8))

            plot_total(energy, self.counts['total'], self.counts_err['total'], 'total', 'k', ax)

            colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
            lab = ['Stainless Steel', 'Copper', 'Cryogenic Copper', 'Lead', 'Polyethylene', 'Titanium', 'Brass']
            for h, k in enumerate(self.counts):
                if k == 'total':
                    continue
                for j in self.counts[k]:
                    if j == 'total':
                        color = colors[h % len(colors)]
                        plot_total(energy, self.counts[k][j], self.counts_err[k][j], lab[h], color, ax)

            ax.set_yscale('log')
            ax.set_xlabel("Energy [keV]", fontsize=20)
            ax.set_ylabel("Counts / (keV·kg·day)", fontsize=20)
            ax.set_title("Internal Background for Octagonal Shielding Geometry", fontsize=20)
            ax.grid()
            
            ax.legend(loc='upper right', fontsize=20)
            plt.show()

        except Exception as e:
            print(f"Error in plotting spectrum totals: {e}")
            return
        
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

class g4_sim_proc_geo:

    def __init__(self, geo, compoment, folder_path, plots= True):
        self.compoment = compoment
        #self.t = thickness
        self.geometry = geo         
        self.folder_path = folder_path
        print(f"Processing files in {self.folder_path}")
        #To choose between the rock or the shielding 
         # either "rock" or "internals"
        
        # ======== constants ========
        self.detMass = 0.797336 
        self.SecPerDay = float(3600 * 24)
        #self.xLow = 0.001
        self.xLow = 0
        self.xHigh = 5000
        #self.xHigh = 5000.001
        self.Bins = 20 
        self.BinSize = (self.xHigh - self.xLow) / self.Bins
        self.Bin = float(1.0 / self.BinSize)
        
        # ======== variables ========
        self.data = {}
        self.data_counts = {}
        self.counts = {}
        self.counts_err = {}
        self.energy = {}
        
        # ======== functions calling ========
        self.load_raw_data()
        #self.load_h5py_data()
         
        self.normalize_data()
        self.get_totals()
        if plots == True : 
            self.get_spectrum()
            self.get_spectrum_totals()
        self.print_simulation_summary()
    
    def get_log_bins(self):
        self.bins = np.geomspace(self.xLow, self.xHigh, self.Bins)  
        
    def get_root_tree(self, file_path):
        try:
            root_file = uproot.open(file_path)
            tree = root_file["events"]
            Params = tree.keys()
            
            if not tree:
                #print(f"Warning: No tree found in {file_path}")
                return []
            edep = np.array(tree.arrays(Params[5])[Params[5]])
            return edep
        
        except Exception as e:
            #print(f"Error reading {file_path}: {e}")
            return []
    
    def hist_it(self, X):
        #logbins=np.geomspace(self.xLow, self.xHigh, self.Bins)
        counts, edges = np.histogram(X, self.Bins,[self.xLow,self.xHigh]) #, bins = logbins)
        centers = (edges[:-1] + edges[1:]) / 2.0
        return centers, counts, edges
    
    def load_raw_data(self):

        def get_layers_and_isotopes(component):
            if component == "internals":
                layers = materials["Material"].unique().tolist()
                get_isotopes = lambda l: materials[materials["Material"] == l]["Isotope"].unique().tolist()
            else:
                material_type = "Rock" if component == "rock" else "Concrete"
                layers = rock[rock["Material"] == material_type]["Particule"].unique().tolist()
                get_isotopes = lambda l: rock[rock["Particule"] == l]["Isotope"].unique().tolist()
            return layers, get_isotopes

        def build_filepath(component, layer, iso, i):
            if component == "internals":
                #return f"{self.folder_path}/{layer}_{iso}_{i}_{t}_proc.root"
                return f"{self.folder_path}/{layer}_{iso}_{i}_proc.root"

            prefix = "Rock" if component == "rock" else "Concrete"
            #return f"{self.folder_path}/{prefix}_{layer}_{iso}_{i}_{t}_proc.root"
            return f"{self.folder_path}/{prefix}_{layer}_{iso}_{i}_proc.root"

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
                        #file_path = build_filepath(self.compoment, layer, iso, i, self.t)
                        file_path = build_filepath(self.compoment, layer, iso, i)

                        if not os.path.exists(file_path):
                            continue
                        data = self.get_root_tree(file_path)
                        if data is not None and len(data) > 0:
                            self.data[layer][iso].extend(data)
                            self.data_counts[layer][iso] += 1
                        pbar.update(1)

        print("Data loading complete.", self.data_counts)

    def load_h5py_data(self):
        self.h5_file = ghd.load_h5_file(self.folder_path)
        for layer in self.h5_file.keys():
            self.data[layer] = {}
            self.data_counts[layer] = {}
            for iso in self.h5_file[layer].keys():
                self.data[layer][iso] = [self.h5_file[layer][iso]['edep']]
                self.data_counts[layer][iso] = len(self.h5_file[layer][iso]['edep'])
                
        
    def normalize_data(self):
        
        self.shielding = materials["Material"].unique().tolist()
        self.particle = rock["Particule"].unique().tolist()
        
        def compute_normalized_counts(X, Y, norm, norm_err):
            return np.multiply(Y, norm), np.multiply(Y, norm_err), X

        def get_parameters(component, layer, iso):
            if component == "internals":
                m = "Mass_" + str(self.geometry)
                df = materials[(materials["Material"] == layer) & (materials["Isotope"] == iso)]
                return float(df[m]), float(df["Count"]), float(df["Activity"]), float(df["Sigma"])
            else:
                mat = "Rock" if component == "rock" else "Concrete"
                df = rock[(rock["Material"] == mat) & (rock["Particule"] == layer) & (rock["Isotope"] == iso)]
                return float(df["Surface"]), float(df["Count"]), float(df["Flux"]), float(df["Sigma"])

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
                    p1, p2, p3, p4 = get_parameters(self.compoment, layer, iso)

                    normalization_factor = (
                        p1 * self.SecPerDay *
                        (1.0 / (p2 * self.data_counts[layer][iso])) *
                        self.Bin * (1.0 / self.detMass)
                    )
                    norm, norm_err = normalization_factor * p3, normalization_factor * p4
                    #norm, norm_err = 1.0,1.0
                    self.counts[layer][iso], self.counts_err[layer][iso], self.energy[layer][iso] = compute_normalized_counts(X, Y, norm, norm_err)
                
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
            energy_key = 'Cu' if self.compoment == 'internals' else 'Gammas'

            energy = self.energy[energy_key]['K40']

            for k in self.counts:
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
                ax.set_title(f"{self.geometry} {self.compoment} {k}", fontsize=20)
                ax.grid()
                ax.set_xlim(0,2000)
                ax.legend(loc='upper right', fontsize=20)
                plt.show()
        except Exception as e:
            print(f"Error in plotting spectrum: {e}")
            return


    def get_spectrum_totals(self):
        def plot_total(energy, counts, counts_err, label, color, ax):
            ax.step(energy, counts, where='mid', color=color, label=label)
            ax.fill_between(energy, counts + counts_err, counts - counts_err, step='mid', alpha=0.2, color=color)
            ax.scatter(energy, counts, color=color, marker='.')
            ax.errorbar(energy, counts, yerr=np.sqrt(counts), fmt='.', capsize=2, color=color)
        try:
            energy_key = 'Cu' if self.compoment == 'internals' else 'Gammas'
            energy = self.energy[energy_key]['K40']
            fig, ax = plt.subplots(figsize=(12, 8))

            plot_total(energy, self.counts['total'], self.counts_err['total'], 'total', 'k', ax)

            colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
            for h, k in enumerate(self.counts):
                if k == 'total':
                    continue
                for j in self.counts[k]:
                    if j == 'total':
                        color = colors[h % len(colors)]
                        plot_total(energy, self.counts[k][j], self.counts_err[k][j], k, color, ax)

            ax.set_yscale('log')
            ax.set_xlabel("Energy [keV]", fontsize=20)
            ax.set_ylabel("Counts / (keV.kg.day)", fontsize=20)
            ax.set_title(f"{self.geometry} {self.compoment}", fontsize=20)
            ax.grid()
            ax.legend(loc='upper right', fontsize=20)
            plt.show()
        except Exception as e:
            print(f"Error in plotting spectrum totals: {e}")
            return
    
    def get_spectrum_totals2(self):
        def plot_total(energy, counts, counts_err, label, color, ax):
            # Replace zero values in counts with their corresponding error
            counts_clean = np.where(counts == 0, counts_err, counts)
            ax.step(energy, counts_clean, where='mid', color=color, label=label)
            ax.scatter(energy, counts_clean, color=color, marker='.')

        try:
            energy_key = 'Cu' if self.compoment == 'internals' else 'Gammas'
            energy = self.energy[energy_key]['K40']
            fig, ax = plt.subplots(figsize=(12, 8))

            plot_total(energy, self.counts['total'], self.counts_err['total'], 'total', 'k', ax)

            colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
            lab = ['Stainless Steel', 'Copper', 'Cryogenic Copper', 'Lead', 'Polyethylene', 'Titanium', 'Brass']
            for h, k in enumerate(self.counts):
                if k == 'total':
                    continue
                for j in self.counts[k]:
                    if j == 'total':
                        color = colors[h % len(colors)]
                        plot_total(energy, self.counts[k][j], self.counts_err[k][j], lab[h], color, ax)

            ax.set_yscale('log')
            ax.set_xlabel("Energy [keV]", fontsize=20)
            ax.set_ylabel("Counts / (keV·kg·day)", fontsize=20)
            ax.set_title("Internal Background for Octagonal Shielding Geometry", fontsize=20)
            ax.grid()
            
            ax.legend(loc='upper right', fontsize=20)
            plt.show()

        except Exception as e:
            print(f"Error in plotting spectrum totals: {e}")
            return
        
    def print_simulation_summary(self):
        print("\nSimulation Summary")
        print("===========================================")
        print(f"Shielding Type: {self.compoment.capitalize()} and Geometry: {self.geometry.capitalize()}\n")
        print("Total counts per Layer of Shielding:")
        print("____________________________________\n")
        for layer in self.counts:
            if layer == 'total':
                continue
            print(f"{layer}: R = {np.sum(self.counts[layer]['total']):.2e} ± {np.sum(self.counts_err[layer]['total']):.2e} counts/keV.day.kg")
        print("\nTotal Counts:")
        print("_______________\n")
        print(f"R = {np.sum(self.counts['total']):.2e} ± {np.sum(self.counts_err['total']):.2e} counts/keV.day.kg")  

        
              
        
              
