import uproot as ROOT
from importlib_resources import files
import os
import pandas as pd
import h5py
import glob
from tqdm import tqdm
import numpy as np
import warnings
import tesssa.utils
warnings.simplefilter("ignore")

h5_output = files('sim_data_example')

class RootToH5PY:
    
    def __init__(self , output_file, shield):
        self.output_file = output_file
        self.shield = shield
        self.root_keys = ["file", "ID", "eventID", "clusterIndex", "timeStamp", "edep"]
        self.__message__(True)
        self.get_output()
        self.get_files_h5()
                
    def __message__(self, top):
        if top == True:
            print("Processing in progress! Be patient... :)")
        if top == False:
            print("Processing completed! :)")
    
    def get_output(self):
        self.folder_path = h5_output.joinpath(self.shield)   
        self.iin = os.path.join(self.folder_path,"raw")
        self.out = os.path.join(self.folder_path, self.output_file)

    def get_root_data(self, file_path):
        try:
            self.root_file = ROOT.open(file_path)
            self.tree = self.root_file["events"] 
            if not self.tree:
                print(f"Warning: No tree found in {file_path}")
                return [] 
            
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            
        Params = self.tree.keys()
        data = {
            key: np.array(self.tree.arrays(Params[i])[Params[i]]).tolist()
            for i, key in enumerate(self.root_keys)
        }
        return data
   
    def get_files_h5(self):
        file_list = glob.glob(os.path.join(self.iin, "*_proc.root"))
        self.grouped_data = {}

        try:    
            for file_path in tqdm(file_list):
                file_name = os.path.basename(file_path)
        
                if self.shield == "internals":
                    layer, isotope, _, _ = file_name.replace("._proc.root", "").split("_")
                elif self.shield in ("rock", "concrete"):
                    _, layer, isotope, _ = file_name.replace("_proc.root", "").split("_")

                extracted_data = self.get_root_data(file_path)

                # Initialize the nested dictionary
                if layer not in self.grouped_data:
                    self.grouped_data[layer] = {}
                if isotope not in self.grouped_data[layer]:
                    self.grouped_data[layer][isotope] = {}

                for key, value in extracted_data.items():
                    if key not in self.grouped_data[layer][isotope]:
                        self.grouped_data[layer][isotope][key] = []
                    
                self.grouped_data[layer][isotope][key].extend(value)
                
            
            # Open or create an HDF5 file
            with h5py.File(self.out, "w") as f:
                # Loop through the grouped data to write to the file
                for layer, isotopes in self.grouped_data.items():
                    # Create a group for each layer
                    layer_group = f.create_group(layer)
                    
                    for isotope, data in isotopes.items():
                        # Create a group for each isotope under the corresponding layer
                        isotope_group = layer_group.create_group(isotope)
                        
                        # Write the data (lists) as datasets within each isotope group
                        for key, values in data.items():
                            # Convert lists to numpy arrays for better storage
                            isotope_group.create_dataset(key, data=values)
        except Exception as e:
            print(f"Error processing files: {e}")



