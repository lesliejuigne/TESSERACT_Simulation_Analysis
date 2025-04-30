import h5py
import pandas as pd
import numpy as np

def load_h5_file(h5_file_path):
    """
    Reads an HDF5 (.h5) file and loads its datasets into a nested dictionary.

    Args:
        h5_file_path (str): Path to the .h5 file.

    Returns:
        dict: A nested dictionary where keys are layers, isotopes, and values are Pandas DataFrames.
    """
    data_dict = {}

    with h5py.File(h5_file_path, "r") as hdf:
        for layer in hdf.keys():  
            data_dict[layer] = {}
            
            for isotope in hdf[layer].keys(): 
                
                data_dict[layer][isotope] = {}
                
                for key in hdf[layer][isotope].keys():
                    if key =="edep":
                        data_dict[layer][isotope][key] = pd.DataFrame(hdf[layer][isotope][key][:][:,0])
                    else:
                        data_dict[layer][isotope][key] = hdf[layer][isotope][key][:]
                    
    return data_dict
