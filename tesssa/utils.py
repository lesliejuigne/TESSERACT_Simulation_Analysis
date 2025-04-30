from importlib_resources import files
import pandas as pd

#all resources are located here: 
data_folder = files('cached_data')

def get_cached_data(file_name):
    return pd.read_csv(data_folder.joinpath(file_name) ) 

def test_util():
    print("I'm even more helpful!")
