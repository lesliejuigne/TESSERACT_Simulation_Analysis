from importlib_resources import files
import pandas as pd
import os

#all resources are located here: 
cached_data_folder = files('cached_data')


def get_cached_data(file_name):
    return pd.read_csv(cached_data_folder.joinpath(file_name) ) 

def load_style_file(style_file="SetStyle_mplstyle.txt"):
    return str(cached_data_folder.joinpath(style_file))   # matplotlib needs a file path string

def test_util():
    print("I'm even more helpful!")
