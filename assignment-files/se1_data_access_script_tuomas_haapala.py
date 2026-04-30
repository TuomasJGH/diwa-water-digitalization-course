
#reproducibility: installers for libraries?


import os
import requests
from functools import partial
import time
from concurrent.futures import ThreadPoolExecutor


"""Submission for Standard Element 1"""
"""Author: Tuomas Haapala, Aalto University"""
"""No authentication protocol is required - the data is openly available."""


"Variables"
GRIDDED_DATA_TO_DOWNLOAD = ["RRday/rrday_"]
#["RRday/rrday_"] - Precipitation.
#["Tday/tday_"] - Mean daily temperature.
#["Globrad/globrad_"] - Global radiation sum.
#["Rh/rh_"] - Relative humidity.
#["ET0_FAO/ET0_FAO_YYYY_months_4_to_9"] - Potential evapotranspiration for April-September. 1981-2023 available
GRIDDED_DATA_URL = "https://fmi-gridded-obs-daily-1km.s3-eu-west-1.amazonaws.com/Netcdf/"
CHUNK_SIZE = 8192

"Functions"

def construct_gridded_urls(start_year, end_year):
    urls = []
    if end_year is None:
        end_year = start_year
    for var in GRIDDED_DATA_TO_DOWNLOAD:
        for year in range(start_year, end_year+1):
            urls.append([GRIDDED_DATA_URL + var + str(year) + ".nc"]) #"_months_4_to.nc" for PET
    return urls

def download_file(url, save_path=None):
    try:
        local_filename = url.split('/')[-1]
        if save_path is not None:
            local_filename = fix_path(save_path) + local_filename

        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                    f.write(chunk)
    except Exception as e:
        print(e)

def fix_path(path):
    if path[-1] != "/":
        path = path + "/"
    return path

def download_gridded_data(start_year, end_year):
    print("Downloading gridded data from FMI...")
    urls = construct_gridded_urls(start_year, end_year)
    download_func = partial(download_file, save_path=data_folder)
    start_downloading_time = time.time()
    with ThreadPoolExecutor() as executor:
        for i in urls:
            executor.map(download_func, i)
    end_downloading_time = time.time()
    print("Gridded data downloaded in {} s".format(end_downloading_time-start_downloading_time))

def main():
    global start_year, end_year, data_folder
    
    """raw_data is created to the same folder where the script is saved."""
    file_location = os.getcwd() + "/raw_data" 
    data_folder = fr"{file_location}" 
    try:
        os.makedirs(data_folder)
    except OSError as iex:
        print(f"Creation of directories failed! {iex}")
    
    "Program introduction and year inputs"
    print("This is the precipitation data downloader script.")
    start_year = int(input("Please input the start year of the dataset: "))
    end_year = int(input("Please input the end year of the dataset: "))
    if start_year == end_year:
        print("One year of data to be downloaded.")
    else:
        print(end_year - start_year+1, "years of data to be downloaded.")
    input("Press Enter to begin download:")
    download_gridded_data(start_year, end_year)
main()    
    
