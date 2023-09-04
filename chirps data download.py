import os
import numpy as np
import pandas as pd
import xarray as xr
import tarfile
import patoolib
from osgeo import gdal 
import rasterio
import rasterio.mask
import geopandas as gpd
import scipy.sparse as sparse
import fiona
import scipy
import rioxarray as rxr
from shapely.geometry import mapping
import netCDF4

################################### new trial #############################################

import requests 
import numpy as np

####################### downloading chirps data ##########################
os.chdir(r'your directory') # change this to your local directory

years = np.arange(1981,2022)
months = ['01','02','03','04','05','06','07','08','09','10','11','12']

for year in years:
    for mon in months:
        url = 'https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_monthly/tifs/chirps-v2.0.' + str(year)+'.'+ mon + '.tif.gz'
        r = requests.get(url, allow_redirects = True)
        open(str(year)+'_'+mon+'.gz', 'wb').write(r.content)

###################### unzipping it ########################################################

for zipfiles in os.listdir(r'directory where the downloaded files are saved'):  # change this to your local directory

    os.chdir(r'directory where the downloaded files are saved') # change this to your local directory


    if zipfiles[-3: ] == '.gz':
        patoolib.extract_archive(zipfiles, outdir = r'directory where the downloaded files are saved')  # change this to your local directory where the extracted file is to be saved

################## clip the tiff files using the study area shapefile ###############

files = os.listdir(r'directory where the extracted files are saved')

shapefile_ET = gpd.read_file(r'directory where the shapefile of your study area is located')

for f in files:
    #f=files[0]
    if f[-4:]=='.tif':

        #imagery = rasterio.open(f)

        with fiona.open(r'directory where the shapefile of your study area is located','r') as shapefile:
            for feature in shapefile:
                shapes = [feature['geometry']]

        os.chdir(r'directory where the extracted files are saved')

        with rasterio.open(f) as src:
           out_image,out_transform = rasterio.mask.mask(src,shapes,crop=True)
           out_meta = src.meta

        out_meta.update({
            "driver":"Gtiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            'transform':out_transform
        }) 

        os.chdir(r'directory where the tif file of your study area is to be saved')
        
        with rasterio.open(f, "w", **out_meta) as dest:
            dest.write(out_image)

####################### convert each tif files of the study area to netcdf file ###########################        

os.chdir(r'directory where the tif file of your study area is saved')

files = os.listdir(r'directory where the tif file of your study area is saved')

for f in files:
    #f=files[0]
    if f[-4:]=='.tif':

        
        output = f[12:16]+"_"+f[17:19]+'.nc'
        
        ds = gdal.Translate(output, f, format='NetCDF')


################# now lets combine them using xarray  #####################

os.chdir(r'directory where the netcdf file of your study area is saved')

files = os.listdir(r'directory where the netcdf file of your study area is saved')

ds = xr.open_mfdataset(files,combine = 'nested', concat_dim="time") #nested

os.chdir(r'directory where the combined netcdf file of your study area is to be saved')

ds.to_netcdf('study_area2.nc')   #,format='NETCDF3_64BIT'
