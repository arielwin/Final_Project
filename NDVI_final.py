import os
import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
import glob
import rasterio.mask
import fiona
from shapely.geometry import Point, LineString, Polygon
from rasterio.plot import show
from matplotlib import pyplot
import numpy as np
import rasterio
from rasterio import Affine as A
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask

'''////////////////////////////////////////////////////////////'''
            
def ndviReclass(npArray):
    '''
    Reclassify NDVI array to Low (1), Medium (2) and High (3)
    Accepts nummpy array
    
    '''
    return np.where((npArray > -1) & (npArray <= .2), 1,    
        np.where((npArray > .2) & (npArray <= .5), 2,   
        np.where((npArray > .5) & (npArray <= 1), 3, np.nan)))

'''/////////////////////////////////////////////////////////////'''
            
def group_clip(tifs, boundary_layer, out_list):
    '''
    Function that clips a group of tifs to a boundary layer and outputs a list.
    tifs: list of tifs
    boundary_layer: single boundary layer
    out_list: empty list to be populated by arrays clipped to bounding layer
    '''
    print('The layers are being clipped.')
    for tif in tifs:    
        #print('Tif [', tif, '] is being processed')
        array = rasterio.open(tif)
        crs = array.crs
        
        boundary = gpd.read_file(boundary_layer)
        boundary = boundary.to_crs(crs)
    

        outras, outtrans = mask(array,boundary['geometry'], crop=True)
        outmet = array.meta
            
        out_list.append(np.array(outras))
    print('Group Clip Complete')

'''///////////////////////////////////////////////////////////////////////////////////////////'''
    
    
def getNDVI(nirs, reds, out_list):
    '''
    This function simply takes a list of Near Infrared Arrays and Red Band Arrays and gets the NDVI.
    The output is single list of arrays
    
    nirs: list of near infrared band arrays
    reds: list of red band arrays
    out_list: empty list for appending NDVIs
    
    Note: must have an equal amount of NIR and Red arrays
    
    '''
    print(len(nirs), 'NDVI(s) being calculated')
    e = 0
    while e < len(nirs):
        out_list.append(np.array((nirs[e] - reds[e]) / (nirs[e] + reds[e])))
        del nirs[e], reds[e]
    print('getNDVI is complete.')
    
print('Functions initialized.')

'''///////////////////////////////////////////////////////////////////////////////'''

#import layers
reds = glob.glob(r'Final/GRSM_IMAGERY(3)/*RED.TIF')
reds.sort()
nirs = glob.glob(r'Final/GRSM_IMAGERY(3)/*NIR.TIF')
nirs.sort()

#import boundary
boundary_layer = r'Final/Boundary/GRSM_Boundary.shp'
boundary = gpd.read_file(boundary_layer)
boundary.crs

#empty lists for functions to populate
clipped_reds = []
clipped_nirs = []
ndvi_list = []
reclass_list = []

#running the functions
group_clip(reds, boundary_layer, clipped_reds)
group_clip(nirs, boundary_layer, clipped_nirs)

getNDVI(clipped_nirs, clipped_reds, ndvi_list)

for x in ndvi_list:
    reclass_list.append(ndviReclass(x))

#empty lists for spatial stats
high = []
medium = []
low = []
total = []
mean = []
years = [1999, 2001, 2003, 2005, 2007, 2009, 2011, 2015, 2017, 2019]

for x in reclass_list:
    
    ones = np.where(x==1, 1 , np.nan)
    twos = np.where(x==2, 1, np.nan)
    threes = np.where(x==3, 1, np.nan)
    tot = np.where(x>0, 1, np.nan)
    
    low.append(np.nansum(ones))
    medium.append(np.nansum(twos))
    high.append(np.nansum(threes))
    total.append(np.nansum(tot))
    mean.append(np.nanmean(x))

#create a data frame and export a csv
stats = {'Year' : years, 'High NDVI' : high, 'Medium NDVI': medium, 'Low NDVI': low, 'Mean NDVI': mean, 'Total Cells': total}
df = pd.DataFrame(stats)
df.to_csv('NDVI.csv') 