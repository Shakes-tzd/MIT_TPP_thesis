#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 13:37:15 2021
Plot ASGM emission over South America in highest available resolution, focus on Chacaltaya
@author: arifeinberg
"""
#%%
import xarray as xr
import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import geopandas as gpd
from shapely.geometry import Point

#%%
#os.chdir('/Users/arifeinberg/target2/fs03/d0/arifein/python/')
#%% Load ASGM emission data
fn1 = '../../../d1/angot/postdoc/GMA_2018/emissions_inventory/GMA_emissions_ASGM_Hg0.0.25x0.25.2015.nc'
ds1 = xr.open_dataset(fn1)

#Other inventories
#EDGAR
fn2 = '../emissions/v2020-07/EDGAR/2010/EDGAR_gold_A_2010_Hg.nc'
ds2 = xr.open_dataset(fn2)

#STREETS
fn3 = '../emissions/v2020-07/Streets/Streets2019_Hg.nc'
ds3 = xr.open_dataset(fn3)

#%% load variables
lat_GMA = ds1.lat
lon_GMA = ds1.lon
Hg0_emiss_ASGM_GMA = ds1.emi_hg_0

lat_EDGAR = ds2.lat
lon_EDGAR = ds2.lon
Hg0_emiss_ASGM_EDGAR = ds2.emi_hg_g

lat_STREETS = ds3.lat
lon_STREETS = ds3.lon
Hg0_emiss_STREETS = ds3.Hg0.isel(time=15)
#%% load grid box areas
fn_g1 = '../gbox_areas/GMA_025_025_gboxarea.nc'
fn_g2 = '../gbox_areas/EDGAR_gboxarea.nc'
fn_g3 = '../gbox_areas/Streets_gboxarea.nc'

ds_g1 = xr.open_dataset(fn_g1)
ds_g2 = xr.open_dataset(fn_g2)
ds_g3 = xr.open_dataset(fn_g3)

gbox_GMA = ds_g1.cell_area
gbox_EDGAR = ds_g2.cell_area
gbox_STREETS = ds_g3.cell_area

#%% Convert units to kg yr^-1

s_in_yr = 3.154e7 # seconds in a year
unit_conv = s_in_yr

Hg0_emiss_ASGM_GMA = Hg0_emiss_ASGM_GMA * unit_conv * gbox_GMA
Hg0_emiss_ASGM_EDGAR = Hg0_emiss_ASGM_EDGAR * unit_conv * gbox_EDGAR
Hg0_emiss_STREETS = Hg0_emiss_STREETS * unit_conv * gbox_STREETS

#%% Load shapefile for Peru
sf_nm = "../shapefiles/gadm36_PER_shp/gadm36_PER_0.shp"
peru_sf = gpd.read_file(sf_nm)

#%% Map plot

#2015 AMAP/UNEP
f, axes = plt.subplots(1, 3, figsize=[16,8],subplot_kw=dict(projection=ccrs.PlateCarree()),
                       gridspec_kw=dict(hspace=0.3, wspace=0.3))
axes = axes.flatten()
                       
axes[0].coastlines()
h = axes[0].pcolormesh(lon_GMA, lat_GMA, Hg0_emiss_ASGM_GMA, vmin=0, vmax=1000, rasterized = True)
peru_sf.geometry.boundary.plot(ax=axes[0],Color=None, edgecolor='w',linewidth = 1)
axes[0].set_title('2015 AMAP/UNEP ASGM inventory', fontsize = 16, fontweight='bold'); #title
axes[0].set_xlim([-85, -60])
axes[0].set_ylim([-25, -5])
axes[0].add_feature(cf.BORDERS)

axes[1].coastlines()
axes[1].pcolormesh(lon_EDGAR, lat_EDGAR, Hg0_emiss_ASGM_EDGAR, vmin=0, vmax=1000, rasterized = True)
peru_sf.geometry.boundary.plot(ax=axes[1],Color=None, edgecolor='w',linewidth = 1)
axes[1].set_title('2010 EDGAR ASGM inventory', fontsize = 16, fontweight='bold'); #title
axes[1].set_xlim([-85, -60])
axes[1].set_ylim([-25, -5])
axes[1].add_feature(cf.BORDERS)


axes[2].coastlines()
axes[2].pcolormesh(lon_STREETS, lat_STREETS, Hg0_emiss_STREETS, vmin=0, vmax=1000, rasterized = True)
peru_sf.geometry.boundary.plot(ax=axes[2],Color=None, edgecolor='w',linewidth = 1)
axes[2].set_title('2015 STREETS inventory (all emiss)', fontsize = 16, fontweight='bold'); #title
axes[2].set_xlim([-85, -60])
axes[2].set_ylim([-25, -5])
axes[2].add_feature(cf.BORDERS)

#f.tight_layout()
cbar = f.colorbar(h, extend='max',orientation='horizontal',ax=axes.ravel().tolist(), shrink=0.7)
cbar.set_label('Hg$^0$ emissions (kg yr$^{-1}$)', fontsize = 14)

plt.savefig('Figures/temp.pdf',bbox_inches = 'tight')
#%% Calculate total Hg mass within the Madre de Dios province
#Create boolean array for lat-lon mesh
inPeru_GMA = np.zeros((lat_GMA.size, lon_GMA.size), dtype=bool)
inPeru_EDGAR = np.zeros((lat_EDGAR.size, lon_EDGAR.size), dtype=bool)
inPeru_STREETS = np.zeros((lat_STREETS.size, lon_STREETS.size), dtype=bool)

#only run over area where have points
min_lon = peru_sf.bounds.minx.values
max_lon = peru_sf.bounds.maxx.values
min_lat = peru_sf.bounds.miny.values
max_lat = peru_sf.bounds.maxy.values

#%% calculate for GMA
for ilon, ln in enumerate(lon_GMA.values):
    if ln > 180. :
        ln_adj = ln - 360.
    else :
        ln_adj = ln

    if ln_adj < min_lon or ln_adj > max_lon: #skip if far away from Peru
        continue
    
    for ilat, lt in enumerate(lat_GMA.values ):
        if lt < min_lat or lt > max_lat: #skip if far away from Peru
            continue
        
        inPeru_GMA[ilat, ilon] = peru_sf.geometry.contains(Point(ln_adj,lt)).values
#%% calculate for EDGAR
for ilon, ln in enumerate(lon_EDGAR.values):
    if ln > 180. :
        ln_adj = ln - 360.
    else :
        ln_adj = ln

    if ln_adj < min_lon or ln_adj > max_lon: #skip if far away from Peru
        continue
    
    for ilat, lt in enumerate(lat_EDGAR.values ):
        if lt < min_lat or lt > max_lat: #skip if far away from Peru
            continue
        
        inPeru_EDGAR[ilat, ilon] = peru_sf.geometry.contains(Point(ln_adj,lt)).values
#%% calculate for STREETS
for ilon, ln in enumerate(lon_STREETS.values):
    if ln > 180. :
        ln_adj = ln - 360.
    else :
        ln_adj = ln

    if ln_adj < min_lon or ln_adj > max_lon: #skip if far away from Peru
        continue
    
    for ilat, lt in enumerate(lat_STREETS.values ):
        if lt < min_lat or lt > max_lat: #skip if far away from Peru
            continue
        
        inPeru_STREETS[ilat, ilon] = peru_sf.geometry.contains(Point(ln_adj,lt)).values
#%% check if shape plotting worked
f, axes = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()),
                       gridspec_kw=dict(hspace=0.3, wspace=0.3))
                       
axes.coastlines()
h = axes.pcolormesh(lon_STREETS, lat_STREETS, inPeru_STREETS, vmin=0, vmax=1, rasterized = True)
peru_sf.geometry.boundary.plot(ax=axes,Color=None, edgecolor='w',linewidth = 1)
axes.set_xlim([-85, -60])
axes.set_ylim([-25, -5])
axes.add_feature(cf.BORDERS)
#%% Calculate sum in boolean area
sum_Peru_ASGM_GMA = Hg0_emiss_ASGM_GMA.where(inPeru_GMA).sum().values
sum_Peru_ASGM_EDGAR = Hg0_emiss_ASGM_EDGAR.where(inPeru_EDGAR).sum().values
sum_Peru_STREETS = Hg0_emiss_STREETS.where(inPeru_STREETS).sum().values
print("Peru GMA emissions from ASGM:" + str(sum_Peru_ASGM_GMA) + " kg yr^-1") 
print("Peru EDGAR emissions from ASGM:" + str(sum_Peru_ASGM_EDGAR) + " kg yr^-1") 
print("Peru STREETS emissions (total):" + str(sum_Peru_STREETS) + " kg yr^-1") 


