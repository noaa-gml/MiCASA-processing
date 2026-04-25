#!/usr/bin/env python3

import numpy as np
#import pandas as pd
#from numba import njit
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from datetime import datetime

from glob import glob

import xarray
# import pint
import xarray as xr
# import pint_xarray
# import cf_xarray as cfxr # not side-effect free
# import cf_xarray.units
# xr.set_options(keep_attrs=True)
import xesmf as xe
# import xcdat # SO many monkey patches


print('finished loading')

earth_radius = 6371.008
files_2024 = sorted(glob('ERA5/MiCASA_v1.nee.2024????.nc'))

# seem to overlap on end of day, so slice off last piece
era5_DS = xr.open_mfdataset(files_2024, preprocess=lambda ds: ds.isel(time=slice(0, -1)))

print('read ds')

grid_1x1 = xe.util.grid_global(1, 1, cf=True)
grid_1x1_areas = xe.util.cell_area(grid_1x1, earth_radius=earth_radius)
era5_DS['1x1 areas'] = xr.DataArray(grid_1x1_areas.to_numpy(), dims=['lat','lon']) #.pint.quantify('km^2')

print('made grid')

# extensive flux
era5_DS['extensive NEE'] = era5_DS['NEE'] * era5_DS['1x1 areas']
print(era5_DS['extensive NEE'])

# plot monthly averages
month_grouped_data = era5_DS.groupby("time.month").mean()

for month in range(0,11):
    plt.clf()
    month_grouped_data['NEE'][month].plot()
#    print((month_grouped_data['extensive NEE'][month]).sum().values)
    plt.savefig(f'figures/month_avg_2024_{month}.png')

