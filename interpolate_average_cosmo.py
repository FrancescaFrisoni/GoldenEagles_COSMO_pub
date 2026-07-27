#%%
#---------------------------------------------------------
# Import Modules
#---------------------------------------------------------
import sys
sys.path.append('/home/tcarrard/eagle.gravity.waves/utilities/')


import pandas as pd
import numpy as np

import xarray as xr
import xesmf as xe

from time import time

from rotate_grid import rotate_points

t0 = time()


#--------------------------------------------------------
# Settings
#--------------------------------------------------------

#define single levels interpolation data
interp_1l = ['U_10M','V_10M','ASHFL_S','TKE','HSURF','SLO_ASP_x', 'SLO_ASP_y', 'SLO_ANG_x', 'SLO_ANG_y']
interp_2l = ['U','V','W','TKE']



#%%
#--------------------------------------------------------
# Read data
#--------------------------------------------------------
cosmo_path= '/net/thermo/atmosdyn/tcarrard/data/eagles/'
eagle_path= '/home/tcarrard/eagle.gravity.waves/data/'
constpath   = '/net/litho/atmosdyn2/jansingl/KENDA-1/ANA23/'
constfile   = 'laf2023010100-const.nc'


#read constant file
const_ds = xr.open_dataset(constpath+constfile)
z = const_ds['HEIGHT'].values
pole_lon = const_ds['grid_mapping_1'].grid_north_pole_longitude
pole_lat = const_ds['grid_mapping_1'].grid_north_pole_latitude

#define which constant variable need to be interpolated
interp_vars = const_ds[['HSURF','SLO_ASP','SLO_ANG']]
interp_vars[['SLO_ASP','SLO_ANG']] = interp_vars[['SLO_ASP','SLO_ANG']]* 180/np.pi

cosmo_mean = xr.open_dataset(cosmo_path+'mean_cosmo_summer2023.nc')
eagle_points = pd.read_csv(eagle_path+'finaldataset_backgandpresence_timeSteps_forAnnotation.csv')

#Remove 'history' attribute
if "history" in cosmo_mean.attrs:
    del cosmo_mean.attrs["history"]
    
#merge constants to cosmo dataset
cosmo_mean.update(interp_vars)

#%%
#--------------------------------------------------------
# Subset COSMO
#--------------------------------------------------------

#subset COSMO in domain
min_lat = np.floor(eagle_points['y'].min())
max_lat =  np.ceil(eagle_points['y'].max())
min_lon =  np.floor(eagle_points['x'].min())
max_lon =  np.ceil(eagle_points['x'].max())


#rotate coordinates
lim_rlon,lim_rlat = rotate_points(pole_lon,pole_lat,[min_lon,max_lon],[min_lat,max_lat], 'n2r','NH')

#subset cosmo

cosmo_sub = cosmo_mean.sel({'x_1':cosmo_mean['x_1'].where((cosmo_mean['x_1']>=lim_rlon[0]) & (cosmo_mean['x_1'] <= lim_rlon[1]), drop=True)})
cosmo_sub1= cosmo_sub.sel({'y_1':cosmo_sub['y_1'].where((cosmo_sub['y_1']>=lim_rlat[0]) & (cosmo_sub['y_1'] <= lim_rlat[1]), drop=True)})

cosmo_destaggered = cosmo_sub1.rolling({'z_2': 2}, center=False).mean().isel({'z_2': slice(1, None)})

#--------------------------------------------------------
# Convert to single levels
#--------------------------------------------------------

#find which model level is closest to 200m
hbs = z[:,:,:]-cosmo_mean['HSURF'].values
hbs_lev = (hbs[:-1,:,:] + hbs[1:,:,:])/2
hbs_lev = np.mean(hbs_lev,axis=(1,2))

diff_height = abs(np.array([x-200 for x in hbs_lev]))

lev_200 = np.unique(np.where(diff_height == diff_height.min())[0])

#subset cosmo levels between first level and lev 200
cosmo_sub2 = cosmo_destaggered.sel({'z_1':cosmo_destaggered['z_1'].where((cosmo_destaggered['z_1']>=lev_200), drop=True),
                                    'z_2':cosmo_destaggered['z_2'].where((cosmo_destaggered['z_2']>=lev_200), drop=True)})

#average cosmo over model levels
cosmo_sub2[interp_2l] = cosmo_sub2[interp_2l].mean(dim=['z_1','z_2'])

#squeeze dataset 
cosmo_squeezed = cosmo_sub2.squeeze()

#--------------------------------------------------------
# Witch slope angle to cartesian coordinates
#--------------------------------------------------------

#convert to cartesian coordinates
cosmo_squeezed['SLO_ASP_x'] = np.cos(np.deg2rad(90-cosmo_squeezed['SLO_ASP']))
cosmo_squeezed['SLO_ASP_y'] = np.sin(np.deg2rad(90-cosmo_squeezed['SLO_ASP']))
cosmo_squeezed['SLO_ANG_x'] = np.cos(np.deg2rad(cosmo_squeezed['SLO_ANG']))
cosmo_squeezed['SLO_ANG_y'] = np.sin(np.deg2rad(cosmo_squeezed['SLO_ANG']))

#--------------------------------------------------------
# Interpolate
#--------------------------------------------------------

eagle_points['rlon'], eagle_points['rlat'] = rotate_points(pole_lon,pole_lat,eagle_points['x'],eagle_points['y'],'n2r','NH')

len(eagle_points['rlon'].unique())
# Create an xarray DataArray for target points
target_grid = xr.Dataset({
    "lon": (("point1",), eagle_points["rlon"]),
    "lat": (("point2",), eagle_points["rlat"])
})

target_grid = target_grid.set_coords(["lon", "lat"])

cosmo_squeezed = cosmo_squeezed.rename({'x_1':'lon', 'y_1':'lat'})

# Create regridder for bilinear interpolation
regridder = xe.Regridder(cosmo_squeezed, target_grid, method="bilinear")

# Apply regridder to dataset (vectorized for all variables)
ds_interpolated = regridder(cosmo_squeezed[interp_1l+interp_2l])

ds_interpolated['SLO_ANG'] =  np.arctan2(ds_interpolated['SLO_ANG_y'],ds_interpolated['SLO_ANG_x']) #trigonometric in radians
ds_interpolated['SLO_ASP'] =  90-np.rad2deg(np.arctan2(ds_interpolated['SLO_ASP_y'],ds_interpolated['SLO_ASP_x'])) #cardinal in degrees

ds_interpolated.drop_vars(['SLO_ANG_x','SLO_ANG_y','SLO_ASP_x','SLO_ASP_y'])

#store as column in dataframe

for var in ds_interpolated.var():
    
    eagle_points[var] = ds_interpolated[var].values.diagonal(axis1=-2, axis2=-1)
    

#compute orographic lifting proxy
eagle_points['UVspeed_10M'] = np.sqrt(eagle_points['U_10M']**2 + eagle_points['V_10M']**2)
eagle_points['UVdir_10m'] = 270 - (180/np.pi)*np.arctan2(eagle_points['V_10M'],eagle_points['U_10M'])

#calculate orographic uplift
eagle_points['updraft_coeff'] = np.sin(eagle_points['SLO_ANG']) * np.cos(np.deg2rad(eagle_points['UVdir_10m'] - eagle_points['SLO_ASP']))
eagle_points['w_oro'] = eagle_points['UVspeed_10M'] * eagle_points['updraft_coeff']

eagle_points.to_csv(eagle_path+'finaldataset_backgandpresence_timeSteps_'+'annotated.csv')