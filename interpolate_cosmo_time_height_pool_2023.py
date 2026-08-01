#--------------------------------------------------------- 
# Create vertical profiles along eagle trajectories
# Tom Carrard
# 22.11.23
#---------------------------------------------------------

#%%
#---------------------------------------------------------
# Import Modules
#---------------------------------------------------------
import sys
sys.path.append('/home/tcarrard/eagle.gravity.waves/utilities/')


import pandas as pd
import numpy as np
import datetime as dt

import random

from time import time

from netCDF4 import Dataset
import my_dypy.netcdf as netcdf
from rotate_grid import rotate_points

from gridcalc_cosmo import calc_th,calc_uv_vec_new,calc_uv_norm_parallel

from my_dypy.intergrid import Intergrid

from multiprocessing import Pool


t0 = time()
#--------------------------------------------------------
# Settings
#--------------------------------------------------------

#Define variables to extract and to interpolate (Hsurf and z are automatically output)

extractvars = ['P','T','U','V','W','ASHFL_S','T_2M','TKE','U_10M','V_10M'] 

interpvars = ['P','T','U','V','W', 'TH','ASHFL_S','T_2M','TKE','U_10M','V_10M'] 

# define whether event is a full 'traj', a segment 'burst', or a single 'ascent'
event_class = 'individual_local_identifier'

#select specific trajectories (uncomment on the line you want)
cases = "all"
#cases = ['Grosio_3_16349491472','Vrata20_3_16349491472','Sampuoir2_69_14179308152','Dischma2_123_14562392367','Aosta2_20_195_16499840535','Dischma2_208_16654433929','Dischma1_208_16654433929','Dischma2_101_15723964013','Dischma1_61_14362600594','Fahrntal19_57_14359287951','Flüela19_61_14362600594','Almen19_39_14349185515','Dischma1_39_14349185515','Torta19_39_14349185515','Sinestra2_13_14067948153','Sinestra2_32_14113867603','Sinestra1_25_14110618036','Flüela19_122_14547821115','Flüela19_122_14547821115','Fahrntal19_120_14543786250','Aosta1_20_192_16374511443','Matsch19_39_14349185515','Adamello20_199_16502138670','Sinestra1_19_14101325618','Tuors1_203_16654522426','Sinestra2_2_14028254712','Sinestra1_148_14642021476','Sinestra2_6_14044841784','Trimmis20_204_16653337706']
#cases = ['Sinestra1_3915_16447103304','Sinestra2_175_14151651053','Sinestra2_177_14154233332']
#cases  = ['1167139393_31059','1187719633_50841','1187719633_51069','2191907554_26035','2313423539_12172','1589292475_36538']
#cases = 150

#--------------------------------------------------------- 
# auxiliary functions
#---------------------------------------------------------
#***************************************#
def make_datelist(startdate,enddate,hstep):
    nowdate = startdate
    datelist = []
    while nowdate <= enddate:
        datelist.append(nowdate)
        nowdate = nowdate + dt.timedelta(hours=hstep)

    return datelist
#***************************************#

#define function to calculate Brunt-Väisala frequency
def calculate_bruntVaisala(z,TH, g = 9.81):
    
    if z.ndim == 3:
        #calculate vertical change of Theta
        dTHdz_levels = (TH[:-1,:,:] - TH[1:,:,:]) / (z[:-1,:,:] - z[1:,:,:])
        
        #define staggered Theta
        TH_staggered = (TH[:-1,:,:] + TH[1:,:,:])/2
    
    elif z.ndim == 2:
        #calculate vertical change of Theta
        dTHdz_levels = (TH[:-1,:] - TH[1:,:]) / (z[:-1,:] - z[1:,:])
        
        #define staggered Theta
        TH_staggered = (TH[:-1,:] + TH[1:,:])/2
    
    #compute Brunt Väisälä frequency
    N2 = dTHdz_levels / TH_staggered * g
    
    return N2 
#--------------------------------------------------------
# Define data path 
#--------------------------------------------------------

traj_path = '/net/thermo/atmosdyn/tcarrard/data/eagles/data/'

cosmopath   = '/net/litho/atmosdyn2/jansingl/KENDA-1/ANA23/'
constfile   = 'laf2023010100-const.nc'

terrainfile = 'terrain_descriptors.nc'

outpath = '/net/thermo/atmosdyn/tcarrard/data/eagles/data/interpolate/'

#--------------------------------------------------------
# Read trajectories
#--------------------------------------------------------

#uncomment line to choose file
traj_file  = 'ForFrancesca_backgroundPoints_AtmoAvailability_perInd.csv'

df_traj_all = pd.read_csv(traj_path+traj_file)

df_traj_all = df_traj_all.rename({'x': 'location.long','y':'location.lat'},axis=1)

#------------------------------------------------------
# Read constant file (surface and grid information)
#------------------------------------------------------

cfile = cosmopath+constfile
print('import '+cfile)
const = Dataset(cfile,'r')

lon   = const.variables['lon_1'][:].data
lat   = const.variables['lat_1'][:].data
rlon  = const.variables['x_1'][:].data
rlat  = const.variables['y_1'][:].data
level = const.variables['z_1'][:].data
hsurf = const.variables['HSURF'][:].data
z     = const.variables['HEIGHT'][:].data[:-1,:,:]

slo_asp_cosmo = const.variables['SLO_ASP'][:].data * 180/np.pi
slo_ang_cosmo = const.variables['SLO_ANG'][:].data * 180/np.pi

pole_lon = const.variables['grid_mapping_1'].grid_north_pole_longitude
pole_lat = const.variables['grid_mapping_1'].grid_north_pole_latitude

const.close()

#--------------------------------------------------------
# Subset events
#--------------------------------------------------------

if cases=='all':
    df_traj_all = df_traj_all

elif isinstance(cases, int):
    
    #define random seed
    np.random.seed(1997)
    
    #take the first n cases 
    cases_ID = np.random.choice(np.unique(df_traj_all['uplift_ID']),cases,replace=False)
    
    df_traj_all = df_traj_all.loc[np.isin(df_traj_all['uplift_ID'],cases_ID),:]
    
else:
    cases_ind   = np.where(np.isin(df_traj_all['burstID'], cases))[0]
    
    #subset cases
    df_traj_all = df_traj_all.loc[cases_ind,:]
    
#%% 
#--------------------------------------------------------
# Define interpolation function
#--------------------------------------------------------
#record initial time
t0 = time()

#event = df_traj_all['burstID'].iloc[0]

def interp_at_eagle(event):
    
    print(event+': Interpolating along eagle trajectory..')
            
    #subset trajectories single ascent
    df_subset = df_traj_all[df_traj_all[event_class]==event].copy().reset_index()
        
    #define lower and upper time limit (no minutes allowed)
    t_min = dt.datetime.strptime(df_subset['timestamp'].min(), '%Y-%m-%d %H:%M:%S').replace(minute=0, second=0, microsecond=0)
    t_max = dt.datetime.strptime(df_subset['timestamp'].max(), '%Y-%m-%d %H:%M:%S').replace(minute=0, second=0, microsecond=0) + dt.timedelta(hours=1)

    #Define datelist (used to iterate in next loop)
    datelist = np.array(make_datelist(t_min,t_max,1))
    
    #---------------------------------------------------------
    # prepare the interpolation
    #---------------------------------------------------------
    #get eagle coordinates on rotated grid (because we interpolate on the rotated grid)
    df_subset['rlon'],df_subset['rlat'] = rotate_points(pole_lon,pole_lat,df_subset['location.long'],df_subset['location.lat'],'n2r','NH')
    
    #define range for 2d Intergrid
    lo = np.array([rlat[0],rlon[0]])
    hi = np.array([rlat[-1],rlon[-1]])
    
    #create interpolator function for hsurf (function that interpolate hsurf on a new grid)
    f_hsurf = Intergrid(hsurf,lo=lo,hi=hi,verbose=False)
    
    #---------------------------------------------------------
    # Interpolate slope aspect and slope angle
    #---------------------------------------------------------
    #convert aspect to trig coordinate
    slo_asp_cosmo_trig = 90-slo_asp_cosmo
    
    #convert to cartesian coordinates
    x_asp = np.cos(np.deg2rad(slo_asp_cosmo_trig))
    y_asp = np.sin(np.deg2rad(slo_asp_cosmo_trig))
    x_ang = np.cos(np.deg2rad(slo_ang_cosmo))
    y_ang = np.sin(np.deg2rad(slo_ang_cosmo))
    
    #interpolate x and y
    f_asp_x   = Intergrid(x_asp,lo=lo,hi=hi,verbose=False)
    f_asp_y   = Intergrid(y_asp,lo=lo,hi=hi,verbose=False)
    f_ang_x   = Intergrid(x_ang,lo=lo,hi=hi,verbose=False)
    f_ang_y   = Intergrid(y_ang,lo=lo,hi=hi,verbose=False)
    
    #---------------------------------------------------------
    # prepare dictionary for output
    #---------------------------------------------------------
    eagle_outdata = {}
    
    # define dimension of output data
    for var in np.append(interpvars, ['HSURF','z']):
        eagle_outdata[var] = np.zeros([len(z), len(df_subset)])
    
    savedstr = None
    
    #---------------------------------------------------------
    # loop over datelist (necessary because cosmo file are stored for single time steps)
    #---------------------------------------------------------
    timestamp = pd.to_datetime(df_subset['timestamp']).copy()
    
    #create storing array for constants 
    z_interp        = np.zeros([len(z), len(df_subset)])
    hsurf_interp    = np.zeros([len(df_subset)])
    
    slo_ang_interp_x  = np.zeros([len(df_subset)])
    slo_ang_interp_y  = np.zeros([len(df_subset)])

    slo_asp_interp_x  = np.zeros([len(df_subset)])
    slo_asp_interp_y  = np.zeros([len(df_subset)])
        
    for nowdate in datelist:
        
        #---------------------------------------------------------
        # Extract eagle data points within the time window
        #---------------------------------------------------------
        
        #locate points
        min_time_cond = timestamp >= nowdate
        max_time_cond = timestamp < nowdate + dt.timedelta(hours=1)
        now_ind       = np.where(min_time_cond & max_time_cond)[0]
        
        #if some points are within the time window load COSMO
        if now_ind.size > 0:
    
            #---------------------------------------------------------
            # load data of respective time step
            #---------------------------------------------------------
            nowdatestr = nowdate.strftime('%Y%m%d%H')
            #check if we do not have it already
            if savedstr != nowdatestr:
            
                if nowdate.month < 10:
                    monthstr = '0'+str(nowdate.month)
                else:
                    monthstr = str(nowdate.month)
                pfile      = cosmopath+'/laf'+nowdatestr+'.nc'
    
                print('import '+pfile)
                rawdata = Dataset(pfile,'r')
                data_nowdate = {}
                for var in np.append(extractvars, 'z_1'):
    
                    data_nowdate[var] = rawdata[var][:].data.squeeze() #extract data as array
                
                #destagger W
                W_destagger_now  = (data_nowdate['W'][:-1,:,:] + data_nowdate['W'][1:,:,:])/2
                data_nowdate['W']  = W_destagger_now
                
                #calculate potential temperature
                TH_now = calc_th(data_nowdate['T'],data_nowdate['P']/100.)
                data_nowdate['TH'] = TH_now
                
                #sys.exit(1)
    
            #---------------------------------------------------------
            # load data of the consecutive time step
            #---------------------------------------------------------
            nextdate = nowdate + dt.timedelta(hours=1)
            nextdatestr = nextdate.strftime('%Y%m%d%H')
            if nextdate.month < 10:
                monthstr = '0'+str(nextdate.month)
            else:
                monthstr = str(nextdate.month)
            pfile      = cosmopath+'/laf'+nextdatestr+'.nc'
    
            print('import '+pfile)
            rawdata = Dataset(pfile,'r')
            data_nextdate = {}
            for var in np.append(extractvars, 'z_1'):
    
                data_nextdate[var] = rawdata[var][:].data.squeeze()
    
            # destagger W
            W_destagger_next = (data_nextdate['W'][:-1,:,:] + data_nextdate['W'][1:,:,:])/2    
            data_nextdate['W'] = W_destagger_next
            
            #calculate potential temperature
            TH_next = calc_th(data_nextdate['T'],data_nextdate['P']/100.)
            data_nextdate['TH'] = TH_next
                            
            #---------------------------------------------------------
            # interpolate to the eagle positions
            #---------------------------------------------------------
            
            # Step 1: Get eagle positions
            points = np.array([df_subset['rlat'][now_ind],df_subset['rlon'][now_ind]]).T
            
            #Step 2: Interpolate to eagles positions for nowdate and nextdate
            eagle_nowdate  = {}
            eagle_nextdate = {}
            
            for var in interpvars:
                eagle_nowdate[var]  = np.zeros([len(z), len(now_ind)])
                eagle_nextdate[var] = np.zeros([len(z), len(now_ind)])
            
            #compute surface height interpolation
            hsurf_interp[now_ind]  = f_hsurf(points)
            
            slo_ang_interp_x[now_ind] = f_ang_x(points)
            slo_ang_interp_y[now_ind] = f_ang_y(points)
            slo_asp_interp_x[now_ind] = f_asp_x(points)
            slo_asp_interp_y[now_ind] = f_asp_y(points)
            
            #back transform to polar coordinates
            slo_ang_interp = np.arctan2(slo_ang_interp_y,slo_ang_interp_x) #trigonometric in radians
            slo_asp_interp = 90-np.rad2deg(np.arctan2(slo_asp_interp_y,slo_asp_interp_x)) #cardinal in degrees

            #loop over vertical levels
            level = data_nowdate['z_1']
            
            print('Interpolating in 2d for each model level ...' )
            
            for lev in range(len(level)):
        
                f_z = Intergrid(z[lev,:,:],lo=lo,hi=hi,order=1,verbose=False)
                z_interp[lev][now_ind] = f_z(points)
                
                for var in interpvars:
                    
                    if data_nowdate[var].ndim == 2:
                        f_nowdate       = Intergrid(data_nowdate[var],lo=lo,hi=hi,order=1,verbose=False)
                        eagle_nowdate[var][lev,:]  = f_nowdate(points)
                        f_nextdate          = Intergrid(data_nextdate[var],lo=lo,hi=hi,order=1,verbose=False)
                        eagle_nextdate[var][lev,:] = f_nextdate(points)
                        
                    else:   
                        #interpolate
                        f_nowdate           = Intergrid(data_nowdate[var][lev,:,:],lo=lo,hi=hi,order=1,verbose=False)
                        eagle_nowdate[var][lev,:]  = f_nowdate(points)
                        f_nextdate          = Intergrid(data_nextdate[var][lev,:,:],lo=lo,hi=hi,order=1,verbose=False)
                        eagle_nextdate[var][lev,:] = f_nextdate(points)
                    
            #----------------------------------------------------------
            #Step 3: Temporal interpolation
            #----------------------------------------------------------
            #Calculate the temporal weights for interpolation
            #these are zero if we are at nextdate, 1 if we are at nowdate
            timediff = np.array( [( (date - nowdate).days*24*3600 + (date - nowdate).seconds ) for date in timestamp[now_ind]] )
            weights = 1 - (timediff / 3600)
            
            for var in eagle_nowdate.keys():
                
                print('Interpolating ', var, ' to GPS time...')
                
                eagle_outdata[var][:,now_ind] = eagle_nowdate[var]*weights + (1 - weights)*eagle_nextdate[var]
            
            #---------------------------------------------------------
            # save nowdate for the next timestep
            #---------------------------------------------------------
            savedstr     = nextdatestr
            data_nowdate = data_nextdate
    
        #---------------------------------------------------------
        # time steps without eagles
        #---------------------------------------------------------
        else:
            savedstr = 'None'
    
    #store z and hsurf and coordinates (repeat arrays to match output dimensions)
    eagle_outdata['HSURF']     = np.repeat(np.expand_dims(hsurf_interp, axis=0), 80, axis = 0)
    eagle_outdata['z']         = z_interp
    eagle_outdata['half_z']    = (z_interp[:-1,:] + z_interp[1:,:]) / 2
    eagle_outdata['timestamp'] = np.repeat(np.expand_dims(timestamp.to_numpy(), axis=0), 80, axis = 0)
    eagle_outdata['longitude'] = np.repeat(np.expand_dims(df_subset['location.long'].to_numpy(),axis=0), 80, axis = 0)
    eagle_outdata['latitude']  = np.repeat(np.expand_dims(df_subset['location.lat'].to_numpy(),axis=0), 80, axis = 0)
    
    #-------------------------------------------------------------
    # Compute secondary variables 
    #-------------------------------------------------------------
    #compute height above surface
    eagle_outdata['height.above.surf'] = eagle_outdata['z'] - eagle_outdata['HSURF']
    
    #compute Brunt-Vaisala frequency
    eagle_outdata['N2'] = calculate_bruntVaisala(z = eagle_outdata['z'], TH=eagle_outdata['TH'])
    
    #compute orographic uplift
    #calculate horizontal wind speed and direction
    eagle_outdata['vel_10m'] = np.sqrt(eagle_outdata['U_10M']**2 + eagle_outdata['V_10M']**2)
    eagle_outdata['dir_10m'] = 270 - (180/np.pi)*np.arctan2(eagle_outdata['V_10M'],eagle_outdata['U_10M'])

    #calculate orographic uplift
    eagle_outdata['updraft_coeff'] = np.sin(slo_ang_interp) * np.cos(np.deg2rad(eagle_outdata['dir_10m'] - slo_asp_interp))
    
    eagle_outdata['w_oro'] = eagle_outdata['vel_10m'] * eagle_outdata['updraft_coeff']
    eagle_outdata['w_oro'][eagle_outdata['w_oro'] < 0] = 0
        
    print('Successful interpolation for ', event)
    
    #-------------------------------------------------------------
    # Reshape wind vector
    #-------------------------------------------------------------
    
    #define timestep (~seconds)
    timestep = 5
    
    #rotate wind vectors from grid north to true north (because cosmo grid is a rotated lon-lat grid)
    eagle_outdata['U'],eagle_outdata['V'] = calc_uv_vec_new(eagle_outdata['U'],eagle_outdata['V'],eagle_outdata['longitude'],eagle_outdata['latitude'],pollon=pole_lon,pollat=pole_lat)
    
    #create empty array, then create new winds at array location
    eagle_outdata['U_norm']     = np.zeros(eagle_outdata['U'].shape)
    eagle_outdata['U_parallel'] = np.zeros(eagle_outdata['U'].shape)
    
    for i in range(5, len(df_subset)-5):
        
        eagle_outdata['U_norm'][:,i],eagle_outdata['U_parallel'][:,i] = calc_uv_norm_parallel(eagle_outdata['U'][:,i],eagle_outdata['V'][:,i],df_subset.loc[i-5, 'location.long'],df_subset.loc[i-5, 'location.lat'],df_subset.loc[i+5, 'location.long'],df_subset.loc[i+5, 'location.lat'])
    
    #-------------------------------------------------------------
    #save as netcdf
    #-------------------------------------------------------------
    print('Saving of ', event+"time_height_interpolated.nc ...", )
        
    interp_output = Dataset(outpath+event+"_time_height_interpolated.nc", "w", format="NETCDF4")
    
    #create dimensions
    interp_output.createDimension("level", len(z))
    interp_output.createDimension("timestep", len(timestamp))
    interp_output.createDimension("halflevel", len(z)-1)
    
    #create variables
    levels = interp_output.createVariable('level','u8',("level"))
    
    levels[:] = np.arange(1,81)
    
    #loop over variables
    for key in eagle_outdata.keys():
        
        if key == 'timestamp':
            
            #create variable
            interp_output.createVariable(key,'u8',('level','timestep'))
            
            #update variable values
            interp_output[key][:] = eagle_outdata['timestamp']
            
        elif key == "N2" or key == "half_z":
            
            #create variable
            interp_output.createVariable(key,'f8',('halflevel','timestep'))
            
            #update variable values 
            interp_output[key][:] = eagle_outdata[key]
                    
        else:
            #create variable
            interp_output.createVariable(key,'f8',('level','timestep'))
            
            #update variable values
            interp_output[key][:] = eagle_outdata[key]
            
    interp_output.close()
            
    print(event+"time_height_interpolated.nc saved in "+ outpath)
    
#---------------------------------------------------
# Apply function to events
#----------------------------------------------------
    
#applly function to selected events
if __name__ == "__main__":
    with Pool(15) as pool:
      result = pool.map(interp_at_eagle, np.unique(df_traj_all[event_class]))
    print("Program finished!")

tend = time() - t0
print("Task took ", tend/60, " minutes to execute")
# %%
