#--------------------------------------------------------- 
# Interpolate and annotate vertical velocity
# Tom Carrard
# 05.08.24
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

from my_dypy.intergrid import Intergrid

from multiprocessing import Pool


t0 = time()
#--------------------------------------------------------
# Settings
#--------------------------------------------------------

#Define variables to extract and to interpolate (Hsurf and z are automatically output)

extractvars = ['W','U','V'] 

interpvars = ['W','U','V'] 

# define whether event is a full 'traj', a segment 'burst', or a single 'ascent'
event_class = 'burst'

#select specific trajectories (uncomment on the line you want)
cases = "all"
#cases = ['Grosio_3_16349491472','Vrata20_3_16349491472','Sampuoir2_69_14179308152','Dischma2_123_14562392367','Aosta2_20_195_16499840535','Dischma2_208_16654433929','Dischma1_208_16654433929','Dischma2_101_15723964013','Dischma1_61_14362600594','Fahrntal19_57_14359287951','Flüela19_61_14362600594','Almen19_39_14349185515','Dischma1_39_14349185515','Torta19_39_14349185515','Sinestra2_13_14067948153','Sinestra2_32_14113867603','Sinestra1_25_14110618036','Flüela19_122_14547821115','Flüela19_122_14547821115','Fahrntal19_120_14543786250','Aosta1_20_192_16374511443','Matsch19_39_14349185515','Adamello20_199_16502138670','Sinestra1_19_14101325618','Tuors1_203_16654522426','Sinestra2_2_14028254712','Sinestra1_148_14642021476','Sinestra2_6_14044841784','Trimmis20_204_16653337706']
#cases = ['Sinestra1_3915_16447103304','Sinestra2_175_14151651053','Sinestra2_177_14154233332']
#cases = 10

#define upper and lower time limits (based on cosmo availability)
t_lower  = dt.datetime(2020,3,1,00) 
t_higher = dt.datetime(2020,10,29,12)

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

#--------------------------------------------------------
# Define data path 
#--------------------------------------------------------

traj_path = '/net/thermo/atmosdyn/tcarrard/data/eagles/data/'

cosmopath   = '/net/litho/atmosdyn2/jansingl/cosmo.analysis.eagle.thermals/'
constfile   = 'const-cosmo1.nc'

outpath = '/net/thermo/atmosdyn/tcarrard/eagles/data/interpolate/'

#--------------------------------------------------------
# Read trajectories
#--------------------------------------------------------

#uncomment line to choose file

traj_file  = 'GPS_allsoaringpoints_areavol_subset_2020.csv'
#traj_file  = 'uplifts_case_studies_labelled_new_corrected.csv'

df_traj_all = pd.read_csv(traj_path+traj_file)

df_traj_all['timestamp'] = pd.to_datetime(df_traj_all['timestamp'])

#len(np.unique(df_traj_all['unique_segmID']))

#remove out-of-bound timesteps
df_traj_all = df_traj_all.loc[(df_traj_all['timestamp'] < t_higher) & (df_traj_all['timestamp'] > t_lower)]

#------------------------------------------------------
# Subset (optional)
#------------------------------------------------------

##subset df
#case = np.unique(df_traj_all['burstID'])[1250:1350]
##1350 index
#df_traj_all = df_traj_all.loc[np.isin(df_traj_all['burstID'],case),:]

#------------------------------------------------------
# Read constant file (surface and grid information)
#------------------------------------------------------

cfile = cosmopath+'/const/'+constfile
print('import '+cfile)
const = Dataset(cfile,'r')

lon   = const.variables['lon_1'][:].data
lat   = const.variables['lat_1'][:].data
rlon  = const.variables['x_1'][:].data
rlat  = const.variables['y_1'][:].data
level = const.variables['z_1'][:].data
hsurf = const.variables['HSURF'][:].data
z     = (const.variables['HHL'][1:,...].data+const.variables['HHL'][:-1, ...].data)/2

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
    cases_ID = np.random.choice(np.unique(df_traj_all['burstID']),cases,replace=False)
    
    df_traj_all = df_traj_all.loc[np.isin(df_traj_all['burstID'],cases_ID),:]
    
else:
    cases_ind   = np.where(np.isin(df_traj_all['ascent_ID'], cases))[0]
    
    #subset cases
    df_traj_all = df_traj_all.loc[cases_ind,:]

#%% 
#--------------------------------------------------------
# Define interpolation function
#--------------------------------------------------------
#event = df_traj_all['ascent_ID']
def interp_at_eagle(event,event_class=event_class):
    
    print(event+': Interpolating along eagle trajectory..')
    
    if event_class == 'ascent':
        
        #subset trajectories single ascent
        df_subset = df_traj_all[df_traj_all['uplift_ID']==event].copy().reset_index()
        
    elif event_class == 'burst':
        
        #subset
        df_subset = df_traj_all[df_traj_all['burstID']==event].copy().reset_index()

    else:
        #subset trajectories single case study
        df_subset = df_traj_all[df_traj_all['case_studyID']==event].copy().reset_index()
        
    print(df_subset.index.size)

    #define lower and upper time limit (no minutes allowed)
    t_min = dt.datetime.strptime(str(df_subset['timestamp'].min()), '%Y-%m-%d %H:%M:%S').replace(minute=0, second=0, microsecond=0)
    t_max = dt.datetime.strptime(str(df_subset['timestamp'].max()), '%Y-%m-%d %H:%M:%S').replace(minute=0, second=0, microsecond=0) + dt.timedelta(hours=1)

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
                pfile      = cosmopath+'/'+str(nowdate.year)+'/'+monthstr+'/laf'+nowdatestr+'.nc'
    
                print('import '+pfile)
                rawdata = Dataset(pfile,'r')
                data_nowdate = {}
                for var in np.append(extractvars, 'z_1'):
    
                    data_nowdate[var] = rawdata[var][:].data.squeeze() #extract data as array
                
                #destagger W
                W_destagger_now  = (data_nowdate['W'][:-1,:,:] + data_nowdate['W'][1:,:,:])/2
                data_nowdate['W']  = W_destagger_now
                
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
            pfile      = cosmopath+'/'+str(nextdate.year)+'/'+monthstr+'/laf'+nextdatestr+'.nc'
    
            print('import '+pfile)
            rawdata = Dataset(pfile,'r')
            data_nextdate = {}
            for var in np.append(extractvars, 'z_1'):
    
                data_nextdate[var] = rawdata[var][:].data.squeeze()
    
            # destagger W
            W_destagger_next = (data_nextdate['W'][:-1,:,:] + data_nextdate['W'][1:,:,:])/2    
            data_nextdate['W'] = W_destagger_next
            
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
    
    #store z and hsurf and coordinates (repeat arrays to match output dimensions)
    eagle_outdata['HSURF']     = np.repeat(np.expand_dims(hsurf_interp, axis=0), 80, axis = 0)
    eagle_outdata['z']         = z_interp
    eagle_outdata['half_z']    = (z_interp[:-1,:] + z_interp[1:,:]) / 2
    eagle_outdata['timestamp'] = np.repeat(np.expand_dims(timestamp.to_numpy(), axis=0), 80, axis = 0)
    eagle_outdata['longitude'] = np.repeat(np.expand_dims(df_subset['location.long'].to_numpy(),axis=0), 80, axis = 0)
    eagle_outdata['latitude']  = np.repeat(np.expand_dims(df_subset['location.lat'].to_numpy(),axis=0), 80, axis = 0)
    #eagle_outdata['distance']  = np.repeat(np.expand_dims(df_subset['distance'].to_numpy(),axis=0), 80, axis=0)
    
    #-------------------------------------------------------------
    # Compute secondary variables 
    #-------------------------------------------------------------
    #compute height above surface
    eagle_outdata['height.above.surf'] = eagle_outdata['z'] - eagle_outdata['HSURF']
        
    print('Successful interpolation for ', event)
    
    #compute height above surface
    df_subset['height.above.surf'] = df_subset['height.above.sea.level'] - eagle_outdata['HSURF'][1,:]
    
        
        #find closest grid
    loc_lev1 = np.array([np.argmin(np.abs(eagle_outdata['height.above.surf'][:,k] -df_subset['height.above.surf'].iloc[k])) for k in range(eagle_outdata['height.above.surf'].shape[1])])
    
    #find if level is lower or higher
    sign = np.sign([eagle_outdata['height.above.surf'][loc_lev1[k],k] - df_subset['height.above.surf'].iloc[k] for k in range(eagle_outdata['height.above.surf'].shape[1])])

    #find second level to interpolate
    loc_lev2 = (loc_lev1 + sign).astype(int)

    #locate here occurences of point where the eagle is below lowest level

    df_subset['W'] = float('NaN')

    cond_low = (loc_lev1==79) & (sign == 1)

    cond_na = df_subset['height.above.surf'].isna() #some heights are NaN

    #set everything to NA if height.above.surf is NA somewhere

    if any(cond_na):
        
        for var in interpvars:
            df_subset[var] = 'NaN'
        
        df_subset['interp'] = 'NaN'

    elif any(cond_low):

        ind_low = np.where(cond_low)

        #fill loc_lev with lower level to avoid out-of-bound error in weight calculation
        loc_lev2[ind_low] = 78
        
        #df_subset['W'] = W_val1*weights_1/(weights_1+weights_2) + W_val2*weights_2/(weights_1+weights_2)

        #find weights (distance to closest model level)
        weights_1 = np.abs([eagle_outdata['height.above.surf'][loc_lev1[k],k] for k in range(len(loc_lev1))]-df_subset['height.above.surf'])
        
        weights_2 = np.abs([eagle_outdata['height.above.surf'][loc_lev2[k],k] for k in range(len(loc_lev2))]-df_subset['height.above.surf'])

        for var in interpvars:
        
            W_val1 = [eagle_outdata[var][loc_lev1[k],k] for k in range(len(loc_lev1))]
            W_val2 = [eagle_outdata[var][loc_lev2[k],k] for k in range(len(loc_lev2))]
    
            df_subset.loc[:,var] =  W_val1 * (1-(weights_1/(weights_1+weights_2))) + W_val2*(1-(weights_2/(weights_1+weights_2)))
            
            #attribute boundary cases with indexing
            df_subset.loc[:,var].iloc[ind_low] = eagle_outdata[var][79,ind_low]

        df_subset['interp'] = ~cond_low

    else:

        #find weights (distance to closest model level)
        weights_1 = np.abs([eagle_outdata['height.above.surf'][loc_lev1[k],k] for k in range(len(loc_lev1))]-df_subset['height.above.surf'])
        
        weights_2 = np.abs([eagle_outdata['height.above.surf'][loc_lev2[k],k] for k in range(len(loc_lev2))]-df_subset['height.above.surf'])

        for var in interpvars:
        
            W_val1 = [eagle_outdata[var][loc_lev1[k],k] for k in range(len(loc_lev1))]
            W_val2 = [eagle_outdata[var][loc_lev2[k],k] for k in range(len(loc_lev2))]
            
            #interpolate at eagle location
            #should be: 
            df_subset.loc[:,var] =  W_val1 * (1-(weights_1/(weights_1+weights_2))) + W_val2*(1-(weights_2/(weights_1+weights_2)))
    
        df_subset['interp'] = True

    return(df_subset)
    
#uplift_df.to_csv(traj_path+traj_file[-4]+'annotated.csv')
    
#---------------------------------------------------
# Apply function to events
#----------------------------------------------------


#apply to ascent or full trajectory segment
if event_class == 'ascent':
    
    #applly function to selected events
    if __name__ == "__main__":
        with Pool() as pool:
          result = pool.map(interp_at_eagle, np.unique(df_traj_all['uplift_ID']))
        print("Program finished!")
        
elif event_class == 'burst':
    
    #applly function to selected events
    if __name__ == "__main__":
        with Pool(20) as pool:
          result = pool.map(interp_at_eagle, np.unique(df_traj_all['burstID']))
        print("Program finished!")
    
else:   
    #applly function to selected events
    if __name__ == "__main__":
        with Pool() as pool:
          result = pool.map(interp_at_eagle, np.unique(df_traj_all['case_studyID']))
        print("Program finished!")
        
#---------------------------------------------------
# Merge and save
#----------------------------------------------------
df_all = pd.concat(result)        

df_all.to_csv(traj_path+traj_file[:-4]+'_annotated_fullwind.csv')

tend = time() - t0
print("Task took ", tend/60, " minutes to execute")
# %%
#
#len(np.unique(df13['unique_segmID']))
#
#len(np.unique(df_all['unique_segmID']))
#
df0 =pd.read_csv(traj_path+traj_file[:-4]+'_annotated_W.csv')
df1 =pd.read_csv(traj_path+traj_file[:-4]+'_annotated_fullwind.csv')

df0['W'].isna().sum()
df1['W'].isna().sum()

import matplotlib.pyplot as plt

plt.scatter(df0['W'],df1['W'])

#
#segm_ann = np.unique(df_all['unique_segmID'])
#segm_0   = np.unique(df0['unique_segmID'])
#
#segm_pb = segm_0[np.isin(segm_0,segm_ann, invert=True)]
#
#np.where(np.unique(df0['unique_segmID'])=='Aosta2_20 (eobs 7558)_1193205346_3088_soar_1')
#dtype_dic = {'unique_segmID':str}
#df_all = pd.read_csv(traj_path+traj_file[:-8]+'annotated1.csv',dtype=dtype_dic)
#len(np.unique(df_all['unique_segmID']))
#
#any(df_all['unique_segmID'] == float('NaN'))

