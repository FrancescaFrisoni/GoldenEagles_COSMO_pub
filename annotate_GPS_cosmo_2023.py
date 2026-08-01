#--------------------------------------------------------- 
# Annotate GPS data by interpolating COSMO (for 2023 data)
# Works with files interpolated with a different 
# script (interpolate_cosmo_time_height_pool_2023.py)
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

from time import time

from netCDF4 import Dataset
from multiprocessing import Pool

t0 = time()
#--------------------------------------------------------
# Settings
#--------------------------------------------------------

#shut down chain assigment warning
pd.options.mode.chained_assignment = None 

colid = 'polyID'
colid = 'individual_local_identifier'

#--------------------------------------------------------
# Define data path 
#--------------------------------------------------------

traj_path = '/net/thermo/atmosdyn/tcarrard/data/eagles/data/'

cosmopath   = '/net/thermo/atmosdyn/tcarrard/data/eagles/data/interpolate/'

constfile   = 'const-cosmo1.nc'

terrainfile = 'terrain_descriptors.nc'

outpath = '/net/thermo/atmosdyn/tcarrard/data/eagles/data/interpolate/'    

#--------------------------------------------------------
# Read trajectories
#--------------------------------------------------------

traj_file  = 'ForFrancesca_backgroundPoints_AtmoAvailability_perInd.csv'

df_traj_all = pd.read_csv(traj_path+traj_file)

#rename column height.above.sea.level
df_traj_all = df_traj_all.rename(columns={'medHeight_asl':'height.above.sea.level'})

#len(np.unique(df_traj_all['unique_segmID']))

event = df_traj_all[colid].iloc[1]

#%%
#--------------------------------------------------------
# Define annotation function
#--------------------------------------------------------

#case = np.unique(df_traj_all['burstID'])[1724:1725]
#df_traj_all = df_traj_all.loc[np.isin(df_traj_all['burstID'],case),:]
#
#np.where(np.unique(df_traj_all['burstID']) == '2191907554_24463')
#event = np.unique(df_traj_all['burstID'])[1724]
#event = 'Adamello20 (eobs 7548)_2023-03-09'

def annotate_at_eagle(event):   
    
    print('Annotating ', event)
    
    #read corresponding COSMO data
    ncfile = cosmopath+event+'_time_height_interpolated.nc'

    interp_data = Dataset(ncfile,'r')

    data = {}                      
    for var in interp_data.variables:
        if var == 'level':
            data[var] = interp_data[var][:]
        else:
            data[var] = interp_data[var][:] #select only ascents
    
    #close netcdf file
    interp_data.close()

    #subset
    df_subset = df_traj_all[df_traj_all[colid]==event].copy().reset_index() 
    
    #compute height above surface
    #df_subset['height.above.surf'] = df_subset['height.above.sea.level'] - data['HSURF'][1,:]
    df_subset['height.above.surf'] = df_subset['medHeight_agl'] + data['HSURF'][1,:]
    
    #compute max N2 below eagle height
    df_subset['N2_max_h'] = 0
    df_subset['N2_max_value'] = 0
    
    for var in ['w_oro', 'ASHFL_S']:
        df_subset[var] = data[var][0,:]
    
    for i in range(len(df_subset)):
                
        eagle_h = df_subset['height.above.surf'].iloc[i]
        
        below_eagle = np.where(data['half_z'][:,i] < eagle_h)
        
        half_z_h = data['half_z'][:,i].squeeze() - data['HSURF'][0,i]
        
        below_3000 = np.where(half_z_h < 2000)
        
        #height of max N2 in the lower 2000m
        df_subset['N2_max_h'].iloc[i] = half_z_h[below_3000][np.argmax(data['N2'][below_3000,i])]
        
        #value of max N2 in the lower 2000m
        df_subset['N2_max_value'].iloc[i] = data['N2'][below_3000,i].max()
        
        #find closest grid
        loc_lev1 = np.array([np.argmin(np.abs(data['height.above.surf'][:,k] -df_subset['height.above.surf'].iloc[k])) for k in range(data['height.above.surf'].shape[1])])

        #find if level is lower or higher
        sign = np.sign([data['height.above.surf'][loc_lev1[k],k] -df_subset['height.above.surf'].iloc[k] for k in range(data['height.above.surf'].shape[1])])
        
        #find second level to interpolate
        loc_lev2 = (loc_lev1 + sign).astype(int)

        cond_low = (loc_lev1==79) & (sign == 1)

        cond_na = df_subset['height.above.surf'].isna() #some heights are NaN

        #set everything to NA if height.above.surf is NA somewhere

        if any(cond_na):
            
            for var in ['U','V','W']:
                df_subset[var] = 'NaN'
            
            df_subset['interp'] = 'NaN'
    
        elif any(cond_low):
    
            ind_low = np.where(cond_low)
    
            #fill loc_lev with lower level to avoid out-of-bound error in weight calculation
            loc_lev2[ind_low] = 78
            
            #df_subset['W'] = W_val1*weights_1/(weights_1+weights_2) + W_val2*weights_2/(weights_1+weights_2)
    
            #find weights (distance to closest model level)
            weights_1 = np.abs([data['height.above.surf'][loc_lev1[k],k] for k in range(len(loc_lev1))]-df_subset['height.above.surf'])
            
            weights_2 = np.abs([data['height.above.surf'][loc_lev2[k],k] for k in range(len(loc_lev2))]-df_subset['height.above.surf'])

            for var in ['U','V','W']:
            
                W_val1 = [data[var][loc_lev1[k],k] for k in range(len(loc_lev1))]
                W_val2 = [data[var][loc_lev2[k],k] for k in range(len(loc_lev2))]
        
                df_subset.loc[:,var] =  W_val1 * (1-(weights_1/(weights_1+weights_2))) + W_val2*(1-(weights_2/(weights_1+weights_2)))
                
                #attribute boundary cases with indexing
                df_subset.loc[:,var].iloc[ind_low] = data[var][79,ind_low]
            
            df_subset['interp'] = ~cond_low

        else:

            #find weights (distance to closest model level)
            weights_1 = np.abs([data['height.above.surf'][loc_lev1[k],k] for k in range(len(loc_lev1))]-df_subset['height.above.surf'])
            
            weights_2 = np.abs([data['height.above.surf'][loc_lev2[k],k] for k in range(len(loc_lev2))]-df_subset['height.above.surf'])
    
            for var in ['U','V','W']:
            
                W_val1 = [data[var][loc_lev1[k],k] for k in range(len(loc_lev1))]
                W_val2 = [data[var][loc_lev2[k],k] for k in range(len(loc_lev2))]
                
                #interpolate at eagle location
                #should be: 
                df_subset.loc[:,var] =  W_val1 * (1-(weights_1/(weights_1+weights_2))) + W_val2*(1-(weights_2/(weights_1+weights_2)))
        
            df_subset['interp'] = True
        
        ##create condition to output NA if the eagle is below the lowest model level
        #if any(np.append(loc_lev1,loc_lev2)>79) or any(np.append(loc_lev1,loc_lev2) ==0):
        #    
        #    #adapt here based on other script
        #    df_subset['U'] = float('NaN')
        #    df_subset['V'] = float('NaN')
        #    df_subset['W'] = float('NaN')
        #else:   
        #    
        #    #loc_lev2 == 80 means that the eagle is above the highest model level which does not make sense
        #    #find where the error lies
        #    
        #    #find weights (distance to closest model level)
        #    weights_1 = np.abs([data['height.above.surf'][loc_lev1[k],k] for k in range(len(loc_lev1))]-df_subset['height.above.surf'])
        #    
        #    weights_2 = np.abs([data['height.above.surf'][loc_lev2[k],k] for k in range(len(loc_lev2))]-df_subset['height.above.surf'])
        #    
        #    U_val1 = [data['U'][loc_lev1[k],k] for k in range(len(loc_lev1))]
        #    U_val2 = [data['U'][loc_lev2[k],k] for k in range(len(loc_lev2))]
        #    V_val1 = [data['V'][loc_lev1[k],k] for k in range(len(loc_lev1))]
        #    V_val2 = [data['V'][loc_lev2[k],k] for k in range(len(loc_lev2))]
        #    W_val1 = [data['W'][loc_lev1[k],k] for k in range(len(loc_lev1))]
        #    W_val2 = [data['W'][loc_lev2[k],k] for k in range(len(loc_lev2))]
        #    
        #    #interpolate at eagle location
        #    df_subset['U'] = U_val1*weights_1/(weights_1+weights_2) + U_val2*weights_2/(weights_1+weights_2)
        #    df_subset['V'] = V_val1*weights_1/(weights_1+weights_2) + V_val2*weights_2/(weights_1+weights_2)
        #    df_subset['W'] = W_val1*weights_1/(weights_1+weights_2) + W_val2*weights_2/(weights_1+weights_2)
    
    return(df_subset)

#%%
#applly function to selected events
if __name__ == "__main__":
    with Pool(25) as pool:
      result = pool.map(annotate_at_eagle, np.unique(df_traj_all[colid]))
    print("Program finished!")
    
#---------------------------------------------------
# Merge and save
#----------------------------------------------------
df_all = pd.concat(result)        

df_all.to_csv(traj_path+traj_file[:-4]+'_annotated_corrected_interp.csv')

df = pd.read_csv(traj_path+traj_file[:-4]+'_annotated_corrected_interp.csv')

df_old = pd.read_csv(traj_path+traj_file[:-4]+'_annotated_new.csv')

import matplotlib.pyplot as plt
plt.scatter(df['W'],df_old['W'])

tend = time() - t0
print("Task took ", tend/60, " minutes to execute")
##inspect NAs
#df['height.above.surf'].hist()
#
#len(df['polyID'].unique())