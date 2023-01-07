import iris
import Utils.iris_utils as iris_utils
import pandas as pd
import numpy as np
import os
import logging
import esmvalcore.preprocessor
import glob
import warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
import xarray as xr
from xmip.preprocessing import rename_cmip6

exps = ["ssp370", "ssp126"]
#exps = ["historical"]
#exps = ['ssp245']
#exps = ['ssp126']
dirs = []
var_path = "Amon/tas"
for experiment in exps:
    if experiment == "historical":
        exp_set = "CMIP"
    else:
        exp_set = "ScenarioMIP"
    for x in glob.glob('/badc/cmip6/data/CMIP6/{es}/*/*/{e}/*/{v}/*/latest/'.format(es=exp_set, e=experiment, v=var_path)):
        dirs.append(x)
dirs.reverse()
print(len(dirs))

def preprocess_multi_ens(folder, arctic_cut_off=66):
    
    """ makes a df contining global mean and {global except arctic} mean
        temp by year, for the first ensemble member of each model """
    
    outpath = 'int_outputs/temperature_multi_ens/{M}_{Exp}_{Ens}.csv'.format(
                M=folder.split('/')[7], Exp=folder.split('/')[8], Ens=folder.split('/')[9])
    
    if os.path.exists(outpath):
        return
    else:                                                              
        try:
            data = rename_cmip6(xr.open_mfdataset(folder + "*.nc", use_cftime=True))
            winter_mask = data.time.dt.month.isin([12,1,2])
            jan_mask = data.time.dt.month.isin([1]) 
            name = str(folder.split('/')[7] + '_' + folder.split('/')[8] + '_' + folder.split('/')[9])
            
            at_data = data['tas']#.isel(member_id=0)
            
            #at_data = at_data[winter_mask]
            
            years = data.time.dt.year[jan_mask].compute()
            
            # month_length = data.time.dt.days_in_month
            # weights = ( month_length.groupby("time.season") / month_length.groupby("time.season").mean())
            world_annual = (at_data).groupby("time.year").mean(dim="time")
            arctic_annual = world_annual.sel(y=slice(-90,arctic_cut_off))
            
            world_w = world_annual.weighted(weights=np.cos(np.deg2rad(world_annual.y)))
            arctic_w = arctic_annual.weighted(weights=np.cos(np.deg2rad(arctic_annual.y)))
            
            df = pd.DataFrame({'no_arctic_tas':arctic_w.mean(("x","y")).compute().values,
                               'world_tas':world_w.mean(("x","y")).compute().values,
                               'year':years.values})
            
            df.set_index('year',inplace=True)
            df.sort_index(inplace=True)
            df['Model'] = folder.split('/')[7]
            df['Experiment'] = folder.split('/')[8]
            df['Ensemble_member'] = folder.split('/')[9]
            df.to_csv(outpath) 
        except:
            print(name)

for dir in tqdm(dirs):
    preprocess_multi_ens(dir)
