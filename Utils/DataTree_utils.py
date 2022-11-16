#import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os
import glob
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
import pandas as pd
#import iris
#import iris.quickplot as qplt
#import iris.plot as iplt
#from iris.experimental.regrid import regrid_weighted_curvilinear_to_rectilinear
import json
import cftime
from itertools import product
from cftime import DatetimeNoLeap
import iris_utils
import Gridding
#from nc_processing import *
#from JASMIN_utils import *
#from analysis import * 
#from plotting import *
#from plotting_utils import *
#import esmvalcore.preprocessor
#import xesmf as xe
#import warnings
#%matplotlib inline
#import seaborn as sns
#sns.set()

from datatree import DataTree
#from xmip.preprocessing import rename_cmip6
#from DataTree_utils import *





def add_ensemble_mean(dt, models, exps, ensemble_member_field_str = "ensemble_member"):
    """
    takes a dt, and returns a new dt with additional 'ensemble_mean' ensemble member
    currently only works for a dt with the following structures: 
    /model/experiment/ensemble_member
    or 
    /experiment/ensemble_member
    
    set model=None to use the second dt format. 
    """
    if models:
        for model in models:
            for exp in exps:        
                dt_sub = dt['/{m}/{e}'.format(m=model, e=exp)]
                groups = dt_sub.groups
                leaves = []
                for x in groups:
                    if dt[x].is_leaf:
                        leaves.append(x)
                da_list = []
                for x in leaves:
                    da_list.append(dt[x].ds)
                ds = xr.concat(da_list, ensemble_member_field_str).mean(dim=ensemble_member_field_str)
                #dt_ens_mean['/{m}/{e}'.format(m=model, e=exp)].children = ds
                dt['/{m}/{e}/ensemble_mean'.format(m=model, e=exp)] = DataTree()
                dt['/{m}/{e}/ensemble_mean'.format(m=model, e=exp)].ds = ds
        return dt
    else:
        if exps:
            for exp in exps:        
                dt_sub = dt['/{e}'.format(e=exp)]
                groups = dt_sub.groups
                leaves = []
                for x in groups:
                    if dt[x].is_leaf:
                        leaves.append(x)
                da_list = []
                for x in leaves:
                    da_list.append(dt[x].ds)
                ds = xr.concat(da_list, ensemble_member_field_str).mean(dim=ensemble_member_field_str)
                #dt_ens_mean['/{m}/{e}'.format(m=model, e=exp)].children = ds
                dt['/{e}/ensemble_mean'.format(e=exp)] = DataTree()
                dt['/{e}/ensemble_mean'.format(e=exp)].ds = ds
            return dt 
        else:
            dt_sub = dt
            groups = dt_sub.groups
            leaves = []
            for x in groups:
                if dt[x].is_leaf:
                    leaves.append(x)
            da_list = []
            for x in leaves:
                da_list.append(dt[x].ds)
            ds = xr.concat(da_list, ensemble_member_field_str).mean(dim=ensemble_member_field_str)
            #dt_ens_mean['/{m}/{e}'.format(m=model, e=exp)].children = ds
            dt['/ensemble_mean'] = DataTree()
            dt['/ensemble_mean'].ds = ds 

