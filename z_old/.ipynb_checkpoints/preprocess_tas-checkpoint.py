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

def temperatures(folder, var_long = 'air_temperature', arctic_cut_off=66):
    """ 
    calculate global, arctic and antarctic monthly temperatures, and return as a DF
    
    """
    
    ### read in any files in folder under each variable name, and concatenate on time (in case there is more than one)
    files = []
    #hack to take only the 2015-2100 file in the CESM ssp585 folder:
    if folder.split('/')[7] == "CESM2-WACCM" and folder.split('/')[8] == "ssp585":
        files.append(folder + os.listdir(folder)[0])
    else:
        for x in os.listdir(folder):
            #hack to prevent post-2100 projections for ssp585 being included (these break the concatenation)
            if not "2300" in x:
                files.append(folder + x)
    cubes = iris.load(files, var_long)
    cube = esmvalcore.preprocessor.concatenate(cubes)
    
    
    ### convert back to cube and sort out coordinates + add some handy coords to use later
    #SI_cube = xrds_out[var].to_iris()  
    SI_cube = cube.copy()
    #SI_cube.coord(loncoord).rename('latitude')
    #SI_cube.coord(lon_coord).rename('longitude')
    SI_cube.coord('latitude').units = 'degrees'
    SI_cube.coord('longitude').units = 'degrees'
    area_weights = iris_utils.get_grid_areas(SI_cube, normalize=False)
    area_weights_normalised = iris.analysis.cartography.area_weights(SI_cube, normalize=False)
    grid_areas = iris.analysis.cartography.area_weights(SI_cube)

    iris_utils.addyearcoords(SI_cube)
    iris_utils.addmonthcoords(SI_cube)
    ens_coord = iris.coords.AuxCoord(folder.split('/')[9], long_name='Ensemble_member', units='no_unit')
    SI_cube.add_aux_coord(ens_coord)
    model_coord = iris.coords.AuxCoord(folder.split('/')[7], long_name='Model', units='no_unit')
    SI_cube.add_aux_coord(model_coord)
    exp_coord = iris.coords.AuxCoord(folder.split('/')[8], long_name='Experiment', units='no_unit')
    SI_cube.add_aux_coord(exp_coord)  
    
    
    ### calculate arctic temp, global temp, deltas, and AA_ratio
    No_arctic_cube = SI_cube.extract(iris.Constraint(latitude=lambda cell: -90 < cell < arctic_cut_off))
    grid_areas = iris.analysis.cartography.area_weights(No_arctic_cube)
    No_arctic_cube = No_arctic_cube.collapsed(['longitude','latitude'], iris.analysis.MEAN, weights=grid_areas)
    #Arctic_cube_70N = SI_cube.extract(iris.Constraint(latitude=lambda cell: 70 < cell < 90))
    #grid_areas = iris.analysis.cartography.area_weights(Arctic_cube_70N)
    #Arctic_cube_70N = Arctic_cube.collapsed(['longitude','latitude'], iris.analysis.MEAN, weights=grid_areas)
    grid_areas = iris.analysis.cartography.area_weights(SI_cube)
    Global_cube = SI_cube.collapsed(['longitude','latitude'], iris.analysis.MEAN, weights=grid_areas)
    
    
    ### calculate arctic temp, global temp, deltas, and AA_ratio
    
    
    
    ### create a DF with full monthly time resolved data
    df = pd.DataFrame(columns = ["Month", "Year", "No_arctic_temp", "Global_temp", "Units", 
                                                           "Experiment", "Model", "Ensemble_member"
                                ])
    l=0
    for x in No_arctic_cube.data:
        df.at[l, 'No_arctic_temp'] = x
        df.at[l, "Month"] = No_arctic_cube.coord('month').points[l]
        df.at[l, "Year"] = No_arctic_cube.coord('year').points[l]
        df.at[l, "Units"] = No_arctic_cube.units
        df.at[l, "Experiment"] = No_arctic_cube.coord('Experiment').points[0]
        df.at[l, "Model"] = No_arctic_cube.coord('Model').points[0]
        df.at[l, "Ensemble_member"] = No_arctic_cube.coord('Ensemble_member').points[0]
        l=l+1
    
    l=0
    for x in Global_cube.data:
        df.at[l, 'Global_temp'] = x
        l=l+1

    model = No_arctic_cube.coord('Model').points[0]
    experiment = No_arctic_cube.coord('Experiment').points[0]
    ensemble_member = No_arctic_cube.coord('Ensemble_member').points[0]
            
    df.to_csv('int_outputs/temperature/No_arctic_and_global_temp_monthly_{M}_{Exp}_{Ens}.csv'.format(M=model, Exp=experiment, Ens=ensemble_member))    
    #return df


exps = ["ssp370", "historical"]
dirs = []
var_path = "Amon/tas"
for experiment in exps:
    if experiment == "historical":
        exp_set = "CMIP"
    else:
        exp_set = "ScenarioMIP"
    for x in glob.glob('/badc/cmip6/data/CMIP6/{es}/*/*/{e}/*/{v}/*/latest/'.format(es=exp_set, e=experiment, v=var_path)):
        dirs.append(x)
        
#dirs = dirs[10:13]
    
for path in tqdm(dirs):
    try:
        temperatures(path)
        print("done: " + path)
    except:
        print("error: " + path)
