### contains functions to quickly generate lists of file paths to CMIP6 runs inside the JASMIN system
import glob
    
def get_all_G6_ssp585_ssp245_runs(var_path='Amon/tas', models=['IPSL-CM6A-LR', 'UKESM1-0-LL', 'MPI-ESM1-2-LR', 
                                                               'CESM2-WACCM', 'MPI-ESM1-2-HR'],
                                 quiet=False):
    dirs_g6sulfur = []
    dirs_g6solar = []
    dirs_ssp245 = []
    dirs_ssp585 = []

    #var_path =  'Amon/tas' 
    #var_path = 'SImon/siconc'

    for x in glob.glob('/badc/cmip6/data/CMIP6/GeoMIP/*/*/G6sulfur/*/{}/*/latest/'.format(var_path)):
        dirs_g6sulfur.append(x)

    #print(dirs_g6sulfur)

    for x in glob.glob('/badc/cmip6/data/CMIP6/GeoMIP/*/*/G6solar/*/{}/*/latest/'.format(var_path)):
        dirs_g6solar.append(x)

    dirs_geo = dirs_g6solar + dirs_g6sulfur


    ### need a path to csv to read in the names of lat lon coords for different models.
    ### THIS CSV IS A LIKELY BREAK POINT - if models used etc. are varied, these names will...
    ### ...need to be updated manually. 
    #lat_lon_names_df = pd.read_csv('OceanModel_latlon_names.csv')
    ### 
    ###

    if quiet == False:
        for model in models:
            print(model)

    for model in models:
        for x in glob.glob('/badc/cmip6/data/CMIP6/ScenarioMIP/*/{m}/ssp245/*/{v}/*/latest/'.format(m=model, v=var_path)):
            dirs_ssp245.append(x)

    #print(dirs_ssp245)

    for model in models:
        for x in glob.glob('/badc/cmip6/data/CMIP6/ScenarioMIP/*/{m}/ssp585/*/{v}/*/latest/'.format(m=model, v=var_path)):
            dirs_ssp585.append(x)

    dirs_nogeo = dirs_ssp245 + dirs_ssp585    
    dirs_all = dirs_geo + dirs_nogeo
    
    if quiet == False:
        print(dirs_all)
    
    return dirs_all, dirs_geo, dirs_nogeo


