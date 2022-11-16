import numpy as np
import iris 
import pandas as pd
import iris.coord_categorisation as coord_cat
import cf_units
import os
import iris.coords as coords
import logging 

logger = logging.getLogger(__name__)

def make_cube_from_all_files_in_folder(folder):
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
    iris_utils.addyearcoords(cube)
    iris_utils.addmonthcoords(cube)
    ens_coord = iris.coords.AuxCoord(folder.split('/')[9], long_name='Ensemble_member', units='no_unit')
    cube.add_aux_coord(ens_coord)
    model_coord = iris.coords.AuxCoord(folder.split('/')[7], long_name='Model', units='no_unit')
    cube.add_aux_coord(model_coord)
    exp_coord = iris.coords.AuxCoord(folder.split('/')[8], long_name='Experiment', units='no_unit')
    cube.add_aux_coord(exp_coord)  
    return cube



def _fix_cube_attributes(cubes):
    """Unify attributes of different cubes to allow concatenation."""
    attributes = {}
    for cube in cubes:
        for (attr, val) in cube.attributes.items():
            if attr not in attributes:
                attributes[attr] = val
            else:
                if not np.array_equal(val, attributes[attr]):
                    attributes[attr] = '{};{}'.format(str(attributes[attr]),
                                                      str(val))
    for cube in cubes:
        cube.attributes = attributes
        
def _by_two_concatenation(cubes):
    """Perform a by-2 concatenation to avoid gaps."""
    concatenated = iris.cube.CubeList(cubes).concatenate()
    if len(concatenated) == 1:
        return concatenated[0]

    concatenated = _concatenate_overlapping_cubes(concatenated)
    if len(concatenated) == 2:
        _get_concatenation_error(concatenated)
    else:
        return concatenated[0]

def _get_concatenation_error(cubes):
    """Raise an error for concatenation."""
    # Concatenation not successful -> retrieve exact error message
    try:
        iris.cube.CubeList(cubes).concatenate_cube()
    except iris.exceptions.ConcatenateError as exc:
        msg = str(exc)
    logger.error('Can not concatenate cubes into a single one: %s', msg)
    logger.error('Resulting cubes:')
    for cube in cubes:
        logger.error(cube)
        time = cube.coord("time")
        logger.error('From %s to %s', time.cell(0), time.cell(-1))

    raise ValueError(f'Can not concatenate cubes: {msg}')

def _concatenate_overlapping_cubes(cubes):
    """Concatenate time-overlapping cubes (two cubes only)."""
    # we arrange [cube1, cube2] so that cube1.start <= cube2.start
    if cubes[0].coord('time').points[0] <= cubes[1].coord('time').points[0]:
        cubes = [cubes[0], cubes[1]]
        logger.debug(
            "Will attempt to concatenate cubes %s "
            "and %s in this order", cubes[0], cubes[1])
    else:
        cubes = [cubes[1], cubes[0]]
        logger.debug(
            "Will attempt to concatenate cubes %s "
            "and %s in this order", cubes[1], cubes[0])

    # get time end points
    time_1 = cubes[0].coord('time')
    time_2 = cubes[1].coord('time')
    if time_1.units != time_2.units:
        raise ValueError(
            f"Cubes\n{cubes[0]}\nand\n{cubes[1]}\ncan not be concatenated: "
            f"time units {time_1.units}, calendar {time_1.units.calendar} "
            f"and {time_2.units}, calendar {time_2.units.calendar} differ")
    data_start_1 = time_1.cell(0).point
    data_start_2 = time_2.cell(0).point
    data_end_1 = time_1.cell(-1).point
    data_end_2 = time_2.cell(-1).point

    # case 1: both cubes start at the same time -> return longer cube
    if data_start_1 == data_start_2:
        if data_end_1 <= data_end_2:
            logger.debug(
                "Both cubes start at the same time but cube %s "
                "ends before %s", cubes[0], cubes[1])
            logger.debug("Cube %s contains all needed data so using it fully",
                         cubes[1])
            cubes = [cubes[1]]
        else:
            logger.debug(
                "Both cubes start at the same time but cube %s "
                "ends before %s", cubes[1], cubes[0])
            logger.debug("Cube %s contains all needed data so using it fully",
                         cubes[0])
            cubes = [cubes[0]]

    # case 2: cube1 starts before cube2
    else:
        # find time overlap, if any
        start_overlap = next((time_1.units.num2date(t)
                              for t in time_1.points if t in time_2.points),
                             None)
        # case 2.0: no overlap (new iris implementation does allow
        # concatenation of cubes with no overlap)
        if not start_overlap:
            logger.debug(
                "Unable to concatenate non-overlapping cubes\n%s\nand\n%s"
                "separated in time.", cubes[0], cubes[1])
        # case 2.1: cube1 ends after cube2 -> return cube1
        elif data_end_1 > data_end_2:
            cubes = [cubes[0]]
            logger.debug("Using only data from %s", cubes[0])
        # case 2.2: cube1 ends before cube2 -> use full cube2 and shorten cube1
        else:
            logger.debug(
                "Extracting time slice between %s and %s from cube %s to use "
                "it for concatenation with cube %s", "-".join([
                    str(data_start_1.year),
                    str(data_start_1.month),
                    str(data_start_1.day)
                ]), "-".join([
                    str(start_overlap.year),
                    str(start_overlap.month),
                    str(start_overlap.day)
                ]), cubes[0], cubes[1])
            c1_delta = extract_time(cubes[0], data_start_1.year,
                                    data_start_1.month, data_start_1.day,
                                    start_overlap.year, start_overlap.month,
                                    start_overlap.day)
            # convert c1_delta scalar cube to vector cube, if needed
            if c1_delta.data.shape == ():
                c1_delta = iris.util.new_axis(c1_delta, scalar_coord="time")
            cubes = iris.cube.CubeList([c1_delta, cubes[1]])
            logger.debug("Attempting concatenatenation of %s with %s",
                         c1_delta, cubes[1])
            try:
                cubes = [iris.cube.CubeList(cubes).concatenate_cube()]
            except iris.exceptions.ConcatenateError as ex:
                logger.error('Can not concatenate cubes: %s', ex)
                logger.error('Cubes:')
                for cube in cubes:
                    logger.error(cube)
                raise ex

    return cubes
        
def concatenate(cubes):
    """Concatenate all cubes after fixing metadata."""
    if not cubes:
        return cubes
    if len(cubes) == 1:
        return cubes[0]

    _fix_cube_attributes(cubes)

    if len(cubes) > 1:
        # order cubes by first time point
        try:
            cubes = sorted(cubes, key=lambda c: c.coord("time").cell(0).point)
        except iris.exceptions.CoordinateNotFoundError as exc:
            msg = "One or more cubes {} are missing".format(cubes) + \
                  " time coordinate: {}".format(str(exc))
            raise ValueError(msg)

        # iteratively concatenate starting with first cube
        result = cubes[0]
        for cube in cubes[1:]:
            result = _by_two_concatenation([result, cube])

    _fix_aux_factories(result)

    return result        


def bbox_extract_2Dcoords(cube, bbox):
    """
    Extract a sub-set of a cube inside a lat range
    bbox=[lat_min lat_max].
    This is a work around too subset an iris cube that has
    2D lon, lat coords.
    
    """   
    minmax = lambda x: (np.min(x), np.max(x))
    lats = cube.coord('latitude').points
    
    inregion = np.logical_and(lats > bbox[0], lats < bbox[1])
    region_inds = np.where(inregion)
    imin, imax = minmax(region_inds[0])
    jmin, jmax = minmax(region_inds[1])
    return cube[..., imin:imax+1, jmin:jmax+1]




def ensemblecollapse(ans):
    if ans.coords('ensemble_member'):
        ens_mean = ans.collapsed(['ensemble_member'], iris.analysis.MEAN)
        return ens_mean      
                  

def cop_metadata_callback(cube, field, filename):
    """ A function which adds an "Experiment" coordinate which comes from the filename. """

    # Extract the experiment name (such as a1b or e1) from the filename (in this case it is just the parent folder's name)
    containing_folder = os.path.dirname(filename)
    experiment_label = os.path.basename(containing_folder)

    # Create a coordinate with the experiment label in it
    exp_coord = coords.AuxCoord(experiment_label, long_name='Experiment', units='no_unit')

    # and add it to the cube
    cube.add_aux_coord(exp_coord)

def subset(x, time_bounds, lat_bounds, lon_bounds, season):
        x = iris_utils.time_select(x, time_bounds)
        x = iris_utils.select_season(x, season)
        x = iris_utils.lat_lon_select(x, lat_bounds, lon_bounds)
        return x
    
def time_select(y, tbnds):
    tcoord = y.coord('time')
    tcoord.units = cf_units.Unit(tcoord.units.origin, calendar='gregorian')
    if not y.coords('year'):
        iris.coord_categorisation.add_year(y, 'time', name='year')
    x = y.extract(iris.Constraint(year=lambda cell: tbnds[0] < cell < tbnds[1]))
    return x

def lat_lon_select(y, lat_bnds, lon_bnds=None):
    """ this assumes contiguity over 0 longitude is desired, rather than contiguity at           180. Use a second constraint instead of intersection in the other case. """
    lats = iris.Constraint(latitude=lambda cell: lat_bnds[0] < cell < lat_bnds[1])
    x = y.extract(lats)
    if lon_bnds:
        x = x.intersection(longitude=(lon_bnds[0], lon_bnds[1]), ignore_bounds=True)
    return x

def select_season(ans, season):
    if not ans.coords('clim_season'):
        iris.coord_categorisation.add_season(ans, 'time', name='clim_season')  
    x = ans.extract(iris.Constraint(clim_season=season))
    return x


def set_eof_sign_convention(lat_lon_coord_with_positive_EOF_value, cube_list):
    lat_lon = [('latitude', lat_lon_coord_with_positive_EOF_value[0]), ('longitude', lat_lon_coord_with_positive_EOF_value[1])]
    for x in cube_list:
        if x.interpolate(lat_lon, iris.analysis.Linear()).data < 0:
            x = -x
    return cube_list
  
def seasonal_mean(ans): 
    if not ans.coords('clim_season'):
        iris.coord_categorisation.add_season(ans, 'time', name='clim_season')  
    mean=ans.collapsed(['clim_season'], iris.analysis.MEAN)
    return mean
    
def addmonthcoords(ans):
    if not ans.coords('month'):
        coord_cat.add_month(ans, 'time', name='month')  
        
def addyearcoords(ans):
    if not ans.coords('year'):
        coord_cat.add_year(ans, 'time', name='year')
             
    
def addseasoncoords(ans):
    if not ans.coords('clim_season'):
        coord_cat.add_season(ans, 'time', name='clim_season')  
        
def yearcollapse(ans): 
    mean=ans.collapsed(['year'], iris.analysis.MEAN)
    return mean             

def MAM(ans):
    MAM = ans.extract(iris.Constraint(clim_season='mam'))
    return MAM

def get_grid_areas(ans, normalize=False):
    if (not ans.coord('latitude').has_bounds()):
        ans.coord('latitude').guess_bounds()
    if (not ans.coord('longitude').has_bounds()):
        ans.coord('longitude').guess_bounds()
    grid_areas=iris.analysis.cartography.area_weights(ans, normalize=normalize)
    return grid_areas

def latlonmean(ans): 
    if (not ans.coord('latitude').has_bounds()):
        ans.coord('latitude').guess_bounds()
    if (not ans.coord('longitude').has_bounds()):
        ans.coord('longitude').guess_bounds()
    grid_areas=iris.analysis.cartography.area_weights(ans)
    mean=ans.collapsed(['longitude','latitude'], iris.analysis.MEAN, weights=grid_areas)
    return mean 

def latlonsum(ans):
    if (not ans.coord('latitude').has_bounds()):
        ans.coord('latitude').guess_bounds()
    if (not ans.coord('longitude').has_bounds()):
        ans.coord('longitude').guess_bounds()
    grid_areas=iris.analysis.cartography.area_weights(ans)
    mean=ans.collapsed(['longitude','latitude'], iris.analysis.sum, weights=grid_areas)
    return mean


def latmean(ans): 
    if (not ans.coord('latitude').has_bounds()):
        ans.coord('latitude').guess_bounds()
    grid_areas=iris.analysis.cartography.area_weights(ans)
    mean=ans.collapsed(['latitude'], iris.analysis.MEAN, weights=grid_areas)
    return mean 


def lonmean(ans): 
    mean=ans.collapsed(['longitude'], iris.analysis.MEAN)
    return mean 



def Station_NAO_index(cube, Iceland_station="Akureyri"):
    Ponta_Delgada = [37.7, 25.7]
    Akureyri_Iceland = [65.7, 18.1]
    Stykkisholmur_Iceland = [65.0, 22.8]

    South_lat_lon = [('latitude', Ponta_Delgada[0]), ('longitude', Ponta_Delgada[1])]
    if Iceland_station == "Akureyri":
        North_lat_lon = [('latitude', Akureyri_Iceland[0]), ('longitude', Akureyri_Iceland[1])]
    elif Iceland_station == "Stykkisholmur":
        North_lat_lon = [('latitude', Stykkisholmur_Iceland[0]), ('longitude', Stykkisholmur_Iceland[1])]
    else:
        print("error: please select either Akureyri or Stykkisholmur as Iceland location")

    NAO = cube.interpolate(South_lat_lon, iris.analysis.Linear()) - cube.interpolate(North_lat_lon, iris.analysis.Linear())
    return NAO


def Area_NAO_index_Stephenson(cube):
    southern_region_lons = [-90, 60]
    southern_region_lats = [20, 55]
    northern_region_lons = southern_region_lons
    northern_region_lats = [55, 90]
    
    northern_region_psl = lat_lon_select(cube, northern_region_lats, northern_region_lons)
    southern_region_psl = lat_lon_select(cube, southern_region_lats, southern_region_lons)
    
    northern_region_psl = latlonmean(northern_region_psl)
    southern_region_psl = latlonmean(southern_region_psl)
    
    NAO =  southern_region_psl - northern_region_psl
    return NAO

def Convert_Pr_units_from_flux_to_rate(pr_cube, conversion_factor=86400):
    pr_cube.rename('precipitation_rate')
    pr_cube.units = 'mm day-1'
    pr_cube = pr_cube*conversion_factor
    #86400 is water density times the second to day conversion...
    #...to go from kgm-2s-1 to mm day-1
    
    return pr_cube

def _fix_aux_factories(cube):
    """Fix :class:`iris.aux_factory.AuxCoordFactory` after concatenation.

    Necessary because of bug in :mod:`iris` (see issue #2478).
    """
    coord_names = [coord.name() for coord in cube.coords()]

    # Hybrid sigma pressure coordinate
    # TODO possibly add support for other hybrid coordinates
    if 'atmosphere_hybrid_sigma_pressure_coordinate' in coord_names:
        new_aux_factory = iris.aux_factory.HybridPressureFactory(
            delta=cube.coord(var_name='ap'),
            sigma=cube.coord(var_name='b'),
            surface_air_pressure=cube.coord(var_name='ps'),
        )
        for aux_factory in cube.aux_factories:
            if isinstance(aux_factory, iris.aux_factory.HybridPressureFactory):
                break
        else:
            cube.add_aux_factory(new_aux_factory)

    # Hybrid sigma height coordinate
    if 'atmosphere_hybrid_height_coordinate' in coord_names:
        new_aux_factory = iris.aux_factory.HybridHeightFactory(
            delta=cube.coord(var_name='lev'),
            sigma=cube.coord(var_name='b'),
            orography=cube.coord(var_name='orog'),
        )
        for aux_factory in cube.aux_factories:
            if isinstance(aux_factory, iris.aux_factory.HybridHeightFactory):
                break
        else:
            cube.add_aux_factory(new_aux_factory)

    # Atmosphere sigma coordinate
    if 'atmosphere_sigma_coordinate' in coord_names:
        new_aux_factory = iris.aux_factory.AtmosphereSigmaFactory(
            pressure_at_top=cube.coord(var_name='ptop'),
            sigma=cube.coord(var_name='lev'),
            surface_air_pressure=cube.coord(var_name='ps'),
        )
        for aux_factory in cube.aux_factories:
            if isinstance(aux_factory,
                          iris.aux_factory.AtmosphereSigmaFactory):
                break
        else:
            cube.add_aux_factory(new_aux_factory)