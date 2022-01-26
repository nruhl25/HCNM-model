# Author: Nathaniel Ruhl
# Date: 1/26/2022
# This is a script that I sent to a collaborator to demonstrate my method for defining the latitude and longitude
# that are used as inputs to extract a density profile from the MSIS atmospheric density model.

# Script to identify the latitude and longitude at the "graze point", the point on
# NICER's line of sight that grazes Earth's surface

import numpy as np
import numbers
from numpy.linalg import norm
import datetime
import pymap3d as pm
# Pymap3d: https://geospace-code.github.io/pymap3d/ecef.html


# The fields below are calculated in another script, but I've copied them here
# starECI is the unit direction vector to the source in ECI coords
# t0_mkf is a time in seconds from Jan 1, 2014
# r0_mkf is  the position vector of the ISS (km) at t0 (ellipsoid denotes that the WGS84 ellipsoid is used)
feb3_dict = {} # V4641 Sgr
crab03_dict = {}  # Crab Nebula

feb3_dict['starECI'] = np.array([0.07619743, -0.90006328, -0.42904549])
feb3_dict['t0_mkf_ellipsoid'] = 192224367.25
feb3_dict['r0_mkf_ellipsoid'] = np.array([-4512.39672815, 3844.34234759, -3326.13647349])

crab03_dict['r0_mkf_ellipsoid'] = np.array([4154.73325568, -897.93001878, -5310.1136966])
crab03_dict['starECI'] = np.array([0.10280821, 0.92137081, 0.37484171])
crab03_dict['t0_mkf_ellipsoid'] = 240165079.01


# This function converts seconds from Jan 1 2014 into a datetime in UTC
def convert_time(time):
    timezero = datetime.datetime(year=2014, month=1,
                                 day=1, hour=0, minute=0, second=0)
    if isinstance(time, numbers.Real):
        new_time = timezero+datetime.timedelta(seconds=time)
        return new_time
    elif isinstance(time, list) or isinstance(time, np.ndarray):
        new_time_list = []
        for index, t in enumerate(time):
            new_time = timezero + datetime.timedelta(seconds=t)
            new_time_list.append(new_time)
        return np.array(new_time_list)
    else:
        raise RuntimeError('time must be MET number or array/list of times')

# This function takes in an array of eci points (in km) and returns latitude, longitude, and altitude using pymap3d
#  note that time seconds is a single time, float number
def eci2geodetic_pymap(los_point_array, time_seconds):
    if los_point_array.ndim == 1:
        los_point_m = los_point_array * 1000   # convert to [m]
        lat, lon, alt = pm.eci2geodetic(
            los_point_m[0], los_point_m[1], los_point_m[2], convert_time(time_seconds))
        # return altitude in km, lat and lon in degrees
        return lat, lon, alt / 1000
    else:
        # multiple los_points, time_seconds is either a number or array
        los_point_array_m = los_point_array * 1000   # convert to [m]
        lats, lons, alts = pm.eci2geodetic(los_point_array_m[:, 0], los_point_array_m[:, 1], los_point_array_m[:, 2], convert_time(time_seconds))
        return lats, lons, alts / 1000


# Line of sight from the predicted satellite position r(t)
def los_line(r0_guess, s_vec, n_list):
    if isinstance(n_list, numbers.Real):
        # n_list is not a list, but a single number
        n = n_list
        return r0_guess + n * s_vec
    else:
        n_column_vec = n_list.reshape((len(n_list), 1))
        starArray = np.ones((len(n_list), 3)) * s_vec
        return r0_guess + n_column_vec * starArray

def find_lat_lon_graze_point(obs_dict):
    s = obs_dict['starECI']
    t0 = obs_dict['t0_mkf_ellipsoid']
    r0_vec = obs_dict['r0_mkf_ellipsoid']  # not a unit vector

    # Define LOS, find closest approach to earth
    # Note that we're defining lats, lons, alts for a single LOS
    n_step_size = 0.5  # km step along the LOS
    n_list = np.arange(0, 4000, n_step_size)
    los_points = los_line(r0_vec, s, n_list)
    lats, lons, alts = eci2geodetic_pymap(los_points, t0)
    graze_index = np.argmin(alts)
    return lats[graze_index], lons[graze_index], alts[graze_index]

def main():

    crab_lat, crab_lon, crab_alt = find_lat_lon_graze_point(crab03_dict)
    print(f'crab_lat={crab_lat}')
    print(f'crab_lon={crab_lon}')
    print(f'crab_alt={crab_alt}')
    print('--------------------')
    v4641_lat, v4641_lon, v4641_alt = find_lat_lon_graze_point(feb3_dict)
    print(f'v4641_lat={v4641_lat}')
    print(f'v4641_lon={v4641_lon}')
    print(f'v4641_alt={v4641_alt}')
    return 0

if __name__ == '__main__':
    main()
