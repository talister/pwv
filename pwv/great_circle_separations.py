# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 21:35:38 2018

@author: rstreet
"""

import numpy as np

def locate_point_on_map(latitudes,longitudes,map_data,location):
    """Function to identify the grid square containing the location given"""

    dsigma = calc_great_circle_distance(latitudes,longitudes,location)

    idx = np.where(dsigma == dsigma.min())

    return idx

def calc_great_circle_distance(latitudes,longitudes,location):
    """Function to calculate the Great Circle distance between location and
    array of positions in (latitudes, longitudes)

    Inputs:
        latitudes  array  1D or 2D array Earth latitude positions
        longitudes array  1D or 2D array Earth longitude positions
        location   EarthLocation  Single position on Earth

    Outputs:
        dsigma     array  1D or 2D array of separations of each input position
                            from location
    """

    dlong = np.abs(np.radians(longitudes) - location.lon.rad)

    dsigma_pre = np.sin(np.radians(latitudes))*np.sin(location.lat.rad) + \
                        (np.cos(np.radians(longitudes))*np.cos(location.lon.rad)*np.cos(dlong))
    dsigma = np.arccos(dsigma_pre)
    return dsigma
