from datetime import datetime, timedelta

from urllib.parse import urljoin
from astropy.table import QTable
import astropy.units as u

def convert_decimal_day(decimal_day):

    dt = datetime(datetime.utcnow().year-1, 12, 31) + timedelta(decimal_day)
    if dt.microsecond >= 500000:
        dt += timedelta(seconds=1)
        dt = dt.replace(microsecond = 0)
    else:
        dt= dt.replace(microsecond = 0)
    return dt

def fetch_pwv(site, year=datetime.utcnow().year):

    url = "http://www.suominet.ucar.edu/data/staYrDayGlob/"

    gps_file = "{}_{:4d}global.plt".format(site, year)
    dl_link = urljoin(url, gps_file)
    mbar = u.def_unit(['millibar', 'millibars'], 1000.0*u.bar)
    units={ 'DayOfYear' : u.day,
            'PWV' : u.mm,
            'PWVerr' : u.mm,
            'ZenithDelay' : u.mm,
            'SurfacePressure' : mbar,
            'SurfaceTemp' : u.deg_C,
            'SurfaceRH' : u.percent
           }
    table = QTable.read(dl_link, format='ascii', names=["DayOfYear", "PWV", "PWVerr", "ZenithDelay", "SurfacePressure", "SurfaceTemp", "SurfaceRH"])
    for column in table.columns:
        table[column].unit = units[column]
    return table

def map_LCO_to_GPS_sites(site_code):
    """Maps LCO site codes (e.g. 'LSC') to SuomiNet sites (e.g. 'CTIO'). None
    is returned if the LCO site code was not found or there isn't a
    corresponding GPS site"""
    
    mapping = { 'LSC' : 'CTIO', 
                'CPT' : 'SUTM', 
                'TFN' : None, 
                'COJ' : None, 
                'ELP' : 'MDO1', 
                'OGG' : 'MAUI' }

    return mapping.get(site_code.upper(), None)
