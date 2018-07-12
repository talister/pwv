from datetime import datetime, timedelta
import os

from urllib.parse import urljoin
from astropy.table import QTable
import astropy.units as u
import netCDF4 as nc

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

def read_merra2(hdf_path, datafile, columns=['PS', 'T2M', 'QV2M', 'TO3', 'TQV']):

    # Open file for reading, map into memory
    data = nc.Dataset(os.path.join(hdf_path,datafile), mode='r', diskless=True)

    # Extract longitude and latitude grid
    lons = data.variables['lon'][:]
    lats = data.variables['lat'][:]
    table = { 'lons' : lons, 'lats' : lats }

    # Read columns
    for column in columns:
        coldata = data.variables[column][:,:,:]
        if column == 'TQV':
             # (convert total preciptable water vapour in kg m^-2 to mm by dividing be density of water)
            coldata = coldata/0.997
            table['PWV'] = coldata
        else:
            table[column] = coldata

    return table

def plot_merra2_pwv(hdf_path, datafile):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from mpl_toolkits.basemap import Basemap

    table = read_merra2(hdf_path, datafile, columns=['TQV'])

    lons = table['lons']
    lats = table['lats']
    pwv = table['PWV']

    if pwv.shape[0] > 1:
        # Not a monthly mean, extract only 1st time slice
        pwv = pwv[0,:,:]
    plt.subplot(1,1,1)
    map = Basemap(resolution='l', projection='eck4', lat_0=0, lon_0=0)
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = map(lon, lat)
    cs = map.pcolor(xi,yi,np.squeeze(pwv), vmin=np.min(pwv), vmax=np.max(pwv), cmap=cm.jet)
    cs.set_edgecolor('face')
    map.drawparallels(np.arange(-90., 90., 15.), labels=[1,0,0,0], fontsize=5)
    map.drawmeridians(np.arange(-180., 180., 30.), labels=[0,0,0,1], fontsize=4)

    # Add Coastlines, States, and Country Boundaries
    map.drawcoastlines()
    map.drawstates()
    map.drawcountries()

    cbar = map.colorbar(cs, location='bottom', pad="10%")
    cbar.set_label('mm')
    cbar.ax.tick_params(labelsize=10)

    date_chunk = datafile.split('.')[2]
    if date_chunk.isdigit() is False:
        date_chunk = ''
    plt.title('MERRA-2 PWV ({})'.format(date_chunk))
    filename = 'MERRA2_PWV_{:s}_TEST.pdf'.format(date_chunk)
    plt.savefig(filename, format='pdf', dpi=360)
    plt.close()

    return filename
