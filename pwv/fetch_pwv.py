from datetime import datetime, timedelta
import os

try:
    from urllib.parse import urljoin
except ImportError:
    from urlparse import urljoin
import urllib.request
from subprocess import check_output

from astropy.table import QTable, Column
import astropy.units as u
import numpy as np
import netCDF4 as nc

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator, AutoMinorLocator

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
    mbar = u.def_unit(['millibar', 'millibars'], 1.0*u.hPa)
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

    # Create new datetime column
    dt = []
    for day in table['DayOfYear']:
        new_dt = convert_decimal_day(day.value)
        dt.append(new_dt)
    aa = Column(dt, name='UTC Datetime')
    table.add_column(aa)

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

def read_ascii(filepath):
    """Read single parameter ASCII format files extracted from the GrADS server e.g.
     https://goldsmr4.gesdisc.eosdis.nasa.gov/dods/M2T1NXSLV.ascii?tqv"""

    foo_fh = open(filepath, 'r')

    in_data = False
    data = {}
    for line in foo_fh:
        line = line.rstrip()
        if len(line) == 0:
            continue
        if line.count(',') == 1:
            if line[0] != '[':
                if in_data is True:
                    data[array_name] = array
                    array = []
                    array_name = line.split(',')[0]
                else:
                    in_data = True
                    array = []
                    array_name = line.split(',')[0]
            else:
                value = float(line.split(',')[1].strip())
                array.append(value)
        elif line.count(',') >= 1 and in_data is True:
            values = [float(x.strip()) for x in line.split(',')]
            data[array_name] = values
            in_data = False
    foo_fh.close()

    return data

def determine_cell(location, lats, longs):
    """Determines the latitude and longitude array indices corresponding to
    the latitude and longitude of the passed <location> (this should either
    be a `astropy.coordinates.earth.EarthLocation` or a dictionary with 'lat'
    and 'long' keys with the values in degrees"""

    long_idx = None
    lat_idx = None

    site_long = None
    try:
        site_long = location.lon.deg
    except AttributeError:
        site_long = location.get('long', None)

    site_lat = None
    try:
        site_lat = location.lat.deg
    except AttributeError:
        site_lat = location.get('lat', None)

    if site_long is not None and site_lat is not None:
        if lats.ndim == 1 and longs.ndim == 1:
            # 1d arrays
            long_idx = (np.abs(longs-site_long)).argmin()
            lat_idx  = (np.abs(lats-site_lat)).argmin()
        elif lats.ndim == 2 and longs.ndim == 2:
            # 2d arrays
            X = np.abs(longs-site_long)
            long_idx = np.where(X == X.min())
            Y = np.abs(lats-site_lat)
            lat_idx = np.where(Y == Y.min())
        else:
            print("Unknown number of dimensions: %d x %d" % (lats.ndim, longs.ndim))

    return lat_idx, long_idx

def determine_index(location, lon0=-180, lat0=-90, lonres=0.625, latres=0.5):
    long_idx = int((location.lon.degree - lon0)/lonres) + 1
    lat_idx  = int((location.lat.degree - lat0)/latres) + 1

    return lat_idx, long_idx

def extract_timeseries(data, long_idx, lat_idx):
    """Subset the passed MERRA-2 dataset in <data> for the
    data[t,lats,longs]"""

    return data[:,lat_idx,long_idx]

def find_modis_data(location, date=datetime.utcnow(), products=['MODATML2','MYDATML2'], collection='61', dbg=False):
    import modapsclient as m
    a = m.ModapsClient()

    start_date = date.date()
    end_date = start_date+timedelta(days=1)
    lat = location.lat.deg
    lon = location.lon.deg
    delta = 0.1
    north_val = max(lat+delta, lat-delta)
    west_val  = min(lon+delta, lon-delta)
    east_val  = max(lon+delta, lon-delta)
    south_val = min(lat+delta, lat-delta)
    dl_files = {}
    for product in products:
        if dbg: print(product)
        collections = a.getCollections(product)
        if collection not in collections:
            print("Did not find specified collection (%s) in %s".format(collection, collections))
            return None
        files = a.searchForFiles(product, start_date, end_date, north=north_val, west=west_val, east=east_val, south=south_val, collection=collection)
        if dbg: print(files)
        for hdf_file in files:
            url = a.getFileUrls(hdf_file)
            props = a.getFileProperties(hdf_file)
            if len(props) == 1:
                props = props[0]
            else:
                print("Unexpected number of file properties found")
                continue
            if props.get('online', 'false') == 'true':
               dl_files[hdf_file] = { 'url' : url,
                                      'checksum' : props.get('checksum', ''),
                                      'size' : props.get('fileSizeBytes', ''),
                                      'name' : props.get('fileName', '')
                                    }
            else:
                print("File id %d (%s) is not online".format(hdf_file, props.get('fileName', 'UNKNOWN')))
    return dl_files

def download_files(files, outpath):
    try:
        os.makedirs(outpath)
    except FileExistsError:
        pass

    num_dl = 0
    for file_id in files:
        url = files[file_id]['url'][0]
        dl_file = files[file_id]['name']
        file_name = os.path.join(outpath, os.path.basename(dl_file))
        if os.path.exists(file_name) is False:
            urllib.request.urlretrieve(url, file_name)
            num_dl += 1
            output = check_output(['cksum', file_name])
            chunks = output.decode('ascii', 'ignore').split()
            if chunks[0] == files[file_id]['checksum'] and chunks[1] == files[file_id]['size']:
                print("Downloaded {}".format(dl_file))
            else:
                print("Checksum/file size check failed for {}".format(dl_file))
        else:
            print("{} already downloaded".format(dl_file))
    return num_dl

def dataset_mapping(product):
    mappings = { 'MODATML2' : { 'PWV' : 'Precipitable_Water_Infrared_ClearSky',},
                 'MYD05_L2' : { 'PWV' : 'Water_Vapor_Infrared', },
                 'MOD05_L2' : { 'PWV' : 'Water_Vapor_Infrared', },
               }
    return mappings.get(product, {})

def read_modis_pwv(datafile, DATAFIELD_NAME = 'Water_Vapor_Infrared'):
    from pyhdf.SD import SD, SDC

    hdf = SD(datafile, SDC.READ)

    # Read dataset.
    data2D = hdf.select(DATAFIELD_NAME)
    data = data2D[:,:].astype(np.double)

    hdf_geo = SD(datafile, SDC.READ)

    # Read geolocation dataset from MOD03 product.
    lat = hdf_geo.select('Latitude')
    latitude = lat[:,:]
    lon = hdf_geo.select('Longitude')
    longitude = lon[:,:]

    # Retrieve attributes.
    attrs = data2D.attributes(full=1)
    lna=attrs["long_name"]
    long_name = lna[0]
    aoa=attrs["add_offset"]
    add_offset = aoa[0]
    fva=attrs["_FillValue"]
    _FillValue = fva[0]
    sfa=attrs["scale_factor"]
    scale_factor = sfa[0]
    vra=attrs["valid_range"]
    valid_min = vra[0][0]
    valid_max = vra[0][1]
    try:
        ua=attrs["unit"]
    except KeyError:
        ua=attrs["units"]
    units = ua[0]

    invalid = np.logical_or(data > valid_max,
                            data < valid_min)
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = (data - add_offset) * scale_factor
    data = np.ma.masked_array(data, np.isnan(data))

    return data, latitude, longitude

def plot_merra2_pwv(hdf_path, datafile):
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

def plot_pwv_timeseries(CTIO_table, times, CTIO_pwv, filename='CTIO_GPS_MERRA2_comp.png'):
    fig, ax = plt.subplots()
    ax.plot(CTIO_table['UTC Datetime'], CTIO_table['PWV'], label='CTIO GPS')
    ax.plot(times, CTIO_pwv, label='CTIO MERRA-2')
    ylims = ax.get_ylim()
    ax.set_ylim(0, ylims[1])
    ax.set_title('Precipitable Water Vapour')
    ax.set_ylabel('PWV (mm)')

    hours = mdates.HourLocator(range(0,25,3))
    hoursFmt = mdates.DateFormatter('%m-%d %H:%M')
    ax.xaxis.set_major_locator(hours)
    ax.xaxis.set_major_formatter(hoursFmt)
    hours = mdates.HourLocator(range(0,25,1))
    ax.xaxis.set_minor_locator(hours)
    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d %H:%M')

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    minorLocator = AutoMinorLocator()
    ax.yaxis.set_minor_locator(minorLocator)

    plt.legend()
    fig.autofmt_xdate()
    plt.savefig(filename, format='png')

    return filename
