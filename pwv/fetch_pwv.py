from datetime import datetime, timedelta, time
import os
from glob import glob

try:
    from urllib.parse import urljoin
except ImportError:
    from urlparse import urljoin
try:
    import urllib2 as ur
except ImportError:
    import urllib.request as ur
try:
    from cookielib import CookieJar
except ImportError:
    from http.cookiejar import CookieJar
from netrc import netrc
from subprocess import check_output

from astropy.table import QTable, Column
import astropy.units as u
from astropy.coordinates import EarthLocation
import numpy as np
import netCDF4 as nc

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from pwv.utils import determine_time_index, time_index_to_dt, make_bounding_box

def convert_decimal_day(decimal_day, year=datetime.utcnow().year):

    if type(year) == datetime:
        year = year.year
    dt = datetime(year-1, 12, 31, 0, 0, 0, 0) + timedelta(decimal_day)
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
        new_dt = convert_decimal_day(day.value, year)
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

def map_LCO_to_location(site_code):

    mapping = { 'LSC' : EarthLocation(lon=-70.806885, lat= -30.169165, height= 2218.45),
                'FTN' : EarthLocation(lon=-156.257029, lat=20.706657, height=3046.52),
                'FTS' : EarthLocation(lon=149.070277778, lat=-31.2731666667, height=1111.8),
                'ELP' : EarthLocation(lon=-104.02199,  lat=30.68049, height=2057.185),
                'CPT' : EarthLocation(lon=20.8101083333, lat=-32.3806083333, height=1807),
                'TFN' : EarthLocation(lon=-16.5117027778, lat=28.3003083333, height=2390.0),
                'COJ' : EarthLocation(lon=149.070705556, lat=-31.2729791667, height=1168.0)
              }
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

def earthdata_login(username=None, password=None):
    """Login into NASA Earthdata system with specified username and password.
    Adapted from https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
    Returns urllib opener..."""

    host = 'urs.earthdata.nasa.gov'
    if username is None or password is None:
        user_netrc = netrc()
        auth = user_netrc.authenticators(host)
        if auth is not None:
            username = auth[0]
            password = auth[2]
        else:
            print("Could not find authentication for '{}' in $HOME/.netrc")
            return None

    # Create a password manager to deal with the 401 reponse that is returned from
    # Earthdata Login
    password_manager = ur.HTTPPasswordMgrWithDefaultRealm()
    password_manager.add_password(None, "https://{}".format(host), username, password)
    cookie_jar = CookieJar()

    opener = ur.build_opener(ur.HTTPBasicAuthHandler(password_manager),
        #ur.HTTPHandler(debuglevel=1),    # Uncomment these two lines to see
        #ur.HTTPSHandler(debuglevel=1),   # details of the requests/responses
        ur.HTTPCookieProcessor(cookie_jar))

    return opener

def determine_gds_url(location, start=None, end=None, product='M2T1NXSLV', quantity='tqv'):
    """Fetches a time series of a quantity from MERRA-2 products at a single point"""

    GDS_URL = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/dods'
    start = start or datetime.utcnow()
    end = end or start + timedelta(days=1)

    # Determine spatial index
    iy, ix = determine_index(location)

    # Determine time index values
    start_index = determine_time_index(start)
    end_index = determine_time_index(end)

    url = "{:s}/{:s}.ascii?{}[{}:{}][{:03d}][{:03d}]".format(GDS_URL, product, \
        quantity, start_index, end_index, iy, ix)
    return url

def fetch_merra2_ascii_timeseries(location, filename=None, start=None, end=None, product='M2T1NXSLV', quantity='tqv'):
    """Fetches a time series from [start] to [end] of a quantity from MERRA-2
    products at a single point <location> and saves it in [filename].

    Returns a tuple of the HTML status code and the filename"""

    start = start or datetime.utcnow()
    end = end or start + timedelta(days=1)

    url = determine_gds_url(location, start, end, product, quantity)

    status_code = -1
    opener = earthdata_login()
    if opener:
        ur.install_opener(opener)
        request = ur.Request(url)
        r = ur.urlopen(request)

        date_fmt = "%Y%m%dT%H:%M:%S"
        if start.time() == time(0,0) and end.time() == time(0,0):
            # No time part, in start or end, remove from filename
            date_fmt = "%Y%m%d"

        filename = filename or "{}_{}_{}-{}.asc".format(product, quantity, start.strftime(date_fmt), end.strftime(date_fmt))
        with open(filename, 'wb') as f:
            f.write(r.read())
        status_code = r.status

    return status_code, filename

def read_ascii(filepath):
    """Read single parameter ASCII format files extracted from the GrADS server e.g.
     https://goldsmr4.gesdisc.eosdis.nasa.gov/dods/M2T1NXSLV.ascii?tqv 
     or from the output of fetch_merra2_ascii_timeseries()"""

    with open(filepath, 'r') as foo_fh:

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

    if data.get('time', None) is not None:
        times = time_index_to_dt(data['time'])
        data['datetime'] = times

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

def find_modis_data(location, date=datetime.utcnow(), ndays=1, products=['MODATML2','MYDATML2'], collection='61', dbg=False):
    """Finds near-realtime MODIS products from the Aqua/Terra satellites for a
    specific location <location>, starting at a specifc ([date]; defaults to
    'now') for a number of days ([ndays]; defaults to 1). By default the MODIS
    Joint Atmosphere products, MODATML2/MYDATML2 in Collection 61:
    https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/l2-joint-atmosphere/MYDATML2/
    are searched but this can be overridden by passing in [products] and
    [collection]. 
    A dictionary with the name of HDF file products as keys is returned; the
    values for each HDF file is a dictionary of the following form:
        url: URL link to retrieve the HDF file,
        checksum: CRC checksum of the file (`cksum` can be used to generate
                  the checksum and size,
        size: size of the HDF file in bytes,
        name: name of the HDF file,
    """
    import modapsclient as m
    import logging
    logging.basicConfig(level=logging.WARNING)
    a = m.ModapsClient()

    start_date = date.date()
    end_date = start_date+timedelta(days=ndays)

    west_val, south_val, east_val, north_val = make_bounding_box(location)

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
            ur.urlretrieve(url, file_name)
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
    """Define mappings for particular MODIS products <product> to dataset names.
    If a single dataset name is given, this uses the standard latitude, longitude
    grid; otherwise a tuple of the latitude, longitude grid names to use is given"""

    joint_atm = { 'PWV' : 'Precipitable_Water_Infrared_ClearSky',
                  'O3'  : 'Total_Ozone',
                  'AOD' : ('AOD_550_Dark_Target_Deep_Blue_Combined', ('Longitude_10km', 'Latitude_10km')),
                  'AAE' : ('Aerosol_Angstrom_Exponent_Ocean', ('Longitude_10km', 'Latitude_10km')),
                }
    mappings = { 'MODATML2' : joint_atm,
                 'MYDATML2' : joint_atm,
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

    # Read geolocation dataset.
    lat = hdf.select('Latitude')
    latitude = lat[:,:]
    lon = hdf.select('Longitude')
    longitude = lon[:,:]

    # Retrieve attributes for dataset
    attrs = data2D.attributes(full=1)
    lna = attrs["long_name"]
    long_name = lna[0]
    aoa = attrs["add_offset"]
    add_offset = aoa[0]
    fva = attrs["_FillValue"]
    _FillValue = fva[0]
    sfa = attrs["scale_factor"]
    scale_factor = sfa[0]
    vra = attrs["valid_range"]
    valid_min = vra[0][0]
    valid_max = vra[0][1]
    try:
        ua = attrs["unit"]
    except KeyError:
        ua = attrs["units"]
    units = ua[0]

    # Retrieve attributes for lat/lon, correct for scale factor and offset
    # if needed
    attrs = lat.attributes(full=1)
    aoa = attrs.get("add_offset", [0.0])
    lat_add_offset = aoa[0]
    sfa = attrs.get("scale_factor", [1.0])
    lat_scale_factor = sfa[0]

    attrs = lon.attributes(full=1)
    aoa = attrs.get("add_offset", [0.0])
    lon_add_offset = aoa[0]
    sfa = attrs.get("scale_factor", [1.0])
    lon_scale_factor = sfa[0]

    # NetCDF does the following automatically but PyHDF doesn't
    latitude = (latitude - lat_add_offset) * lat_scale_factor
    longitude = (longitude - lon_add_offset) * lon_scale_factor

    invalid = np.logical_or(data > valid_max,
                            data < valid_min)
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = (data - add_offset) * scale_factor
    data = np.ma.masked_array(data, np.isnan(data))

    return data, latitude, longitude, units, long_name

def split_MODIS_filename(filename):
    """Extracts the product, datetime and collection from a MODIS filename"""

    filename = os.path.basename(filename)
    chunks = filename.split('.')
    dt = datetime.strptime(chunks[1]+chunks[2], 'A%Y%j%H%M')
    return chunks[0], dt, str(int(chunks[3]))

def extract_MODIS_pwv_timeseries(hdf_path, location):

    pwv_values = []
    times = []
    for hdf_file in glob(os.path.join(hdf_path, '*.hdf')):
        product, date, collection = split_MODIS_filename(hdf_file)
        mapping = dataset_mapping(product)
        if len(mapping) == 0:
            print("Unable to find product mapping for {}".format(product))
            continue
        dataset_name = mapping.get('PWV', None)
        if dataset_name is not None:
            data, latitude, longitude, name, units = read_modis_pwv(hdf_file, dataset_name)
            lat_idx, long_idx = determine_cell(location, latitude, longitude)
            print(os.path.basename(hdf_file), lat_idx, long_idx, latitude.shape, longitude.shape, data.shape)
            pwv_value = data[lat_idx[0], long_idx[1]]
            times.append(date)
            pwv_values.append(pwv_value[0])
            print(pwv_value[0])
        else:
            print("Unable to find dataset mapping for PWV in {}".format(mapping.keys()))
    return times, pwv_values

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

def plot_modis(basename, location, data, latitude, longitude, units, long_name):

    from mpl_toolkits.basemap import Basemap

    delta = 15.0
    lat = location.lat.deg
    lon = location.lon.deg
    north_val = max(lat+delta, lat-delta)
    west_val  = min(lon+delta, lon-delta)
    east_val  = max(lon+delta, lon-delta)
    south_val = min(lat+delta, lat-delta)
    plt.clf()
    m = Basemap(projection='stere', resolution='l',
                    llcrnrlon=west_val,llcrnrlat=south_val,urcrnrlon=east_val,urcrnrlat=north_val,lat_0=lat,lon_0=lon)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 90., 10.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 15), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True)

    cb = m.colorbar()
    cb.set_label(units)
    basename = os.path.basename(basename)
    plt.title('{0}\n{1}\n'.format(basename, long_name))
    fig = plt.gcf()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)
    return pngfile
