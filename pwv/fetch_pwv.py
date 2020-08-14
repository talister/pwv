from datetime import datetime, timedelta, time
import os
from glob import glob
from operator import itemgetter

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
import requests
import xml.etree.ElementTree as etree
from subprocess import check_output

from astropy.table import QTable, Column, join
from astropy.io.ascii import InconsistentTableError
import astropy.units as u
from astropy.coordinates import EarthLocation
import numpy as np
import netCDF4 as nc
from pyhdf.SD import SD, SDC
from pydap.client import open_url
from pydap.cas.urs import setup_session

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from pwv.utils import convert_decimal_day, determine_time_index, time_index_to_dt, make_bounding_box, compute_pwv
from pwv.telemetry import map_quantity_to_LCO_datum, query_LCO_telemetry, interpolate_LCO_telemetry

def fetch_pwv(site_code, start=None, end=None):
    """Download the GPS-based precipitable water vapor for the specified LCO
    <site_code> (e.g. 'FTN', 'COJ') between the specified [start] and [end]
    values (defaults to Jan 1 of current year and midnight on the current date
    if not specified).
    If the PWV is bad (mean < 0), it will attempt to repopulate it from the
    Total Zenith Distance by obtaining the pressure and temperature from telemetry
    and computing the Zenith Hydrostatic Delay.
    Returns an AstroPy Q(uantity)Table with columns:
        DayOfYear: decimal day of the year
        PWV: precipitable water vapour (in mm; -9.9 for missing value)
        PWVerr: error on precipitable water vapour (in mm)
        TotalZenithDelay: total zenith delay (in mm)
        SurfacePressure: surface pressure (in millibars)
        SurfaceTemp: surface temperature (in degrees C)
        SurfaceRH: surface relative humidty (%)
        UTC Datetime: day of the year converted to a datetime (computed, not in original)
    Masking of bad/missing values is not done.
    Reference: https://www.suominet.ucar.edu/data.html"""
    table = None
    start = start or datetime(datetime.utcnow().year, 1, 1)
    end = end or datetime.utcnow().replace(hour=0, minute=0, second=0, microsecond=0)

    GPS_site_code = map_obs_to_GPS_sites(site_code)
    if GPS_site_code:
        print("Fetching for", GPS_site_code)
        if start.year != end.year:
            print("Need to have start and end in the same year")
        if GPS_site_code != 'ORM_PWV':
            # Query SuomiNet
            GPS_table = fetch_GPS_pwv(GPS_site_code, start.year)
        else:
            # Query ORM Sky Quality Group
            GPS_table = fetch_ORM_pwv(None)

        # Check if the mean of the PWV is <0, indicating it's mostly bad values (-9.9)
        if GPS_table['PWV'].mean() < 0.0:
            print("PWV estimate bad, trying to repopulate")
            weather_table = fetch_LCO_weather(site_code, start, end, 300)
            location = map_LCO_to_location(site_code)
            if location and 'pressure' in weather_table.colnames and 'temperature' in weather_table.colnames:
                combined_table = join(GPS_table, weather_table)
                table = populate_PWV_column(location, combined_table)
            else:
                if not location:
                    print("Couldn't determine location for site code {}".format(site_code))
                else:
                    print("Couldn't find weather data for site code {}".format(site_code))
                table = GPS_table
        else:
            table = GPS_table
    else:
        print("No site code found")
    return table

def fetch_GPS_pwv(site, year=datetime.utcnow().year):
    """Download the SuomiNet GPS-based precipitable water vapor for the specified <site>
    and the particular [year] (defaults to current one if not specified).
    Returns an AstroPy Q(uantity)Table with columns:
        DayOfYear: decimal day of the year
        PWV: precipitable water vapour (in mm; -9.9 for missing value)
        PWVerr: error on precipitable water vapour (in mm)
        TotalZenithDelay: total zenith delay (in mm)
        SurfacePressure: surface pressure (in millibars)
        SurfaceTemp: surface temperature (in degrees C)
        SurfaceRH: surface relative humidty (%)
        UTC Datetime: day of the year converted to a datetime (computed, not in original)
    Masking of bad/missing values is not done.
    Reference: https://www.suominet.ucar.edu/data.html"""

    url = "https://www.suominet.ucar.edu/data/staYrDayGlob/"

    gps_file = "{}_{:4d}global.plt".format(site, year)
    dl_link = urljoin(url, gps_file)
    mbar = u.def_unit(['millibar', 'millibars'], 1.0*u.hPa)
    units={ 'DayOfYear' : u.day,
            'PWV' : u.mm,
            'PWVerr' : u.mm,
            'TotalZenithDelay' : u.mm,
            'SurfacePressure' : mbar,
            'SurfaceTemp' : u.deg_C,
            'SurfaceRH' : u.percent
           }

    try:
        table = QTable.read(dl_link, format='ascii', names=["DayOfYear", "PWV", "PWVerr", "TotalZenithDelay", "SurfacePressure", "SurfaceTemp", "SurfaceRH"])
    except InconsistentTableError:
        table = QTable.read(dl_link, format='ascii', names=["DayOfYear", "PWV", "PWVerr", "TotalZenithDelay", "SurfacePressure", "SurfaceTemp", "SurfaceRH", "Unknown1", "Unknown2", "Unknown3"])
    except ur.HTTPError:
        # CONUS not global site
        url = "https://www.suominet.ucar.edu/data/staYrDay/"
        gps_file = "{}pp_{:4d}.plt".format(site, year)
        dl_link = urljoin(url, gps_file)
        table = QTable.read(dl_link, format='ascii', names=["DayOfYear", "PWV", "PWVerr", "TotalZenithDelay", "SurfacePressure", "SurfaceTemp", "SurfaceRH", "Unknown1", "Unknown2", "Unknown3"])
    table_columns = table.columns.copy()
    for column in table_columns:
        if "Unknown" in column:
            table.remove_column(column)
        else:
            table[column].unit = units[column]

    # Create new datetime column
    dt = np.zeros(len(table['DayOfYear']), dtype='datetime64[s]')
    for i, day in enumerate(table['DayOfYear']):
        new_dt = convert_decimal_day(day.value, year)
        dt[i] = new_dt
    aa = Column(dt, name='UTC Datetime')
    table.add_column(aa)

    return table

def fetch_ORM_pwv(url_or_datafile=None):
    """Fetch the precipitable water vapor from the Roque de los Muchachos
    Observatory (ORM) (if <url_or_datafile> is None) or from that file if
    not None.
    (In principle, different amounts of data can be downloaded but this doesn't
    seem to work for num_days != 30)
    Returns an AstroPy Q(uantity)Table with columns:
        UTC Datetime: UTC datetime (computed from 'DATE' and 'UT' columns, not in original)
        PWV: precipitable water vapour (in mm; -9.9 for missing value)
        PWVerr: error on precipitable water vapour (in mm)
        TotalZenithDelay: total zenith delay (in mm)
        SurfacePressure: surface pressure (in hPa)
        SurfaceTemp: surface temperature (in degrees C)
        DataSource: Source of the meteorological data (model or sensor)
    Masking of bad/missing values is not done.
    Reference: https:/vivaldi.ll.iac.es/proyecto/site-testing/index.php?option=com_wrapper&Itemid=148
    """

    if url_or_datafile is not None:
        with open(url_or_datafile, 'rb') as f:
            data = f.read()
    else:
        iac_url = 'http://gaulli.ll.iac.es/data_file/download/'
        params = { "num_days" : 30, "type" : "post" }
        resp = requests.get(iac_url, params)
        if resp.status_code in [200, 201]:
            data = resp.content
    table_lines = data.decode('latin1').split('\n')

    units={ 'UTC Datetime' : None,
            'PWV' : u.mm,
            'PWVerr' : u.mm,
            'TotalZenithDelay' : u.mm,
            'SurfacePressure' : u.hPa,
            'SurfaceTemp' : u.deg_C,
            'DataSource' : None
           }
    table = QTable.read(table_lines, format="csv",header_start=7, data_start=8, delimiter=",", names=["UTC Date", "UTC Time", "PWV", "PWVerr", "TotalZenithDelay", "SurfacePressure", "SurfaceTemp", "DataSource"])
    # Create new datetime column
    dt = np.zeros(len(table['UTC Date']), dtype='datetime64[s]')
    for i, day in enumerate(table['UTC Date']):
        time = table['UTC Time'][i]
        new_dt = datetime.strptime(day+time, "%Y/%m/%d%H:%M")
        dt[i] = new_dt
    aa = Column(dt, name='UTC Datetime')
    table.add_column(aa, index=0)
    table.remove_columns(['UTC Date', 'UTC Time'])

    # Add units
    for column in table.columns:
        table[column].unit = units[column]

    return table

def map_obs_to_GPS_sites(site_code):
    """Maps LCO or other site codes (e.g. 'LSC' or 'NOT') to SuomiNet sites
    (e.g. 'CTIO') or the ORM PWV source on La Palma.
    None is returned if the LCO site code was not found or there isn't a
    corresponding GPS site"""

    mapping = { 'LSC' : 'CTIO',
                'CPT' : 'SUTM',
                'TFN' : 'ORM_PWV',
                'COJ' : None,
                'ELP' : 'MDO1',
                'OGG' : 'MAUI',
                'NOT' : 'ORM_PWV',
                'TNG' : 'ORM_PWV',
                'WHT' : 'ORM_PWV' }

    return mapping.get(site_code.upper(), None)

def map_LCO_to_location(site_code):
    """Map LCO site codes/telescope names (https://lco.global/observatory/sites/) to
    AstroPy EarthLocations.
    Returns None if no match found."""

    mapping = { 'LSC' : EarthLocation(lon=-70.806885, lat= -30.169165, height= 2218.45),
                'FTN' : EarthLocation(lon=-156.257029, lat=20.706657, height=3046.52),
                'OGG' : EarthLocation(lon=-156.257029, lat=20.706657, height=3046.52),
                'FTS' : EarthLocation(lon=149.070277778, lat=-31.2731666667, height=1111.8),
                'ELP' : EarthLocation(lon=-104.02199,  lat=30.68049, height=2057.185),
                'CPT' : EarthLocation(lon=20.8101083333, lat=-32.3806083333, height=1807),
                'TFN' : EarthLocation(lon=-16.5117027778, lat=28.3003083333, height=2390.0),
                'COJ' : EarthLocation(lon=149.070705556, lat=-31.2729791667, height=1168.0),
                'SBA' : EarthLocation(lon=-119.862856667, lat=34.4325633333, height=15.9)
              }
    return mapping.get(site_code.upper(), None)

def read_merra2(hdf_path, datafile, columns=['PS', 'T2M', 'QV2M', 'TO3', 'TQV']):
    """Read a MERRA-2 NetCDF file <datafile> located in <hdf_path>, extracting
    the quantities in [columns] (defaults to:
    PS: surface pressure,
    T2M: temperature at 2m level,
    QV2M: specific humidity at 2m level,
    TO3: total ozone,
    TQV: total precipitable water vapor (converted from kg m^-2 to mm by diving by 0.997),
    Returns a dictionary of the columns along with the latitude (`lats`) and
    longitude (`lons`)"""

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

def get_netrc_credentials(host='urs.earthdata.nasa.gov'):
    """Obtains the username and password for logging into the specified [host]
    from the ~/.netrc file"""

    username = None
    password = None

    user_netrc = netrc()
    auth = user_netrc.authenticators(host)
    if auth is not None:
        username = auth[0]
        password = auth[2]
    else:
        print("Could not find authentication for '{}' in $HOME/.netrc")

    return username, password

def earthdata_login(username=None, password=None):
    """Login into NASA Earthdata system with specified username and password.
    Adapted from https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
    Returns urllib opener..."""

    host = 'urs.earthdata.nasa.gov'
    if username is None or password is None:
        username, password = get_netrc_credentials(host)

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
    """Build a GDS URL to enable fetching of a time series of a [quantity] from
    MERRA-2 products ([product]) at a single point"""

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
    quantity = quantity.lower()

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

def read_ascii(filepath, dbg=False):
    """Read single parameter ASCII format files extracted from the GrADS server e.g.
     https://goldsmr4.gesdisc.eosdis.nasa.gov/dods/M2T1NXSLV.ascii?tqv
     or from the output of fetch_merra2_ascii_timeseries()"""

    with open(filepath, 'r') as foo_fh:

        first_line = foo_fh.readline()
        if 'Dataset:' in first_line.strip():
            # JD for 2002-01-01, starting point of aggregation values
            t0 = 2452275.5
            # OpenDAP format file
            data = {}
            for line in foo_fh:
                line = line.rstrip()
                if len(line) == 0:
                    continue
                if dbg: print(line)
                if '[' not in line:
                    chunks = line.split(',')
                    if len(chunks) == 2:
                        quantity = chunks[0]
                        value = chunks[1]
                        try:
                            value = float(value)
                        except ValueError:
                            pass
                        value = [value,]
                    else:
                        quantity = chunks[0]
                        value = [float(x.strip()) for x in chunks[1:]]
                    data[quantity] = value
                else:
                    chunks = line.split(',')
                    array_name = chunks[0]
                    value = chunks[1]
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                    index = array_name.find('[')
                    if index > 0:
                        array_name = array_name[0:index]
                    if array_name not in data:
                        # New quantity
                        data[array_name] = []
                    data[array_name].append(value)
        else:
            # GDS format file
            t0=1721423.5
            status = foo_fh.seek(0)

            in_data = False
            data = {}
            for line in foo_fh:
                line = line.rstrip()
                if len(line) == 0:
                    continue
                if dbg: print(line, in_data)
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
                            if dbg: print("Created array_name", array_name)
                    else:
                        value = float(line.split(',')[1].strip())
                        array.append(value)
                elif line.count(',') >= 1 and in_data is True:
                    values = [float(x.strip()) for x in line.split(',')]
                    data[array_name] = values
                    in_data = False
        foo_fh.close()

    if data.get('time', None) is not None:
        times = np.array(time_index_to_dt(data['time'], t0), dtype='datetime64[s]')
        data['datetime'] = times

    table = QTable()
    for key in data.keys():
        if dbg: print(key,len(data[key]))
        if type(data[key]) != float and len(data[key]) > 1:
            # Skip 1D values such as latitude or longitude
            aa = Column(data[key], name=key)
            table.add_column(aa)
    return table

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
            coordinates = np.unravel_index((np.abs(lats - site_lat) + np.abs(longs - site_long)).argmin(), lats.shape)
            lat_idx = coordinates[0]
            long_idx = coordinates[1]
        else:
            print("Unknown number of dimensions: %d x %d" % (lats.ndim, longs.ndim))

    return lat_idx, long_idx

def determine_index(location, lon0=-180, lat0=-90, lonres=0.625, latres=0.5):
    if lon0 > 0:
        long_idx = int((lon0 - location.lon.degree)/lonres) + 1
    else:
        long_idx = int((location.lon.degree - lon0)/lonres) + 1
    if lat0 > 0:
        lat_idx  = int((lat0- location.lat.degree)/latres) + 1
    else:
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
    """Downloads a list of HDF files <files> (produced by `find_modis_data()`) to
    the specified <outpath> (which is created if necessary). The file size and
    checksum (calculated by `cksum`) is compared with the expected values in
    the dictionary of the original file. Existing files are skipped and not
    re-downloaded.
    The number of files successfully downloaded is returned."""

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
    """Define mappings for particular EOS products <product> to dataset names.
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
                 # "Realtime" products
                 'O3_RT'    : { 'server' : 'https://acdisc.gesdisc.eosdis.nasa.gov/',
                                'level'  : 'Aqua_AIRS_Level3',
                                'product': 'AIRS3STD.006',
                                'O3'     : 'TotO3_D' },
                 'PWV_RT'   : { 'server' : 'https://acdisc.gesdisc.eosdis.nasa.gov/',
                                'level'  : 'Aqua_AIRS_Level3',
                                'product': 'AIRS3STD.006',
                                'PWV'    : 'TotH2OVap_D' },
               }
    return mappings.get(product, {})

def read_modis_pwv(datafile, DATAFIELD_NAME = 'Water_Vapor_Infrared'):
    """Open the Aqua/Terra MODIS HDF4-EOS file <datafile> and extract the
    data for the specified data field ([DATAFIELD_NAME]).
    This is returned in <data> along with the latitude and longitude
    geospatial arrays, the units and data field name"""

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

def extract_MODIS_pwv_timeseries(hdf_path, location, dbg=False):
    """Loop through a series of Aqua/Terra MODIS HDF files in <hdf_path> and
    extract the PWV at specified EatthLocation <location>"""

    pwv_values = []
    times = []
    for hdf_file in glob(os.path.join(hdf_path, 'M?D*.hdf')):
        product, date, collection = split_MODIS_filename(hdf_file)
        mapping = dataset_mapping(product)
        if len(mapping) == 0:
            print("Unable to find product mapping for {}".format(product))
            continue
        dataset_name = mapping.get('PWV', None)
        if dataset_name is not None:
            data, latitude, longitude, units, name = read_modis_pwv(hdf_file, dataset_name)
            lat_idx, long_idx = determine_cell(location, latitude, longitude)
            if dbg: print(os.path.basename(hdf_file), lat_idx, long_idx, latitude.shape, longitude.shape, data.shape)
            pwv_value = data[lat_idx, long_idx]
            times.append(date)
            pwv_values.append(pwv_value)
            if dbg: print(pwv_value)
        else:
            print("Unable to find dataset mapping for PWV in {}".format(mapping.keys()))
    if len(times) > 0 and len(pwv_values) > 0:
        # Sort the datetimes and PWV values by the time of the datasets
        sorted_times, sorted_pwv_values = [list(x) for x in zip(*sorted(zip(times, pwv_values), key=itemgetter(0)))]
        times = sorted_times
        pwv_values = sorted_pwv_values

    return times, pwv_values

def determine_opendap_agg_base_url(day, opendap_server, process_level, product):
    """Return the base URL for the AIRS Aggregated Data Service in OPeNDAP
    This takes <day> (a `datetime`), the <opendap_server> (e.g. https://acdisc.gesdisc.eosdis.nasa.gov/),
    the <process_level> (e.g. 'Aqua_AIRS_Level3') and the product (e.g. 'AIRS3STD.006')
    """

    pieces = [ 'opendap', 'ncml/aggregation', product, product+"_Aggregation_"+str(day.year)+".ncml.ascii"]
    path = "/".join(s for s in pieces)
    url = urljoin(opendap_server, path)

    return url

def determine_opendap_agg_url(location, start, end, server, level, product, variables):
    """Build a OpenDAP aggregation URL to enable fetching of a time series of a
    set of <variables> from the specified <server> for the particular
    products ([product]) at a single point (<location>) between the
    timespan of <start> -> <end>"""

    url = None

    url = determine_opendap_agg_base_url(start, server, level, product)
    lat_index, long_index = determine_index(location, lat0=90, lonres=1, latres=1)
    lat_index = max(0,lat_index-1)
    long_index = max(0,long_index-1)
    start_time_index = int(start.strftime("%j"))-1
    end_time_index = int(end.strftime("%j"))-1

    url += "?"
    for quantity in variables:
        url += "{:s}[{:d}:{:d}][{:d}:{:d}][{:d}:{:d}],".format(quantity, start_time_index, end_time_index, lat_index, lat_index, long_index,long_index)
    # Add trailing spatial and temporal indexes
    url += "Latitude[{:d}:{:d}],Longitude[{:d}:{:d}],time[{:d}:{:d}]".format(lat_index, lat_index, long_index,long_index, start_time_index, end_time_index)
    return url

def find_opendap_catalog(day, opendap_server, process_level, product):
    """Return the path to the XML OpenDap catalog of products.
    This takes <day> (a `datetime`), the <opendap_server> (e.g. https://acdisc.gesdisc.eosdis.nasa.gov/),
    the <process_level> (e.g. 'Aqua_AIRS_Level3') and the product (e.g. 'AIRS3STD.006')
    """

    pieces = [ 'opendap', process_level, product, str(day.year), 'catalog.xml']
    path = "/".join(s for s in pieces)
    url = urljoin(opendap_server, path)

    return url

def find_opendap_products(catalog_url, max_files=400):
    """Parses the passed THREDDS Catalog XML file at <catalog_url> and returns
    the names of the HDF files in the catalog as a list"""
    hdf_files = []

    with ur.urlopen(catalog_url) as response:
        tree = etree.parse(response)
        for item in tree.iter():
            if item.tag.endswith('dataset'):
                if item.attrib['name'].endswith('.hdf'):
                    hdf_files.append(item.attrib['name'])
            if len(hdf_files) >= max_files:
                break
        hdf_files.sort()

    return hdf_files

def determine_opendap_url(day, hdf_files_list, catalog_url):
    """Look for HDF files containing the <day> in the list of <hdf_files_list> (from
    find_opendap_products()).
    Returns a list of URLs to obtain the products or for remote query with PyDAP"""

    hdf_files = [hdf for hdf in hdf_files_list if day.strftime("%Y.%m.%d") in hdf]
    hdf_urls = [catalog_url.replace('catalog.xml', hdf) for hdf in hdf_files]

    return hdf_urls

def extract_opendap_data(opendap_url, DATAFIELD_NAME):
    username, password = get_netrc_credentials()
    if username and password:
        session = setup_session(username, password, check_url=opendap_url)
        dataset = open_url(opendap_url, session=session)
        data3D = dataset[DATAFIELD_NAME]
        data = data3D.data[:,:]

        # Read geolocation dataset.
        # For...reasons... the Latitude and Longitude arrays come back as 1D
        # via pydap rather than the correct 2D ones, via pyhdf. So we need
        # to grow them back to the right shape
        lat = dataset['Latitude'].data
        lon = dataset['Longitude'].data
        latitude = np.array([lat[:,:],]*lon.shape[0]).transpose()
        longitude = np.array([lon[:,:],]*lat.shape[0])

        # Handle fill value.
        attrs = data3D.attributes
        fillvalue=attrs["_FillValue"]
        data[data == fillvalue] = np.nan
        data = np.ma.masked_array(data, np.isnan(data))

    return data, latitude, longitude

def fetch_realtime(day, products=['O3_RT', 'PWV_RT'], dbg=False):

    prior_catalog_url = None
    prior_hdf_file = None
    datasets = {}
    for product in products:
        mapping = dataset_mapping(product)
        if len(mapping) == 4:
            catalog_url = find_opendap_catalog(day, mapping['server'], mapping['level'], mapping['product'])
            hdf_files = find_opendap_products(catalog_url)
            if dbg:
                print("Found {} products at {}".format(len(hdf_files), catalog_url))
            data_url = determine_opendap_url(day, hdf_files, catalog_url)
            quantity = product.split('_')[0]
            if dbg:
                print("Fetching {}".format(data_url))
            data, latitude, longitude = extract_opendap_data(data_url[0], mapping[quantity])
            datasets[quantity] = {'data' : data, 'lat' : latitude, 'long' : longitude}
    return datasets

def fetch_airs_ascii_timeseries(location, filename=None, start=None, end=None, products=['O3_RT', 'PWV_RT']):

    mapping = dataset_mapping(products[0])
    variables = [dataset_mapping(product)[product.split('_')[0]] for product in products]

    agg_url = determine_opendap_agg_url(location, start, end, mapping['server'], mapping['level'], mapping['product'], variables)

    quantities = "_".join(variables)
    status_code = -1
    opener = earthdata_login()
    if opener:
        ur.install_opener(opener)
        request = ur.Request(agg_url)
        r = ur.urlopen(request)

        date_fmt = "%Y%m%dT%H%M%S"
        if start.time() == time(0,0) and end.time() == time(0,0):
            # No time part, in start or end, remove from filename
            date_fmt = "%Y%m%d"

        filename = filename or "{}_{}_{}-{}.asc".format(mapping['level'], quantities, start.strftime(date_fmt), end.strftime(date_fmt))
        with open(filename, 'wb') as f:
            f.write(r.read())
        status_code = r.status

    return status_code, filename

def fetch_LCO_weather(site_code, start=None, end=None, interval=600, dbg=False):
    """Fetch LCO weather (temperature and pressure) for LCO <site_code> (e.g. 'ogg')
    between [start] and [end] (defaults to start of current year and now) interpolated
    to a spacing of [interval] seconds (defaults to 600s).
    Returns an AstroPy QTable of UTC datetime and temperature and pressure"""

    site_code = site_code.lower()
    start = start or datetime(datetime.utcnow().year, 1, 1)
    end = end or datetime.utcnow().replace(hour=0, minute=0, second=0, microsecond=0)
    nrows = int((end-start)/ timedelta(seconds=interval))

    # Construct a Q(uantity)Table and a column of datetime64's from <start> to
    # <end> with [interval] spacing
    table = QTable()
    dt = np.arange(start, end, step=interval, dtype='datetime64[s]')
    aa = Column(dt, name='UTC Datetime')
    table.add_column(aa)

    for quantity in ['temperature', 'pressure']:
        datum = map_quantity_to_LCO_datum(quantity)
        data = query_LCO_telemetry(site_code, start, end, datum)
        if len(data) > 0:
            interp_timestamps, interp_values = interpolate_LCO_telemetry(data, interval)
            unit = interp_values[0].unit
            if dbg: print(quantity, interp_timestamps[0], interp_timestamps[-1], len(interp_timestamps), len(interp_values))

            # Check if the interpolated dataset starts late or ends early and pad
            # accordingly
            num_before = max(int((interp_timestamps[0]-start) / timedelta(seconds=interval)), 0)
            num_after  = max(int((end-interp_timestamps[-1]) / timedelta(seconds=interval)), 0)
            pad_values = np.pad(interp_values, (num_before, num_after), 'edge')
            # Put units back
            pad_values = pad_values * unit
            if dbg: print("Padding by {} before, {} after, new length={}".format(num_before, num_after, len(pad_values)))
            # Trim padded array to right length, turn into a column and add to table
            col = Column(pad_values[0:nrows], name=quantity)
            table.add_column(col)
        else:
            print("Found no data for {} at {} between {}->{}".format(datum, site_code, start, end))
    return table

def populate_PWV_column(location, combined_table, dbg=False):
    """Repopulates the PWV column of the passed <combined_table> by calculating
    the PWV from the TotalZenithDelay, temperature and pressure. The value is not
    updated if already positive (>0) as -9.9 is normally used to indicate missing
    data. The PWVerr column is set to -9.9 mm"""

    for i, row in enumerate(combined_table):
        if row['PWV'] < 0.0:
            pwv = compute_pwv(row['TotalZenithDelay'], location, row['pressure'], row['temperature'])
            if dbg: print(row['TotalZenithDelay'], row['pressure'], row['temperature'], pwv)
            combined_table[i]['PWV'] = pwv
            combined_table[i]['PWVerr'] = -9.9 * u.mm
    return combined_table

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

def plot_pwv_timeseries(pwv_table, times=None, merra2_pwv=None, site='CTIO', filename=None):
    """Plot a PWV timeseries table (from fetch_pwv()) and optionally a MERRA-2
    comparison from extract_MODIS_pwv_timeseries() ([times] and [merra2_pwv]). [site]
    is used as the prefix to the labels and to the filename (if not specified)"""

    fig, ax = plt.subplots()
    ax.plot(pwv_table['UTC Datetime'], pwv_table['PWV'], label=site+' GPS')
    if times and merra2_pwv:
        ax.plot(times, merra2_pwv, label=site+' MERRA-2')
    else:
        times = [datetime.max, datetime.min]
    ylims = ax.get_ylim()
    ax.set_ylim(0, ylims[1])
    ax.set_title('Precipitable Water Vapour')
    ax.set_ylabel('PWV (mm)')

    timespan = max(pwv_table['UTC Datetime'][-1], times[1]) - min(pwv_table['UTC Datetime'][0], times[0])
    if timespan.days <= 2:
        hours = mdates.HourLocator(range(0,25,3))
        hoursFmt = mdates.DateFormatter('%m-%d %H:%M')
        ax.xaxis.set_major_locator(hours)
        ax.xaxis.set_major_formatter(hoursFmt)
        hours = mdates.HourLocator(range(0,25,1))
        ax.xaxis.set_minor_locator(hours)
        ax.format_xdata = mdates.DateFormatter('%Y-%m-%d %H:%M')
    else:
        dates = mdates.AutoDateLocator()
        ax.xaxis.set_major_locator(dates)
        datesFmt = mdates.AutoDateFormatter(dates)
        ax.xaxis.set_major_formatter(datesFmt)

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    minorLocator = AutoMinorLocator()
    ax.yaxis.set_minor_locator(minorLocator)

    plt.legend()
    fig.autofmt_xdate()
    if not filename:
        filename = site+'_GPS_MERRA2_comp.png'
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
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)
    return pngfile
