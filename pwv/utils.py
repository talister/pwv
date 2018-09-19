from math import cos, acos, degrees
from datetime import datetime, timedelta
import time
import warnings

import astropy.units as u
from astropy.time import Time


PascalPerMillibar = 100.0

def compute_zhd(pres, location):
    """Compute Zenith Hydrostatic Delay (ZHD) from total pressure at the
    site <location>. The model is that of Saastamoinen (1972)"""

    top = 0.0022767 * pres.to(u.hPa)
    bottom = 1.0 - (0.00266 * cos(2.0*location.lat.radian)) - (0.00028 * location.height.to(u.km).value)
    bottom = bottom
    zhd = (top.value/bottom) * u.m
    return zhd

def compute_kappa(T_2m):
    """Compute dimensionless constant kappa' for zenith wet delay (ZWD) from
    the temperature at 2m altitude (<T_2m>; in K)"""

    rho = 998.0 *(u.kg/u.m**3) # density of water
    R_v = 461.5 * (u.J/(u.kg*u.K)) # specific gas constant for water vapour
    c_1 = 77.604 * (u.K/u.Pa)
    c_2 = 17.0  * (u.K/u.Pa)
    c_3 = 3.776*10**5 * (u.K**2/u.Pa)
    m = 0.62198
    # Compute mean temperature. Need equivalencies as astropy thinks degrees C
    # and degrees K aren't the same and can't be converted <doh>
    T_mean = (70.2 * u.K) + (0.72 * T_2m.to(u.K, equivalencies=u.temperature()))
    kappa = 10**8/(rho * R_v * ((c_3/T_mean) + c_2 - (m*c_1)))
    return kappa.decompose()

def compute_zwd(ztd, location, pres, T_2m):
    """Compute the Zenith Wet Delay (ZWD) from the Zenith Total Delay (ZTD; <ztd>)
    minus the Zenith Hydrostatic Delay (ZHD). This also requires the location
    which should be a `astropy.coordinates.earth.EarthLocation`, the pressure
    at the location (<pres> as an astropy pressure Quantity) and temperature
    at the location (<pres> as an astropy temperature Quantity)
    The ZWD is returned as an astropy.Quantity with units of meters"""

    zhd = compute_zhd(pres, location)
    zwd = ztd.to(u.mm) - zhd.to(u.mm)

    return zwd

def compute_pwv(ztd, location, pres, T_2m):
    """Compute the Precipitable Water Vapourfrom the Zenith Total Delay (ZTD; <ztd>)
    minus the Zenith Hydrostatic Delay (ZHD). This also requires the location
    which should be a `astropy.coordinates.earth.EarthLocation`, the pressure
    at the location (<pres> as an astropy pressure Quantity) and temperature
    at the location (<pres> as an astropy temperature Quantity)
    The ZWD is returned as an astropy.Quantity with units of meters"""

    zwd = compute_zwd(ztd, location, pres, T_2m)
    kappa = compute_kappa(T_2m)
    pwv = kappa * zwd.to(u.mm)

    return pwv

def determine_time_index(date_or_str, t0='00:30Z01Jan1980'):

    if type(t0) != datetime:
        t0 = datetime.strptime(t0, '%H:%MZ%d%b%Y')

    if type(date_or_str) != datetime:
        t = datetime.strptime(date_or_str, '%H:%MZ%d%b%Y')
    else:
        t = date_or_str
    delta_t = t-t0
    index = int(delta_t.total_seconds()/3600) + 1

    return index

def time_index_to_dt(date, t0=1721423.5):
    """Convert decimal date number from GrADS server downloads from read_ascii()
    to datetime. The default offset [t0] is the Julian Date of 0001-01-01 00:00
    calculated by the USNO Julian Date Converter:
    http://aa.usno.navy.mil/jdconverter?ID=AA&year=1&month=1&day=1&era=1&hr=0&min=0&sec=0.0
    This is because astropy and SLALIB produce the wrong answers for these early dates.
    """

    if type(date) == float or type(date) == int:
        dates = [date,]
    else:
        dates = date
    times = Time(t0, dates, format='jd', scale='utc')
    times_dt = []
    for t in  times.datetime:
        if t.microsecond >= 500000:
            t += timedelta(seconds=1)
            t = t.replace(microsecond = 0)
        else:
            t = t.replace(microsecond = 0)
        times_dt.append(t)

    if type(date) == float:
        times_dt = times_dt[0]

    return times_dt

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def pascal_to_mbar(pascals):
    """Convert pressure from Pascals to millibars
    """
    try:
        return pascals.value/PascalPerMillibar
    except AttributeError:
        return pascals/PascalPerMillibar

def mbar_to_pascal(mbar):
    """Convert pressure from millibars to Pascals
    """
    return mbar*PascalPerMillibar*u.Pa

def co2_ppmv(date=datetime.utcnow()):
    """Returns the CO2 concentration in ppmv (or micromol/mol) as a function
    of the year. This is from a quadratic fit to the Mauna Loa dataset:
    ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_annmean_mlo.txt
    The fit uses 1976.0 as a pivot point as this is the year that the US
    standard atmosphere models are based.
    """

    if type(date) == int:
        date = datetime(date, 1, 1)
    year = toYearFraction(date)
    if year < 1950 or year > 2018:
        warnings.warn("Warning, year is beyond the end of the fit (1950..2017)")
    x = year-1976.0
    a0 = 333.038137957249
    a1 = 1.25213771117096
    a2 = 0.0124809170121469

    co2 = a0 + a1 * x + a2 * x**2

    return co2

def make_bounding_box(location, delta=0.1, num_dp=-1):
    """Calculate the West, South, East, North values of a bounding box around
    <location> with [delta] degrees (defaults to 0.1)"""

    lat = location.lat.deg
    lon = location.lon.deg

    north_val = max(lat+delta, lat-delta)
    if num_dp >= 0: north_val = round(north_val, num_dp)
    west_val  = min(lon+delta, lon-delta)
    if num_dp >= 0: west_val = round(west_val, num_dp)
    east_val  = max(lon+delta, lon-delta)
    if num_dp >= 0: east_val = round(east_val, num_dp)
    south_val = min(lat+delta, lat-delta)
    if num_dp >= 0: south_val = round(south_val, num_dp)

    return west_val, south_val, east_val, north_val

def airmass_to_sza(airmass):
    """Convert airmass to solar zenith angle (for e.g. libRadtran)
    The simple 1/cos(X) formula is used although is not very good at large
    zenith distances"""

    sza = degrees(acos(1.0/airmass))

    return sza
