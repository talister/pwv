from math import cos
from datetime import datetime, timedelta

import astropy.units as u
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
