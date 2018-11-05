# pwv
Code to fetch Precipitable Water Vapo(u)r measurements and compute atmospheric transmission.

## Installation

A Python Virtual Environment (virtualenv) is recommended. This can be constructed by executing:
```bash
python3 -m venv <path to virtualenv>
source <path to virtualenv>/bin/activate # for bash-shells
```

or:

`source <path to virtualenv>/bin/activate.csh # for (t)csh-shells`

then:
`pip3 install numpy`
`pip3 install -r requirements.txt`

Numpy needs to be installated separately or the latter command will
fail.
Depending on your distribution you may also have to install development packages for libhdf4 and python.
If you want to compute atmospheric transmissions, you will need to get and install libRadtran from http://www.libradtran.org/doku.php?id=download
Only the basic program is needed at this time, none of the additional modules are needed.

## Usage

### Atmospheric Parameters

There are a couple of ways of obtaining the Precipitable Water Vapour (PWV) and other parameters needed to model atmospheric transmission. The most direct way if there is a SuomiNet (https://www.suominet.ucar.edu/) site nearby is to use `pwv.fetch_pwv.fetch_GPS_pwv(<suominet site code>)`. It may be that there are SuomiNet data files for a site which don't contain valid PWV measurements but do contain a Total Zenith Delay (ZTD). In this case if you can obtain a source of pressure and temperature from another source, it is possible to recreate the PWV using `pwv.fetch_pwv.populate_PWV_column()`; `pwv.fetch_pwv.fetch_pwv()` illustrates this by using the LCO site telemetry to backfill the PWV value.

There is also support for obtaining the PWV for the telescopes on La Palma, Canary Islands (possibly also valid on Tenerife as well) through the data of the Sky Quality Group of the Roque de los Muchachos Observatory (ORM) (http://vivaldi.ll.iac.es/proyecto/site-testing/index.php?option=com_wrapper&Itemid=122) through the `pwv.fetch_ORM_pwv()` method (or through `pwv.fetch_pwv.fetch_pwv(<site_code>)` where `<site_code>` is one of `{'NOT', 'TNG', 'WHT', 'TFN'}`)

Other sources of data are the MERRA-2 reprocessed data or the MODIS and AIRS instruments on the NASA Aqua/Terra satellites; methods to query these are provided. These require that you have an EarthData login by registering at https://urs.earthdata.nasa.gov/ Once this is done, you should create a file in your home directory called `.netrc` which should contain:
`machine urs.earthdata.nasa.gov login <earthdata login> password <earthdata password>`
(You should do a `chmod 600 ~/.netrc` on it and not use the same password anywhere else as it will be stored in plaintext on your computer)

### Atmospheric Transmission

Once you have an idea of the typical values and range of the following parameters:
* atmospheric pressure
* precipitable water vapor
* amount of ozone
* amount of aerosols
from GPS, satellite or local data for your site, you can make use of libRadtran to compute a grid of atmospheric transmission profiles as a function of wavelength.
`pwv.compute_atm` will produce a grid of atmospheric profiles as a function of airmass (1..3 by default) and for ranges of the PWV and ozone amounts.

## Notes

This experimental code under development... Feel free to contact me (Tim Lister) through github if there are issues or for improvements in the code or the documentation (not difficult)...
