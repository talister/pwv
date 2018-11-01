import os
from subprocess import Popen, PIPE, STDOUT, call

import astropy.units as u
import numpy as np

from utils import airmass_to_sza

class UVspec:
    """Class to wrap setup and invocation of uvspec from libRadtran"""

    def __init__(self, path='', atmosphere='us'):
        self.Prog = 'RT'
        self.Ver = '2.0.2'
        self.Rte = 'pp'     # pp for parallel plane or ps for pseudo-spherical
        self.Proc = 'as'    # as for aerosol special
        if path == '':
            self.path = os.getenv('LIBRADTRANDIR', os.getenv('HOME'))
        else:
            self.path = path
        self.atmosphere = atmosphere
        self.molmodel = 'rt'
        self.molresol = 'coarse'
        data_path = os.path.join(self.path, 'data')
        self.inp = { 'data_files_path'  : data_path,
                     'atmosphere_file'  : os.path.join(data_path, 'atmmod', self.atmosphere + '.dat'),
                     'source'           : 'solar ' + os.path.join(data_path, 'solar_flux', 'kurudz_1.0nm.dat'), # This typo of Kurucz is correct...
                     'albedo'           : 0.2,
                     'mol_abs_param'    : self.molmodel + ' ' + self.molresol,
                     'output_quantity'  : 'reflectivity',
                     'phi0'             : '0',
                     'rte_solver'       : 'disort',
                     'output_user'      : 'lambda edir',
                     'quiet'            : ''
                    }

    @property
    def atmosphere(self):
        return self._atmosphere

    @atmosphere.setter
    def atmosphere(self, atm='us'):
        atmosphere_map = {  'us' : 'afglus',
                            'ms' : 'afglms',
                            'mw' : 'afglmw',
                            'tp' : 'afglt',
                            'ss' : 'afglss',
                            'sw' : 'afglsw'
                          }
        if atm not in atmosphere_map:
            raise ValueError("Invalid atmosphere choice. Valid choices are: " + ", ".join(atmosphere_map.keys()))
        self._atmosphere = atmosphere_map[atm]
        self.atmkey = atm

    @property
    def molmodel(self):
        return self._molmodel

    @molmodel.setter
    def molmodel(self, mol='rt'):
        molmodel_map = { 'rt' : 'reptran',
                         'lt' : 'lowtran',
                         'kt' : 'kato',
                         'k2' : 'kato2',
                         'fu' : 'fu',
                         'cr' : 'crs'
                       }
        if mol not in molmodel_map:
            raise ValueError("Invalid molecular model choice. Valid choices are: " + ", ".join(molmodel_map.keys()))
        self._molmodel = molmodel_map[mol]
        self.molkey = mol

    def write_input(self, fn):
        with open(fn, 'w') as f:
            for key in sorted(self.inp):
                if key=="mol_modify2":
                    f.write( "mol_modify" + ' ' + str(self.inp[key]) + '\n')
                else:
                    f.write( key + ' ' + str(self.inp[key]) + '\n')

            f.close()

    def run(self, inp, out, verbose):
        if verbose:
            print("Running uvspec with input file: ", inp)
            print("Output to file                : ", out)

        cmd = os.path.join(self.path, 'bin', 'uvspec') +  ' < ' + inp  +  ' > ' + out
        if verbose:
            print("uvspec cmd: ", cmd)
        p = call(cmd,shell=True,stdin=PIPE,stdout=PIPE)
        return p

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(f):
        os.makedirs(f)

def Sim_Atm(site, airmass, pwv, oz, wl0, tau0, press, altitude):
    """Compute simulated atmospheric transmission using libRadTran's uvspec"""

    verbose = True
    uvspec = UVspec()

    rootdir = os.path.join('simulations', uvspec.Prog, uvspec.Ver, site)
    ensure_dir(rootdir)
    BaseFilename_part1 = uvspec.Prog + '_' + site + '_' +uvspec.Rte+'_'

    # manage input and output directories and vary the ozone
    rootdir2 = os.path.join(rootdir, uvspec.Rte, uvspec.atmkey, uvspec.Proc, uvspec.molkey)
    ensure_dir(rootdir2)
    input_dir = os.path.join(rootdir2, 'in')
    ensure_dir(input_dir)
    output_dir = os.path.join(rootdir2, 'out')
    ensure_dir(output_dir)

    #water vapor
    pwv_str = 'H2O ' + str(pwv) + ' MM'
    wvfileindex = pwv

    #aerosols
    aerosol_str = str(wl0) + ' ' + str(tau0)
    aer_index = int(tau0*1000.)

    # airmass
    airmass = airmass
    amfileindex = int(airmass*10)

    # Ozone
    oz_str = 'O3 ' + str(oz) + ' DU'
    ozfileindex = int(oz/1.)

    BaseFilename = BaseFilename_part1 + uvspec.atmkey + '_' + uvspec.Proc +\
        '_' + uvspec.molkey + '_z' + str(amfileindex) + '_wv' + str(wvfileindex) +\
        '_oz' + str(ozfileindex) + '_aer' + str(aer_index)

    uvspec.inp["aerosol_default"] = ''
    uvspec.inp["aerosol_set_tau_at_wvl"] = aerosol_str

    # set up the ozone value
    uvspec.inp["mol_modify"] = pwv_str
    uvspec.inp["mol_modify2"] = oz_str

    # rescale pressure   if reasonable pressure values are provided
    if press > 600.0 and press < 1015.0:
        uvspec.inp["pressure"] = press
    else:
        raise ValueError("crazy pressure p=" + press + ' hPa')

    uvspec.inp["altitude"] = altitude
    sza = airmass_to_sza(airmass)
    uvspec.inp["sza"] = str(sza)
    uvspec.inp["wavelength"] = '300.0 1200.0'

    input_filename = BaseFilename + '.INP'
    output_filename = BaseFilename + '.OUT'
    inp = os.path.join(input_dir, input_filename)
    out = os.path.join(output_dir, output_filename)

    uvspec.write_input(inp)
    uvspec.run(inp, out, verbose)

    return output_dir, output_filename

def get_site_params(site):
    """Return the default parameters for the passed site"""

    site = site.lower()

    site_mapping = { 'FTN' :  { 'altitude' :  3.055 , # km
                                'avg_pressure' : 710,   # millibars
                                'aerosol_tau0' : 0.075, # Mean value for 2016
                                'pwv_range' : (5, 50),
                                'ozone_range' : (240, 300),
                              },
                     'FTS' :  { 'altitude' :  1.116 , # km
                                'avg_pressure' : 890,   # millibars
                                'aerosol_tau0' : 0.067,
                                'pwv_range' : (10, 60),
                                'ozone_range' : (230, 310),
                              },
                    'CTIO' :  { 'altitude' :   2.200 , # km
                                'avg_pressure' : 780,   # millibars
                                'aerosol_tau0' : 0.018,
                                'pwv_range' : (10, 60),
                                'ozone_range' : (230, 310),
                              }
                   }

    return site_mapping.get(site, None)

if __name__ == "__main__":

    if os.getenv('LIBRADTRANDIR', None) is None:
        raise EnvironmentError("Need to set LIBRADTRANDIR to the root of the libRadTran installation")

    wl0 = 550.0 # Wavelength for Total Aerosol Extinction (from MODIS)

    site = 'NOT'
    site_params = get_site_params(site)
    if site_params is None:
        site_params = { 'altitude' :   2.382,  # km,
                        'avg_pressure' : 770,   # millibars
                        'aerosol_tau0' : 0.018,
                        'pwv_range' : (10, 60),
                        'ozone_range' : (230, 310),
                      }

    print('*****************************************************')
    print('Computing for {} (altitude= {}m)'.format(site, site_params['altitude']))
    print('*****************************************************')
    print()

    for airmass in np.linspace(1,3,21):
        airmass_nb = round(airmass, 1)

        for pwv_nb in np.arange(site_params['pwv_range'][0], site_params['pwv_range'][1]+0.01, 5):
            for oz_nb in np.arange(site_params['ozone_range'][0], site_params['ozone_range'][1]+0.01, 10):
                path, outputfile = Sim_Atm(site,airmass_nb,pwv_nb,oz_nb,wl0,site_params['aerosol_tau0'],site_params['avg_pressure'], site_params['altitude'])
                print('*****************************************************')
                print(' path       = ', path)
                print(' outputfile =  ', outputfile)
                print('*****************************************************')
