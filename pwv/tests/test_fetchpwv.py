from datetime import datetime

from astropy.coordinates import EarthLocation
import astropy.units as u

from pwv.fetch_pwv import *


class TestDetermineIndex(object):
    def setup_method(self):
        self.FTN = EarthLocation(lon=-156.257029, lat=20.706657, height=3046.52)
        self.GES_example = EarthLocation(lon=-75.625, lat=37.375, height=0)

    def test_merra2_ftn(self):
        expected_long_index = 38
        expected_lat_index = 222

        lat_index, long_index = determine_index(self.FTN)

        assert expected_long_index == long_index
        assert expected_lat_index == lat_index

    def test_gldas_example(self):
        expected_long_index = 418
        expected_lat_index = 390

        lat_index, long_index = determine_index(self.GES_example, lon0=-179.875, lat0=-59.875, lonres=0.25, latres=0.25)

        assert expected_long_index == long_index
        assert expected_lat_index == lat_index

    def test_airs_ftn(self):
        expected_long_index = 23
        expected_lat_index = 69

        lat_index, long_index = determine_index(self.FTN, lon0=-179.030, lat0=89.03, lonres=359/360.0, latres=179/180.0)

        assert expected_long_index == long_index
        assert expected_lat_index == lat_index

class TestSplitMODISfilename(object):

    def test_joint_atm1(self):
        expected_product = 'MODATML2'
        expected_dt = datetime(2018, 7, 11, 8, 35)
        expected_collection = '61'

        filename = 'MODATML2.A2018192.0835.061.2018192195509.hdf'

        product, dt, collection = split_MODIS_filename(filename)

        assert expected_product == product
        assert expected_dt == dt
        assert expected_collection == collection

    def test_joint_atm_path(self):
        expected_product = 'MODATML2'
        expected_dt = datetime(2018, 7, 11, 8, 35)
        expected_collection = '61'

        filename = '/foo/bar/MODATML2.A2018192.0835.061.2018192195509.hdf'

        product, dt, collection = split_MODIS_filename(filename)

        assert expected_product == product
        assert expected_dt == dt
        assert expected_collection == collection

class TestDetermineGDSUrl(object):

    def test_ges_tutorial_example(self):
        expected_url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/dods/M2T1NXSLV.ascii?u50m[184105:192864][273][124]"

        loc = EarthLocation(lon=-103.125*u.deg, lat=46.05*u.deg)
        start = datetime(2001, 1, 1, 0, 30, 0)
        end = datetime(2001, 12, 31, 23, 30, 0)

        url = determine_gds_url(loc, start, end, quantity='u50m')

        assert expected_url == url

class TestDatasetMapping(object):

    def test_airs_ozone(self):
        expected_mapping = {'server' : 'https://acdisc.gesdisc.eosdis.nasa.gov/',
                            'level' : 'Aqua_AIRS_Level3',
                            'product' : 'AIRS3STD.006',
                            'O3' : 'TotO3_D'
                           }

        product = 'O3_RT'

        mapping = dataset_mapping(product)

        assert expected_mapping == mapping

    def test_airs_water(self):
        expected_mapping = {'server' : 'https://acdisc.gesdisc.eosdis.nasa.gov/',
                            'level' : 'Aqua_AIRS_Level3',
                            'product' : 'AIRS3STD.006',
                            'PWV' : 'TotH2OVap_D'
                           }

        product = 'PWV_RT'

        mapping = dataset_mapping(product)

        assert expected_mapping == mapping

class TestFindOpenDapcatalog(object):

    def test_airs(self):
        expected_url = 'https://acdisc.gesdisc.eosdis.nasa.gov/opendap/Aqua_AIRS_Level3/AIRS3STD.006/2018/catalog.xml'

        day = datetime(2018, 9, 15, 10, 15)
        server = 'https://acdisc.gesdisc.eosdis.nasa.gov/'
        level = 'Aqua_AIRS_Level3'
        product = 'AIRS3STD.006'

        url = find_opendap_catalog(day, server, level, product)

        assert expected_url == url

class TestDetermineOpenDAPAggBaseUrl(object):

    def test_airs_2009(self):
        expected_url = 'https://acdisc.gesdisc.eosdis.nasa.gov/opendap/ncml/aggregation/AIRS3STD.006/AIRS3STD.006_Aggregation_2009.ncml.ascii'

        day = datetime(2009, 9, 15, 10, 15)
        server = 'https://acdisc.gesdisc.eosdis.nasa.gov/'
        level = 'Aqua_AIRS_Level3'
        product = 'AIRS3STD.006'

        url = determine_opendap_agg_base_url(day, server, level, product)

        assert expected_url == url

    def test_airs_2018(self):
        expected_url = 'https://acdisc.gesdisc.eosdis.nasa.gov/opendap/ncml/aggregation/AIRS3STD.006/AIRS3STD.006_Aggregation_2018.ncml.ascii'

        day = datetime(2018, 9, 15, 10, 15)
        server = 'https://acdisc.gesdisc.eosdis.nasa.gov/'
        level = 'Aqua_AIRS_Level3'
        product = 'AIRS3STD.006'

        url = determine_opendap_agg_base_url(day, server, level, product)

        assert expected_url == url

class TestDetermineOpenDAPAggUrl(object):

    def test_airs_2009_O3(self):
        expected_url = 'https://acdisc.gesdisc.eosdis.nasa.gov/opendap/ncml/aggregation/AIRS3STD.006/AIRS3STD.006_Aggregation_2009.ncml.ascii?TotO3_D[0:364][50:50][103:103],Latitude[50:50],Longitude[103:103],time[0:364]'

        start = datetime(2009, 1, 1, 10, 15)
        end = datetime(2009, 12, 31, 23, 59)
        server = 'https://acdisc.gesdisc.eosdis.nasa.gov/'
        level = 'Aqua_AIRS_Level3'
        product = 'AIRS3STD.006'
        location =  EarthLocation(lon=-76.877, lat=39.1495, height=42)
        variables = ['TotO3_D']

        url = determine_opendap_agg_url(location, start, end, server, level, product, variables)

        assert expected_url == url

    def test_airs_2009_O3_H2O(self):
        expected_url = 'https://acdisc.gesdisc.eosdis.nasa.gov/opendap/ncml/aggregation/AIRS3STD.006/AIRS3STD.006_Aggregation_2009.ncml.ascii?TotO3_D[0:364][50:50][103:103],TotH2OVap_D[0:364][50:50][103:103],Latitude[50:50],Longitude[103:103],time[0:364]'

        start = datetime(2009, 1, 1, 10, 15)
        end = datetime(2009, 12, 31, 23, 59)
        server = 'https://acdisc.gesdisc.eosdis.nasa.gov/'
        level = 'Aqua_AIRS_Level3'
        product = 'AIRS3STD.006'
        location =  EarthLocation(lon=-76.877, lat=39.1495, height=42)
        variables = ['TotO3_D', 'TotH2OVap_D']

        url = determine_opendap_agg_url(location, start, end, server, level, product, variables)

        assert expected_url == url

class TestReadAscii(object):

    def setup_method(self):
        self.DAP_location =  EarthLocation(lon=-76.877, lat=39.1495, height=42)
        self.merra2_ascii = os.path.join('pwv', 'tests', 'data', 'M2T1NXSLV_tqv_20090101-20090104.asc')
        self.aqua_ascii = os.path.join('pwv', 'tests', 'data', 'Aqua_TotO3_D_TotH2OVap_D_test.asc')

    def test_merra2(self):
        expected_keys = ['tqv', 'time', 'lat', 'datetime']
        expected_num = 73
        expected_tqv_0 = 3.0582914
        expected_tqv_last = 8.244712
        expected_dt_0 = datetime(2009, 1, 1, 0, 30)
        expected_dt_last = datetime(2009, 1, 4, 0, 30)

        data = read_ascii(self.merra2_ascii)

        assert expected_keys == list(data.keys())
        for key in expected_keys:
            if key not in ['lat', 'lon']:
                assert expected_num == len(data[key]), "Check on {} failed".format(key)
        assert expected_tqv_0 == data['tqv'][0]
        assert expected_tqv_last == data['tqv'][-1]
        assert expected_dt_0 == data['datetime'][0]
        assert expected_dt_last == data['datetime'][-1]

    def test_aqua(self):
        expected_keys = ['Longitude', 'Latitude', 'TotH2OVap_D', 'TotO3_D', 'time', 'datetime']
        expected_num = 3
        expected_pwv_0 =  2.91797
        expected_pwv_last =  6.09766
        expected_o3_0 =  353.0
        expected_o3_last =  355.25
        expected_dt_0 = datetime(2009, 1, 1, 0, 0)
        expected_dt_last = datetime(2009, 1, 3, 0, 0)

        data = read_ascii(self.aqua_ascii)

        assert expected_keys == list(data.keys())
        for key in expected_keys:
            if key not in ['Latitude', 'Longitude']:
                assert expected_num == len(data[key]), "Check on {} failed".format(key)
        assert expected_pwv_0 == data['TotH2OVap_D'][0]
        assert expected_pwv_last == data['TotH2OVap_D'][-1]
        assert expected_o3_0 == data['TotO3_D'][0]
        assert expected_o3_last == data['TotO3_D'][-1]
        assert expected_dt_0 == data['datetime'][0]
        assert expected_dt_last == data['datetime'][-1]
