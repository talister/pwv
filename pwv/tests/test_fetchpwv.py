from datetime import datetime

from astropy.coordinates import EarthLocation
import astropy.units as u

from pwv.fetch_pwv import *

class TestConvertDecimalDay(object):

    def test_Jan1(self):
        expected_dt = datetime(datetime.utcnow().year, 1, 1, 4, 15, 0)

        decimal_day = 1.17708

        dt = convert_decimal_day(decimal_day)

        assert expected_dt == dt

    def test_Jun12(self):
        expected_dt = datetime(datetime.utcnow().year, 6, 12, 18, 45, 0)

        decimal_day = 163.78125

        dt = convert_decimal_day(decimal_day)

        assert expected_dt == dt

    def test_2016Jun11(self):
        # 2016 is a leap year so this is 1 day earlier
        expected_dt = datetime(2016, 6, 11, 18, 45, 0)

        decimal_day = 163.78125

        dt = convert_decimal_day(decimal_day, 2016)

        assert expected_dt == dt

    def test_2016Jun11_dt(self):
        expected_dt = datetime(2016, 6, 11, 18, 45, 0)

        decimal_day = 163.78125

        dt = convert_decimal_day(decimal_day, year=datetime(2016,1,1))

        assert expected_dt == dt

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
