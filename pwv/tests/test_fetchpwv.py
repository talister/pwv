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
    def test_airs(self):
        expected_mapping = {'server' : 'https://acdisc.gesdisc.eosdis.nasa.gov/',
                            'level' : 'Aqua_AIRS_Level3',
                            'product' : 'AIRS3STD.006',
                            'O3' : 'TotO3_D'
                           }

        product = 'O3_RT'

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
