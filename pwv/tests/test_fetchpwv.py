from datetime import datetime

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
