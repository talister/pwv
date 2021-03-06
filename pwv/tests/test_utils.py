from astropy.coordinates import EarthLocation
import astropy.units as u
from astropy.units.cds import mmHg
from astropy.tests.helper import assert_quantity_allclose

import pytest

from pwv.utils import *

class TestComputeKappa(object):

    def test_in_K(self):
        expected_kappa = 0.16046773 * u.dimensionless_unscaled
        temp = 281.35*u.K

        kappa = compute_kappa(temp)

        assert_quantity_allclose(expected_kappa, kappa)

    def test_in_C(self):
        expected_kappa = 0.16046773 * u.dimensionless_unscaled
        temp = 8.2*u.deg_C

        kappa = compute_kappa(temp)

        assert_quantity_allclose(expected_kappa, kappa)

class TestZHD(object):

    def setup_method(self):
        self.FTN = EarthLocation(lon=-156.257029, lat=20.706657, height=3046.52)
        self.ELP = EarthLocation(lon=-104.02199,  lat=30.68049, height=2057.185)

    def test_FTN_hPa(self):
        expected_zhd = 1.71240176993675 * u.m

        pressure = 750.0 * u.hPa

        zhd = compute_zhd(pressure, self.FTN)

        assert_quantity_allclose(expected_zhd, zhd)

    def test_FTN_mmHg(self):
        expected_zhd = 1.71240176993675 * u.m

        pressure = 562.546171 * mmHg

        zhd = compute_zhd(pressure, self.FTN)

        assert_quantity_allclose(expected_zhd, zhd)

    def test_ELP_inHg(self):
        expected_zhd = 1.82597434743557 * u.m

        pressure = (23.64 * u.imperial.inch).to(u.mm).value * mmHg

        zhd = compute_zhd(pressure, self.ELP)

        assert_quantity_allclose(expected_zhd, zhd)

class TestComputePWV(object):

    def setup_method(self):
        self.CTIO = EarthLocation(lon=-70.806885, lat= -30.169165, height= 2218.45)

    def test_CTIO(self):
        # From SuomiNet CTIO_2018global.plt
        # DayOfYear   PWV PWVerr ZTD    Press Temp(C) RH(%)
        # 158.98958   3.3   1.1 1801.3  780.4  13.7 100.0

        expected_pwv = 3.44 * u.mm

        ztd = 1801.3 * u.mm
        pressure = 780.4 * u.hPa
        temp = 13.7 * u.deg_C

        pwv = compute_pwv(ztd, self.CTIO, pressure, temp)

        assert_quantity_allclose(expected_pwv, pwv, rtol=1e-3)


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


class TestDetermineTimeIndex(object):

    def test_t1(self):
        expected_value = 184105

        t_value = determine_time_index('00:30Z01JAN2001')

        assert expected_value == t_value

    def test_t2(self):
        expected_value = 192864

        t_value = determine_time_index('23:30Z31DEC2001')

        assert expected_value == t_value

    def test_t1_dt(self):
        expected_value = 184105

        t_value = determine_time_index(datetime(2001,1,1,0,30,0))

        assert expected_value == t_value

    def test_t2_dt(self):
        expected_value = 192864

        t_value = determine_time_index(datetime(2001,12,31,23,30))

        assert expected_value == t_value


class TestTimeIndexToDt(object):
    def test_scalar(self):
        expected_dt = datetime(2018, 1, 1, 0, 30, 0, 0)

        times = 736696.0208333334

        dt = time_index_to_dt(times)

        assert expected_dt == dt

    def test_array(self):
        expected_dt = [ datetime(2018, 4, 30, 21, 30),
                        datetime(2018, 4, 30, 22, 30),
                        datetime(2018, 4, 30, 23, 30)
                      ]

        times = [736815.8958333334, 736815.9375, 736815.9791666666]

        dts = time_index_to_dt(times)

        assert expected_dt == dts

class TestPascalMbarConversion(object):

    def test_pascal_to_mbar(self):
        expected_mbar = 0.01
        assert expected_mbar == pascal_to_mbar(1)

    def test_mbar_to_pascal(self):
        expected_pascals = 1.0
        assert expected_pascals == mbar_to_pascal(0.01)


    def test_pascal_to_mbar_Quantity(self):
        expected_mbar = 0.01
        assert expected_mbar == pascal_to_mbar(1 * u.Pa)

    def test_mbar_to_pascal(self):
        expected_pascals = 1.0 * u.Pa
        assert expected_pascals == mbar_to_pascal(0.01)
