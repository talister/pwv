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
