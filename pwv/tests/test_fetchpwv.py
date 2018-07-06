from datetime import datetime

import pytest

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
