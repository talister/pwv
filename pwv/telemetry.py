from datetime import datetime, timedelta

from elasticsearch import Elasticsearch
from elasticsearch_dsl import Search
import numpy as np
from scipy.interpolate import interp1d

from pwv.utils import round_datetime

def map_quantity_to_LCO_datum(quantity):
    """Simple mapper from concepts (<quantity>) to LCO telemetry datum name"""

    mapping = { 'pressure' : 'Weather Barometric Pressure Value',
                'temperature' : 'Weather Air Temperature Value'
              }

    return mapping.get(quantity.lower(), None)

def query_LCO_telemetry(site, start, end, datum='Weather Barometric Pressure Value', dbg=False):
    """Query the ElasticSearch service at LCO for telemetry for the specific <site>
    between <start> and <end> (datetimes) for the particular [datum] (defaults
    to 'Weather Barometric Pressure Value')

    A list of dictionaries, containing {'UTC Datetime', <datum>} is returned, sorted into chronological order"""

    client = Elasticsearch(hosts='elasticsearch.lco.gtn', retry_on_timeout=True)
    s = Search(using=client, index='mysql-telemetry-*')
    s = s.query('match', site=site)
    s = s.query('match', datumname=datum)
    s = s.filter('range', timestamp={'gte': start.strftime("%Y-%m-%dT%H:%M:%S"),
                                     'lte': end.strftime("%Y-%m-%dT%H:%M:%S")})
    s = s.filter('range', value_float = {'gt': 0.0})

    data = []
    for hit in s.scan():
        result = {  'UTC Datetime' : datetime.strptime(hit.timestampmeasured, "%Y-%m-%dT%H:%M:%S.%fZ"),
                    datum : hit.value_float}
        data.append(result)
        if dbg: print('HIT: site={site} datumname={datumname} timestampmeasured={timestampmeasured} value={value_float}'.format(**hit.to_dict()))

    sorted_data = sorted(data, key=lambda k:  k['UTC Datetime'])

    return sorted_data

def interpolate_LCO_telemetry(data, interval=300, kind='linear'):
    """Interpolates the LCO telemetry in <data> (assumed to be list of dicts
    containing 'UTC Datetime' and a datum - from query_LCO_telemetry()) to
    a uniform spacing of [interval] seconds via interpolation method [kind] (as
    defined in `scipy.interpolate.interp1d`
    (https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html)
    Returns numpy arrays of the interpolated timestamps (datetimes) and values"""

    times = [i['UTC Datetime'] for i in data]
    timestamps = [t.timestamp() for t in times]
    quantity = [key for key in data[0].keys() if key != 'UTC Datetime'][0]
    values = [x[quantity] for x in data]

    # Determine start and end times (rounded to the interval) and number of steps
    start = round_datetime(times[0], interval/60.0, round_up=True)
    end = round_datetime(times[-1], interval/60.0, round_up=False)
    delta = timedelta(seconds=interval)
    steps = int((end - start) / delta)
    increments = range(0, steps) * np.array([delta]*steps)
    # Generate new datetimes and timestamps
    interp_times = start + increments
    interp_timestamps = [x.timestamp() for x in interp_times]

    # Interpolate the data and then evaluate it at each of the
    # interpolated timestamps
    interp_func = interp1d(timestamps, values, kind)
    interp_values = interp_func(interp_timestamps)

    return interp_times, interp_values
