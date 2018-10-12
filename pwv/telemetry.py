from datetime import datetime, timedelta

from elasticsearch import Elasticsearch
from elasticsearch_dsl import Search

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

    client = Elasticsearch(hosts='elasticsearch.lco.gtn')
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
