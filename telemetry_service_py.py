import requests
from datetime import datetime
from bokeh.plotting import figure, show
from bokeh.palettes import Dark2_8 as palette

BASE_URL = 'http://telemetryservice.lco.gtn/get'


def get_telemetry(datum_name, start, end, site=None, observatory=None, telescope=None, cadence=10):
    data = {'datum': datum_name,
            'start': start.strftime('%Y-%m-%dT%H:%M:%S.0'),
            'end': end.strftime('%Y-%m-%dT%H:%M:%S.0'),
            'cadence': cadence}
    if site:
        data['site'] = site
    if observatory:
        data['observatory'] = observatory
    if telescope:
        data['telescope'] = telescope

    result = requests.get(BASE_URL, params=data)

    results = []
    if result in [200,201]:
        results = result.json()
    else:
        print("Error retrieving results")
    return results


def clean_telemetry(raw_telemetry):
    values = []
    timestamps = []
    previous_timestamp = datetime(2000,1,1,1)
    for rt in reversed(raw_telemetry):
        if rt['value'] != 'NaN':
            converted_timestamp = datetime.fromtimestamp(rt['timestamp'] / 1000)
            if converted_timestamp < previous_timestamp:
                print("Timestamps out of order: {} is not less than {}. Ignoring.".format(previous_timestamp, converted_timestamp))
            else:
                values.append(rt['value'])
                timestamps.append(converted_timestamp)
                previous_timestamp = converted_timestamp

    return {'values': values, 'timestamps': timestamps}


def plot(title, data_dict, values_label='values', timestamps_label='timestamps'):
    p = figure(title=title, x_axis_label=timestamps_label, x_axis_type='datetime', y_axis_label=values_label)
    i = 0
    for key, item in data_dict.items():
        p.line(item['timestamps'], item['values'], legend=key, line_width=2, color=palette[i % 8])
        i += 1
    show(p)


def output_data(data_dict, filename):
    with open(filename, 'w') as out_fh:
        for i in zip(data_dict['timestamps'], data_dict['values']):
            print(i[0],i[1],file=out_fh)
    return

def main():
    data_dict = {}

    site = 'lsc'
    observatory = 'lsc'
    start = datetime(2016,7,1)
    end = datetime(2016,7,31,23,59,59)

    datum = 'Weather Barometric Pressure Value'
    results = get_telemetry(datum, start, end, site=site, observatory=observatory, cadence=60*60)
    print("Fetch complete")
    data_dict['Pressure'] = clean_telemetry(results)
    print("Clean complete")

    filename = "{}_pressure_{}-{}.dat".format(site, start.strftime("%Y%m%d"), end.strftime("%Y%m%d"))
    output_data(data_dict['Pressure'], filename)
#    datum = 'Weather Particulate Concentration 1.0 Micron Value'
#    results = get_telemetry(datum, start, end, site=site, cadence=60*60)
#    print("Fetch complete")
#    data_dict['PM1.0'] = clean_telemetry(results)
#    print("Clean complete")

#    plot('Pressure', data_dict)


if __name__ == "__main__":
    main()
