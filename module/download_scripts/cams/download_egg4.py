import cdsapi

c = cdsapi.Client()

c.retrieve(
    'cams-global-ghg-reanalysis-egg4',
    {
        'variable': [
            'carbon_dioxide', 'methane',
        ],
        'date': '2020-08-06/2020-08-06',
        'step': [
            '12', '15',
        ],
        'area': [
            54.98, 5.99, 47.3,
            15.25,
        ],
        'format': 'grib',
        'pressure_level': [
            '1', '2', '3',
            '5', '7', '10',
            '20', '30', '50',
            '70', '100', '150',
            '200', '250', '300',
            '400', '500', '600',
            '700', '800', '850',
            '900', '925', '950',
            '1000',
        ],
    },
    'download.grib')
