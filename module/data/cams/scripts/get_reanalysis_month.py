 #!/usr/bin/python
 #*** Script to download CAMS reanalysis data (https://confluence.ecmwf.int/display/CKB/CAMS%3A+Reanalysis+data+documentation#CAMS:Reanalysisdatadocumentation-Dataorganisationandaccess)
 #*** as daily files
 #*** Run it by entering "python get_reanalysis_day.py startdate stopdate"
 #*** Dates in format yyyymmdd
 #*** Andre Butz
 #*** Nov. 2020

import cdsapi
import yaml
from datetime import *
import calendar
import argparse
import sys
import os
### The new ECMWF-ADS/CDS batch acces requires ~/.cdsapirc

def main(startdate,stopdate):
    
    # For retrieving a month at once
    yi=int(startdate[0:4])
    mi=int(startdate[4:6])

    yf=int(stopdate[0:4])
    mf=int(stopdate[4:6])
    
    for i in month_year_iter(mi,yi,mf,yf):
        year = str(i[0])
        month = str('{:02d}'.format(i[1]))
        day = str(calendar.monthrange(i[0],i[1])[1])
        date_string = year +'-'+ month +'-'+ '01/' + year +'-'+ month +'-'+ day
        file_string = 'cams_rea'+year + month + '.grb'
        file_string_sfc = 'cams_reasfc'+year + month + '.grb'
        print('retrieving '+year+'-'+month)
        retrieve_data(date_string,file_string,file_string_sfc)

    sys.exit()

##########################################################
#
def retrieve_data(date_string, file_string, file_string_sfc):

    os.system('cp ~/.cdsapirc_ads ~/.cdsapirc')
    c = cdsapi.Client()

    c.retrieve(
    'cams-global-ghg-reanalysis-egg4',
    {
        'format': 'grib',
        'variable': [
            'carbon_dioxide', 'methane', 
        ],
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
        'date': date_string,
        'step': [
            '0', '6', '12', '18',
        ],
        'area': [
            54.98, 5.99, 47.3, 15.25
        ],
    },
    'download/'+file_string)
    
    # c.retrieve(
    # 'cams-global-ghg-reanalysis-egg4',
    # {
    #     'format': 'grib',
    #     'variable': 'surface_pressure',
    #     'date': date_string,
    #     'step': [
    #         '0', '6', '12', '18',
    #     ],
    #     'area': [
    #         54.98, 5.99, 47.3, 15.25
    #     ],
    # },
    # 'download/'+file_string_sfc)    
    
##########################################################
#
#*** Returns tuples (year,month) for given range
def month_year_iter( start_month, start_year, end_month, end_year ):
    ym_start = 12 * start_year + start_month - 1
    ym_end   = 12 * end_year   + end_month   
    for ym in range( ym_start, ym_end ):
        y, m = divmod( ym, 12 )
        yield y, m+1

#########################################################################################
# Argument Parser
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to download daily CAMS files with variables CO2, CH4, surf p.")
    parser.add_argument("startdate", help="Start date in format yyyymmdd.")
    parser.add_argument("stopdate", help="Stop date in format yyyymmdd.")    
    args = parser.parse_args()
    return args
    
#########################################################################################

if __name__ == '__main__':
    args = parse_arguments()
    main(**vars(args))
