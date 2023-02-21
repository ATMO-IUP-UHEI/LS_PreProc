#!/usr/bin/python
# Script to download ERA5 data from CDS (https://apps.ecmwf.int/data-catalogues/era5/?type=an&class=ea&stream=oper&expver=1) as monthly files
# Run it by entering "python get_cds_month.py startdate enddate"
# Dates in format yyyymm

import cdsapi
import calendar
import argparse
import os

### The CDS batch access requires ~/.cdsapirc with user credentials
# data download specifications: (see also https://confluence.ecmwf.int/display/UDOC/Keywords+in+MARS+and+Dissemination+requests)
cls     = "ea"              # do not change
expver  = "1"               # do not change
grid    = "0.25/0.25"       # 0.25degx0.25deg grid
levlistall = "1/to/137",    # every level
levlistsfc = "1",           # surface model level
levtype = "ml"              # do not change
levtypesfc = "sfc"          # do not change
paramsall  = "130/133"      # parameters to download on all levels: temperature, specific humidity
paramssfc  = "129.128/134.128/165.128/166.128"   # parameters to download on sfc level: u-, v-component of wind, geopotential, logarithm of surface pressure
stream  = "oper"            # do not change
tp      = "an"              # type: Use "an" (analysis) unless you have a particular reason to use "fc" (forecast). (https://confluence.ecmwf.int/pages/viewpage.action?pageId=85402030)
area = [54.98, 5.99, 47.3, 15.25]
time    = '14/15'     # time: ERA5 data is hourly. Specify a single time as "00:00:00", or a range as "00:00:00/01:00:00/02:00:00" or "00:00:00/to/23:00:00/by/1".


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
        date_string = year + month + '01/to/' + year + month + day
        file_string = year + month + '.grb'
        print('retrieving '+year+'-'+month)
        retrieve_data(date_string,file_string)

##########################################################
#
def retrieve_data(date_string,file_string):

    os.system('cp ~/.cdsapirc_cds ~/.cdsapirc')
    c = cdsapi.Client()

	# all levels
    c.retrieve('reanalysis-era5-complete', {
    'class': cls,
    'date': date_string,
    'expver': expver,
    'grid': grid,
    'levelist': levlistall,
    'levtype': levtype,
    'param': paramsall,
    'stream': stream,
    'time': time,
	'area' : area,
    'type': tp,
    }, './ml' + file_string
    )
    # surface only
    c.retrieve('reanalysis-era5-complete', {
    'class': cls,
    'date': date_string,
    'expver': expver,
    'grid': grid,
    'levelist': levlistsfc,
    'levtype': levtypesfc,
    'param': paramssfc,
    'stream': stream,
    'time': time,
	'area' : area,
    'type': tp,
    }, './sfc' + file_string
    )


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
    parser = argparse.ArgumentParser(description="Script to download monthly CDS ERA5 files with variables temperature, specific humidity, geopotential height, ln(surf p).")
    parser.add_argument("startdate",    help="Start date in format yyyymm.")
    parser.add_argument("stopdate",     help="Stop date in format yyyymm.")
    args = parser.parse_args()
    return args

#########################################################################################

if __name__ == '__main__':
    args = parse_arguments()
    main(**vars(args))
