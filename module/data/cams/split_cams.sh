# Bash script to split 1-month CAMS grib files into 6-hour files.
# Call
# bash split_cams.sh cams73_latest_co2_conc_surface_ 200901 201912
# Input
# arg1 - skeleton filename
# arg2 - startdate (yyyy or yyyymm)
# arg3 - enddate (yyyy or yyyymm)

filetrunk=$1
startdate=$2
enddate=$3

# yyyy input format
if [ ${#startdate} = 4 ]; then
    startyear=$startdate
    startmonth=01
    endyear=$enddate
    endmonth=12
# yyyymm input format
elif [ ${#startdate} = 6 ]; then
    startyear=${startdate:0:4}
    startmonth=${startdate:4:2}
    endyear=${enddate:0:4}
    endmonth=${enddate:4:2}
# wrong input format
else
    echo "Wrong date format. Either yyyy or yyyymm."
fi

# split loop
for ((year=startyear; year<=endyear; year++)); do
    tmpstartmonth=01
    tmpendmonth=12
    if [ $year = $startyear ]; then
        tmpstartmonth=${startmonth#0}
    fi
    if [ $year = $endyear ]; then
        tmpendmonth=${endmonth#0}
    fi
    for ((tmpmonth=tmpstartmonth; tmpmonth<=tmpendmonth; tmpmonth++)); do
        printf -v month "%02d" $tmpmonth
        # split to days        
        echo "Splitting "${filetrunk}${year}${month}".grb into days (takes ca. 8min/month)..."
        cdo splitday ./${filetrunk}${year}${month}.grb ./${filetrunk}${year}${month}
        
        # remove 1month files
        echo "Removing "${filetrunk}${year}${month}".grb ..."
        #rm ./${filetrunk}${year}${month}.grb
        
        # split to hours
        #echo "Splitting rea"${year}${month}"xx.grb into hours and removing 1-day files..."
        #for file in ./rea${year}${month}*.grb; do
        #    cdo splithour $file ${file:0:13}
        #    # remove 1day files
        #   rm $file
        #done
    done
done
