# Bash script to split 1-month ERA5 GRIB2 files into 1-hour files.
# Input
# arg1 - startdate (yyyy or yyyymm)
# arg2 - enddate (yyyy or yyyymm)
# February 2020
# Markus Haun

startdate=$1
enddate=$2

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
        echo "Splitting ml"${year}${month}".grb, sfc"${year}${month}".grb into days (takes ca. 8min/month)..."
        cdo splitday ./ml${year}${month}.grb ./ml${year}${month}
        cdo splitday ./sfc${year}${month}.grb ./sfc${year}${month}

        # remove 1month files
        echo "Removing ml"${year}${month}".grb, sfc"${year}${month}".grb..."
        rm ./ml${year}${month}.grb ./sfc${year}${month}.grb

        # split to hours
        echo "Splitting ml"${year}${month}"xx.grb, sfc"${year}${month}"xx.grb into hours and removing 1-day files..."
        for file in ./ml${year}${month}*.grb; do
            cdo splithour $file ${file:0:12}
            # remove 1day files
            rm $file
        done
        for file in ./sfc${year}${month}*.grb; do
	    echo $file
            cdo splithour $file ${file:0:13}
            # remove 1day files
            rm $file
        done
    done
done
