# This example shell script calls a Python script and passes in all arguments
# To locate your interpreter or compiler, run "whereis [interpreter]"
# Edit as necessary for different interpreters and/or file locations
# Note: "$@" is functionally-equivalent to "$1 $2 $3 ..."
# Come to hours or post on Piazza if you have questions about scripting!


# arguments: 1 - reads file, 2 - vector file, 3 - k for contamination, 4 - output file, 5 - k for correction, 
# 6 - t, 7 - d, 8 - output file 2
python3 contamination.py $1 $2 $3 $4
python3 correction.py $4 $5 $6 $7 $8
python3 debruijn.py $8