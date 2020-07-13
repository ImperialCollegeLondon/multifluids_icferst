#!/bin/sh -f
rem OutputFile: $2/$1.log
rm "$2/$1.log"
cd ..
rm 1.tmp
echo "Generating output file ......" >>$2/$1.log
"$3/GID_B3D_Linux" "$2/$1.dat" "$3/B3D.par"
if exist 1.tmp goto end
:end 
rm 1.tmp
rm "$2/$1.log"

