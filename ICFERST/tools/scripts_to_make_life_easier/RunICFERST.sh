#!/bin/bash

myvar=$1
chrlen=${#myvar}

if [ "$chrlen" -ge 4 ]; then
    cp /<ABSOLUTE_PATH_TO_INSTALLATION_FOLDER>/bin/icferst ./  && ./icferst ./$myvar
else
    cp /<ABSOLUTE_PATH_TO_INSTALLATION_FOLDER>/bin/icferst ./  && ./icferst ./*.mpml
fi






