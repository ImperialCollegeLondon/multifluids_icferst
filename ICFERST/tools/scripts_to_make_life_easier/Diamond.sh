#!/bin/bash

myvar=$1
chrlen=${#myvar}

if [ "$chrlen" -ge 4 ]; then
    diamond -s /<ABSOLUTE_PATH_TO_INSTALLATION_FOLDER>/ICFERST/schemas/multiphase.rng ./$myvar
else
    diamond -s /<ABSOLUTE_PATH_TO_INSTALLATION_FOLDER>/ICFERST/schemas/multiphase.rng ./*.mpml
fi
