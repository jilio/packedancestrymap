#!/bin/bash

command -v curl >/dev/null 2>&1 || { 
    echo >&2 "I require curl but it's not installed. Aborting.";
    exit 1;
}

mkdir -p test_data
curl -L https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/Harney_et_al_2019_ancient_genotypes_0.zip | tar -xf - -C test_data
curl -L https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/Data_S1_EA_aDNA.tar.gz | tar -xf - -C test_data
curl -L https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/Olalde_et_al_genotypes.zip | tar -xf - -C test_data

rm -rf test_data/__MACOSX