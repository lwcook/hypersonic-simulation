#!/bin/bash

touch ./doc/source/index.rst

sphinx-build -b html ./doc/source ./doc/build

#cp -r ./doc/build/. ~/aerotools/hypersonicsimulation/documentation/
scp -r ./doc/build/. edclinux:/export/www/edc/aerotools/hypersonicsimulation/documentation
