#!/bin/csh
set echo
module load intel impi netcdf-c netcdf-f
cp intel.mk Makefile.fms_spectral_solo Makefile obj
cd obj
make
cd ../../src/mppnccombine
./mppnccombine_compile
