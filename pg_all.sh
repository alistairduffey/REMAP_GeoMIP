#!/bin/bash
module load pgi
module load netcdf/4.2.1-pgf90



rm bands_ECHAM6 write_WACCM_p write_WACCM_earth write_WACCM_sun *.o




pgf90 -I${NETCDF}/include -o bands_ECHAM6  bands_ECHAM6.f selib.f -L${NETCDF}/lib -lcdfio -lcdfplus -ltra -lipo -lgm2em -lcalvar -lnetcdf -lnetcdff

./bands_ECHAM6

rm *.o 



pgf90 -I${NETCDF}/include -o write_WACCM_p  write_WACCM_p.f selib.f -L${NETCDF}/lib -lcdfio -lcdfplus -ltra -lipo -lgm2em -lcalvar -lnetcdf -lnetcdff

./write_WACCM_p

rm *.o 



pgf90 -I${NETCDF}/include -o write_WACCM_earth  write_WACCM_earth.f selib.f -L${NETCDF}/lib -lcdfio -lcdfplus -ltra -lipo -lgm2em -lcalvar -lnetcdf -lnetcdff

./write_WACCM_earth

rm *.o 



pgf90 -I${NETCDF}/include -o write_WACCM_sun  write_WACCM_sun.f selib.f -L${NETCDF}/lib -lcdfio -lcdfplus -ltra -lipo -lgm2em -lcalvar -lnetcdf -lnetcdff

./write_WACCM_sun

rm *.o 


echo
echo "The full Mie calculation is completed, output files are ready!"

