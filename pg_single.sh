#!/bin/bash
module load pgi
module load netcdf/4.2.1-pgf90

#set name=$1
name=${1}


# DO NOT DELETE THIS IF-LOOP, DATA LOSS MAY OCCUR
if [[ ${name} == *".f"* ]]; then
	echo "Do not input source files, otherwise you will IRREVOCABLY LOSE them. Instead enter the name without the file extension (.f), i.e. the resulting executable."
	echo "(That was close, man!)"
	
	kill -INT $$
else
	rm ${name}
fi


pgf90 -I${NETCDF}/include -o ${name}  ${name}.f selib.f -L${NETCDF}/lib -lcdfio -lcdfplus -ltra -lipo -lgm2em -lcalvar -lnetcdf -lnetcdff

./${name}

rm *.o 

