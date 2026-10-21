###### This file contains the instructions to use the Mie tool of REMAP ######





### What does it do? ###

The Mie tool reads data produced by WACCM at NCAR for the CCMI experiments and translates aerosol properties onto a different, specified spectral grid. This produces fields of aerosol variables typically ingestable by climate models. With this functionality, the Mie tool is supposed to serve as a data preparation tool for modelling groups to perform CCMI experiments that require a uniform aerosol forcing. The code is based on the REMAP method by Jörimann et al. (2026).





### How is it run? ###

SIMPLE VERSION:

1. Change the bands.dat file, so it contains your desired band number, output wavelength band upper and lower boundaries (in microns), as well as the number of both short-wave and long-wave bands (first line).

2. Source the bash-script with the compilation commands with the following command:

source pg_all.sh

3. Wait up to a few hours. There should be continuous output to the command line and at the end a success message.

4. The newly produced ccmi_WACCM_to_*.nc files contain your output data. Check them thorougly!





MORE IN-DEPTH VERSION

1. Change the bands.dat file, so it contains your desired band number, output wavelength band upper and lower boundaries (in microns), as well as the number of both short-wave and long-wave bands (first line). The layouting is not that strict, it works as long as you start a new line for every band and separate all numbers by at least one space.

2. Source the bash-script with the compilation commands with the following command:

source pg_all.sh

This requires the pgf90 fortran compiler and several libraries to be installed (see compilation command). Another compiler like ifort can be used with an adapted compilation command, of course. If none of that works, you could just run the executables, but, of course, you won't be able to change any of the hard-coded variables and file names. If you just want to change the bands, this should still work.

There is also a file that lets you compile and run the single FORTRAN scripts: pg_single.sh. In the following part, I will lay out the different file contents and how they depend on each other.



------------------

--- Source code ---

All the .f files contain FORTRAN source code.

bands_ECHAM6.f: This file contains a program to calculate a lookup table that provides extinction coefficients as a function of environmental variables (temperature, weight percentage of H2O,...). It must be compiled with selib.f, which it depends on. Line 89 creates the file, where the lookup table is be written into. You could theoretically change the name of the file, if you wanted to. If you do, remember to also change the name of the input file in selib.f (see below).

write_WACCM_p/sun/earth.f: These files contain programs to write the pressure, surface area density & mean radius and short-wave and long-wave optical properties (extinction coefficients, single-scattering albedo, asymmetry factor) respectively. They all also depend on selib.f. You could change the netCDF attributes and the file names of the final output, if you wanted to.

selib.f: This file contains all the subroutines and functions called by the aforementioned files. It is split in two parts; the first part is exclusively used by the bands_ECHAM6.f code and the second part is exclusively used by the write_WACCM_*.f code. While the other source code files are neatly formatted, commented out and do not contain any implicit variables, this file is the Wild West, so to speak. There should only ever be one line that you may want to edit, which is line 69, where the name of the input file (lookup_table.dat) for the subroutine "lookup_radiation" is specified. This would typically be the case, if you changed the name of the output file created by bands_ECHAM6.f, which is the same file you read here. So the two are interdependent and must be identical if you run all the programs consecutively and don't change the name of the file manually in between.

There are two flow diagrams that detail all of the code units and how they go together. This should give an overview over the inner workings of the two archetypal programs - bands and write. Calls are always being made from the higher to the lower connected function in the diagrams. The subroutines are also colored according to their "rank", as in first-, second-, third-,... order-dependent from the main program. Enjoy the beautiful color scheme.


--- Resources ---

There are two repositories:

data_base: The files in this directory are accessed by various subroutines within selib.f. They contain parameters for Mie calculations and are better left alone.

WACCM_data: As the name suggests, the base data for the aerosol properties live here. They represent variables in a geoengineering scenario and were derived in the Geoengineering Large Ensemble Project (GLENS). The experiment is based on a feeback control mechanism that maintains 2020-2030 decadal mean global surface temperatures by injecting SO2 into the stratosphere at four distinct location, essentially forming an aerosol layer that cools the Earth's surface. In principle, other input files could be used, provided that the grids and time axis stayed the same.



--- Output files ---

lookup_table.dat: This file contains a 4D table of extinction coefficients.

ccmi_WACCM_to_*.nc: These files contain the derived optical properties on the specified spectral grid. They are meant to serve as input data for climate models and should be checked for validity and completeness before being used.

-----------------



3. Wait up to a few hours. There should be continuous output to the command line and at the end a success message. The output of the lookup table calculation gives the four parameters that the optical properties depend on for each cycle of the calculation, i.e. the values of the four dimensions of the table, for which the optical properties are currently being calculated. For the write_* programs, the output shows example values for each month (outermost loop iteration). If any instabilities occur, you can see it in these values.


4. The newly produced ccmi_WACCM_to_*.nc files contain your output data. Check them thorougly! Usually, I use ncview for a quick visual check and panoply for quantitative, but still quick analysis.





### What could go wrong? ###

Nothing, of course...

But just in case, a few comments: If you don't change anything in the source code or the layout of the input file bands.dat, there should be no problems. The only obvious thing that we've encountered so far is the netCDF file limit of ~4 GB. If you use more than 20 bands for a single program (sun/earth), it might fail, because the supposed size of the file that the netCDF library will try to create is too large. In that case, it is easiest to split the output into the single variables, e.g. make one file with extinction coefficients and single-scattering albedo and one file with the asymmetry factor. You can still combine them after, especially if you change the netCDF standard to 4, where a 64-bit limit applies.

Have fun!










Much of the original code was produced by Beiping Luo at ETH Zürich.

READ_ME, flow diagrams, and code polishing, reformatting & comments by Andrin Jörimann, ETH Zürich





