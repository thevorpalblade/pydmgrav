# pydmgrav
This library is for reading and analysing the IGETS data files (superconducting gravimeter records) from:
http://isdc.gfz-potsdam.de/igets-data-base/

The library depends on matplotlib, numpy, scipy, and ipyparallel for parallelization.

The dependency on ipyparallel can be easily lifted by removing the appropriate lines and
replacing the parallel call to map() with the native python call to map()

If using ipyparallel, it is important to first start a ipyparallel cluster with:
$ ipcluster start
on the command line outside of python.

