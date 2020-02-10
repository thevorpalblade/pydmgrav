# pydmgrav
This library is for reading and analysing the IGETS data files (superconducting gravimeter records) from:
http://isdc.gfz-potsdam.de/igets-data-base/

The library depends on matplotlib, numpy, and scipy.

To get started quickly, try (in ipython or jupyter or similar)

import pydmgrav
freqs, ASD, ASD_err = pydmgrav.main_raw(basedir="/path/to/igets/data")



