import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (InsetPosition, inset_axes,
                                                  mark_inset)

f = np.load("level2_fft_nobaseline.npz")
freqs = f['arr_0']
intensities = f['arr_1']

fig, ax1 = plt.subplots()
ax1.plot(freqs, intensities)

ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel('Acceleration/Hz (nm/sec)')
# ax1.legend(loc=0)

ax1.annotate(
    '',
    xy=(.0003, 3200),
    xytext=(.0003, 6000),
    arrowprops=dict(width=1, headwidth=6, facecolor='black'),
)
# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
left, bottom, width, height = [0.3, 0.3, 0.4, 0.4]
ax2 = fig.add_axes([left, bottom, width, height])
## Manually set the position and relative size of the inset axes within ax1
##ip = InsetPosition(ax1, [0.4, 0.2, 0.5, 0.5])
##ax2.set_axes_locator(ip)
## Mark the region corresponding to the inset axes on ax1 and draw lines
## in grey linking the two axes.
##mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

ax2.plot(freqs[1000:2500], intensities[1000:2500])
ax2.annotate(
    '',
    xy=(.0003, 3200),
    xytext=(.0003, 4000),
    arrowprops=dict(facecolor='black', width=1, headwidth=6, shrink=0.05),
)

        # Some ad hoc tweaks.
#ax1.set_ylim(0, 26)
#ax2.set_yticks(np.arange(0, 2, 0.4))
#ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
#ax2.tick_params(axis='x', which='major', pad=8)

plt.show()
