import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (InsetPosition, inset_axes,
                                                  mark_inset)

l2 = np.load("level2_final.npz")
l3 = np.load("level3_final.npz")
l2_freqs, l2_accs = l2['arr_0']
l3_freqs, l3_accs = l3['arr_0']


annotation_x = 1/3300

fig1, ax11 = plt.subplots()
ax11.plot(l2_freqs, l2_accs)

ax11.set_xlabel("Frequency (Hz)")
ax11.set_ylabel('Acceleration/Hz (nm/sec)')
# ax1.legend(loc=0)

ax11.annotate(
    '',
    xy=(annotation_x, 3200),
    xytext=(annotation_x, 5000),
    arrowprops=dict(width=1, headwidth=6, facecolor='black'),
)
# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
left, bottom, width, height = [0.4, 0.4, 0.4, 0.4]
ax12 = fig1.add_axes([left, bottom, width, height])
## Manually set the position and relative size of the inset axes within ax1
##ip = InsetPosition(ax1, [0.4, 0.2, 0.5, 0.5])
##ax2.set_axes_locator(ip)
## Mark the region corresponding to the inset axes on ax1 and draw lines
## in grey linking the two axes.
##mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

ax12.plot(l2_freqs[1000:2500], l2_accs[1000:2500])
ax12.annotate(
    '',
    xy=(annotation_x, 3200),
    xytext=(annotation_x, 3600),
    arrowprops=dict(facecolor='black', width=1, headwidth=6, shrink=0.05),
)

        # Some ad hoc tweaks.
#ax1.set_ylim(0, 26)
#ax2.set_yticks(np.arange(0, 2, 0.4))
#ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
#ax2.tick_params(axis='x', which='major', pad=8)

fig2, ax21 = plt.subplots()
ax21.plot(l3_freqs, l3_accs)
ax21.set_ylim((0, 400))

ax21.set_xlabel("Frequency (Hz)")
ax21.set_ylabel('Acceleration/Hz (nm/sec)')
# ax1.legend(loc=0)

ax21.annotate(
    '',
    xy=(annotation_x, 170),
    xytext=(annotation_x, 300),
    arrowprops=dict(width=1, headwidth=6, facecolor='black'),
)
# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
left, bottom, width, height = [0.4, 0.4, 0.4, 0.4]
ax22 = fig2.add_axes([left, bottom, width, height])
## Manually set the position and relative size of the inset axes within ax1
##ip = InsetPosition(ax1, [0.4, 0.2, 0.5, 0.5])
##ax2.set_axes_locator(ip)
## Mark the region corresponding to the inset axes on ax1 and draw lines
## in grey linking the two axes.
##mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

ax22.plot(l3_freqs[1000:2500], l3_accs[1000:2500])
ax22.annotate(
    '',
    xy=(annotation_x, 160),
    xytext=(annotation_x, 300),
    arrowprops=dict(facecolor='black', width=1, headwidth=6, shrink=0.05),
)

plt.show()
