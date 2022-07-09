import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from os.path import dirname as dirname

xaxis = np.arange(-1,11,0.05)
yaxis = np.sin(xaxis+np.pi/2) + 1
elipsesize = [0.1*16, 0.1*4]
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(xaxis, yaxis, color = 'mediumblue')
ax.set_ylabel(r'$P_{\mathrm{out}}$ [a.u.]', usetex = True, fontsize = 11)
ax.set_xlabel(r'$\Delta L$ [a.u.]', usetex = True, fontsize = 11)
ax.add_patch(Ellipse((2*np.pi,2), *elipsesize, color='r', fill=False))
ax.add_patch(Ellipse((np.pi,0), *elipsesize, color='r', fill=False))
ax.set_xlim([0, 8])
# powerpoint is 13.333 inches wide by 7.5 inches high
# a4 paper is 8.25 by 11.75 in
# geometry package boundary -0.875in*2 horizontally
# figures set to 0.85 * textwidth
fig.set_size_inches(0.8*(8.25-0.875*2), 2) 
fig.tight_layout()
# plt.show()
dir = dirname(__file__)
fig.savefig(__file__[:-3] + '.png', format= 'png', dpi = 300, pad_inches = 0.01)
