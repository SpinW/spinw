__author__ = 'github.com/wardsimon'
__version__ = '0.0.1'

import numpy as np
from pyspinw import Matlab

try:
    from matplotlib import pyplot as plt
except ImportError:
    plt = None

m = Matlab()
# An example of specifying the MATLAB version and path
# m = Matlab(matlab_version='R2023a', matlab_path='/usr/local/MATLAB/R2023a/')

# Suppress output to make it less verbose in CI output
m.swpref().fid = 0

# Create a spinw model, in this case a triangular antiferromagnet
s = m.sw_model('triAF', 1)
print(s)

# Specify the start and end points of the q grid and the number of points
q_start = [0, 0, 0]
q_end = [1, 1, 0]
pts = 501

# Calculate the spin wave spectrum, apply an energy grid and convolute with a Gaussian
spec = m.sw_egrid(m.spinwave(s, [q_start, q_end, pts]))
spec2 = m.sw_instrument(spec, dE=0.3)

# Plot the result if matplotlib is available
if plt is not None:
    ax = plt.imshow(np.flipud(spec2['swConv']),
                    aspect='auto',
                    extent=[q_start[-1], q_end[0], spec2["Evect"][0][0], spec2["Evect"][0][-1]])
    ax.set_clim(0, 0.15)
    plt.xlabel('Q [q, q, 0] (r.l.u)')
    plt.ylabel('Energy (meV)')
    plt.title('Spectra of a triangular antiferromagnet')
    plt.savefig('pyspinw.png')
    plt.show()
