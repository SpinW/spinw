# [v3.2.0](https://github.com/pace-neutrons/libpymcr/compare/0.0.1...v3.2.0)

## Initial public beta of PySpinW

This is an initial public beta version of PySpinW released on PyPI.

Please install using:

```bash
pip install spinw
```

This will install a module called `pyspinw` (note the `py` at the start).

You can then run SpinW with:

```python
import numpy as np
import matplotlib.pyplot as plt
from pyspinw import Matlab
m = Matlab()
swobj = m.spinw()
swobj.genlattice('lat_const', [3, 3, 6], 'angled', [90, 90, 120], 'sym', 'P 1');
swobj.addatom('r', [0, 0, 0], 'S', 1/2, 'label', 'MCu2')
swobj.gencoupling('maxDistance', 5)
swobj.addmatrix('label', 'J1', 'value', 1.00, 'color', 'g')
swobj.addcoupling('mat', 'J1', 'bond', 1)
swobj.genmagstr('mode', 'helical', 'k', [-1/3, -1/3, 0], 'n',[0, 0, 1], 'unit', 'lu', 'S', [[1], [0], [0]])
spec = swobj.spinwave([[-1/2, 0, 0], [0, 0, 0], [1/2, 1/2, 0], 100], 'hermit', False)
spec = m.sw_egrid(spec, 'component', 'Sxx+Syy', 'imagChk', False, 'Evect', np.linspace(0, 3, 100))
ax = plt.imshow(np.real(np.flipud(spec['swConv'])), aspect='auto', vmax=1)
plt.show()
```

On Windows and Linux systems, as long as you're running PySpinW locally, Matlab plotting commands like `m.plot(swobj)` will work. This is not the case on MacOS (a known bug) and on remote systems (e.g. via JupyterHub).


# [v0.0.1](https://github.com/spinw/spinw/compare/v3.1.2...0.0.1)

## pySpinW

This is an initial release of pySpinW as a `pip` installable wheel for python >= 3.8 and MATLAB >= R2021a

### Installation

Please install with

```bash
pip install pyspinw*.whl
```

This package can now be used in python if you have a version of MATLAB or MCR available on the machine. 
The package will try to automatically detect your installation, however if it is in a non-standard location, the path and version will have to be specified.

```python
from pyspinw import Matlab
m = Matlab(matlab_version='R2023a', matlab_path='/usr/local/MATLAB/R2023a/')
```

### Example

An example would be:

```python
import numpy as np
from pyspinw import Matlab

m = Matlab()

# Create a spinw model, in this case a triangular antiferromagnet
s = m.sw_model('triAF', 1)

# Specify the start and end points of the q grid and the number of points
q_start = [0, 0, 0]
q_end = [1, 1, 0]
pts = 501

# Calculate the spin wave spectrum
spec = m.spinwave(s, [q_start, q_end, pts])
```

### Known limitations

At the moment graphics will not work on macOS systems and is disabled.
