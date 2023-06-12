# pySpinW

This is an intial release of pySpinW as a `pip` installable wheel for python >= 3.8 and MATLAB >= R2021a

## Installation 

Please install with

```bash
pip install pyspinw*.whl
```

This package can now be used in python if you have a version of MATLAB or MCR available on the machine. 
The package will try to automatically detect your installation, however if it is in a non-standard location, the path and version will have to be specified.

```python
m = Matlab(matlab_version='R2023a', matlab_path='/usr/local/MATLAB/R2023a/')
```

## Example

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

## Known limitations

At the moment graphics will not work on macOS systems and is disabled.
