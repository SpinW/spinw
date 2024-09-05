# [v4.0.1](https://github.com/spinw/spinw/compare/v4.0.0...v4.0.1)

Bugfix release

# [v4.0.0](https://github.com/spinw/spinw/compare/v3.2.0...v4.0.0)

## New Features

- Add a function to output Mantid MDHistogramWorkspaces (`sw_spec2MDHisto`)
- Add Python plotting of magnetic structure using [vispy](https://vispy.org/)
- Add mex files to compute main loop in `spinwave()` enabling a 2x - 4x speed up depending on system size
- Add python wrapper for all Matlab functions so e.g. `sw_plotspec` etc can be called without the `m.` prefix in pyspinw:q

## Improvements

- Replace `spinwavefast()` method with a new `fastmode` option in the main `spinwave()` method to reduce confusion
- Adds a `neutron_output` option to `spinwave()` to compute only the neutron (`Sperp`) cross-section and not the full spin-spin correlation tensor, saving memory

## Bug Fixes

- Corrects equation for Q-range from `thetaMin` in `sw_instrument` and add `thetaMax` option
- Fixes `sw_issymspec` to recognise powder spectra
- Fixes a parsing error in the `spinw.fourier` method if no sublat option given.
- Fixes several bugs in `sw_plotspec` where it ignores user options in `'auto'` mode, and where it inverts user supplied colormaps.
- Fixes several bugs in `.fitspec()` for handling twins and where it only outputs the final chi^2 values.


# [v3.2.0](https://github.com/spinw/spinw/compare/0.0.1...v3.2.0)

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


# [v3.1.2](https://github.com/spinw/spinw/compare/v3.1.0...v3.1.2)

## Improvements

- Change to preallocation of output energies in `spinwave` to reduce memory usage and improve calculation speed
- Use a mex function (if mex is enabled) for matrix multiplication in `spinwave` with `hermit=false` that reduces memory usage and improves calculation speed for large magnetic cells (in an example with 216 magnetic atoms the execution time was reduced by ~65%)


## Bug Fixes

- Fix generation of lattice from basis vectors in `genlattice`, see issue [#28](https://github.com/SpinW/spinw/issues/28)
- `sortMode` in `spinwave` now correctly sorts the spin wave modes within each twin
- A `spinw` object can now be correctly created from a structure figure
- `.cif` files with a mixture of tabs and spaces or containing a `?` in the comments can now be read correctly
- Rotation matrix `rotC`  in `addtwin` is now required to be a valid rotation or reflection matrix.
- Spin of atom in `addatom` must have `S>=0`.
- Anisotropic g-tensor in `addg` must be physically valid - i.e. :math:`g^\dagger.g` must be a symmetric positive definite matrix.
- Fix bug in addcoupling that did not allow user to supply 'atom' with numeric array of atom indices (previously only worked for string or cell of strings corresponding to atom labels).
- Renamed undocumented `gencoupling` parameter `tol` to `tolMaxDist` (see doc string of `gencoupling` for more details).
- Added validation to `gencoupling` to ensure `maxDistance > dMin`.
- Fixed uncaught error in `gencoupling` by checking if any bonds have length < `maxSym`
- A warning will now be emitted if `saveSabp` is requested in `spinwave` for a commensurate structure
- Fix bug in definition of rotation matrix transforming to spinw coordinate system when left-handed set of basis vectors supplied to `genlattice`, see issue [#57](https://github.com/SpinW/spinw/issues/57)
- Validation added for `perm` and `origin` arguments supplied to `genlattice` (and warn users that these will be ignored if no symmetry/spacegroup is supplied in the same function call).
- Deprecated `spgr` argument to `genlattice` (users should use `sym` instead).
- Fix `MATLAB:nonLogicalConditional` error raised when using multiple k in `genmagstr`  with `helical` mode
- Raise error if invalid shape `S` or `k` is provided to `genmagstr`, previously they would be silently set to zero
- Raise error if wrong number of spins `S` is provided to `genmagstr` in `helical` mode. Previously the structure would be silently initialised to a random structure.
- Raise error if a complex spin `S` is provided to `genmagstr` in `helical` mode. Previously this meant it would silently ignore the `n` option, and behave exactly like `fourier` mode.
- Raise error if `rotate` mode is used without first initialising a magnetic structure
- Emit deprecation warning if the undocumented `extend` mode is used in `genmagstr`
- Raise error if the first spin is parallel to `n` and no rotation angle is provided in `rotate` mode in `genmagstr`. Previously this would silently result in `NaN`
- Raise error if `phi` or `phid` is not real in `rotate` mode in `genmagstr`. This was an undocumented feature which has been removed.
- Emit warning that the spin amplitude will be moderated if components of `S` are parallel to `n` in `helical` mode in `genmagstr`
- Emit warning if  `nExt` is unnecessarily large compared to `k` in `helical` and `fourier` modes in `genmagstr`
- Emit warning if arguments that will be ignored are passed to a particular mode in `genmagstr` (e.g. `S` is passed to `random`)
- Raise error if complex values is provided for `n` in `genmagstr`.  Previously this would've caused a crash.
- Fix error when plotting progress of `optmagsteep` without existing figure
- Correctly report magnetic moments in each iteration of `optmagsteep`.
- Fix errors when calling `intmatrix` with dipolar bonds and symbolic spinw object with fitmode true and false
- Ensure biquadratic exchange interactions are isotropic in `addcoupling` (previously checked in `intmatrix`)
- Raise error if invalid shape `kbase` is provided to `optmagk`, previously it would be silently set to empty
- Ensure varargin is correctly passed through to `ndbase.pso` from `optmagk`. Previously user provided `TolFun`, `TolX` and `MaxIter` would be overwritten by the defaults.
- Warn users that that the results of `spinwave` have not been scientifically validated for supercell structures with an incommensurate modulation.
- Emit warning if wrong length `xmin`, `xmax` or `x0` is passed to `optmagstr`. Previously they would be silently ignored.
- No longer require a magnetic structure be initialised with `genmagstr` before using `optmagstr`. If not intialised, a default `nExt` of `[1 1 1]` is used. This has also been clarified in the docstring.
- Fix bug where powder spectra was not recognised in `sw_plotspec`, introduced by a previous update to provide more helpful error messages.
- `sw_instrument` now calculates the limits for thetaMax, before it was using the continuation of the thetaMin line to high Q which is incorrect.
- Fixes a parsing error in the `spinw.fourier` method if no sublat option given.

