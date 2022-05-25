`Unreleased <https://github.com/SpinW/spinw/compare/v3.1.2...HEAD>`_
--------------------------------------------------------------------

New Features
############

Improvements
############
- Change to preallocation of output energies in ``spinwave`` to reduce
  memory usage and improve calculation speed

Bug Fixes
#########
- Fix generation of lattice from basis vectors in ``genlattice``, see issue
  `#28 <https://github.com/SpinW/spinw/issues/28>`_
- ``sortMode`` in ``spinwave`` now correctly sorts the spin wave modes
  within each twin
- A ``spinw`` object can now be correctly created from a structure figure
- ``.cif`` files with a mixture of tabs and spaces or containing a ``?``
  in the comments can now be read correctly
