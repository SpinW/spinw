`Unreleased <https://github.com/SpinW/spinw/compare/v3.1.2...HEAD>`_
----------

- Bug fixes:

  - Fix generation of lattice from basis vectors in ``genlattice``, see issue
    #28
  - ``sortMode`` in ``spinwave`` now correctly sorts the spin wave modes
    within each twin
  - A ``spinw`` object can now be correctly created from a structure figure
  - ``.cif`` files with a mixture of tabs and spaces or containing a ``?``
    in the comments can now be read correctly