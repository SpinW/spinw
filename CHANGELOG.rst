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
- Rotation matrix ``rotC``  in ``addtwin`` is now required to be a valid
  rotation or reflection matrix.
- Spin of atom in ``addatom`` must have ``S>=0``.
- Anisotropic g-tensor in ``addg`` must be physically valid - i.e.
  :math:`g^\dagger.g` must be a symmetric positive definite matrix.
- Fix bug in addcoupling that did not allow user to supply 'atom' with
  numeric array of atom indices (previously only worked for string or
  cell of strings corresponding to atom labels).
- Renamed undocumented ``gencoupling`` parameter ``tol`` to ``tolMaxDist``
  (see doc string of ``gencoupling`` for more details).
- Added validation to ``gencoupling`` to ensure ``maxDistance > dMin``.
- Fixed uncaught error in ``gencoupling`` by checking if any bonds have
  length < ``maxSym``
- A warning will now be emitted if ``saveSabp`` is requested in ``spinwave``
  for a commensurate structure
- Fix bug in definition of rotation matrix transforming to spinw coordinate system when left-handed set of
  basis vectors supplied to ``genlattice``, see issue `#57 <https://github.com/SpinW/spinw/issues/57>`_
- Validation added for ``perm` and ``origin`` arguments supplied to ``genlattice`` (and warn users that these will be
  ignored if no symmetry/spacegroup is supplied in the same function call).
- Deprecated ``spgr`` argument to ``genlattice`` (users should use ``sym`` instead).
- Fix ``MATLAB:nonLogicalConditional`` error raised when using multiple
  k in ``genmagstr``  with ``helical`` mode
- Raise error if invalid shape ``S`` or ``k`` is provided to ``genmagstr``,
  previously they would be silently set to zero
- Raise error if wrong number of spins ``S`` is provided to ``genmagstr`` in
  ``helical`` mode. Previously the structure would be silently initialised
  to a random structure.
- Raise error if a complex spin ``S`` is provided to ``genmagstr`` in
  ``helical`` mode. Previously this meant it would silently ignore the
  ``n`` option, and behave exactly like ``fourier`` mode.
- Raise error if ``rotate`` mode is used without first initialising
  a magnetic structure
- Emit deprecation warning if the undocumented ``extend`` mode is used
  in ``genmagstr``
- Raise error if the first spin is parallel to ``n`` and no rotation
  angle is provided in ``rotate`` mode in ``genmagstr``. Previously
  this would silently result in ``NaN``
- Raise error if ``phi`` or ``phid`` is not real in ``rotate`` mode in
  ``genmagstr``. This was an undocumented feature which has been removed.
- Emit warning that the spin amplitude will be moderated if components
  of ``S`` are parallel to ``n`` in ``helical`` mode in ``genmagstr``
- Emit warning if  ``nExt`` is unnecessarily large compared to ``k`` in
  ``helical`` and ``fourier`` modes in ``genmagstr``
- Emit warning if arguments that will be ignored are passed to a particular
  mode in ``genmagstr`` (e.g. ``S`` is passed to ``random``)
- Raise error if complex values is provided for ``n`` in ``genmagstr``.
  Previously this would've caused a crash.
- Fix error when plotting progress of ``optmagsteep`` without existing figure
- Correctly report magnetic moments in each iteration of ``optmagsteep``.
