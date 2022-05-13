classdef spinw < handle & matlab.mixin.SetGet
    % class to store and solve magnetic Hamiltonians
    %
    % ### Syntax
    %
    % `obj = spinw`
    %
    % `obj = spinw(obj)`
    %
    % `obj = spinw(source)`
    %
    % `obj = spinw(figure_handle)`
    %
    % ### Description
    %
    % `obj = spinw` constructs an empty spinw class object.
    %
    % `obj = spinw(obj)` constructs a spinw class object from the
    % parameters defined in `obj`. If `obj` is spinw class, it only checks
    % the integrity of its internal data structure. If `obj` is struct
    % type, it creates new spinw object and checks data integrity.
    %
    % `obj = spinw(source)` construct new spinw class object, where
    % `source` is either a file path pointing to a local `cif` or `fst`
    % file or a link to an online file.
    %
    % `obj = spinw(figure_handle)` copy the spinw object stored in a
    % previous structural 3D plot figure, referenced by `figure_handle`.
    %
    %
    % The data structure within the spinw object can be accessed by using
    % [spinw.struct] method. All fields of the struct type data behind the
    % spinw object are accessible through the main field names of the `obj`
    % object. For example the lattice parameters can be accessed using:
    %
    % ```
    % abc = obj.unit_cell.lat_const
    % ```
    %
    % spinw is a handle class, which means that only the handle of the
    % object is copied in an assinment command `swobj1 = swobj2`. To create
    % a copy (clone) of an spinw object use:
    %
    % ```
    % swobj1 = swobj2.copy
    % ```
    %
    % ### Properties
    %
    % The data within the `spinw` object is organized into a tree structure
    % with the main groups and the type of data they store are the
    % following:
    %
    % * [spinw.lattice] unit cell parameters
    % * [spinw.unit_cell] atoms in the crystallographic unit cell
    % * [spinw.twin] crystal twin parameters
    % * [spinw.matrix] 3x3 matrices for using them in the Hailtonian
    % * [spinw.single_ion] single ion terms of the Hamiltonian
    % * [spinw.coupling] list of bonds
    % * [spinw.mag_str] magnetic structure
    % * [spinw.unit] physical units for the Hamiltonian
    % * [spinw.cache] temporary values
    %
    % ### Methods
    %
    % Methods are the different commands that require a `spinw` object as a
    % first input, thus they can be called as `method1(obj,...)`,
    % alternatively the equivalent command is `obj.method1(...)`. The list
    % of public methods is below.
    %
    % #### Lattice operations
    %
    %   spinw.genlattice
    %   spinw.basisvector
    %   spinw.rl
    %   spinw.nosym
    %   spinw.newcell
    %   spinw.addatom
    %   spinw.unitcell
    %   spinw.abc
    %   spinw.atom
    %   spinw.matom
    %   spinw.natom
    %   spinw.formula
    %   spinw.disp
    %   spinw.symmetry
    %   
    % #### Plotting
    %
    %   spinw.plot
    %
    % #### Crystallographic twin operations
    %
    %   spinw.addtwin
    %   spinw.twinq
    %   spinw.notwin
    %
    % #### Magnetic structure operations
    %
    %   spinw.genmagstr
    %   spinw.magstr
    %   spinw.magtable
    %   spinw.nmagext
    %   spinw.optmagstr
    %   spinw.optmagk
    %   spinw.optmagsteep
    %   spinw.anneal
    %   spinw.annealloop
    %   spinw.structfact
    %   
    % #### Matrix operations
    %
    %   spinw.addmatrix
    %   spinw.getmatrix
    %   spinw.setmatrix
    %   
    % #### Spin Hamiltonian generations
    %
    %   spinw.quickham
    %   spinw.gencoupling
    %   spinw.addcoupling
    %   spinw.addaniso
    %   spinw.addg
    %   spinw.field
    %   spinw.temperature
    %   spinw.intmatrix
    %   spinw.symop
    %   spinw.setunit
    %   
    % #### Solvers
    %
    %   spinw.spinwave
    %   spinw.powspec
    %   spinw.energy
    %   spinw.moment
    %   spinw.spinwavesym
    %   spinw.symbolic
    %   spinw.fourier
    %   spinw.fouriersym
    %
    % #### Fitting spin wave spectrum
    %
    %   spinw.fitspec
    %   spinw.matparser
    %   spinw.horace
    %   
    % #### Miscellaneous
    %
    %   spinw.copy
    %   spinw.export
    %   spinw.table
    %   spinw.validate
    %   spinw.version
    %   spinw.struct
    %   spinw.clearcache
    %   spinw.spinw
    %
    % ### See also
    %
    % [spinw.copy], [spinw.struct], [Comparing handle and value classes](https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html)
    %

    properties (SetObservable)
        % stores the unit cell parameters
        %
        % ### Sub fields
        %
        % `lat_const`
        % : Lattice constants in a $[1\times 3]$ vector in units defined in
        %   [spinw.unit] (default value is \\ang).
        %
        % `angle`
        % : `[\\alpha,\\beta,\\gamma]` angles in a $[1\times 3]$ vector in
        %   radian units.
        %
        % `sym`
        % : Symmetry operators stored in matrix with dimensions of
        %   $[3\times 4 \times n_{op}]$.
        %
        % `origin`
        % : Origin of the cell in lattice units.
        % 
        % `label`
        % : Label of the space group.
        %
        lattice
        % stores the atoms in the crystallographic unit cell
        %
        % ### Sub fields
        %
        % `r`
        % : Positions of the atoms in the unit cell, stored in a
        %   matrix with dimensions of $[3\times n_{atom}]$, values are
        %   in lattice units.
        %
        % `S`
        % : Spin quantum number of the atoms, stored in a row vector with
        %   $n_{atom}$ number of elements, non-magnetic atoms have `S=0`.
        %
        % `label`
        % : Label of the atom, strings stored in a $[1\times n_{atom}]$
        %   cell.
        %
        % `color`
        % : Color of the atom stored in a matrix with dimensions of $[3\times n_{atom}]$, where every
        %   column defines an RGB color with values between 0 and 255.
        %
        % `ox`
        % : Oxidation number of the atom, stored in a $[1\times n_{atom}]$
        %   matrix.
        %
        % `occ`
        % : Site occupancy in a $[1\times n_{atom}]$ matrix.
        %
        % `b`
        % : Scattering length of the atoms for neutron and x-ray
        %   stored in a $[2\times n_{atom}]$ matrix, first row is neutron,
        %   second row is for x-ray.
        %
        % `ff`
        % : Form factor of the site stored in a $[2\times 9\times
        %   n_{atom}]$ matrix, first row is the magnetic form factor for
        %   neutrons, the second row is the charge form factor for x-ray
        %   cross section.
        %
        % `Z`
        % : Atomic number in a row vector.
        %
        % `A`
        % : Atomic mass (N+Z) for isotopes and -1 for natural
        %   distribution of isotopes stored in a row vector.
        %
        % `biso`
        % : Isotropic displacement factors in units of \\ang$^2$.
        %   Definition is the same as in
        %   [FullProf](https://www.ill.eu/sites/fullprof/), defining the
        %   Debye-Waller factor as $W(d) = 1/8*b_{iso}/d^2$ which is
        %   included in the structure factor as $\exp(-2W(d))$.
        %
        unit_cell
        % stores the crystal twin parameters
        %
        % ### Sub fields
        %
        % `rotc`
        % : Rotation matrices in the $xyz$ coordinate system for
        %   every twin, stored in a matrix with dimensions of $[3\times
        %   3\times n_{twin}]$.
        %
        % `vol`
        % : Volume ratio of the different twins, stored in a
        %    row vector with $n_{twin}$ elements.
        %
        twin
        % stores 3x3 matrices for using them in the Hailtonian
        %
        % ### Sub fields
        %
        % `mat`
        % : Stores the actual values of 3x3 matrices, in a matrix with
        % dimensions of $[3\times 3\times n_{matrix}]$, if assigned for a 
        % bond, the unit of energy is stored in [spinw.unit] (default value 
        % is meV).
        %
        % `color`
        % : Color assigned for every matrix, stored in a
        %   matrix with dimensions of $[3\times n_{matrix}]$, with each
        %   column defining an RGB value.
        %
        % `label`
        % : Label for every matrix, stored as string in a cell with
        %   dimensions of $[1\times n_{matrix}]$.
        %
        matrix
        % stores single ion terms of the Hamiltonian
        %
        % ### Sub fields
        %
        % `aniso`
        % : Row vector that contains $n_{magatom}$ integers, each integer
        %   assignes one of the matrices from the [spinw.matrix] property
        %   to a magnetic atom in the generated [spinw.matom] list as a single
        %   ion anisotropy. Zero value of `aniso` means no single ion
        %   anisotropy for the corresponding magnetic atom.
        %
        % `g`
        % : Row vector with $n_{magatom}$ integers, each integer
        %   assignes one of the matrices from the [spinw.matrix] property
        %   to a magnetic atom in the spinw.matom list as a
        %   g-tensor. Zero value of `g` means a default g-value of 2 for
        %   the corresponding atoms.
        %
        % `field`
        % : External magnetic field stored in a row vector with 3 elements,
        %   unit is defined in [spinw.unit] (default unit is Tesla).
        %
        % `T`
        % : Temperature, scalar, unit is defined in [spinw.unit] (default
        %   unit is Kelvin).
        %
        single_ion
        % stores the list of bonds
        %
        % ### Sub fields
        %
        % `dl`
        % : Distance between the unit cells of two interacting
        %   spins, stored in a $[3\times n_{coupling}]$ matrix.
        %
        % `atom1`
        % : First magnetic atom, pointing to the list of
        %   magnetic atoms in the list generated by `spinw.matom`, stored in a
        %   row vector with $n_{coupling}$ element.
        %
        % `atom2`
        % : Second magnetic atom, same as `atom1`.
        %
        % `mat_idx`
        % : Stores pointers to matrices for every coupling in a
        %   $[3\times n_{coupling}]$ matrix, maximum three matrix per
        %   coupling (zeros for no coupling) is allowed.
        %
        % `idx`
        % : Neighbor index, increasing indices for the equivalent
        %   couplings, starting with 1,2,... which means first and second
        %   neighbor bonds, respectively.
        %
        % `type`
        % : Type of coupling corresponding to `mat_idx` matrices.
        %   Default is 0 for quadratic exchange, `type = 1` for
        %   biquadratic exchange.
        % 
        % `sym`
        % : If `true` the bond symmetry operators will be applied
        %   on the assigned matrix.
        % 
        % `rdip`
        % : Maximum distance until the dipolar interaction is
        %   calculated. Zero value means no dipolar interactions
        %   are considered.
        % 
        % `nsym`
        % : The largest bond `idx` value that is generated
        %   using the space group operators. Typically very long bonds for
        %   dipolar interactions won't be calculated using space group
        %   symmetry.
        %
        coupling
        % stores the magnetic structure
        %
        % ### Sub fields
        %
        % `F`
        % : Complex magnetization (strictly speaking complex
        %   spin expectation value) for every spin in the magnetic
        %   cell, represented by a matrix with dimensions of $[3\times
        %   n_{magext}\times n_k]$,
        %   where `nMagExt = nMagAtom*prod(N_ext)` and $n_k$ is the number
        %   of the magnetic propagation vectors.
        %
        % `k`
        % : Magnetic propagation vectors stored in a matrix with dimensions
        %   of $[3\times n_k]$.
        %
        % `N_ext`
        % : Size of the magnetic supercell in lattice units, default value
        %   is `[1 1 1]` emaning that the magnetic cell is identical to the
        %   crystallographic cell. The $[1\times 3]$ vector extends the cell
        %   along the $a$, $b$ and $c$ axes.
        %
        mag_str
        % stores the physical units for the Hamiltonian
        %
        % Default values are meV, T, \\ang and K for energy, magnetic
        % field, length and temperature, respectively.
        %
        % ### Sub fields
        %
        % `kB`
        % : Boltzmann constant, default value is 0.0862 meV/K.
        %
        % `muB`
        % : Bohr magneton, default values is 0.0579 meV/T.
        %
        % `mu0`
        % : Vacuum permeability, default value is 201.335431 T$^2$\\ang$^3$/meV.
        %
        % `label`
        % : Labels for distance, energy, magnetic field and temperature
        % stored in a cell with dimensions $[1\times 4]$.
        %
        % `nformula`
        % : Number of formula units in the unit cell.
        %
        % `qmat`
        % : Transformation matrix that converts the given $Q$ values into
        % the internal reciprocal lattice. The matrix has dimensions of
        % $[3\times 3]$.
        %
        unit
        % stores temporary values
        %
        % This property should be only used to check consistency of the code.
        %
        % {% include warning.html content="Changing these values is strongly 
        % discouraged as it can break the code!" %}
        % 
        % ### Sub fields
        %
        % `matom`
        % : Generated data of the magnetic unit cell.
        %
        % `symop`
        % : Generated symmetry operators for each bond.
        %
        cache = struct('matom',[],'symop',[]);
        
    end
    
    properties (Access = private)
        % stores the property change listener handles
        propl = event.proplistener.empty;
        % stores whether the couplings are generated under symmetry constraints
        sym   = false;
        % stores whether the calculation are done symbolically
        symb  = false;
        % use the version property as contant, this will be executed only
        % once
        ver   = sw_version;
    end
    
    methods (Static)
        % static methods
        validate(varargin)
    end
        
    methods
        function obj = spinw(varargin)
            % spinw constructor
            %
            % ### Syntax
            %
            % `obj = spinw`
            %
            % `obj = spinw(struct)` 
            %
            % `obj = spinw(hFigure)`
            %
            % `obj = spinw(fName)`
            %
            % `obj = spinw(obj)`
            %
            % ### Description
            %
            % `obj = spinw` creates an empty SpinW object with default
            % values.
            %
            % `obj = spinw(struct)` creates a SpinW object from a structure
            % which has fields that are compatible with the SpinW property
            % structure.
            %
            % `obj = spinw(hFigure)` clones SpinW object from an swplot
            % figure or spectral plot figure.
            %
            % `obj = spinw(fName)` imports the file referenced by `fName`.
            % SpinW is able to import .cif/.fts files for crystal or
            % magnetic structure from a local file or a web address.
            %
            % `obj = spinw(obj)` checks the input SpinW object for
            % consistency.
            %

            if nargin==0
                objS = initfield(struct);
                fNames = fieldnames(objS);
                for ii = 1:length(fNames)
                   obj.(fNames{ii}) = objS.(fNames{ii});
                end
                return
            end
            
            firstArg = varargin{1};
            
            if ishandle(firstArg)
                % get spinw object from graphics handle
                figDat = getappdata(firstArg);
                if isfield(figDat, 'spectra')
                    figDat = figDat.spectra;
                end
                if isfield(figDat, 'obj') && isa(figDat.obj, 'spinw')
                    obj = copy(figDat.obj);
                    return
                else
                    error('spinw:spinw:WrongInput', ...
                          'No spinw object could be found in the input figure');
                end
            end
            
            if isa(firstArg, 'spinw')
                %  it is used when objects are passed as arguments.
                obj = copy(firstArg);
                return
            end
            
            if isstruct(firstArg)
                objS = initfield(firstArg);
                
                % change lattice object
                if size(objS.lattice.lat_const,1)==3
                    objS.lattice.lat_const = objS.lattice.lat_const';
                end
                if size(objS.lattice.angle,1)==3
                    objS.lattice.angle = objS.lattice.angle';
                end
                
                spinw.validate(objS);
                fNames = fieldnames(objS);
                for ii = 1:length(fNames)
                    obj.(fNames{ii}) = objS.(fNames{ii});
                end
                return;
            end
            if ischar(firstArg)
                % import data from file (cif/fst are supported)
                
                objS = initfield(struct);
                fNames = fieldnames(objS);
                for ii = 1:length(fNames)
                    obj.(fNames{ii}) = objS.(fNames{ii});
                end
                
                obj = sw_import(firstArg,false,obj);
                
            end
            
        end % .spinw
        
        
        function nMagExt = nmagext(obj)
            % number of magnetic sites
            %
            % ### Syntax
            %
            % `nMagExt = nmagext(obj)`
            %
            % ### Description
            %
            % `nMagExt = nmagext(obj)` returns the number of magnetic sites
            % in the magnetic supercell. If the magnetic supercell (stored
            % in `spinw.mag_str.nExt` is identical to the crystal lattice)
            % the number of magnetic sites is equal to the number of
            % magnetic atoms in the unit cell. Where the number of magnetic
            % atoms in the unit cell can be calculated using [spinw.matom].
            %
            % ### See Also
            %
            % [spinw.matom] \| [spinw.natom]
            %
            
            nMagExt = size(obj.mag_str.F,2);
        end
        function nAtom = natom(obj)
            % number of symmetry unrelated atoms
            %
            % ### Syntax
            %
            % `nAtom = natom(obj)`
            %
            % ### Description
            %
            % `nAtom = natom(obj)` return the number of symmetry unrelated
            % atoms stored in `obj`.
            %
            % ### See Also
            %
            % [spinw.nmagext] \| [spinw.atom]
            %
            
            nAtom = size(obj.unit_cell.r,2);
        end
                
        function clearcache(obj, chgField)
            % clears the cache
            %
            % ### Syntax
            %
            % `clearcache(obj)`
            %
            % ### Description
            %
            % `clearcache(obj)` clears the cache that contains
            % precalculated magnetic structure and bond symmetry operators.
            % It is not necessary to clear the cache at any point as SpinW
            % clears it whenever necessary. 
            %
            % ### See Also
            %
            % [spinw.cache]
            %
            
            % listening to changes of the spinw object to clear cache is
            % necessary
            if nargin<2
                % delete the existing listener handles
                delete(obj.propl(ishandle(obj.propl)));
                % remove cache
                obj.cache.matom = [];
                obj.cache.symop = [];
                return
            end
            
            switch chgField
                case 1
                    % magnetic atoms: delete the stored magnetic atom positions
                    obj.cache.matom = [];
                    % remove the listeners
                    delete(obj.propl(1:2));
                case 2
                    % bond symmetry operators: delete the stored operators
                    obj.cache.symop = [];
                    % remove the listeners
                    delete(obj.propl(3:5));
            end
        end

    end
    
    methods(Hidden=true,Static=true)
            function obj = loadobj(obj)
                % restore property listeners
                % add new listeners to the new object
                if ~isempty(obj.cache.matom)
                    % add listener to lattice and unit_cell fields
                    obj.addlistenermulti(1);
                end
                if ~isempty(obj.cache.symop)
                    % add listener to lattice, unit_cell and coupling fields
                    obj.addlistenermulti(2);
                end
            end
    end
        
    methods(Hidden=true)
        
        function addlistenermulti(obj, chgField)
            % create the corresponding listeners to each cache subfield
            
            switch chgField
                case 1
                    % add listener to lattice and unit_cell fields
                    obj.propl(1) = addlistener(obj,'lattice',  'PostSet',@(evnt,src)obj.clearcache(1));
                    obj.propl(2) = addlistener(obj,'unit_cell','PostSet',@(evnt,src)obj.clearcache(1));
                case 2
                    % add listener to lattice, unit_cell and coupling fields
                    obj.propl(3) = addlistener(obj,'lattice',  'PostSet',@(evnt,src)obj.clearcache(2));
                    obj.propl(4) = addlistener(obj,'unit_cell','PostSet',@(evnt,src)obj.clearcache(2));
                    obj.propl(5) = addlistener(obj,'coupling', 'PostSet',@(evnt,src)obj.clearcache(2));
            end
        end
        
        function obj = saveobj(obj)
            % remove property change listeners
            delete(obj.propl);
            % empty pointers
            obj.propl = event.proplistener.empty;
        end

        function lh = addlistener(varargin)
            lh = addlistener@handle(varargin{:});
        end
        function notify(varargin)
            notify@handle(varargin{:});
        end
        function Hmatch = findobj(varargin)
            Hmatch = findobj@handle(varargin{:});
        end
        function p = findprop(varargin)
            p = findprop@handle(varargin{:});
        end
        function TF = eq(varargin)
            TF = eq@handle(varargin{:});
        end
        function TF = ne(varargin)
            TF = ne@handle(varargin{:});
        end
        function TF = lt(varargin)
            TF = lt@handle(varargin{:});
        end
        function TF = le(varargin)
            TF = le@handle(varargin{:});
        end
        function TF = gt(varargin)
            TF = gt@handle(varargin{:});
        end
        function TF = ge(varargin)
            TF = ge@handle(varargin{:});
        end
        
        function varargout = set(varargin)
            varargout = set@handle(varargin{:});
        end
        function varargout = setdisp(varargin)
            varargout = setdisp@handle(varargin{:});
        end
        function varargout = getdisp(varargin)
            varargout = getdisp@handle(varargin{:});
        end
        function varargout = get(varargin)
            varargout = get@handle(varargin{:});
        end
    end % classdef
    
end