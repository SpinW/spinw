classdef unittest_spinw_gencoupling < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        default_coupling = struct('dl', int32([1, 0; 0 1; 0 0]), ...
            'atom1', ones(1,2, 'int32'), 'atom2', ones(1, 2, 'int32'), ...
            'mat_idx', zeros(3, 2, 'int32'), 'idx', int32([1 1]), ...
            'type', zeros(3, 2, 'int32'), 'sym', zeros(3, 2, 'int32'), ...
            'rdip', 0.0, 'nsym',int32(0))
    end
    properties (TestParameter)
        dist_params = {'maxDistance', 'tolMaxDist', 'tolDist', 'dMin', 'maxSym'}
    end
    
    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw(); % default init
            testCase.swobj.genlattice('lat_const',[3 3 5])
            testCase.swobj.addatom('r',[0 0 0],'S',1)
        end
    end

    methods (Test)
            
        function test_gencoupling_requires_magnetic_atom(testCase)
            % change previously added atom to have S=0 (non magnetic)
            testCase.swobj.addatom('r',[0 0 0],'S',0, ...
                'label', 'atom_1', 'update', true)
            testCase.verifyError(...
                @() testCase.swobj.gencoupling(), ...
                'spinw:gencoupling:NoMagAtom')
        end
        
        function test_gencoupling_with_maxDistance_less_than_dMin(testCase)
            testCase.verifyError(...
                @() testCase.swobj.gencoupling('maxDistance', 0.1, 'dMin', 0.5), ...
                'spinw:gencoupling:MaxDLessThanMinD')
        end
        
        function test_gencoupling_with_tolMaxDist_grtr_than_maxDist(testCase)
            testCase.verifyError(...
                @() testCase.swobj.gencoupling('maxDistance', 5, 'tolMaxDist', 10), ...
                'spinw:gencoupling:TolMaxDist')
        end
        
        function test_gencoupling_with_tolDist_grtr_than_maxDist(testCase)
            testCase.verifyError(...
                @() testCase.swobj.gencoupling('maxDistance', 5, 'tolDist', 10), ...
                'spinw:gencoupling:TolDist')
        end
        
        function test_gencoupling_with_atom_dist_less_than_dMin(testCase)
            testCase.verifyError(...
                @() testCase.swobj.gencoupling('dMin', 5), ...
                'spinw:gencoupling:AtomPos')
        end
        
        function test_gencoupling_with_maxDistance_less_than_lattice_param(testCase)
            testCase.verifyError(...
                @() testCase.swobj.gencoupling('maxDistance', ...
                testCase.swobj.lattice.lat_const(1)/2), ...
                'spinw:gencoupling:maxDistance')
        end
        
        function test_gencoupling_with_negative_distances(testCase, dist_params)
            testCase.verifyError(...
                @() testCase.swobj.gencoupling(dist_params, -1), ...
                'spinw:gencoupling:NegativeDistances')
        end
        
        function test_gencoupling_with_tol_when_maxDistance_equal_lattice_param(testCase)
            % check bonds are created for atoms separated by latt. param.
            % when tol > delta
            delta = 1e-2;
            testCase.swobj.gencoupling('tolMaxDist', 10*delta, 'maxDistance', ...
                testCase.swobj.lattice.lat_const(1) - delta)
            testCase.verify_val(testCase.swobj.coupling, ...
                                testCase.default_coupling)
        end
        
        function test_gencoupling_with_tolDist_for_symm_equiv_bonds(testCase)
            % change lattice param to slightly break tetragonal sym.
            delta = 1e-3;
            testCase.swobj.genlattice('lat_const',[3 3+delta 5])
            % check that when tolDist > delta the bonds are equiv.
            testCase.swobj.gencoupling('maxDistance', 4, 'tolDist', 10*delta)
            testCase.verify_val(testCase.swobj.coupling, ...
                                testCase.default_coupling)
            % check that when tolDist > delta the bonds are inequiv.
            testCase.swobj.gencoupling('maxDistance', 4, 'tolDist', 0.1*delta)
            expected_coupling = testCase.default_coupling;
            expected_coupling.idx = int32([1, 2]); % i.e. not sym equiv
            testCase.verify_val(testCase.swobj.coupling, ...
                                expected_coupling)
        end
        
        function test_gencoupling_with_non_P0_spacegroup(testCase)
            % set ortho spacegroup even though lattice is tetragonal
            testCase.swobj.genlattice('sym', 'P 2')  % not overwrite abc
            % test forceNoSym = true (just uses bond length for idx)
            testCase.swobj.gencoupling('maxDistance', 4, 'forceNoSym', true)
            testCase.verify_val(testCase.swobj.coupling, ...
                                testCase.default_coupling)
            % test forceNoSym = false (default) - checks spacegroup
            testCase.swobj.gencoupling('maxDistance', 4)
            expected_coupling = testCase.default_coupling;
            expected_coupling.idx = int32([1, 2]); % i.e. not sym equiv
            expected_coupling.nsym = int32(2);
            testCase.verify_val(testCase.swobj.coupling, ...
                                expected_coupling)
        end
        
        function test_gencoupling_with_nonzero_fid(testCase)
            mock_fprintf = sw_tests.utilities.mock_function('fprintf0');
            fid = 3;
            testCase.swobj.gencoupling('maxDistance', 4, 'fid', fid)
            testCase.assertEqual(mock_fprintf.n_calls, 2);
            % check fid used to write file
            for irow = 1:mock_fprintf.n_calls
                testCase.assertEqual(mock_fprintf.arguments{irow}{1}, fid)
            end
            testCase.verify_val(testCase.swobj.coupling, ...
                                testCase.default_coupling)
        end
        
        function test_gencoupling_overwrites_previous_call(testCase)
            testCase.swobj.gencoupling('maxDistance', 8)
            testCase.swobj.gencoupling('maxDistance', 4)
            testCase.verify_val(testCase.swobj.coupling, ...
                                testCase.default_coupling)
        end
        
        function test_gencoupling_with_maxSym_multiple_bonds(testCase)
            testCase.swobj.genlattice('sym', 'P 4')  % not overwrite abc
            % set maxSym distance to only include bond idx = 1
            testCase.swobj.gencoupling('maxDistance', 4.5 , 'maxSym', 4)
            expected_coupling = testCase.default_coupling;
            expected_coupling.dl = int32([1 0 1 1; 0 1 -1 1; 0 0 0 0]);
            expected_coupling.atom1 = ones(1, 4, 'int32');
            expected_coupling.atom2 = ones(1, 4, 'int32');
            expected_coupling.mat_idx = zeros(3, 4, 'int32');
            expected_coupling.idx = int32([1 1 2 2]);
            expected_coupling.type = zeros(3, 4, 'int32');
            expected_coupling.sym = zeros(3, 4, 'int32');
            expected_coupling.nsym = int32(1);
            testCase.verify_val(testCase.swobj.coupling, ...
                                expected_coupling)
            % increase maxSym to include bond idx = 2 (length sqrt(2)*3)
            testCase.swobj.gencoupling('maxDistance', 4.5 , 'maxSym', 4.25);
            expected_coupling.nsym = int32(2);
            testCase.verify_val(testCase.swobj.coupling, ...
                                expected_coupling)
        end
        
        function test_gencoupling_maxSym_less_than_any_bond_length(testCase)
            testCase.swobj.genlattice('sym', 'P 4')  % not overwrite abc
            testCase.verifyWarning(...
                @() testCase.swobj.gencoupling('maxDistance', 4, ...
                'maxSym', 2), 'spinw:gencoupling:maxSym')
            % check nsym not set (bonds reduced to P0)
            testCase.verify_val(testCase.swobj.coupling, ...
                                testCase.default_coupling)
        end
        
        function test_gencoupling_when_tol_violates_symmetry(testCase)
            testCase.swobj.genlattice('lat_const',[3 3.5 5], 'sym', 'P 4')
            % bond along a and b are sym. equiv. according to P 4
            % but lengths differ by 0.5 Ang
            testCase.verifyError(...
                @() testCase.swobj.gencoupling('maxDistance', 3.25, ...
                    'tolDist', 0.6), 'spinw:gencoupling:SymProblem')
        end
        
        function test_gencoupling_with_two_mag_atoms(testCase)
            testCase.swobj.addatom('r',[0.25 0.25 0.25],'S',2)
            testCase.swobj.gencoupling('maxDistance', 3);
            expected_coupling = testCase.default_coupling;
            expected_coupling.dl = int32([0, 1, 0, 1, 0, 1, 0;
                                          0, 0, 1, 0, 1, 0, 1;
                                          0, 0, 0, 0, 0, 0, 0]);
            expected_coupling.atom1 = int32([1 2 2 1 1 2 2]);
            expected_coupling.atom2 = int32([2 1 1 1 1 2 2]);
            expected_coupling.mat_idx = zeros(3, 7, 'int32');
            expected_coupling.idx = int32([1 2 2 3 3 3 3]);
            expected_coupling.type = zeros(3, 7, 'int32');
            expected_coupling.sym = zeros(3, 7, 'int32');
            testCase.verify_val(testCase.swobj.coupling, ...
                                expected_coupling)
        end
        
        function test_gencoupling_ignores_non_mag_atoms(testCase)
            testCase.swobj.addatom('r',[0.25 0.25 0.25],'S',0)
            testCase.swobj.gencoupling('maxDistance', 4);
            testCase.verify_val(testCase.swobj.coupling, ...
                                testCase.default_coupling)
        end

     end

end
