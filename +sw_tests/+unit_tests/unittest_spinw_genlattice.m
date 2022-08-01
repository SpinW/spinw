classdef unittest_spinw_genlattice < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        default_latt = struct('angle', repmat(pi/2, 1, 3), ...
            'lat_const', repmat(3, 1, 3),  ...
            'sym', zeros(3, 4, 0), ...
            'origin', zeros(1, 3), 'label', 'P 0');
        P2_sym = cat(3, [1 0 0 0; 0 1 0 0; 0 0 1 0], ...
                        [-1 0 0 0; 0 1 0 0; 0 0 -1 0])
    end
    properties (TestParameter)
        param_name = {'angle', 'lat_const'};
        spgr = {'P 2', 3};  % spacegroup and index in symmetry.dat
        sym_param_name = {'sym', 'spgr'};
        basis_vecs =  {[1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2], ... % RH
            [1/2 1/2 0; 1/2 0 1/2; 0 1/2 1/2]}; % LH
    end
    
    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw(); % default init
        end
    end

    methods (Test)     

        function test_params_alter_lattice_fields_same_name(testCase, param_name)
            value = [0.5, 0.5, 0.5];
            testCase.swobj.genlattice(param_name, value);
            expected_latt = testCase.default_latt;
            expected_latt.(param_name) = value;
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_angled_degree_conversion(testCase)
            % check bonds are created for atoms separated by latt. param.
            % when tol > delta
            testCase.swobj.genlattice('angled', [90, 90, 90]);
            testCase.verify_val(testCase.default_latt, ...
                testCase.swobj.lattice)
        end
        
        function test_angle_and_angled_provided(testCase)
            value = [1,1,1];
            testCase.swobj.genlattice('angled',[1,1,1], 'angle', [1,1,1]);
            expected_latt = testCase.default_latt;
            expected_latt.angle = deg2rad(value);  % 'angled' used
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_spacegroup_property(testCase, sym_param_name, spgr)
            testCase.swobj.genlattice(sym_param_name, spgr);
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.label = 'P 2';
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_spacegroup_with_sym_operation_matrix(testCase, sym_param_name)
            testCase.swobj.genlattice(sym_param_name, testCase.P2_sym);
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.label = '';
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_spacegroup_with_sym_operation_string(testCase)
            spgr_str = 'P 2';
            testCase.swobj.genlattice('spgr', spgr_str, 'perm', 'bac');
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.sym(:, :, end) = [1 0 0 0; 0 -1 0 0; 0 0 -1 0];
            expected_latt.label = spgr_str;
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_spacegroup_with_axes_permutation(testCase)
            sym_str = '-x,y,-z';
            testCase.swobj.genlattice('spgr', sym_str);
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.label = sym_str;
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_valid_label(testCase)
            label = 'P 4';
            testCase.swobj.genlattice('label', label);
            expected_latt = testCase.default_latt;
            expected_latt.label = label;
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
                
        function test_origin_set_only_when_spgr_provided(testCase)
            origin = [0.5, 0, 0];
            testCase.swobj.genlattice('origin', origin);
            testCase.verify_val(testCase.default_latt, ...
                testCase.swobj.lattice); % origin unchanged without spgr
            % dewfine spacegroup
            testCase.swobj.genlattice('spgr', testCase.P2_sym, ...
                'origin', origin);
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.label = '';
            expected_latt.origin = origin;
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_nformula_unit(testCase)
            nformula = int32(2);
            testCase.swobj.genlattice('nformula', nformula);
            testCase.verify_val(testCase.swobj.unit.nformula, nformula);
        end
        
        function test_basis_vector_and_rotation_matrix(testCase, basis_vecs)
            args = {'lat_const',[4.5 4.5 4.5], 'spgr','F 2 3'};
            % if no basis vectors are supplied then rot matrix always I
            R = testCase.swobj.genlattice(args{:});
            testCase.verify_val(eye(3), R);
            % add basis vectors for primitive cell
            R = testCase.swobj.genlattice(args{:}, 'bv', basis_vecs); 
            sw_basis_vecs = R*basis_vecs;
            % check first spinwave basis vec is along x
            testCase.verify_val([sqrt(2)/2; 0; 0], sw_basis_vecs(:,1));
        end
        
     end

end
