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
        invalid_spgr = {'P2', 'P 7', '-x,y,-x', '-x,y', eye(2)};
        sym_param_name = {'sym', 'spgr'};
        basis_vecs =  {[1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2], ... % RH
            [1/2 1/2 0; 1/2 0 1/2; 0 1/2 1/2]}; % LH
        nelem = {1,3}; % length of cell input for spgr
        invalid_perm = {[1,4,2],[0,1,2], [1,1,1], 'bad', 'zzz', 'aaa', ...
            {1,2,3}}
        invalid_origin = {[-0.5,0,0], [0,2,0]};
        invalid_label = {1, {'label'}}
        % test user provided label always used for all types of sym input
        spgr_type = {'P 2', 3, '-x,y,-z', [eye(3), zeros(3,1)]};
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
            testCase.swobj.genlattice('angled', [90, 90, 90]);
            testCase.verify_val(testCase.default_latt, ...
                testCase.swobj.lattice)
        end
        
        function test_angle_and_angled_provided(testCase)
            value = [1,1,1];
            testCase.verifyWarning(...
                @() testCase.swobj.genlattice('angled', value, ...
                    'angle', value), ...
                'spinw:genlattice:WrongInput');
            expected_latt = testCase.default_latt;
            expected_latt.angle = deg2rad(value);  % 'angled' used
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_spgr_throws_deprecation_warning(testCase)
            testCase.verifyWarning(...
                @() testCase.swobj.genlattice('spgr', 'P 2'), ...
                'spinw:genlattice:DeprecationWarning');
        end
        
        function test_spgr_and_sym_throws_error(testCase)
            testCase.verifyError(...
                @() testCase.swobj.genlattice('spgr', 3, 'sym', 3), ...
                'spinw:genlattice:WrongInput');
        end
        
        function test_spacegroup_property(testCase, sym_param_name, spgr)
            testCase.swobj.genlattice(sym_param_name, spgr);
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.label = 'P 2';
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_label_always_used(testCase, sym_param_name, spgr_type)
            label = 'label';
            testCase.swobj.genlattice(sym_param_name, spgr_type, ...
                'label', label);
            testCase.verify_val(label, testCase.swobj.lattice.label)
        end
        
        function test_spacegroup_with_sym_operation_matrix(testCase, sym_param_name)
            testCase.swobj.genlattice(sym_param_name, testCase.P2_sym);
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.label = '';
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_spacegroup_with_sym_operation_string(testCase)
            % test perm supplied without symmetry throws warning
            perm = 'bac';
            testCase.verifyWarning(...
                @() testCase.swobj.genlattice('perm', perm), ...
                'spinw:genlattice:WrongInput');
            testCase.verify_val(testCase.default_latt, ...
                testCase.swobj.lattice) % object unchanged
            % supply spacegroup and check a and b swapped
            spgr_str = 'P 2';
            testCase.swobj.genlattice('sym', spgr_str, 'perm', 'bac');
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.sym(:, :, end) = [1 0 0 0; 0 -1 0 0; 0 0 -1 0];
            expected_latt.label = spgr_str;
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
        end
        
        function test_spacegroup_with_axes_permutation(testCase)
            sym_str = '-x,y,-z';
            testCase.swobj.genlattice('sym', sym_str);
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
            testCase.verifyWarning(...
                @() testCase.swobj.genlattice('origin', origin), ...
                'spinw:genlattice:WrongInput');
            testCase.verify_val(testCase.default_latt, ...
                testCase.swobj.lattice); % origin unchanged without spgr
            % define spacegroup
            testCase.swobj.genlattice('sym', testCase.P2_sym, ...
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
            args = {'lat_const',[4.5 4.5 4.5], 'sym','F 2 3'};
            % if no basis vectors are supplied then rot matrix always I
            R = testCase.swobj.genlattice(args{:});
            testCase.verify_val(eye(3), R);
            % add basis vectors for primitive cell
            R = testCase.swobj.genlattice(args{:}, 'bv', basis_vecs);
            sw_basis_vecs = R*basis_vecs;
            % check first spinwave basis vec is along x
            testCase.verify_val([sqrt(2)/2; 0; 0], sw_basis_vecs(:,1));
        end
        
        function test_spacegroup_with_cell_input(testCase, sym_param_name)
            spgr_str = '-x,y,-z';
            label = 'label';
            testCase.swobj.genlattice(sym_param_name, {spgr_str, label});
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.label = label;
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
            % provide label in cell and as separate argument
            new_label = 'new label';
            testCase.verifyWarning(...
                @() testCase.swobj.genlattice(sym_param_name, ...
                    {spgr_str, label},'label', new_label), ...
                'spinw:genlattice:WrongInput');
            expected_latt.label = new_label;
            testCase.verify_val(expected_latt, testCase.swobj.lattice);
        end
        
        function test_spacegroup_cell_input_invalid_numel(testCase, nelem)
            testCase.verifyError(...
                @() testCase.swobj.genlattice('sym', cell(1, nelem)), ...
                'spinw:genlattice:WrongInput');
        end
        
        function test_invalid_permutation(testCase, invalid_origin)
            testCase.verifyError(...
                @() testCase.swobj.genlattice('sym', 'P 2', ...
                    'origin', invalid_origin), ...
                'spinw:genlattice:WrongInput');
        end
        
        function test_invalid_origin(testCase, invalid_perm)
            testCase.verifyError(...
                @() testCase.swobj.genlattice('sym', 'P 2', ...
                    'perm', invalid_perm), ...
                'spinw:genlattice:WrongInput');
        end
        
        function test_invalid_spgr(testCase, invalid_spgr)
            testCase.verifyError(...
                @() testCase.swobj.genlattice('sym', invalid_spgr), ...
                'generator:WrongInput');
        end
        
        function test_non_default_spacegroup_overwritten(testCase)
            testCase.swobj.genlattice('sym','F 2 3');
            testCase.swobj.genlattice('sym', 'P 2');
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.label = 'P 2';
            testCase.verify_val(expected_latt, testCase.swobj.lattice);
        end
        
        function test_zero_spacegroup(testCase)
            testCase.swobj.genlattice('sym','P 2');
            testCase.swobj.genlattice('sym', 0);
            expected_latt = testCase.default_latt;
            expected_latt.label = 'No sym';
            expected_latt.sym = [eye(3) zeros(3,1)];% actually equiv. to P 1
            testCase.verify_val(expected_latt, testCase.swobj.lattice);
        end
        
        function test_invalid_label(testCase, invalid_label)
            testCase.swobj.genlattice('label', invalid_label)
            testCase.verify_val(testCase.default_latt, ...
                testCase.swobj.lattice); % label not overwritten
        end
        
        function test_lookup_new_lines_in_symmetry_dat_file(testCase)
            % copy file as backup (would need to do this anyway)
            spinw_dir = sw_rootdir();
            filesep = spinw_dir(end);
            dat_dir = [spinw_dir 'dat_files' filesep];
            dat_path = [dat_dir 'symmetry.dat'];
            backup_path = [dat_dir 'symmetry_backup.dat'];
            copyfile(dat_path, backup_path);
            % add line to file (same P2 sym op but new number and label)
            extra_line = ' 231  P P        : -x,y,-z\n';
            fid = fopen(dat_path, 'a');
            fprintf(fid, extra_line);
            fclose(fid);
            % run test
            testCase.swobj.genlattice('sym', 231);
            expected_latt = testCase.default_latt;
            expected_latt.sym = testCase.P2_sym;
            expected_latt.label = 'P P';
            testCase.verify_val(expected_latt, testCase.swobj.lattice)
            % restore backup - note movefile errored on windows server
            copyfile(backup_path, dat_path);
            delete(backup_path);
        end
        
     end

end
