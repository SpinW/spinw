classdef unittest_spinw_intmatrix < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        %               icoupling = 1 2 3 4
        default_SS = struct('all', [1 1 0 0; % dl_a
                                    0 0 1 1; % dl_b
                                    0 0 0 0; % dl_c
                                    1 2 1 2; % matom1
                                    1 2 1 2; % matom2
                                    1 1 -2 -2; % J11
                                    0 0 0 0; % J12
                                    0 0 0 0; % J13
                                    0 0 0 0; % J21
                                    1 1 -2 -2; % J22
                                    0 0 0 0; % J23
                                    0 0 0 0; % J31
                                    0 0 0 0; % J32
                                    1 1 -2 -2; % J33
                                    0 0 1 1], ... % 0=quad, 1=biquad
                            'dip', [1       1;
                                    0       0;
                                    0       0;
                                    1       2;
                                    1       2;
                                   -0.0537 -0.0537;
                                    0       0;
                                    0       0;
                                    0       0;
                                    0.0067  0.0067;
                                    0       0;
                                    0       0;
                                    0       0;
                                    0.0067  0.0067;
                                    0      0]);
        default_SI = struct('aniso', repmat([0 0 0; 0 0 0; 0 0 -0.1],1,1,2), ...
                            'g', repmat([2 0 0; 0 1 0; 0 0 1],1,1,2), ...
                            'field', [0 0 0.5]);
        default_RR = [0 0.5; 0 0.5; 0 0.5];
    end

    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw();
            testCase.swobj.genlattice('lat_const', [2 3 5], 'sym', 'I m m m')
            testCase.swobj.addatom('r',[0; 0; 0],'S',1)
            testCase.swobj.addmatrix('label','g1','value',diag([2 1 1]))
            testCase.swobj.addmatrix('label', 'A', ...
                                     'value', diag([0 0 -0.1])) % c easy
            testCase.swobj.addmatrix('label', 'J1', 'value', 1)
            testCase.swobj.addmatrix('label', 'J2', 'value', -2)
            testCase.swobj.addmatrix('label','D','value',[0 -1 0])
            testCase.swobj.addmatrix('label','gen','value', ...
                                     reshape(2:2:18, [3, 3]))
            testCase.swobj.gencoupling('maxDistance', 5);
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1); % bond // a
            testCase.swobj.addcoupling('mat', 'J2', 'bond', 2, 'type', 'biquadratic'); % bond // b
            testCase.swobj.addaniso('A');
            testCase.swobj.addg('g1');
            testCase.swobj.field([0 0 0.5]);
            testCase.swobj.coupling.rdip = 3;  % set max dist for dipolar
        end
    end
    
    methods
        function expected_SS = get_expected_SS_fitmode_false(testCase)
            expected_SS = testCase.default_SS;
            expected_SS.iso = expected_SS.all(1:6,1:2);
            expected_SS.bq = expected_SS.all(1:6,3:4);
            expected_SS.ani = [1; 0; 0; 2; 1; 1; 2; 3];
            expected_SS.dm = [1; 0; 0; 2; 1; 0; -1; 0];
            expected_SS.gen = [1 0 0 2 1 2 4 6 8 10 12 14 16 18]';
            expected_SS.all = [expected_SS.all, zeros(15, 3)];
            expected_SS.all(1:5, 5:end) = repmat(expected_SS.ani(1:5), ...
                                                 1, 3);
            expected_SS.all(1:end-1, 6) = expected_SS.gen;
            expected_SS.all([8, 12], 5) = [-1 1]; % DM
            expected_SS.all([6, 10, 14], 7) = [1, 2, 3]; % aniso
        end
    end
    
    methods (Test)

        function test_intmatrix_no_matoms(testCase)
            [SS, SI, RR] = spinw().intmatrix('fitmode', true);
            expected_SS = struct('all', zeros(15,0), 'dip', zeros(15,0));
            expected_SI = struct('aniso', zeros(3,3,0), 'g', zeros(3,3,0), ...
                                 'field', zeros(1,3));
            expected_RR = zeros(3,0);
            testCase.verify_val(expected_SS, SS)
            testCase.verify_val(expected_SI, SI)
            testCase.verify_val(expected_RR, RR)
        end
        
        function test_intmatrix_no_couplings(testCase)
            sw = spinw();
            sw.addatom('r',[0; 0; 0],'S',1);
            [SS, SI, RR] = sw.intmatrix('fitmode', true);
            expected_SS = struct('all', zeros(15,0), 'dip', zeros(15,0));
            expected_SI = struct('aniso', zeros(3,3), 'g', 2*eye(3), ...
                                 'field', zeros(1,3));
            expected_RR = zeros(3,1);
            testCase.verify_val(expected_SS, SS)
            testCase.verify_val(expected_SI, SI)
            testCase.verify_val(expected_RR, RR)
        end
        
        function test_intmatrix_fitmode_true(testCase)
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', true);
            testCase.verify_val(testCase.default_SS, SS, 'abs_tol', 1e-4)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
        
                
        function test_intmatrix_fitmode_true_plotmode_true(testCase)
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', true, ...
                                                    'plotmode', true);
            expected_SS = testCase.default_SS;
            expected_SS.all = [expected_SS.all; zeros(3, 4)];
            expected_SS.all(15:end,:) = [3 3 4 4;
                                         1 1 2 2;
                                         0 0 1 1;
                                         1 2 3 4];
            testCase.verify_val(expected_SS, SS, 'abs_tol', 1e-4)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
        
        function test_intmatrix_fitmode_true_DM_interaction(testCase)
            testCase.verifyWarning( ...
                @() testCase.swobj.addcoupling('mat', 'D', 'bond', 3, 'subIdx', 1), ...
                'spinw:addcoupling:SymetryLowered');
            
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', true);
            
            dm_elems = [1; 0; 0; 2; 1; 0; 0; -1; 0; 0; 0; 1; 0; 0; 0];
            expected_SS = testCase.default_SS;
            expected_SS.all = [expected_SS.all, dm_elems];
            testCase.verify_val(expected_SS, SS, 'abs_tol', 1e-4)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
        
        function test_intmatrix_fitmode_false(testCase)
            % add all different types of interaction
            testCase.disable_warnings('spinw:addcoupling:SymetryLowered');
            testCase.swobj.addmatrix('label','Janiso','value', ...
                                     diag([1,2,3]))
            for mat_name = {'D', 'gen', 'Janiso'}
                testCase.swobj.addcoupling('mat', mat_name, 'bond', 3, 'subIdx', 1);
            end
             
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', false);
            
            expected_SS = testCase.get_expected_SS_fitmode_false();
            testCase.verify_val(expected_SS, SS, 'abs_tol', 1e-4)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
        
        function test_intmatrix_zeroC_false_removes_zero_matrices(testCase)
            testCase.swobj.addmatrix('label', 'J1', 'value', 0)
            
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', true, ...
                'zeroC', false);
            
            expected_SS = testCase.default_SS;
            expected_SS.all = expected_SS.all(:,3:end);
            testCase.verify_val(expected_SS, SS, 'abs_tol', 1e-4)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
        
        function test_intmatrix_conjugate(testCase)
            % add coupling between two different atoms
            testCase.verifyWarning( ...
                @() testCase.swobj.addcoupling('mat', 'gen', 'bond', 3, 'subIdx', 1), ...
                'spinw:addcoupling:SymetryLowered');
            % zero other couplings for brevity (will be omitted by zeroC)
            testCase.swobj.addmatrix('label', 'J1', 'value', 0)
            testCase.swobj.addmatrix('label', 'J2', 'value', 0)
            
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', true, ...
                'conjugate', true);
     
            expected_SS = testCase.default_SS;
            % two bonds 2->1 and 1->2 with half interaction
            bond1 = [1, 0, 0, 2, 1, 1:9, 0]';
            bond2 = [-1, 0, 0, 1, 2, 1, 4, 7, 2, 5, 8, 3, 6, 9, 0]';
            expected_SS.all = [bond1, bond2];
            expected_SS.dip = repmat(expected_SS.dip, 1, 2);
            expected_SS.dip(1:3,3:end) = -expected_SS.dip(1:3,3:end);
            expected_SS.dip(6:end-1,:) = 0.5*expected_SS.dip(6:end-1,:);
            testCase.verify_val(expected_SS, SS, 'abs_tol', 1e-4)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
               
        function test_extend_false_with_supercell(testCase)
           % make a supercell
           testCase.swobj.genmagstr('mode', 'random', 'nExt', [2 1 1]);
           
           [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', true, ...
                'extend', false);
            
            testCase.verify_val(testCase.default_SS, SS, 'abs_tol', 1e-4)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
        
        function test_extend_true_with_supercell(testCase)
           % make a supercell
           testCase.swobj.genmagstr('mode', 'random', 'nExt', [2 1 1]);
           
           [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', true, ...
                'extend', true);
            
            % SS
            expected_SS = testCase.default_SS;
            expected_SS.all = repmat(expected_SS.all, 1, 2);
            expected_SS.all(1,1:2) = 0;
            expected_SS.all(5,1:2) = [3, 4];
            expected_SS.all(4,5:8) = repmat([3, 4], 1, 2);
            expected_SS.all(5,7:8) = [3, 4];
            expected_SS.dip = repmat(expected_SS.dip, 1, 2);
            expected_SS.dip(1:5,:) = expected_SS.all(1:5,[1:2 5:6]);
            testCase.verify_val(expected_SS, SS, 'abs_tol', 1e-4)
            % SI
            expected_SI = testCase.default_SI;
            expected_SI.aniso = repmat(expected_SI.aniso, 1,1,2);
            expected_SI.g = repmat(expected_SI.g, 1,1,2);
            testCase.verify_val(expected_SI, SI)
            % RR
            expected_RR = [0 0.25 0.5 0.75;
                           0 0.5  0   0.5;
                           0 0.5  0   0.5];
            testCase.verify_val(expected_RR, RR)
        end
        
        function test_sortDM_reorders_bonds(testCase)
            testCase.disable_warnings('spinw:addcoupling:SymetryLowered');
            testCase.swobj.addmatrix('label', 'J2', 'value', 0);
            % make face-centred to have bond order depend on  sortDM
            testCase.swobj.genlattice('sym', 'F m m m');
            testCase.swobj.gencoupling();
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, 'subIdx',4:6);
            
            [SS_sort, ~, ~] = testCase.swobj.intmatrix('fitmode', true, ...
                                                       'sortDM', true);
            [SS_unsort, ~, ~] = testCase.swobj.intmatrix('fitmode', true, ...
                                                         'sortDM', false);
            % first bond same
            testCase.verify_val(SS_sort.all(:,1), SS_unsort.all(:,1));
            % just check the atom indices of the subsequent bonds
            testCase.verify_val([2 1; 1 2], SS_unsort.all(4:5,2:end));
            testCase.verify_val([1 2; 2 1], SS_sort.all(4:5,2:end));
        end
        
        function test_2_atoms_in_unit_cell_P1_sym(testCase)
            % Revert P1 sym and add atom at body-center
            testCase.swobj.genlattice('sym', 'P 1')
            testCase.swobj.addatom('r',[0.5; 0.5; 0.5],'S',1)
            testCase.swobj.gencoupling('maxDistance', 5);
            % need to add each matrix twice (once for each atom)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1:2); % bond // a
            testCase.swobj.addcoupling('mat', 'J2', 'bond', 3:4, 'type', ...
                                       'biquadratic'); % bond // b
            testCase.swobj.addaniso('A');
            testCase.swobj.addg('g1');
            testCase.swobj.field([0 0 0.5]);
            
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', true);
            testCase.verify_val(testCase.default_SS, SS, 'abs_tol', 1e-4)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
    end
     
     methods (Test, TestTags = {'Symbolic'})
        function symbolic_obj_with_fitmode_true(testCase)
            testCase.swobj.symbolic(true);
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', true);
            % replace symbolic variables with 1
            SS = structfun(@sw_sub1, SS, 'UniformOutput', false);
            SI = structfun(@sw_sub1, SI, 'UniformOutput', false);
            expected_SS = testCase.default_SS;
            expected_SS.dip(6:end,:) = 298.3569*expected_SS.dip(6:end,:);
            testCase.verify_val(expected_SS, SS, 'abs_tol', 5e-3)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
        
        function symbolic_obj_with_fitmode_false(testCase)
            % add all different types of interaction
            testCase.disable_warnings('spinw:addcoupling:SymetryLowered');
            testCase.swobj.addmatrix('label','Janiso','value', ...
                                     diag([1,2,3]))
            for mat_name = {'D', 'gen', 'Janiso'}
                testCase.swobj.addcoupling('mat', mat_name, 'bond', 3, ...
                                           'subIdx', 1);
            end
            testCase.swobj.symbolic(true);
            
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode', false);
            % replace symbolic variables with 1
            SS = structfun(@sw_sub1, SS, 'UniformOutput', false);
            SI = structfun(@sw_sub1, SI, 'UniformOutput', false);
            
            expected_SS = testCase.get_expected_SS_fitmode_false();
            expected_SS.dip(6:end,:) = 298.3569*expected_SS.dip(6:end,:);
            
            testCase.verify_val(expected_SS, SS, 'abs_tol', 5e-3)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
        
     
    end

end
