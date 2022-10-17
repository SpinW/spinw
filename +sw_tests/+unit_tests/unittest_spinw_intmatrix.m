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
                                    'dip', zeros(15,0));
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
            testCase.swobj.gencoupling();
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1); % bond // a
            testCase.swobj.addcoupling('mat', 'J2', 'bond', 2, 'type', 'biquadratic'); % bond // b
            testCase.swobj.addaniso('A');
            testCase.swobj.addg('g1')
            testCase.swobj.field([0 0 0.5])
            testCase.swobj.plot('range', [2,1,1])
        end
    end

    methods (Test)
            
%         function test_addg_requires_matrix(testCase)
%             testCase.verifyError(@() testCase.swobj.addg(), ...
%                 'MATLAB:minrhs') % better if sw_readparam:MissingParameter
%         end
        
        function test_intmatrix_fitmode_true(testCase)
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode',1);
            testCase.verify_val(testCase.default_SS, SS)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
        
        function test_intmatrix_fitmode_true_DM_interaction(testCase)
            testCase.swobj.addmatrix('label','D','value',[0 -1 0])
            testCase.swobj.addcoupling('mat', 'D', 'bond', 3, 'subIdx', 1);
            
            [SS, SI, RR] = testCase.swobj.intmatrix('fitmode',1);
            
            dm_elems = [1; 0; 0; 2; 1; 0; 0; -1; 0; 0; 0; 1; 0; 0; 0];
            expected_SS = testCase.default_SS;
            expected_SS.all = [expected_SS.all, dm_elems];
            testCase.verify_val(expected_SS, SS)
            testCase.verify_val(testCase.default_SI, SI)
            testCase.verify_val(testCase.default_RR, RR)
        end
        
        
     end

end