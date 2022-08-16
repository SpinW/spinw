classdef unittest_sw_neutron < sw_tests.unit_tests.unittest_super
    % Runs through unit test for sw_neutron.m

    properties
        swobj = [];
    end

    methods (TestClassSetup)
        function setup_spinw_model(testCase)
            % Just create a very simple FM 1D chain model
            testCase.swobj = spinw;
            testCase.swobj.genlattice('lat_const', [3 8 8], 'angled', [90 90 90]);
            testCase.swobj.addatom('r', [0 0 0],'S', 1, 'label', 'MNi2');
            testCase.swobj.gencoupling('maxDistance', 7);
            testCase.swobj.addmatrix('value', -eye(3), 'label', 'Ja');
            testCase.swobj.addcoupling('mat', 'Ja', 'bond', 1);
            testCase.swobj.genmagstr('mode', 'direct', 'k', [0 0 0], 'n', [1 0 0], 'S', [0; 1; 0]);
        end
    end
    methods (Test)
        function test_formfact(testCase)
            % Tests that the form factor calculation is applied correctly
            hkl = {[0 0 0] [10 0 0] 100};
            % Runs calculation with/without formfactor
            spec_no_ff = sw_neutron(testCase.swobj.spinwave(hkl, 'formfact', false));
            spec_ff = sw_neutron(testCase.swobj.spinwave(hkl, 'formfact', true));
            % The form factor is calculated using sw_mff, and the scaling is F(Q)^2 not F(Q).
            implied_ff = spec_ff.Sperp ./ spec_no_ff.Sperp;
            ff = sw_mff(testCase.swobj.unit_cell.label{1}, spec_ff.hklA);
            testCase.verify_val(ff.^2, implied_ff(1,:), ...
                                'rel_tol', 0.01, 'abs_tol', 1e-6);
        end
    end

end