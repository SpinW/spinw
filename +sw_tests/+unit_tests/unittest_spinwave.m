classdef unittest_spinwave < sw_tests.unit_tests.unittest_super
    % Runs through unit test for @spinw/spinwave.m

    properties
        swobj = [];
        default_spinwave = struct('formfact', boolean(0), ...
                                  'incomm', boolean(0), ...
                                  'helical', boolean(0), ...
                                  'norm', boolean(0), ...
                                  'nformula', 0, ...
                                  'param', struct('notwin', boolean(1), ...
                                                  'sortMode', boolean(1), ...
                                                  'tol', 1e-4, ...
                                                  'omega_tol', 1e-5, ...
                                                  'hermit', boolean(1)), ...
                                  'title', 'Numerical LSWT spectrum', ...
                                  'gtensor', boolean(0), ...
                                  'datestart', '', ...
                                  'dateend', '');
    end

    properties (TestParameter)
        % Test directions and literal qpts work
        qpts_h5 = {{[0 0 0], [1 0 0], 5}, ...
                   [0:0.25:1; zeros(2,5)]}
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
        function test_noInput(testCase)
            % Tests that if call spinwave with no input, it calls the help
            % First mock the help call
            help_function = sw_tests.utilities.mock_function('swhelp');
            testCase.swobj.spinwave();
            testCase.assertEqual(help_function.n_calls, 1);
            testCase.assertEqual(help_function.arguments, {{'spinw.spinwave'}});
        end
        function test_sw_output(testCase, qpts_h5)
            sw_out = testCase.swobj.spinwave(qpts_h5);

            expected_hkl = [0:0.25:1; zeros(2,5)];
            expected_Sab = zeros(3, 3, 2, 5);
            expected_Sab([1 9 10 18 19 27 28 36 37 45 46 ...
                          54 55 63 64 72 73 81 82 90]) = 0.5;
            expected_Sab([7 12 25 30 43 48 61 66 75 88]) =  0.5i;
            expected_Sab([3 16 21 34 39 52 57 70 79 84]) = -0.5i;

            expected_sw_out = testCase.default_spinwave;
            expected_sw_out.omega = [ 1e-5  2.  4.  2. -1e-5; ...
                                     -1e-5 -2. -4. -2.  1e-5];
            expected_sw_out.Sab = expected_Sab;
            expected_sw_out.hkl = expected_hkl;
            expected_sw_out.hklA = expected_hkl*2.0943951023932;
            expected_sw_out.obj = testCase.swobj;
            testCase.verify_spinwave(expected_sw_out, sw_out, ...
                                     'rel_tol', 1e-10);
        end
        function test_symbolic(testCase)
            % Test that spinw.spinwavesym() is called if spinw.symbolic==true and spinw.spinwave() is called
            % Mock the spinwavesym and symbolic methods
            [mocksw, bh] = testCase.createMock(?spinw, 'MockedMethods', {'symbolic', 'spinwavesym'});
            % Make sure that spinw thinks it is symbolic
            testCase.assignOutputsWhen(withAnyInputs(bh.symbolic), true);
            % Make sure that spinwavesym outputs a sample spectra struct
            hkl = [1 2; 3 4; 5 6];
            out_spec = struct('hkl', hkl, 'swConv', [1; 2]);
            testCase.assignOutputsWhen(withAnyInputs(bh.spinwavesym), out_spec);
            % Call spinwave
            spec = mocksw.spinwave(hkl);
            testCase.assertCalled(withExactInputs(bh.spinwavesym()));
            testCase.assertEqual(spec, out_spec);
        end
        function test_incommensurate(testCase)
            % Tests that incommensurate calculation is ok
            hkl = {[0 0 0] [0 1 0] [1 0 0] 5};
            commensurate_spec = testCase.swobj.spinwave(hkl);
            testCase.swobj.genmagstr('mode', 'direct', 'k', [0.123 0 0], 'n', [1 0 0], 'S', [0; 1; 0]);
            % Forcing incomm struct means exchange parameters don't agree, so we need 'hermit' false
            incomm_spec = testCase.swobj.spinwave(hkl, 'hermit', false);
            testCase.assertEqual(size(incomm_spec.omega, 1), size(commensurate_spec.omega, 1) * 3);
            % Rests the structure for the following tests
            testCase.swobj.genmagstr('mode', 'direct', 'k', [0 0 0], 'n', [1 0 0], 'S', [0; 1; 0]);
        end
        function test_twin(testCase)
            % Tests that setting twins gives correct outputs
            testCase.swobj.addtwin('axis', [0 0 1], 'phid', [60 120], 'vol', [1 1]);
            rotc = testCase.swobj.twin.rotc;
            hkl = [1 2; 3 4; 5 6];
            sw_out = testCase.swobj.spinwave(hkl);

            expected_sw = testCase.default_spinwave;
            expected_sw.param.notwin = boolean(0);
            expected_sw.omega = {};
            expected_sw.Sab = {};
            expected_sw.obj = copy(testCase.swobj);
            expected_sw.hkl = hkl;
            expected_sw.hklA = [2/3 4/3; 0.75, 1; 1.25, 1.5]*pi;
            % Recalculate without twins for each set of hkl's and compare
            [~, rotQ] = testCase.swobj.twinq([0;0;0]);
            testCase.swobj.twin = struct('vol', 1, 'rotc', eye(3));
            for ii = 1:3
                sw_single = testCase.swobj.spinwave((hkl' * rotQ(:, :, ii))');
                expected_sw.omega = [expected_sw.omega sw_single.omega];
                rot = rotc(:, :, ii);
                rot_Sab = zeros(3, 3, 2, 2);
                % Twin Sab is a rotation of single Sab
                for jj = 1:2
                    for kk = 1:2
                        rot_Sab(:, :, jj, kk) = rot*sw_single.Sab(:, :, jj, kk)*rot';
                    end
                end
                expected_sw.Sab = [expected_sw.Sab rot_Sab];
            end
            testCase.verify_spinwave(expected_sw, sw_out);
        end
        function test_formfact(testCase)
            % Tests that the form factor calculation is applied correctly
            hkl = {[0 0 0] [10 0 0] 100};
            % Runs calculation with/without formfactor
            spec_no_ff = sw_neutron(testCase.swobj.spinwave(hkl, 'formfact', false));
            spec_ff = sw_neutron(testCase.swobj.spinwave(hkl, 'formfact', true));
            % The form factor is calculated using sw_mff, and the scaling is F(Q)^2 not F(Q).
            implied_ff = spec_ff.Sperp ./ spec_no_ff.Sperp;
            ff = sw_mff(testCase.swobj.unit_cell.label{1}, spec_ff.hklA);
            testCase.verify_val(implied_ff(1,:), ff.^2, 'rel_tol', 0.01, 'abs_tol', 1e-6);
        end
        function test_hermit(testCase)
            % Tests that the 'hermit' option to switch to a non-hermitian calculation works
            % First make the model non-Hermitian by adding a large axial SIA perpendicular to the spins
            testCase.swobj.addmatrix('label', 'K', 'value', diag([-1 0 0]));
            testCase.swobj.addaniso('K');
            hkl = {[0 0 0] [0 1 0] [1 0 0] 50};
            % Check that calling it with 'hermit' on gives an error
            testCase.assertError(@() testCase.swobj.spinwave(hkl, 'hermit', true), ?MException);
            % Now check that there are imaginary eigenvalues/energies in the output with 'hermit' off
            spec = testCase.swobj.spinwave(hkl, 'hermit', false);
            testCase.assertGreaterThan(sum(abs(imag(spec.omega(:)))), 0);
        end
    end

end
