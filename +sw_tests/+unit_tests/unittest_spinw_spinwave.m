classdef unittest_spinw_spinwave < sw_tests.unit_tests.unittest_super
    % Runs through unit test for @spinw/spinwave.m

    properties
        swobj = [];
        swobj_tri = [];
        default_spinwave = struct('formfact', false, ...
                                  'incomm', false, ...
                                  'helical', false, ...
                                  'norm', false, ...
                                  'nformula', 0, ...
                                  'param', struct('notwin', true, ...
                                                  'sortMode', true, ...
                                                  'tol', 1e-4, ...
                                                  'omega_tol', 1e-5, ...
                                                  'hermit', true), ...
                                  'title', 'Numerical LSWT spectrum', ...
                                  'gtensor', false, ...
                                  'datestart', '', ...
                                  'dateend', '');
    end

    properties (TestParameter)
        % Test directions and literal qpts work
        qpts_h5 = {{[0 0 0], [1 0 0], 5}, ...
                   [0:0.25:1; zeros(2,5)]};
    end

    methods (TestClassSetup)
        function setup_chain_model(testCase)
            % Just create a very simple FM 1D chain model
            testCase.swobj = spinw;
            testCase.swobj.genlattice('lat_const', [3 8 8], 'angled', [90 90 90]);
            testCase.swobj.addatom('r', [0 0 0],'S', 1, 'label', 'MNi2');
            testCase.swobj.gencoupling('maxDistance', 7);
            testCase.swobj.addmatrix('value', -eye(3), 'label', 'Ja');
            testCase.swobj.addcoupling('mat', 'Ja', 'bond', 1);
            testCase.swobj.genmagstr('mode', 'direct', 'k', [0 0 0], 'n', [1 0 0], 'S', [0; 1; 0]);
        end
        function setup_tri_model(testCase)
            % Just create a very simple FM 1D chain model
            testCase.swobj_tri = spinw;
            testCase.swobj_tri.genlattice('lat_const', [4 4 6], 'angled', [90 90 120]);
            testCase.swobj_tri.addatom('r', [0 0 0], 'S', 3/2, 'label', 'MCr3');
            testCase.swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], 'n', [0 0 1], 'k', [1/3 1/3 0]);
            J1 = 1;
            testCase.swobj_tri.addmatrix('label','J1','value',J1);
            testCase.swobj_tri.gencoupling;
            testCase.swobj_tri.addcoupling('mat','J1','bond',1);
        end
    end
    methods
        function sw = get_expected_sw_qh5(testCase)
            expected_hkl = [0:0.25:1; zeros(2,5)];
            expected_Sab = zeros(3, 3, 2, 5);
            expected_Sab([1 9 10 18 19 27 28 36 37 45 46 ...
                          54 55 63 64 72 73 81 82 90]) = 0.5;
            expected_Sab([7 12 25 30 43 48 61 66 75 88]) =  0.5i;
            expected_Sab([3 16 21 34 39 52 57 70 79 84]) = -0.5i;

            sw = testCase.default_spinwave;
            sw.omega = [ 1e-5  2.  4.  2. -1e-5; ...
                        -1e-5 -2. -4. -2.  1e-5];
            sw.Sab = expected_Sab;
            sw.hkl = expected_hkl;
            sw.hklA = expected_hkl*2/3*pi;
            sw.obj = copy(testCase.swobj);
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
        function test_sw_qh5(testCase, qpts_h5)
            sw_out = testCase.swobj.spinwave(qpts_h5);
            expected_sw = testCase.get_expected_sw_qh5();
            testCase.verify_spinwave(expected_sw, sw_out);
        end
        function test_sw_qh5_periodic(testCase)
            % Test qpts in different BZ give same omega, Sab
            qpts = [1:0.25:2; ones(2,5)];
            sw_out = testCase.swobj.spinwave(qpts);
            expected_sw = testCase.get_expected_sw_qh5();
            expected_sw.hkl = qpts;
            expected_sw.hklA =  [qpts(1, :)*2/3; qpts(2:end, :)*0.25 ]*pi;
            testCase.verify_spinwave(expected_sw, sw_out);
        end
        function test_sw_qh5_saveH_saveV(testCase)
            sw_out = testCase.swobj.spinwave([0:0.25:1; zeros(2,5)], ...
                                             'saveV', true, 'saveH', true);
            expected_V = repmat(eye(2), 1, 1, 5);
            expected_H = zeros(2, 2, 5);
            expected_H(:, :, [2 4]) = 2*repmat(eye(2), 1, 1, 2);
            expected_H(:, :, 3) = 4*eye(2);

            expected_sw = testCase.get_expected_sw_qh5();
            expected_sw.V = expected_V;
            expected_sw.H = expected_H;
            testCase.verify_spinwave(expected_sw, sw_out);
        end
        function test_sw_qh5_optmem(testCase)
            % At least test that optmem gives the same result
            qpts = [0:0.25:1; zeros(2,5)];
            sw_out = testCase.swobj.spinwave(qpts, 'optmem', 2);
            expected_sw = testCase.get_expected_sw_qh5();
            testCase.verify_spinwave(expected_sw, sw_out);
        end
        function test_sw_qh5_fitmode(testCase)
            qpts = [0:0.25:1; zeros(2,5)];
            sw_out = testCase.swobj.spinwave(qpts, 'fitmode', true);
            % fitmode automatically turns off sortMode
            expected_sw = testCase.swobj.spinwave(qpts, 'sortMode', false);
            expected_sw = rmfield(expected_sw, {'obj', 'datestart', 'dateend'});
            testCase.verify_spinwave(expected_sw, sw_out);
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
            % Create copy to avoid changing obj for other tests
            swobj = copy(testCase.swobj)
            commensurate_spec = swobj.spinwave(hkl);
            swobj.genmagstr('mode', 'direct', 'k', [0.123 0 0], 'n', [1 0 0], 'S', [0; 1; 0]);
            % Forcing incomm struct means exchange parameters don't agree, so we need 'hermit' false
            incomm_spec = swobj.spinwave(hkl, 'hermit', false);
            testCase.assertEqual(size(incomm_spec.omega, 1), size(commensurate_spec.omega, 1) * 3);
        end
        function test_twin(testCase)
            % Tests that setting twins gives correct outputs
            % Create copy to avoid changing obj for other tests
            swobj_twin = copy(testCase.swobj);
            swobj_twin.addtwin('axis', [0 0 1], 'phid', [60 120], 'vol', [1 2]);
            rotc = swobj_twin.twin.rotc;
            hkl = [1 2; 3 4; 5 6];
            sw_out = swobj_twin.spinwave(hkl);

            expected_sw = testCase.default_spinwave;
            expected_sw.param.notwin = false;
            expected_sw.omega = {};
            expected_sw.Sab = {};
            expected_sw.obj = copy(swobj_twin);
            expected_sw.hkl = hkl;
            expected_sw.hklA = [2/3 4/3; 0.75, 1; 1.25, 1.5]*pi;
            % Recalculate without twins for each set of hkl's and compare
            [~, rotQ] = swobj_twin.twinq([0;0;0]);
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
        function test_notwin(testCase)
            % Create copy to avoid changing obj for other tests
            swobj_twin = copy(testCase.swobj);
            swobj_twin.addtwin('axis', [0 0 1], 'phid', [60 120], 'vol', [1 2]);
            qpts = [0:0.25:1; zeros(2,5)];
            sw_out = swobj_twin.spinwave(qpts, 'notwin', true);

            expected_sw = testCase.get_expected_sw_qh5;
            expected_sw.obj.twin = swobj_twin.twin;
            testCase.verify_spinwave(expected_sw, sw_out);
        end
        function test_cmplxBase_equivalent_with_tri(testCase)
            qpts = {[0 0 0], [1 0 0], 5};
            sw_out = testCase.swobj_tri.spinwave(qpts);
            sw_out_cmplx = testCase.swobj_tri.spinwave(qpts, 'cmplxBase', true);
            testCase.verify_spinwave(sw_out, sw_out_cmplx);
        end
        function test_cmplxBase_fails_with_chain(testCase)
            qpts = {[0 0 0], [1 0 0], 5};
            testCase.verifyError(...
                @() testCase.swobj.spinwave(qpts, 'cmplxBase', true), ...
                'spinw:spinwave:NonPosDefHamiltonian');
        end
        function test_formfact(testCase)
            qpts = {[0 0 0] [10 5 1] 19};
            sw_ff = testCase.swobj.spinwave(qpts, 'formfact', true);

            % Test that Sab with the form factor (ff) is explicitly the
            % same as Sab with no ff multiplied by ff
            % ff calculated with sw_mff and the scaling is F(Q)^2.
            expected_sw = testCase.swobj.spinwave(qpts, 'formfact', false);
            ff = sw_mff(testCase.swobj.unit_cell.label{1}, sw_ff.hklA);
            expected_sw.Sab = expected_sw.Sab.*permute(ff.^2, [1 3 4 2]);
            expected_sw.formfact = true;

            testCase.verify_spinwave(expected_sw, sw_ff);
        end
        function test_formfactfun(testCase)
            function F = formfactfun(atom_label, Q)
                F = sum(Q, 1);
            end
            qpts = {[0 0 0] [10 5 1] 19};
            sw_ff = testCase.swobj.spinwave(qpts, 'formfact', true, ...
                                            'formfactfun', @formfactfun);

            % Test that Sab with the form factor (ff) is explicitly the
            % same as Sab with no ff multiplied by ff
            expected_sw = testCase.swobj.spinwave(qpts, 'formfact', false);
            ff = formfactfun(testCase.swobj.unit_cell.label{1}, sw_ff.hklA);
            expected_sw.Sab = expected_sw.Sab.*permute(ff.^2, [1 3 4 2]);
            expected_sw.formfact = true;

            testCase.verify_spinwave(expected_sw, sw_ff, 'rel_tol', 1e-15);
        end
        function test_gtensor(testCase)
            qpts = {[0 0 0], [1 1 1], 5};
            gmat = diag([1, 2, 3]);
            % Create copy to avoid changing obj for other tests
            swobj_g = copy(testCase.swobj);
            swobj_g.addmatrix('label','g_1','value', gmat)
            swobj_g.addg('g_1')
            sw_g = swobj_g.spinwave(qpts, 'gtensor', true);

            % Also check that it warns that gtensor is not being used
            expected_sw = testCase.verifyWarning(...
                @() swobj_g.spinwave(qpts, 'gtensor', false), ...
                'spinw:spinwave:NonZerogTensor');
            expected_sw.Sab = expected_sw.Sab.*[1 2 3; 2 4 6; 3 6 9];
            expected_sw.gtensor = true;
            testCase.verify_spinwave(expected_sw, sw_g);
        end
        function test_gtensor_incomm(testCase)
            qpts = {[0 0 0], [1 1 1], 5};
            gmat = diag([1, 2, 3]);
            % Create copy to avoid changing obj for other tests
            swobj_g = copy(testCase.swobj_tri);
            swobj_g.addmatrix('label','g_1','value', gmat)
            swobj_g.addg('g_1')
            sw_g = swobj_g.spinwave(qpts, 'gtensor', true);

            expected_sw = swobj_g.spinwave(qpts);
            expected_sw.Sab = expected_sw.Sab.*[2.25 2.25 4.5; ...
                                                2.25 2.25 4.5; ...
                                                4.5  4.5  9];
            expected_sw.gtensor = true;
            testCase.verify_spinwave(expected_sw, sw_g);
        end
        function test_hermit(testCase)
            % Create copy to avoid changing obj for other tests
            swobj_h = copy(testCase.swobj);
            % Tests that the 'hermit' option to switch to a non-hermitian calculation works
            % First make the model non-Hermitian by adding a large axial SIA perpendicular to the spins
            swobj_h.addmatrix('label', 'K', 'value', diag([-1 0 0]));
            swobj_h.addaniso('K');
            hkl = {[0 0 0] [0 1 0] [1 0 0] 50};
            % Check that calling it with 'hermit' on gives an error
            testCase.assertError(@() swobj_h.spinwave(hkl, 'hermit', true), ?MException);
            % Now check that there are imaginary eigenvalues/energies in the output with 'hermit' off
            spec = swobj_h.spinwave(hkl, 'hermit', false);
            testCase.assertGreaterThan(sum(abs(imag(spec.omega(:)))), 0);
        end
    end

end
