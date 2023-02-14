classdef unittest_sw_egrid < sw_tests.unit_tests.unittest_super
    % Runs through unit test for @spinw/spinwave.m

    properties
        swobj = [];
        swobj_tri = [];
        spectrum = struct();
        sw_egrid_out = struct();
        qh5 = [0:0.25:1; zeros(2,5)];
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
            testCase.swobj.genmagstr('mode', 'direct', 'k', [0 0 0], ...
                                     'S', [0; 1; 0]);
        end
    end

    methods (TestMethodSetup)
        function setup_default_spectrum_and_egrid(testCase)
            testCase.spectrum.obj = testCase.swobj;
            testCase.spectrum.formfact = false;
            testCase.spectrum.incomm = false;
            testCase.spectrum.helical = false;
            testCase.spectrum.norm = false;
            testCase.spectrum.nformula = int32(0);
            testCase.spectrum.param = struct('notwin', true, ...
                                             'sortMode', true, ...
                                             'tol', 1e-4, ...
                                             'omega_tol', 1e-5, ...
                                             'hermit', true);
            testCase.spectrum.title = 'Numerical LSWT spectrum';
            testCase.spectrum.gtensor = false;
            testCase.spectrum.datestart = '01-Jan-2023 00:00:01';
            testCase.spectrum.dateend = '01-Jan-2023 00:00:02';
            testCase.spectrum.hkl = [0:0.25:1; zeros(2,5)];
            testCase.spectrum.hklA = testCase.spectrum.hkl*2/3*pi;
            testCase.spectrum.omega = [ 1e-5  2.  4.  2. -1e-5; ...
                                       -1e-5 -2. -4. -2.  1e-5];
            
            testCase.spectrum.Sab = zeros(3, 3, 2, 5);
            Sab1 = [0.5 0  0.5j; 0 0 0; -0.5j 0 0.5];
            Sab2 = [0.5 0 -0.5j; 0 0 0;  0.5j 0 0.5];
            testCase.spectrum.Sab(:, :, 1, 1:4) = repmat(Sab1, 1, 1, 1, 4);
            testCase.spectrum.Sab(:, :, 2, 5) = Sab1;
            testCase.spectrum.Sab(:, :, 2, 1:4) = repmat(Sab2, 1, 1, 1, 4);
            testCase.spectrum.Sab(:, :, 1, 5) = Sab2;

            testCase.sw_egrid_out = testCase.spectrum;
            testCase.sw_egrid_out.param.n = [0 0 1];
            testCase.sw_egrid_out.param.pol = false;
            testCase.sw_egrid_out.param.uv = {};
            testCase.sw_egrid_out.param.sumtwin = true;
            testCase.sw_egrid_out.intP = [];
            testCase.sw_egrid_out.Pab = [];
            testCase.sw_egrid_out.Mab = [];
            testCase.sw_egrid_out.Sperp = 0.5*ones(2, 5);
            testCase.sw_egrid_out.T = 0;
            testCase.sw_egrid_out.Evect = linspace(0, 4.4, 501);
            testCase.sw_egrid_out.swConv = zeros(500, 5);
            testCase.sw_egrid_out.swConv([727, 1455, 1727]) = 0.5;
            testCase.sw_egrid_out.swInt = 0.5*ones(2, 5);
            testCase.sw_egrid_out.component = 'Sperp';

        end
    end

    methods (Test)
        function test_noInput(testCase)
            % Tests that if sw_egrid is called with no input, it calls the help
            % First mock the help call
            help_function = sw_tests.utilities.mock_function('swhelp');
            sw_egrid();
            testCase.assertEqual(help_function.n_calls, 1);
            testCase.assertEqual(help_function.arguments, {{'sw_egrid'}});
        end
        function test_defaults(testCase)
            out = sw_egrid(testCase.spectrum);
            testCase.verify_obj(out, testCase.sw_egrid_out);
        end
        function test_Sperp(testCase)
            out = sw_egrid(testCase.spectrum, 'component', 'Sperp');
            testCase.verify_obj(out, testCase.sw_egrid_out);
        end
        function test_Sxx(testCase)
            component = 'Sxx';
            expected_out = rmfield(testCase.sw_egrid_out, ...
                                   {'intP', 'Pab', 'Mab', 'Sperp'});
            expected_out.param = rmfield(expected_out.param, ...
                                         {'n', 'pol', 'uv'});
            expected_out.component = component;
            out = sw_egrid(testCase.spectrum, 'component', component);
            testCase.verify_obj(out, expected_out);
        end
        function test_Sxy(testCase)
            component = 'Sxy';
            expected_out = rmfield(testCase.sw_egrid_out, {'intP', 'Pab', 'Mab', 'Sperp'});
            expected_out.param = rmfield(expected_out.param, {'n', 'pol', 'uv'});
            expected_out.component = component;
            expected_out.swInt = zeros(2, 5);
            expected_out.swConv = zeros(500, 5);
            out = sw_egrid(testCase.spectrum, 'component', component);
            testCase.verify_obj(out, expected_out);
        end
        function test_diff_Sxx_Szz(testCase)
            component = 'Sxx-Szz';
            expected_out = rmfield(testCase.sw_egrid_out, {'intP', 'Pab', 'Mab', 'Sperp'});
            expected_out.param = rmfield(expected_out.param, {'n', 'pol', 'uv'});
            expected_out.component = component;
            expected_out.swInt = zeros(2, 5);
            expected_out.swConv = zeros(500, 5);
            out = sw_egrid(testCase.spectrum, 'component', component);
            testCase.verify_obj(out, expected_out);
        end
        function test_sum_Sxx_Szz_Sxy(testCase)
            component = 'Sxx+Szz+Sxy';
            expected_out = rmfield(testCase.sw_egrid_out, {'intP', 'Pab', 'Mab', 'Sperp'});
            expected_out.param = rmfield(expected_out.param, {'n', 'pol', 'uv'});
            expected_out.component = component;
            expected_out.swInt = 2*expected_out.swInt;
            expected_out.swConv = 2*expected_out.swConv;
            out = sw_egrid(testCase.spectrum, 'component', component);
            testCase.verify_obj(out, expected_out);
        end
    end

end
