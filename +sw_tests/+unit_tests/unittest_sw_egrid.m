classdef unittest_sw_egrid < sw_tests.unit_tests.unittest_super
    % Runs through unit test for @spinw/spinwave.m

    properties
        swobj = [];
        swobj_tri = [];
        spectrum = struct();
        sw_egrid_out = struct();
        sw_egrid_out_pol = struct();
        sw_egrid_out_sperp = struct();
        qh5 = [0:0.25:1; zeros(2,5)];
    end

    properties (TestParameter)
        % Components that require sw_neutron be called first
        pol_components = {'Mxx', 'Pxy', 'Pz'};
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
            % Spectrum input to sw_egrid
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

            % Output from sw_egrid when 'component' is specified
            testCase.sw_egrid_out = testCase.spectrum;
            testCase.sw_egrid_out.param.sumtwin = true;
            testCase.sw_egrid_out.T = 0;
            testCase.sw_egrid_out.Evect = linspace(0, 4.4, 501);
            testCase.sw_egrid_out.swConv = zeros(500, 5);
            testCase.sw_egrid_out.swConv([728, 1455, 1728]) = 0.5;
            testCase.sw_egrid_out.swInt = 0.5*ones(2, 5);
            testCase.sw_egrid_out.component = 'Sperp';

            % Default output from sw_egrid - when Sperp is used there are
            % extra fields
            testCase.sw_egrid_out_sperp = testCase.sw_egrid_out;
            testCase.sw_egrid_out_sperp.intP = [];
            testCase.sw_egrid_out_sperp.Mab = [];
            testCase.sw_egrid_out_sperp.Pab = [];
            testCase.sw_egrid_out_sperp.Sperp = 0.5*ones(2, 5);
            testCase.sw_egrid_out_sperp.param.n = [0 0 1];
            testCase.sw_egrid_out_sperp.param.pol = false;
            testCase.sw_egrid_out_sperp.param.uv = {};
            testCase.sw_egrid_out_sperp.param.sumtwin = true;

            % Output from sw_egrid when polarisation required - when
            % sw_neutron has to be called first
            testCase.sw_egrid_out_pol = testCase.sw_egrid_out_sperp;
            testCase.sw_egrid_out_pol.param.pol = true;
            testCase.sw_egrid_out_pol.intP = 0.5*ones(3, 2, 5);
            testCase.sw_egrid_out_pol.Pab = repmat(diag([-0.5 -0.5 0.5]), 1, 1, 2, 5);
            Mab_mat = [0.5 0 1i*0.5; 0 0 0; -1i*0.5 0 0.5];
            testCase.sw_egrid_out_pol.Mab = repmat(Mab_mat, 1, 1, 2, 5);
            testCase.sw_egrid_out_pol.Mab(:, :, 1, 5) = conj(Mab_mat);
            testCase.sw_egrid_out_pol.Mab(:, :, 2, :) = conj(testCase.sw_egrid_out_pol.Mab(:, :, 1, :));

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
        function test_pol_component_no_sw_neutron_causes_error(testCase, pol_components)
            testCase.verifyError(...
                @() sw_egrid(testCase.spectrum, 'component', pol_components), ...
                'sw_egrid:WrongInput');
        end
        function test_invalid_component(testCase)
            testCase.verifyError(...
                @() sw_egrid(testCase.spectrum, 'component', 'Sab'), ...
                'sw_parstr:WrongString');
        end
        function test_epsilon_deprecated_warning(testCase)
            testCase.verifyWarning(...
                @() sw_egrid(testCase.spectrum, 'epsilon', 1e-5), ...
                'sw_egrid:DeprecationWarning');
        end
        function test_defaults(testCase)
            out = sw_egrid(testCase.spectrum);
            testCase.verify_obj(out, testCase.sw_egrid_out_sperp);
        end
        function test_Sperp(testCase)
            out = sw_egrid(testCase.spectrum, 'component', 'Sperp');
            testCase.verify_obj(out, testCase.sw_egrid_out_sperp);
        end
        function test_Sxx(testCase)
            component = 'Sxx';
            expected_out = testCase.sw_egrid_out;
            expected_out.component = component;
            out = sw_egrid(testCase.spectrum, 'component', component);
            testCase.verify_obj(out, expected_out);
        end
        function test_Sxy(testCase)
            component = 'Sxy';
            expected_out = testCase.sw_egrid_out;
            expected_out.component = component;
            expected_out.swInt = zeros(2, 5);
            expected_out.swConv = zeros(500, 5);
            out = sw_egrid(testCase.spectrum, 'component', component);
            testCase.verify_obj(out, expected_out);
        end
        function test_diff_Sxx_Szz(testCase)
            component = 'Sxx-Szz';
            expected_out = testCase.sw_egrid_out;
            expected_out.component = component;
            expected_out.swInt = zeros(2, 5);
            expected_out.swConv = zeros(500, 5);
            out = sw_egrid(testCase.spectrum, 'component', component);
            testCase.verify_obj(out, expected_out);
        end
        function test_sum_Sxx_Szz_Sxy(testCase)
            component = 'Sxx+Szz+Sxy';
            expected_out = testCase.sw_egrid_out;
            expected_out.component = component;
            expected_out.swInt = 2*expected_out.swInt;
            expected_out.swConv = 2*expected_out.swConv;
            out = sw_egrid(testCase.spectrum, 'component', component);
            testCase.verify_obj(out, expected_out);
        end
        function test_cell_array_component(testCase)
            component = {'Sxx', 'Sxy'};
            expected_out = testCase.sw_egrid_out;
            expected_out.component = component;
            expected_out.swInt = {expected_out.swInt; zeros(2, 5)};
            expected_out.swConv = {expected_out.swConv; zeros(500, 5)};
            out = sw_egrid(testCase.spectrum, 'component', component);
            testCase.verify_obj(out, expected_out);
        end
        function test_Mzz(testCase)
            component = 'Mzz';
            expected_out = testCase.sw_egrid_out_pol;
            expected_out.component = component;

            neutron_out = sw_neutron(testCase.spectrum, 'pol', true);
            out = sw_egrid(neutron_out, 'component', component);
            testCase.verify_obj(out, expected_out);
        end
        function test_sum_Pzz_Mxx(testCase)
            component = 'Pzz+Mxx';
            expected_out = testCase.sw_egrid_out_pol;
            expected_out.component = component;
            expected_out.swInt = 2*expected_out.swInt;
            expected_out.swConv = 2*expected_out.swConv;

            neutron_out = sw_neutron(testCase.spectrum, 'pol', true);
            out = sw_egrid(neutron_out, 'component', component);
            testCase.verify_obj(out, expected_out);
        end

        function test_Px(testCase)
            component = 'Px';
            expected_out = testCase.sw_egrid_out_pol;
            expected_out.component = component;
            neutron_out = sw_neutron(testCase.spectrum, 'pol', true);
            out = sw_egrid(neutron_out, 'component', component);
            testCase.verify_obj(out, expected_out);
        end

        function test_fName_component(testCase)
            % add field to spectrum
            component = 'Sperp';
            spectrum = testCase.spectrum;
            spectrum.(component) = ones(2,5);  % double actual
            
            % note do not add fields normally added to output when Sperp
            expected_out = testCase.sw_egrid_out;
            expected_out.component = component;
            expected_out.(component) = spectrum.(component);
            expected_out.swInt = 2*expected_out.swInt;
            expected_out.swConv = 2*expected_out.swConv;

            out = sw_egrid(spectrum, 'component', component);

            testCase.verify_obj(out, expected_out);
        end
        function test_Evect(testCase)
            Evect = linspace(1, 3, 201);
            expected_out = testCase.sw_egrid_out_sperp;
            expected_out.swConv = zeros(200, 5);
            expected_out.swConv([301, 701]) = 0.5;
            expected_out.Evect = Evect;
            out = sw_egrid(testCase.spectrum, 'Evect', Evect);
            testCase.verify_obj(out, expected_out);
        end
        function test_Evect_cbin(testCase)
            Evect_in = linspace(1.005, 2.995, 200);
            expected_out = testCase.sw_egrid_out_sperp;
            expected_out.swConv = zeros(200, 5);
            expected_out.swConv([301, 701]) = 0.5;
            expected_out.Evect = linspace(1, 3, 201);
            out = sw_egrid(testCase.spectrum, 'Evect', Evect_in, 'binType', 'cbin');
            testCase.verify_obj(out, expected_out);
        end
        function test_temp(testCase)
            temp = 300;
            expected_out = testCase.sw_egrid_out_sperp;
            expected_out.swConv([728, 1728]) = 6.70976913583173;
            expected_out.swConv(1455) = 3.48826657260066;
            expected_out.T = temp;
            out = sw_egrid(testCase.spectrum, 'T', temp);
            testCase.verify_obj(out, expected_out, 'rel_tol', 1e-10);
        end
        function test_single_ion_temp(testCase)
            temp = 300;
            spectrum = testCase.spectrum;
            % Copy swobj here so don't interfere with other tests
            spectrum.obj = copy(testCase.swobj);
            spectrum.obj.single_ion.T = temp;
            expected_out = testCase.sw_egrid_out_sperp;
            expected_out.obj = spectrum.obj;
            expected_out.swConv([728, 1728]) = 6.70976913583173;
            expected_out.swConv(1455) = 3.48826657260066;
            expected_out.T = temp;
            out = sw_egrid(spectrum);
            testCase.verify_obj(out, expected_out, 'rel_tol', 1e-10);
        end
        function test_twin(testCase)
            swobj_twin = copy(testCase.swobj);
            swobj_twin.addtwin('axis', [0 0 1], 'phid', [60 120], 'vol', [1 2]);
            spectrum = swobj_twin.spinwave(testCase.spectrum.hkl);
            spectrum.datestart = testCase.spectrum.datestart;
            spectrum.dateend = testCase.spectrum.dateend;
            out = sw_egrid(spectrum);

            expected_out = testCase.sw_egrid_out_sperp;
            expected_out.obj = spectrum.obj;
            expected_out.omega = spectrum.omega;
            expected_out.Sab = spectrum.Sab;
            expected_out.param.notwin = false;
            expected_out.swConv = zeros(500, 5);
            expected_out.swConv([728, 1455, 1728]) = 0.125;
            expected_out.swConv([567, 1228, 1888, 2455]) = 0.65625;
            expected_out.swInt = 0.78125*ones(2, 5);
            expected_out.intP = cell(1,3);
            expected_out.Pab = cell(1,3);
            expected_out.Mab = cell(1,3);
            expected_out.Sperp = {0.5*ones(2, 5) 0.875*ones(2, 5) 0.875*ones(2, 5)};

            testCase.verify_obj(out, expected_out);
        end
        function test_twin_nosum(testCase)
            swobj_twin = copy(testCase.swobj);
            swobj_twin.addtwin('axis', [0 0 1], 'phid', [60 120], 'vol', [1 2]);
            spectrum = swobj_twin.spinwave(testCase.spectrum.hkl);
            spectrum.datestart = testCase.spectrum.datestart;
            spectrum.dateend = testCase.spectrum.dateend;
            out = sw_egrid(spectrum, 'sumtwin', false);

            expected_out = testCase.sw_egrid_out_sperp;
            expected_out.obj = spectrum.obj;
            expected_out.omega = spectrum.omega;
            expected_out.Sab = spectrum.Sab;
            expected_out.param.notwin = false;
            expected_out.param.sumtwin = false;
            expected_out.component = {'Sperp'};
            expected_out.swConv = cell(1, 3);
            expected_out.swConv{1} = testCase.sw_egrid_out_sperp.swConv;
            expected_out.swInt = cell(1, 3);
            expected_out.swInt{1} = testCase.sw_egrid_out_sperp.swInt;
            expected_out.intP = cell(1,3);
            expected_out.Pab = cell(1,3);
            expected_out.Mab = cell(1,3);
            expected_out.Sperp = {0.5*ones(2, 5) 0.875*ones(2, 5) 0.875*ones(2, 5)};

            testCase.verify_obj(out, expected_out);
        end
        function test_imagChk(testCase)
            dE = 1;
            testCase.spectrum.omega(1) = 0 + 2i*dE; % imag > dE bins 
            testCase.verifyError(...
                @() sw_egrid(testCase.spectrum, ...
                            'Evect', 0:dE:4, 'imagChk', true), ...
                'egrid:BadSolution');
        end
        function test_autoEmin(testCase)
            eps_imag = 1e-8;
            imag_omega = 0 + 1i*eps_imag;
            testCase.spectrum.omega(1) = imag_omega;
            component = 'Sxx';
            expected_out = testCase.sw_egrid_out;
            expected_out.component = component;
            expected_out.Evect(1) = expected_out.Evect(1) + eps_imag;
            expected_out.omega(1) = imag_omega;
            expected_out.swConv(1) = 0;
            out = sw_egrid(testCase.spectrum, 'component', component, 'autoEmin', true);
            testCase.verify_obj(out, expected_out);
        end
        function test_modeIdx(testCase)
            % only consisder -ve  mode (magnon anhillation/ energy gain)
            component = 'Sxx';
            expected_out = testCase.sw_egrid_out;
            expected_out.component = component;
            for modeIdx = 1:2
                out = sw_egrid(testCase.spectrum, 'component', component, ...
                           'modeIdx', modeIdx);
                if modeIdx == 2
                    % only consisder -ve  mode (magnon anhillation/energy gain)
                    % which is not in range of default energy bins
                    expected_out.swConv = zeros(500, 5);
                end
                testCase.verify_obj(out, expected_out);
            end
        end
        function test_zeroEnergyTol(testCase)
            % set zeroEnergyTol > max energy to prodcuce zero intensity
            out = sw_egrid(testCase.spectrum, ...
                           'binType', 'ebin' , 'zeroEnergyTol', 5);
            expected_out = testCase.sw_egrid_out_sperp;
            expected_out.swConv = zeros(size(expected_out.swConv));
            testCase.verify_obj(out, expected_out);
        end
        function test_negative_zeroEnergyTol(testCase)
            out = sw_egrid(testCase.spectrum, 'zeroEnergyTol', -1);
            expected_out = testCase.sw_egrid_out_sperp;
            expected_out.swConv([1,2001]) = 0.5; % intensity at zero energy
            testCase.verify_obj(out, expected_out);
        end
        function test_maxDSF(testCase)
            % set maxDSF low to zero all intensity
            out = sw_egrid(testCase.spectrum, ...
                           'binType', 'ebin' , 'maxDSF', 1e-2);
            expected_out = testCase.sw_egrid_out_sperp;
            expected_out.swConv = zeros(size(expected_out.swConv));
            testCase.verify_obj(out, expected_out);
        end
    end

end
