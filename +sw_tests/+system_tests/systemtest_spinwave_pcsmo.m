classdef systemtest_spinwave_pcsmo < sw_tests.system_tests.systemtest_spinwave

    properties
        reference_data_file = 'systemstest_spinwave_pcsmo.mat';
    end

    properties (TestParameter)
        usehorace = {false true};
    end

    methods (TestMethodSetup)
        function prepareForRun(testCase)
            % Runs the Goodenough model for (Pr,Ca)SrMn2O7, based on PRL, 109, 237202 (2012)
            JF1 = -11.39; JA = 1.5; JF2 = -1.35; JF3 = 1.5; Jperp = 0.88; D = 0.074;
            lat = [5.408 5.4599 19.266]; alf = [90 90 90];
            SM4 = 7/4;   % Spin length for Mn4+
            SM3 = 7/4;   % Spin length for Mn3+
            pcsmo = spinw;
            pcsmo.genlattice('lat_const', lat.*[2 2 1], 'angled', alf, 'sym', 'x,y+1/2,-z');
            [~,ffn3] = sw_mff('MMn3');
            [~,ffn4] = sw_mff('MMn4');
            myaddatom3 = @(x,y,z) pcsmo.addatom('label', x, 'r', y, 'S', SM3, 'color', z, ...
                'formfactn', ffn3, 'formfactx', 'MMn3', 'Z', 25, 'b', sw_nb('MMn3'));
            myaddatom4 = @(x,y,z) pcsmo.addatom('label', x, 'r', y, 'S', SM4, 'color', z, ...
                'formfactn', ffn4, 'formfactx', 'MMn4', 'Z', 25, 'b', sw_nb('MMn4'));
            myaddatom4('Mn4-up', [0 0 0.1], 'gold');
            myaddatom4('Mn4-up', [0.5 0.5 0.1], 'gold');
            myaddatom4('Mn4-dn', [0 0.5 0.1], 'gold');
            myaddatom4('Mn4-dn', [0.5 0 0.1], 'gold');
            myaddatom3('Mn3-up', [0.25 0.75 0.1], 'black');
            myaddatom3('Mn3-up', [0.75 0.75 0.1], 'black');
            myaddatom3('Mn3-dn', [0.25 0.25 0.1], 'black');
            myaddatom3('Mn3-dn', [0.75 0.25 0.1], 'black');
            % Generate the CE magnetic structure
            S0 = [0; 1; 0];
            spin_up = find(~cellfun(@isempty, strfind(pcsmo.table('matom').matom, 'up')));
            spin_dn = find(~cellfun(@isempty, strfind(pcsmo.table('matom').matom, 'dn')));
            SS = zeros(3, 16);
            SS(:, spin_up) = repmat(S0, 1, numel(spin_up));
            SS(:, spin_dn) = repmat(-S0, 1, numel(spin_dn));
            pcsmo.genmagstr('mode', 'direct', 'S', SS)
            % Generate the exchange interactions
            pcsmo.gencoupling('forceNoSym', true)
            pcsmo.addmatrix('label', 'JF1', 'value', JF1, 'color', 'green');
            pcsmo.addmatrix('label', 'JA', 'value', JA, 'color', 'yellow');
            pcsmo.addmatrix('label', 'JF2', 'value', JF2, 'color', 'white');
            pcsmo.addmatrix('label', 'JF3', 'value', JF3, 'color', 'red');
            pcsmo.addmatrix('label', 'Jperp', 'value', Jperp, 'color', 'blue');
            pcsmo.addmatrix('label', 'D', 'value', diag([0 0 D]), 'color', 'white');
            % The zig-zag chains couple Mn3-Mn4 with same spin.
            pcsmo.addcoupling('mat', 'JF1', 'bond', 1, 'atom', {'Mn3-up', 'Mn4-up'})
            pcsmo.addcoupling('mat', 'JF1', 'bond', 1, 'atom', {'Mn3-dn', 'Mn4-dn'})
            % And vice-versa for the inter-chain interaction
            pcsmo.addcoupling('mat', 'JA', 'bond', 1, 'atom', {'Mn3-up', 'Mn4-dn'})
            pcsmo.addcoupling('mat', 'JA', 'bond', 1, 'atom', {'Mn3-dn', 'Mn4-up'})
            pcsmo.addcoupling('mat', 'Jperp', 'bond', 2)
            % JF3 couples Mn3 within the same zig-zag (same spin)
            pcsmo.addcoupling('mat', 'JF3', 'bond', 3, 'atom', 'Mn3-up')
            pcsmo.addcoupling('mat', 'JF3', 'bond', 3, 'atom', 'Mn3-dn')
            % Find indexes of the Mn4+ atoms which have a=0.5:
            idmid = find((~cellfun(@isempty, strfind(pcsmo.table('matom').matom, 'Mn4'))) ...
                .* (pcsmo.table('matom').pos(:,1)==0.5));
            bond8 = pcsmo.table('bond', 8);
            % Finds the bonds which start on one of these atoms and goes along +b
            idstart = find(ismember(bond8.idx1, idmid) .* (bond8.dr(:,2)>0));
            % Finds the bonds which ends on one of these atoms and goes along -b
            idend = find(ismember(bond8.idx2, idmid) .* (bond8.dr(:,2)<0));
            pcsmo.addcoupling('mat', 'JF2', 'bond', 8, 'subIdx', [idstart; idend]')
            pcsmo.addaniso('D')
            % Define twins
            pcsmo.addtwin('rotC', [0 1 0; 1 0 0; 0 0 1]);
            pcsmo.twin.vol = [0.5 0.5];
            pcsmo.unit.qmat = diag([2 2 1]);
            % Assign to property
            testCase.swobj = pcsmo;
            testCase.absToll = 2e-6;
        end
    end

    methods (Test)
        function test_pcsmo(testCase, usehorace)
            % (002) is problematic - Goldstone mode gives indexing error in different Matlab versions
            qln = {[0 0 0] [1.98 0 0] 50};
            if usehorace
                hkl = sw_qscan(qln);
                % spinwavefast doesn't work for twins yet
                testCase.swobj.notwin();
                [w0, s0] = testCase.swobj.horace(hkl(:,1), hkl(:,2), hkl(:,3), 'dE', 2, 'useFast', false);
                [w1, s1] = testCase.swobj.horace(hkl(:,1), hkl(:,2), hkl(:,3), 'dE', 2, 'useFast', true);
                testCase.verify_test_data({w0(1:numel(w1)) s0(1:numel(w1))}, {w1 s1});
                testCase.generate_or_verify_generic({w0 s0}, 'data_horace');
            else
                spec = testCase.swobj.spinwave(qln, 'formfact', true, 'saveV', true, 'saveH', true, 'optmem', 2);
                spec = sw_egrid(spec, 'Evect', linspace(0, 100, 200));
                spec = sw_neutron(spec);
                testCase.generate_or_verify(spec, {}, struct('V', spec.V, 'H', spec.H), 'approxSab', 0.1);
            end
        end
    end

end
