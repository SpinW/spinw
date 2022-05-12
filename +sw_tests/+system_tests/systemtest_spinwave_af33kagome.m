classdef systemtest_spinwave_af33kagome < sw_tests.system_tests.systemtest_spinwave

    properties
        reference_data_file = 'systemstest_spinwave_af33kagome.mat';
    end

    methods (TestMethodSetup)
        function prepareForRun(testCase)
            % From Tutorial 8, a sqrt(3) x sqrt(3) kagome AFM to test incommensurate calculations
            AF33kagome = spinw;
            AF33kagome.genlattice('lat_const',[6 6 40],'angled',[90 90 120],'spgr','P -3');
            AF33kagome.addatom('r',[1/2 0 0],'S', 1,'label','MCu1','color','r');
            AF33kagome.gencoupling('maxDistance',7);
            AF33kagome.addmatrix('label','J1','value',1.00,'color','g');
            AF33kagome.addcoupling('mat','J1','bond',1);
            testCase.swobj = AF33kagome;
            testCase.relToll = 0.027;
            testCase.absToll = 1.2e-5;
        end
    end

    methods (Test)
        function test_af33kagome(testCase)
            AF33kagome = testCase.swobj;
            S0 = [0 0 -1; 1 1 -1; 0 0 0];
            AF33kagome.genmagstr('mode','helical','k',[-1/3 -1/3 0],'n',[0 0 1],'unit','lu','S',S0,'nExt',[1 1 1]);
            kag33Spec = AF33kagome.spinwave({[-1/2 0 0] [0 0 0] [1/2 1/2 0] 100},'hermit',false,'saveSabp',true);
            kag33Spec = sw_egrid(kag33Spec,'component','Sxx+Syy','imagChk',false, 'Evect', linspace(0, 3, 100));
            % Reduce values of S(q,w) so it falls within tolerance (rather than change tolerance for all values)
            kag33Spec.swConv = kag33Spec.swConv / 2e5;
            % Ignores swInt in this case
            kag33Spec.swInt = 0;
            testCase.generate_or_verify(kag33Spec, {}, struct('energy', AF33kagome.energy, 'Sabp', kag33Spec.Sabp), 'approxSab', 0.5);
        end
    end

end
