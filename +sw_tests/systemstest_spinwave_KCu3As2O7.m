classdef systemstest_spinwave_KCu3As2O7 < sw_tests.systemstest_spinwave

    properties
        reference_data_file = 'systemstest_spinwave_KCu3As2O7.mat';
    end

    methods (TestMethodSetup)
        function prepareForRun(testCase)
            % From Tutorial 18, a distorted kagome lattice from PRB 89, 140412 (2014)
            % To test incommensurate spin wave calculations and also the structure optimisation routine
            J   = -2; Jp  = -1; Jab = 0.75; Ja  = -J/.66 - Jab; Jip = 0.01;
            hK = spinw;
            hK.genlattice('lat_const',[10.2 5.94 7.81],'angled',[90 117.7 90],'spgr','C 2/m');
            hK.addatom('r',[0   0   0],'S',1/2,'label','MCu2','color','b');
            hK.addatom('r',[1/4 1/4 0],'S',1/2,'label','MCu2','color','k');
            hK.gencoupling();
            hK.addmatrix('label','J-',  'color','r',   'value',J);
            hK.addmatrix('label','J''','color','g',   'value',Jp);
            hK.addmatrix('label','Ja', 'color','b',   'value',Ja);
            hK.addmatrix('label','Jab','color','cyan','value',Jab);
            hK.addmatrix('label','Jip','color','gray','value',Jip);
            hK.addcoupling('mat','J-','bond',1);
            hK.addcoupling('mat','J''','bond',2);
            hK.addcoupling('mat','Ja','bond',3);
            hK.addcoupling('mat','Jab','bond',5);
            hK.addcoupling('mat','Jip','bond',10);
            testCase.swobj = hK;
        end
    end

    methods (Test)
        function test_KCu3As2O7(testCase)
            hK = testCase.swobj;
            hK.genmagstr('mode','helical','n',[0 0 1],'S',[1 0 0]','k',[0.77 0 0.115],'next',[1 1 1]);
            optpar.func = @gm_planar;
            optpar.nRun = 10;
            optpar.xmin = [    zeros(1,6), 0.5 0 0.0, 0 0];
            optpar.xmax = [2*pi*ones(1,6), 1.0 0 0.5, 0 0];
            magoptOut = hK.optmagstr(optpar);
            opt_energy = hK.energy;
            % Optmised structure with optmagstr not constant enough for check (will vary within a phase factor)
            % Use optmagsteep structure instead, but check its ground state energy.
            hK.genmagstr('mode','helical','n',[0 0 1],'S',[1 0 0]','k',[0.77 0 0.115],'next',[1 1 1]);
            magoptOut = hK.optmagsteep('nRun', 100);
            hkSpec = hK.spinwave({[0 0 0] [1 0 0] 100},'hermit',false);
            hkSpec = sw_neutron(hkSpec);
            hkSpec = sw_egrid(hkSpec,'Evect',linspace(0,5,100),'imagChk',false);
            % Remove problematic indices
            hkSpec.Sab(:,:,7,[97 99]) = 0;
            hkSpec.swInt(7,[97 99]) = 0;
            testCase.generate_or_verify(hkSpec, {}, struct('opt_energy', opt_energy, 'energy', hK.energy), 'approxSab', 0.5);
        end
    end

end
