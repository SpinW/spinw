classdef systemtest_spinwave_incommensurate_and_supercell_consistency < sw_tests.system_tests.systemtest_spinwave

    properties
        reference_data_file = [];
    end
    
    methods (Test)
        function test_AFM_kagome(testCase)
            % setup structure (taken from tutorial 8)
            AF33kagome = spinw;
            AF33kagome.genlattice('lat_const',[6 6 40], ...
                                  'angled',[90 90 120], 'sym','P -3')
            AF33kagome.addatom('r',[1/2 0 0],'S', 1, ...
                               'label','MCu1','color','r')
            AF33kagome.gencoupling('maxDistance',7);
            AF33kagome.addmatrix('label','J1','value',1.00)
            AF33kagome.addcoupling('mat','J1','bond',1);
            % sqrt3 x sqrt(3) magnetic structure 
            k = [-1/3 -1/3 0];
            n = [0, 0, 1];
            S = [0 0 -1; 1 1 -1; 0 0 0];
            % binning for spinwave spectrum
            qarg = {[-1/2 0 0] [0 0 0] [1/2 1/2 0] 3};
            evec = 0:0.5:2;
            
            % use structural unit cell with incommensurate k
            AF33kagome.genmagstr('mode','helical','unit','lu', 'k', k,...
                                 'n',n, 'S', S, 'nExt',[1 1 1]);
            spec_incom = AF33kagome.spinwave(qarg, 'hermit', true);
            spec_incom = sw_egrid(spec_incom, 'component','Sperp', 'Evect',evec);
            % use supercell k=0 structure
            AF33kagome.genmagstr('mode','helical','unit','lu', 'k', k,...
                                 'n',n, 'S', S, 'nExt',[3 3 1]);
            spec_super = AF33kagome.spinwave(qarg, 'hermit', true);
            spec_super = sw_egrid(spec_super, 'component','Sperp', 'Evect',evec);
            
            
            % test cross-section in q,En bins
            testCase.verify_test_data(spec_incom.swConv, ...
                                      spec_super.swConv)
            % check correct number of modes (2*nmag)
            nmode = 18;
            testCase.assertEqual(size(spec_incom.Sperp, 1), nmode)
            testCase.assertEqual(size(spec_super.Sperp, 1), 3*nmode)
        end
    end

end
