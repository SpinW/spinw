classdef systemtest_spinwave_incommensurate_and_supercell_consistency < sw_tests.system_tests.systemtest_spinwave

    properties
        reference_data_file = [];
        tol = 1e-5;
    end
    methods (Static)
        function om = remove_ghosts(spec, tol)
            om = spec.omega(find(abs(spec.Sperp) > tol));
            om = sort(unique(round(om / tol) * tol));
        end
    end
    methods
        function assert_super_and_incom_consistency(testCase, swobj, ...
                                                    spec_super, spec_incom, ...
                                                    ghost_tol)
            if nargin < 5
                ghost_tol = testCase.tol;
            end
            % test cross-section in q,En bins
            testCase.verify_test_data(spec_incom.swConv, ...
                                      spec_super.swConv)
            % check correct number of modes (2*nmag)
            nExt = swobj.mag_str.nExt;
            n_matom = numel(swobj.matom().S);
            testCase.assertEqual(size(spec_incom.Sperp, 1), ...
                                 6*n_matom);
            testCase.assertEqual(size(spec_super.Sperp, 1), ...
                                 2*prod(nExt)*n_matom);
            testCase.assertEqual(testCase.remove_ghosts(spec_super, ghost_tol), ...
                                 testCase.remove_ghosts(spec_incom, ghost_tol));
        end
    end
    
    methods (Test)
        function test_AFM_kagome(testCase)
            % setup structure (taken from tutorial 8)
            AF33kagome = spinw;
            AF33kagome.genlattice('lat_const',[6 6 40], ...
                                  'angled',[90 90 120], 'sym','P -3')
            AF33kagome.addatom('r',[1/2 0 0],'S', 1, 'label','MCu1')
            AF33kagome.gencoupling('maxDistance',7);
            AF33kagome.addmatrix('label','J1','value',1)
            AF33kagome.addcoupling('mat','J1','bond',1);
            % sqrt3 x sqrt(3) magnetic structure 
            k = [-1/3 -1/3 0];
            n = [0, 0, 1];
            S = [0 0 -1; 1 1 -1; 0 0 0];
            % binning for spinwave spectrum
            qarg = {[-1/2 0 0] [0 0 0] [1/2 1/2 0] 50};
            evec = 0:0.1:1.5;
            
            % use structural unit cell with incommensurate k
            AF33kagome.genmagstr('mode','helical','unit','lu', 'k', k,...
                                 'n',n, 'S', S, 'nExt',[1 1 1]);
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            spec_incom = AF33kagome.spinwave(qarg, 'hermit', true);
            spec_incom = sw_egrid(spec_incom, 'component','Sperp', 'Evect', evec, ...
                                  'zeroEnergyTol', 1e-2);
            % use supercell k=0 structure
            AF33kagome.genmagstr('mode','helical','unit','lu', 'k', k,...
                                 'n',n, 'S', S, 'nExt', [3,3,1]);

            spec_super = AF33kagome.spinwave(qarg, 'hermit', true);
            spec_super = sw_egrid(spec_super, 'component','Sperp', 'Evect', evec);
            
            testCase.assert_super_and_incom_consistency(AF33kagome, ...
                                                        spec_super, ...
                                                        spec_incom, 5e-2);
        end
        
        function test_two_matom_per_unit_cell(testCase)
            % setup structure (taken from tutorial 19)
            FeCuChain = spinw;
            FeCuChain.genlattice('lat_const',[3 8 4],'sym','P 1')
            FeCuChain.addatom('label','MCu2','r',[0 0 0])
            FeCuChain.addatom('label','MFe2','r',[0 1/2 0])

            FeCuChain.gencoupling
            FeCuChain.addmatrix('label','J_{Cu-Cu}','value',1)
            FeCuChain.addmatrix('label','J_{Fe-Fe}','value',1)
            FeCuChain.addmatrix('label','J_{Cu-Fe}','value',-0.1)

            FeCuChain.addcoupling('mat','J_{Cu-Cu}','bond',1)
            FeCuChain.addcoupling('mat','J_{Fe-Fe}','bond',2)
            FeCuChain.addcoupling('mat','J_{Cu-Fe}','bond',[4 5])
            % AFM  structure 
            k = [1/2, 0, 0];
            S = [0 0;1 1;0 0];
            % binning for spinwave spectrum
            qarg = {[0 0 0] [0, 0.5, 0] 5};
            evec = 0:0.5:5;
            
            % use structural unit cell with incommensurate k
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian', ...
                                      'spinw:magstr:NotExact', ...
                                      'spinw:spinwave:Twokm');
            FeCuChain.genmagstr('mode','helical','k', k,...
                                'S', S, 'nExt',[1 1 1]);
            spec_incom = FeCuChain.spinwave(qarg, 'hermit', true);
            spec_incom = sw_egrid(spec_incom, 'component','Sperp', 'Evect',evec);
            % use supercell k=0 structure
            FeCuChain.genmagstr('mode','helical','k', k,...
                                'S', S, 'nExt', [2,1,1]);
            spec_super = FeCuChain.spinwave(qarg, 'hermit', true);
            spec_super = sw_egrid(spec_super, 'component','Sperp', 'Evect',evec);
            
            testCase.assert_super_and_incom_consistency(FeCuChain, ...
                                                        spec_super, ...
                                                        spec_incom);
        end
        
        function test_two_sym_equiv_matoms_per_unit_cell(testCase)
            sw = spinw;
            sw.genlattice('lat_const',[4,5,12],'sym','I m m m')
            sw.addatom('S', 1, 'r',[0 0 0]);
            sw.addmatrix('label', 'J', 'value', 1);
            sw.addmatrix('label', 'A', 'value', diag([0 0 -0.1]))
            sw.gencoupling;
            sw.addcoupling('mat','J','bond', 1);
            sw.addaniso('A');
            % AFM structure
            S = [0; 0; 1];
            k = [0.5, 0 0];
            % binning for spinwave spectrum
            qarg = {[0 0 0] [1/2 0 0] 5};
            evec = 0:0.5:1.5;
            
            % use structural unit cell with incommensurate k
            testCase.disable_warnings('spinw:genmagstr:SnParallel');
            sw.genmagstr('mode','helical','k', k,...
                                'S', S, 'nExt',[1 1 1]);
            spec_incom = sw.spinwave(qarg, 'hermit', true);
            spec_incom = sw_egrid(spec_incom, 'component','Sperp', 'Evect', evec);
            % use supercell k=0 structure
            sw.genmagstr('mode','helical','k', k,...
                                'S', S, 'nExt', [2 1 1]);
            spec_super = sw.spinwave(qarg, 'hermit', true);
            spec_super = sw_egrid(spec_super, 'component','Sperp', 'Evect', evec);
            
            testCase.assert_super_and_incom_consistency(sw, ...
                                                        spec_super, ...
                                                        spec_incom);
        end
    end

end
