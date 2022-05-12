classdef systemtest_spinwave_biquadratic < sw_tests.system_tests.systemtest_spinwave

    properties
        reference_data_file = 'systemstest_spinwave_biquadratic.mat';
    end

    methods (TestMethodSetup)
        function prepareForRun(testCase)
            % From Tutorial 28, to test the biquadratic interactions functionality, theory from PHYSICAL REVIEW B 85, 054409 (2012)
            S  = 5/2; J2 = 5.5; J1 = 5.0; Q = 0.01*J2;
            fcc = spinw;
            fcc.genlattice('lat_const',[8 8 8]);
            fcc.addatom('r',[0   0   0],'S',S);
            fcc.addatom('r',[1/2 1/2 0],'S',S);
            fcc.addatom('r',[1/2 0 1/2],'S',S);
            fcc.addatom('r',[0 1/2 1/2],'S',S);
            fcc.gencoupling();
            fcc.addmatrix('label','J1','value',J1,'color','b');
            fcc.addmatrix('label','J2','value',J2,'color','g');
            % there is a factor 2 difference between SpinW and paper
            fcc.addmatrix('label','B','value',-0.5*Q/S^3*2,'color','r');
            fcc.addcoupling('mat','J1','bond',1);
            fcc.addcoupling('mat','J2','bond',2);
            fcc.addcoupling('mat','B','bond',1,'type','biquadratic');
            fcc.genmagstr('mode','helical','S',[1 -1 -1 -1;1 -1 -1 -1;zeros(1,4)],'k',[1/2 1/2 1/2],'next',[2 2 2],'n',[0 0 1]);
            testCase.swobj = fcc;
        end
    end

    methods (Test)
        function test_biquadratic(testCase)
            fcc = testCase.swobj;
            spec = fcc.spinwave({[1 0 0] [0 0 0] [1/2 1/2 0] [1/2 1/2 1/2] [0 0 0] 50});
            spec = sw_egrid(spec);
            spec = sw_omegasum(spec,'zeroint',1e-5,'emptyval',0,'tol',1e-4);
            testCase.generate_or_verify(spec, {}, struct(), 'approxSab', 0.01);
        end
    end

end
