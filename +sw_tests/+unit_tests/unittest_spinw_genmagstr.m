classdef unittest_spinw_genmagstr < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
    end
    
    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw(); % default init
        end
    end

    methods (Test)
        function test_no_magnetic_atoms_raises_error(testCase)
            testCase.verifyError(...
                @() testCase.swobj.genmagstr(), ...
                'spinw:genmagstr:NoMagAtom')
        end
    end

end
