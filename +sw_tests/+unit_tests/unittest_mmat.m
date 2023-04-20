classdef unittest_mmat < sw_tests.unit_tests.unittest_super
    % Runs through unit test for sw_neutron.m

    properties
        swobj = [];
    end

    methods (Test)
        function test_mmat_branch(testCase)
            matA = rand(15,10);
            matB = rand(10,5,100);
            % Run it normally
            mat0 = mmat(matA, matB);
            % Force sw_freemem to return only 100 bytes available
            mock_freemem = sw_tests.utilities.mock_function('sw_freemem', 100);
            % Checks that with low memory the routine actually does not call bsxfun
            mock_bsxfun = sw_tests.utilities.mock_function('bsxfunsym');
            mat1 = mmat(matA, matB);
            testCase.verify_val(mat0, mat1, 'rel_tol', 0.01, 'abs_tol', 1e-6);
            testCase.assertEqual(mock_freemem.n_calls, 1);
            testCase.assertEqual(mock_bsxfun.n_calls, 0);
        end
    end

end

