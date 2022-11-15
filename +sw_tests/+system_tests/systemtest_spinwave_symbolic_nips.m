classdef (TestTags = {'Symbolic'}) systemtest_spinwave_symbolic_nips < sw_tests.system_tests.systemtest_spinwave

    properties
        reference_data_file = '';
        swobj_nn = [];
        symspec = [];
    end

    properties (TestParameter)
        test_spectra_function_name = {'sw_neutron', 'sw_egrid', 'sw_instrument', ...
                                      'sw_omegasum', 'sw_plotspec', 'sw_tofres'};
    end

    methods (TestClassSetup)
        function prepareForRun(testCase)
            % Symbolic calculation, based on "Magnetic dynamics of NiPS3", A.R.Wildes et al., Phys. Rev. B in press
            nips = spinw();
            nips.genlattice('lat_const', [5.812, 10.222, 6.658], 'angled', [90, 107.16, 90], 'sym', 12);
            nips.addatom('r', [0 0.333 0], 'S', 1, 'label', 'MNi2');
            nips.gencoupling();
            nips.addmatrix('label', 'J1', 'mat', 1);
            nips.addcoupling('mat', 'J1', 'bond', 1);
            nips.addcoupling('mat', 'J1', 'bond', 2);
            nips.genmagstr('mode', 'direct', 'k', [0 1 0], 'S', [1 0 0; -1 0 0; -1 0 0; 1 0 0]');
            % The full system is too complicated to determine the general dispersion 
            % It will run for several hours and then run out of memory.
            % Instead for some tests we simplify it by including only nearest neighbour interactions
            % to be able to run through the full calculation.
            testCase.swobj_nn = nips.copy();
            testCase.swobj_nn.symbolic(true);
            testCase.symspec = testCase.swobj_nn.spinwavesym();
            nips.addmatrix('label', 'J2', 'mat', 1);
            nips.addcoupling('mat', 'J2', 'bond', 3);
            nips.addcoupling('mat', 'J2', 'bond', 4);
            nips.addmatrix('label', 'J3', 'mat', 1);
            nips.addcoupling('mat', 'J3', 'bond', 7);
            nips.addcoupling('mat', 'J3', 'bond', 8);
            nips.addmatrix('label', 'Jp', 'mat', 1);
            nips.addcoupling('mat', 'Jp', 'bond', 5);
            nips.addcoupling('mat', 'Jp', 'bond', 6);
            nips.addmatrix('mat', diag([0 0 1]), 'label','D');
            nips.addaniso('D');
            nips.symbolic(true);
            testCase.swobj = nips;
        end
    end

    methods (Test)
        function test_symbolic_hamiltonian(testCase)
            % Calculates the symbolic Hamiltonian and check it is Hermitian
            nips = testCase.swobj;
            % Specify 'eig', false to force not calculate dispersion (see general_hkl test below)
            symSpec = nips.spinwavesym('eig', false);
            ham = symSpec.ham;
            % Checks hamiltonian is hermitian
            import matlab.unittest.constraints.IsEqualTo
            testCase.verifyThat(simplify(ham - conj(transpose(ham))), IsEqualTo(sym(zeros(8))));
        end
        function test_symbolic_hamiltonian_white(testCase)
            % Calculates the symbolic Hamiltonian and check it obeys definition
            % of White et al., PR 139 A450 (1965)
            nips = testCase.swobj;
            symSpec = nips.spinwavesym('hkl', [0; 0; 0]);
            ham = symSpec.ham;
            g = sym(diag([ones(1, 4) -ones(1,4)]));
            Echeck = simplify(eig(g * symSpec.ham));
            import matlab.unittest.constraints.IsEqualTo
            testCase.verifyThat(simplify(symSpec.omega), IsEqualTo(Echeck));
        end
        function test_symbolic_hamiltonian_squared(testCase)
            % Calculates the symbolic Hamiltonian and its squared eigenvalues
            % agree with the block determinant identity
            nips = testCase.swobj;
            symSpec = nips.spinwavesym('hkl', [0.5; 0.5; 0.5]);
            g = sym(diag([ones(1, 4) -ones(1,4)]));
            hamsq = (g * symSpec.ham)^2;
            % Split squared hamiltonian into blocks - hamsq = [A B; C D]
            A = hamsq(1:4, 1:4);
            B = hamsq(1:4, 5:8);
            C = hamsq(5:8, 1:4);
            D = hamsq(5:8, 5:8);
            % Now the normal hamiltonian has the form: ham = [U V; V' U]
            % so when squared we get hamsq = [U^2+V^2 2*U*V; 2*U*V U^2+V^2] - e.g. [A B; B A]
            % This would satisfy the block determinant identity
            % det([A B; B A]) = det(A - B) * det(A + B)
            % So the eigenvalues of the full matrix [A B; B A] are those of [A-B] and [A+B]
            % (From wikipedia: https://en.wikipedia.org/wiki/Block_matrix#Block_matrix_determinant)
            % But actually for spin wave Hamiltonians, a stricter criteria applies, with B = -B'
            % so the eigenvalues of (A+B) is the same as (A-B)
            % These eigenvalues are the squared magnon frequencies here, and the positive (negative)
            % roots represents magnon creation (anihilation) modes.
            import matlab.unittest.constraints.IsEqualTo
            testCase.verifyThat(simplify(A), IsEqualTo(simplify(D)));
            testCase.verifyThat(simplify(B), IsEqualTo(simplify(C)));
            E1 = sort(simplify(eig(A + B)));
            E2 = sort(simplify(eig(A - B)));
            testCase.verifyThat(E1, IsEqualTo(E2));
            omega = sort(simplify(symSpec.omega.^2));
            testCase.verifyThat(unique(omega), IsEqualTo(unique(E1)));
        end
        function test_symbolic_general_hkl(testCase)
            % Test we can run the full spin wave calc outputing the spin-spin correlation matrix Sab
            variables = [sym('J1'), sym('h'), sym('k')]; % No 'l' because no out-of-plane coupling
            import matlab.unittest.constraints.IsEqualTo
            testCase.verifyThat(sort(symvar(testCase.symspec.omega)), IsEqualTo(variables));
        end
        function test_symbolic_spectra(testCase, test_spectra_function_name)
            % Tests that running standard functions with symbolic spectra gives error
            test_fun = eval(['@' test_spectra_function_name]);
            testCase.verifyError( ...
                @() feval(test_spectra_function_name, testCase.symspec), ...
                [test_spectra_function_name ':SymbolicInput']);
        end
    end

end
