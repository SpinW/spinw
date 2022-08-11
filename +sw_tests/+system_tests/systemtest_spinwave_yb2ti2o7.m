classdef systemtest_spinwave_yb2ti2o7 < sw_tests.system_tests.systemtest_spinwave

    properties
        reference_data_file = 'systemstest_spinwave_yb2ti2o7.mat';
    end

    properties (TestParameter)
        B = {2 5};
        Q = {{[-0.5 -0.5 -0.5] [2 2 2]} {[1 1 -2] [1 1 1.5]} {[2 2 -2] [2 2 1.5]} {[-0.5 -0.5 0] [2.5 2.5 0]} {[0 0 1] [2.3 2.3 1]}};
    end

    methods (TestMethodSetup)
        function prepareForRun(testCase)
            % From Tutorial 20, Yb2Ti2O7, based on PRX 1, 021002 (2011)
            % First set up the crystal structure
            symStr = '-z, y+3/4, x+3/4; z+3/4, -y, x+3/4; z+3/4, y+3/4, -x; y+3/4, x+3/4, -z; x+3/4, -z, y+3/4; -z, x+3/4, y+3/4';
            yto = spinw;
            a = 10.0307;
            yto.genlattice('lat_const',[a a a],'angled',[90 90 90],'sym',symStr,'label','F d -3 m Z')
            yto.addatom('label','Yb3+','r',[1/2 1/2 1/2],'S',1/2)
            % We generate the list of bonds.
            yto.gencoupling
            % We create two 3x3 matrix, one for the first neighbor anisotropic exchange
            % and one for the anisotropic g-tensor. And assign them appropriately.
            yto.addmatrix('label', 'J1', 'color', [255 0 0], 'value', 1)
            yto.addmatrix('label', 'g0', 'color', [0 0 255], 'value', -0.84*ones(3)+4.32*eye(3));
            yto.addcoupling('mat', 'J1', 'bond', 1)
            yto.addg('g0')
            % Sets the correct values for the matrix elements of J1
            J1 = -0.09; J2 = -0.22; J3 = -0.29; J4 = 0.01;
            yto.setmatrix('mat','J1','pref',[J1 J3 J2 -J4]);
            testCase.swobj = yto;
        end
    end

    methods (Test)
        function test_yto(testCase, B, Q)
            n = [1 -1 0];
            % set magnetic field
            testCase.swobj.field(n/norm(n)*B);
            % create fully polarised magnetic structure along the field direction
            testCase.swobj.genmagstr('S',n','mode','helical');
            % find best structure using steepest descendend
            testCase.swobj.optmagsteep;
            ytoSpec = testCase.swobj.spinwave([Q {50}],'gtensor',true);
            ytoSpec = sw_neutron(ytoSpec);
            % bin the spectrum in energy
            ytoSpec = sw_egrid(ytoSpec,'Evect',linspace(0,2,100),'component','Sperp');
            %figure; sw_plotspec(ytoSpec,'axLim',[0 0.5],'mode',3,'dE',0.09,'colorbar',false,'legend',false); title(''); caxis([0 60]); colormap(jet);
            testCase.generate_or_verify(ytoSpec, {B Q}, struct(), 'approxSab', 0.5);
        end
        function test_yto_twin(testCase)
            % Adds a twin and runs test with single field/Q
            testCase.swobj.addtwin('axis', [1 -1 0], 'phid', 90);
            testCase.test_yto(4, {[-0.5 -0.5 -0.5] [2 2 2]});
        end
    end

end
