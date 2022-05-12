classdef systemtest_spinwave < matlab.unittest.TestCase
    % Base class for all systems test of spinwave.m based on tutorials

    properties
        generate_reference_data = false;
        reference_data = [];
        reference_data_dir = fullfile('.', 'test_data');
        relToll = 0.01;
        absToll = 1e-6;
        swobj = [];
    end

    methods (TestClassSetup)
        function get_reference_data(testCase)
            fname = fullfile(testCase.reference_data_dir, testCase.reference_data_file);
            if ~exist(testCase.reference_data_dir, 'dir')
                mkdir(testCase.reference_data_dir);
            end
            if ~exist(fname, 'file') || testCase.generate_reference_data
                testCase.generate_reference_data = true;
                tmp = []; save(fname, 'tmp');
                warning('Generating reference data');
            else
                testCase.reference_data = load(fname);
            end
        end
    end

    methods (TestClassTeardown)
        function save_reference_data(testCase)
            if testCase.generate_reference_data
                fname = fullfile(testCase.reference_data_dir, testCase.reference_data_file);
                ref_dat = load(fname);
                ref_dat = rmfield(ref_dat, 'tmp');
                save(fname, '-struct', 'ref_dat');
            end
        end
    end

    methods (Static)
        function out = get_hash(obj)
            % Calculates a hash for an object or struct using undocumented built-ins
            % Based on DataHash (https://uk.mathworks.com/matlabcentral/fileexchange/31272-datahash)
            Engine = java.security.MessageDigest.getInstance('MD5');
            Engine.update(getByteStreamFromArray(obj));
            out = typecast(Engine.digest, 'uint8');
        end
        function out = sanitize_data(in_dat)
            if isstruct(in_dat)
                out = sanitize_struct(in_dat);
            elseif iscell(in_dat)
                out = sanitize_cell(in_dat);
            elseif isnumeric(in_dat)
                out = sanitize(in_dat);
            else
                out = in_dat;
            end
        end
    end

    methods
        function out = approxMatrix(testCase, actual, expected, frac_not_match)
            % Checks if two arrays are approximately the same with most entries equal but a fraction not
            if iscell(actual)
                out = actual;
                for ii = 1:numel(actual)
                    out{ii} = testCase.approxMatrix(actual{ii}, expected{ii}, frac_not_match);
                end
            else
                diff = abs(actual - expected);
                rel_diff = diff ./ expected;
                if (numel(find((diff > testCase.absToll) & (rel_diff > testCase.relToll))) / numel(actual)) < frac_not_match
                    out = expected;
                else
                    out = actual;
                end
            end
        end
        function [actual, expected] = verify_eigval_sort(testCase, actual, expected, nested)
            if nargin < 4
                nested = 0;
            end
            if iscell(actual)
                for ii = 1:numel(actual)
                    [actual{ii}, expected{ii}] = testCase.verify_eigval_sort(actual{ii}, expected{ii}, nested);
                end
            else
                % Checks if actual and expected eigenvalues match, otherwise try a different sorting
                import matlab.unittest.constraints.*
                theseBounds = RelativeTolerance(testCase.relToll) | AbsoluteTolerance(testCase.absToll);
                comparator = IsEqualTo(expected, 'Within', theseBounds);
                if ~comparator.satisfiedBy(actual)
                    if nested > 1
                        actual = sort(abs(actual));
                        expected = sort(abs(expected));
                    else
                        actual = sort(actual, 'ComparisonMethod', 'real');
                        expected = sort(expected, 'ComparisonMethod', 'real');
                    end
                    if nested < 2
                        [actual, expected] = testCase.verify_eigval_sort(actual, expected, nested + 1);
                    end
                end
            end
        end
        function fieldname = get_fieldname(testCase, pars)
            if isempty(pars)
                fieldname = 'data';
            elseif ischar(pars);
                fieldname = pars;
            else
                fieldname = ['d' reshape(dec2hex(testCase.get_hash(pars)),1,[])];
            end
        end
        function save_test_data(testCase, data, pars)
            filename = fullfile(testCase.reference_data_dir, testCase.reference_data_file);
            tmpstr.(testCase.get_fieldname(pars)) = data;
            save(filename, '-append', '-struct', 'tmpstr');
        end
        function verify_test_data(testCase, test_data, ref_data)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            theseBounds = RelativeTolerance(testCase.relToll) | AbsoluteTolerance(testCase.absToll);
            test_data = testCase.sanitize_data(test_data);
            ref_data = testCase.sanitize_data(ref_data);
            testCase.verifyThat(test_data, IsEqualTo(ref_data, 'Within', theseBounds));
        end
        function generate_or_verify(testCase, spec, pars, extrafields, approxSab, tolSab)
            if nargin < 5
                approxSab = false;
            elseif nargin == 5
                tolSab = 0.05;
            end
            if testCase.generate_reference_data
                data.input = struct(testCase.swobj);
                data.spec = {spec.omega spec.Sab};
                if isfield(spec, 'swConv'); data.spec = [data.spec {spec.swConv}]; end
                if isfield(spec, 'swInt');  data.spec = [data.spec {spec.swInt}];  end
                if nargin > 3
                    extras = fieldnames(extrafields);
                    for ii = 1:numel(extras)
                        data.(extras{ii}) = extrafields.(extras{ii});
                    end
                end
                testCase.save_test_data(data, pars);
            else
                ref_data = testCase.reference_data.(testCase.get_fieldname(pars));
                test_data.input = struct(testCase.swobj);
                [spec.omega, ref_data.spec{1}] = testCase.verify_eigval_sort(spec.omega, ref_data.spec{1});
                test_data.spec = {spec.omega spec.Sab};
                if isfield(spec, 'swConv'); test_data.spec = [test_data.spec {spec.swConv}]; end
                if isfield(spec, 'swInt');  test_data.spec = [test_data.spec {spec.swInt}];  end
                if nargin > 3
                    extras = fieldnames(extrafields);
                    for ii = 1:numel(extras)
                        test_data.(extras{ii}) = extrafields.(extras{ii});
                    end
                end
                if any(approxSab)
                    % For the Sab or Sabp tensor, just check that a fraction of entries match
                    test_data.spec{2} = testCase.approxMatrix(spec.Sab, ref_data.spec{2}, tolSab);
                    if numel(test_data.spec) == 4
                        test_data.spec{4} = testCase.approxMatrix(spec.swInt, ref_data.spec{4}, tolSab);
                    end
                    if isfield(test_data, 'Sabp')
                        test_data.Sabp = testCase.approxMatrix(test_data.Sabp, ref_data.Sabp, tolSab);
                    end
                    if isfield(test_data, 'V')
                        test_data.V = testCase.approxMatrix(test_data.V, ref_data.V, tolSab);
                    end
                end
                testCase.verify_test_data(test_data, ref_data);
            end
        end
        function generate_or_verify_generic(testCase, data, fieldname)
            if testCase.generate_reference_data
                testCase.save_test_data(data, fieldname);
            else
                testCase.verify_test_data(data, testCase.reference_data.(fieldname));
            end
        end
    end

end

function out = sanitize(array)
    out = array;
    out(abs(out) > 1e8) = 0;
end

function sanitized = sanitize_struct(in_dat)
    fnam = fieldnames(in_dat);
    for ii = 1:numel(fnam)
        if isnumeric(in_dat.(fnam{ii}))
            sanitized.(fnam{ii}) = sanitize(in_dat.(fnam{ii}));
        elseif isstruct(in_dat.(fnam{ii}))
            sanitized.(fnam{ii}) = sanitize_struct(in_dat.(fnam{ii}));
        elseif iscell(in_dat.(fnam{ii}))
            sanitized.(fnam{ii}) = sanitize_cell(in_dat.(fnam{ii}));
        else
            sanitized.(fnam{ii}) = in_dat.(fnam{ii});
        end
    end
end

function sanitized = sanitize_cell(in_dat)
    sanitized = in_dat;
    for ii = 1:numel(in_dat)
        if isnumeric(in_dat{ii})
            sanitized{ii} = sanitize(in_dat{ii});
        elseif isstruct(in_dat{ii})
            sanitized{ii} = sanitize_struct(in_dat{ii});
        elseif iscell(in_dat{ii})
            sanitized{ii} = sanitize_cell(in_dat{ii});
        end
    end
end
