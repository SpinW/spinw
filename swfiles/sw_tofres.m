function spec = sw_tofres(spec, varargin)
% convolutes the spectrum with a Q bin
% 
% ### Syntax
% 
% `spectra = sw_tofres(spectra,Name,Value)`
% 
% ### Description
% 
% `spectra = sw_tofres(spectra,Name,Value)` simulates the finite bin size
% of the cuts of direct TOF neutron scattering data. It calculates the
% spectrum at multiple points within the given bin volume and sums them up.
% The function is usefull if relatively large bins were used to analyse the
% data due to low signal to noise ratio of the measurement.
% 
% ### Input Arguments
% 
% `spectra`
% : Input structure, contains calculated correlation functions
%   withouth the resolution effect.
% 
% ### Name-Value Pair Arguments
% 
% `'method'`
% : String that determines the method to generate the $Q$ points, options:
%   * `'random'`    The bin volume will be randomly sampled.
%   * `'grid'`      The bin volume will be split up to a finer regular
%                   grid.
% 
% `'dQ'`
% : Row vector with 3 elements. The width of the $Q$ bin
%   along the three reciprocal lattice directions. The spectrum
%   will be integrated in the $Q\pm (\delta Q/2)$ range. Default value is
%   `[0.1 0.1 0.1]`.
% 
% `'nQ'`
% : Row vector with 3 elements when `method` is `grid` and gives the
%   number of $Q$ points along the three reciprocal lattice directions to
%   average over. When `method` is `random` it is a scalar that determines
%   the number of random $Q$ points.
%
% `'fid'`
% : Defines whether to provide text output. The default value is determined
%   by the `fid` preference stored in [swpref]. The possible values are:
%   * `0`   No text output is generated.
%   * `1`   Text output in the MATLAB Command Window.
%   * `fid` File ID provided by the `fopen` command, the output is written
%           into the opened file stream.
%
% `'tid'`
% : Determines if the elapsed and required time for the calculation is
%   displayed. The default value is determined by the `tid` preference
%   stored in [swpref]. The following values are allowed (for more details
%   see [sw_timeit]):
%   * `0` No timing is executed.
%   * `1` Display the timing in the Command Window.
%   * `2` Show the timing in a separat pup-up window.
% 
% ### Output Arguments
% 
% `spectra`
% : Same as the input except that it contains the calculated intensity in
%   the `swConv` field.
% 
% ### See Also
% 
% [sw_egrid] \| [sw_instrument]
%
% *[TOF]: Time Of Flight
%


% help when executed without argument
if nargin==0
    swhelp sw_tofres
    return
end

if sw_issymspec(spec)
    error('sw_tofres:SymbolicInput', 'This function does not handle symbolic spectra');
end

dQ0 = ones(1,3)*0.1;
nQ0 = ones(1,3)*5;

inpForm.fname  = {'method'  'dQ'  'nQ' 'fid' 'tid'};
inpForm.defval = {'grid'    dQ0    nQ0 -1    -1   };
inpForm.size   = {[1 -1] [1 -2] [1 -3] [1 1] [1 1]};

param = sw_readparam(inpForm, varargin{:});
pref = swpref;
obj   = spec.obj;

if param.fid == -1
    % Print output into the following file
    fid = pref.fid;
else
    fid = param.fid;
end

if param.tid == -1
    param.tid = pref.tid;
end

hkl   = spec.hkl;
Evect = spec.Evect;
nE    = numel(Evect)-1;
dQ    = param.dQ;
nQ    = param.nQ;

switch param.method
    case 'random'
        
        if numel(param.dQ) == 1
            dQ = ones(1,3)*dQ;
        end
        
        conv0 = zeros(nE,size(hkl,2));
        
        fprintf0(fid,'Calculating TOF Q-binning using random method...\n')
        
        sw_timeit(0,1,param.tid,'TOF Q-binning resolution calculation');
        
        for ii = 1:prod(nQ)
            hklC = bsxfun(@minus,bsxfun(@times,dQ',rand(size(hkl))),dQ'/2)+hkl;
            spec0 = obj.spinwave(hklC,'formfact',spec.formfact,'fid',0);
            spec0 = sw_egrid(spec0,'Evect',Evect,'component',spec.component);
            conv0 = conv0 + spec0.swConv;
            
            sw_timeit(ii/prod(nQ)*100,0,param.tid);
        end
        
    case 'grid'
        
        
        if numel(param.dQ) == 1
            dQ = ones(1,3)*dQ;
        end
        
        if numel(nQ) == 1
            nQ = ones(1,3)*nQ;
        end
        
        % vector along each rl direction
        Qv = cell(1,3);
        for ii = 1:3
            Qv{ii}= linspace(-dQ(ii),dQ(ii),nQ(ii));
        end
        
        dhkl = cell(1,3);
        [dhkl{1},dhkl{2},dhkl{3}] = ndgrid(Qv{:});
        
        
        conv0 = zeros(nE,size(hkl,2));
        
        fprintf0(fid,'Calculating TOF Q-binning using grid method...\n')
        
        sw_timeit(0,1,param.tid,'TOF Q-binning resolution calculation');
        
        
        for ii = 1:prod(nQ)
            hklC = bsxfun(@plus,hkl,[dhkl{1}(ii); dhkl{1}(ii); dhkl{1}(ii)]);
            spec0 = obj.spinwave(hklC,'formfact',spec.formfact);
            spec0 = sw_egrid(spec0,'Evect',Evect,'component',spec.component);
            conv0 = conv0 + spec0.swConv;
            
            sw_timeit(ii/prod(nQ)*100,0,param.tid);
        end
        
    otherwise
        error('sw_tofres:WrongOption','Wrong method option!')
end

sw_timeit(100,2,param.tid);
obj.fileid(fid);
spec.swConv = conv0/prod(nQ);
fprintf0(fid,'Calculation finished.\n')

end
