function addtwin(obj,varargin)
% adds crystallographic twins
% 
% ### Syntax
% 
% `addtwin(obj,Name,Value)`
% 
% ### Description
% 
% `addtwin(obj,Name,Value)` adds crystallographic twins defined by a
% rotation matrix and its volume fraction. Using crystallographic twins,
% SpinW can simulate imperfect samples and if the relative orientation of
% the crystallographic twins are knows, SpinW simulations can be directly
% compared to the expeiments on the inperfect sample.
% 
% ### Examples
% 
% This example shows how to add two extra crystallographic twins to the
% crystal.  Together with the original orientation there will be three
% twins with equal volumes.
%
% ```
% cryst.addtwin('axis',[0 0 1],'phid',[60 120],'vol',[1 1])
% ```
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
% 
% `'axis'`
% : Defines the axis of rotation to generate twins in the xyz
%   coordinate system, dimensions are $[1\times 3]$.
% 
% `'phi'`
% : Defines the angle of rotation to generate twins in radian
%   units. Several twins can be defined parallel if `phi` is a
%   row vector. Dimensions are $[1\times n_{twin}]$.
% 
% `'phid'`
% : Alternative to `phi` but the unit is \\deg.
% 
% `'rotC'`
% : Rotation matrices, that define crystallographic twins. This is an
%   alternative to the `axis`-`phi` parameter pair. Matrix dimensions are 
%   $[3\times 3\times n_{twin}]$.
% 
% `'vol'`
% : Volume fractions of the twins stored in a row vector with $n_{twin}$
%   elements. Default value is `ones(1,nTwin)`.
% 
% `'overwrite'`
% : If `true`, the last twin will be overwritten, instead of adding a
%   new one. Default is `false`.
% 
% ### Output Arguments
% 
% The function adds extra entries to the [spinw.twin] property.
% 
% ### See Also
% 
% [spinw] \| [spinw.twinq]
%

if nargin == 1
    swhelp spinw.addtwin
    return
end

inpForm.fname  = {'axis'  'phi'  'rotC'   'vol'  'phid' 'overwrite'};
inpForm.defval = {[0 0 0] 0      zeros(3) 1      0      false      };
inpForm.size   = {[1 3]   [1 -1] [3 3 -2] [1 -3] [1 -4] [1 1]      };

param = sw_readparam(inpForm, varargin{:});

if numel(param.phid)>1 || param.phid~=0
    param.phi = param.phid *pi/180;
end

% prefer axis definition over matrix
if any(param.axis)
    nTwin      = size(param.phi,2);
    param.rotC = zeros(3,3,nTwin);
    for ii = 1:nTwin
        [~, param.rotC(:,:,ii)] = sw_rot(param.axis,param.phi(ii));
    end
else
    % rotC matrix explicitly provided
    nTwin = size(param.rotC,3);
    for itwin = 1:nTwin
        %  check valid rotation or reflection (|R|=1 and R.T = R^-1)
        twin_rotc = param.rotC(:,:,itwin);
        if abs(det(twin_rotc))~=1 || ~isequal(twin_rotc*(twin_rotc'), eye(3))
            error('spinw:addtwin:WrongInput', ... 
                'rotC matrix for twin %d must be a valid rotation', itwin);
        end
    end

end

if size(param.vol,2)<nTwin
    param.vol = ones(1,nTwin);
end

if param.overwrite
    % overwrite the last twin
    obj.twin.vol(end)      = param.vol;
    obj.twin.rotc(:,:,end) = param.rotC;
else
    obj.twin.vol  = [obj.twin.vol param.vol];
    obj.twin.rotc = cat(3,obj.twin.rotc, param.rotC);
end

end