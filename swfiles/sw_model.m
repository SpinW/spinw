function obj = sw_model(model, param, fid)
% creates predefined spin models
% 
% ### Syntax
% 
% `obj = sw_model(model, param)`
% 
% ### Description
% 
% `obj = sw_model(model, param)` generates spin models, such as triangular
% lattice antiferromagnet, square lattice, etc. It also generates the
% magnetic ground state. For each lattice arbitrary number of further
% neighbor bonds can be defined using a vector of exchange values.
% 
% ### Input Arguments
% 
% `model`
% : String, name of the model, one of the following:
%   * `'triAF'`     Triangular lattice Heisenberg antiferromagnet
%                   in the $ab$ plane ($a=b=3$ \\ang), with \\gamma =
%                   120\\deg angle and optimised magnetic structure.
%   * `'squareAF'`  Square lattice antiferromagnet.
%   * `'chain'`     Chain with further neighbor interactions.
%   * `swm_*`       Custom models which are in the matlab path can be 
%                 evaluated. Checkout: 
%                 https://www.github.com/spinw/Models for pre-made models.
%
% `param`
% : Input parameters of the model, row vector which gives the values of the
%   Heisenberg exchange for first, second, thirs etc. neighbor bonds stored
%   in `p(1)`, `p(2)`, `p(3)`, etc. respectively.
% 
% ### Output Arguments
% 
% `obj`
% : [spinw] class object with the selected model.
% 
% ### See Also
% 
% [spinw]
%

if nargin == 0
    swhelp sw_model
    return
end

pref = swpref;

if nargin < 3
    fid = pref.fid;
end

fprintf0(fid,'Preparing ''%s'' model ...\n',model);
modelSearch = 'swm_';

obj = spinw;

switch model
    case 'triAF'
        obj.genlattice('lat_const',[3 3 9],'angled',[90 90 120])
        obj.addatom('r',[0 0 0],'S',1,'color','darkmagenta')
        obj.gencoupling('maxDistance',10)
        if nargin > 1
            for ii = 1:numel(param)
                obj.addmatrix('value',param(ii),'label',['J_' num2str(ii)],'color',swplot.color(randi(140),1))
                obj.addcoupling('mat',ii,'bond',ii)
            end
            
            %obj.lattice.lat_const(3) = 5;
            obj.genmagstr('mode','helical','S',[1 0 0]','k',[0 0 0])
            obj.optmagstr('func',@gm_planar,'xmin',[0 0 0 0 0 0],'xmax',[0 1/2 1/2 0 0 0],'nRun',10)
        end
    case 'squareAF'
        obj.genlattice('lat_const',[3 3 9],'angled',[90 90 90])
        obj.addatom('r',[0 0 0],'S',1,'color','darkmagenta')
        obj.gencoupling('maxDistance',10)
        if nargin>1
            for ii = 1:numel(param)
                obj.addmatrix('value',param(ii),'label',['J_' num2str(ii)],'color',swplot.color(randi(140),1))
                obj.addcoupling('mat',ii,'bond',ii)
            end
            
            %obj.lattice.lat_const(3) = 5;
            obj.genmagstr('mode','helical','S',[1 0 0]','k',[0 0 0])
            obj.optmagstr('func',@gm_planar,'xmin',[0 0 0 0 0 0],'xmax',[0 1/2 1/2 0 0 0],'nRun',10)
            tol = 2e-4;
            helical = sum(abs(mod(abs(2*obj.magstr.k)+tol,1)-tol).^2) > tol;
            if ~helical
                obj.genmagstr('mode','helical','nExt',[2 2 1]);
            end
        end
    case 'chain'
        obj.genlattice('lat_const',[3 9 9],'angled',[90 90 90])
        obj.addatom('r',[0 0 0],'S',1,'color','darkmagenta')
        obj.gencoupling('maxDistance',10)
        if nargin>1
            for ii = 1:numel(param)
                obj.addmatrix('value',param(ii),'label',['J_' num2str(ii)],'color',swplot.color(randi(140),1))
                obj.addcoupling('mat',ii,'bond',ii)
            end
            
            %obj.lattice.lat_const(2:3) = 5;
            obj.genmagstr('mode','helical','S',[1 0 0]','k',[0 0 0])
            obj.optmagstr('func',@gm_planar,'xmin',[0 0 0 0 0 0],'xmax',[0 1/2 0 0 0 0],'nRun',10)
            tol = 2e-4;
            helical = sum(abs(mod(abs(2*obj.magstr.k)+tol,1)-tol).^2) > tol;
            if ~helical && any(obj.magstr.k>tol)
                nExt = [1 1 1];
                nExt(obj.magstr.k(1:2)>tol) = 2;
                obj.genmagstr('mode','helical','nExt',nExt);
            end
        end
    otherwise
        % All paths
        allPaths =  strsplit(path,':');
        % Remove the builtin functions
        allPaths =  allPaths(~strncmp(allPaths,fullfile(matlabroot, 'toolbox'), length(fullfile(matlabroot, 'toolbox'))));
        relPaths = allPaths(~strncmp(allPaths,fullfile(matlabroot, 'example'), length(fullfile(matlabroot, 'example'))));
        % Get all files
        allFiles = cellfun(@(x) dir(x), relPaths, 'UniformOutput', false);
        allFiles = cellfun(@(x) {x.name}, allFiles, 'UniformOutput', false);
        allFiles = [allFiles{:}];
        % Search for files which are models
        relFiles = cellfun(@(x) x(1:end-2) ,allFiles(strncmp(allFiles, modelSearch, 3)),'UniformOutput',false);
        if ~any(strcmp(model, relFiles))
            error('sw_model:WrongInput','Model does not exists!')
        end
        % Evaluate the model
        try
            obj = feval(model, param);        
        catch ME
            error('sw_model:ModelError','This model has an error!')
        end
end

fprintf0(fid,'... ready!\n');

end