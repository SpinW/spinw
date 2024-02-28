function sw_spec2MDHisto(spectra,proj, dproj, filename)
% saves spectrum to MDHisto
% 
% ### Syntax
% 
% sw_spec2MDHisto(spectra,proj,dproj,filename)`
% 
% ### Description
% 
% `sw_spec2MDHisto(spectra,proj,dproj,filename)` saves a 
% spectrum that is calculated by sw_egrid
%  
% ### Input Arguments
% spectra: a structure calculated by sw_egrid 
% 
% proj: a 3x3 matrix defining an orthogonal coordinate system 
%       where each column is a vector defining the orientation 
%       of the view. One of the vectors must be along the Q axis 
%       defined by the direction of the calculation. It is also used to define the units along the x axis.  
% 
% dproj: is a 3 vector that is the bin size in each of the 
%        directions defined in proj. For the direction of the 
%        calculation, the value used is internally calcualted from the spectrum.
%        It is wise to enter the step size for clarity.  
%       
% filename: is the name of the nexus file.  It will overwrite the existing
%           file if one already exists
%
% Example:
% q0 = [0 0 0];
% qmax = [2 0 0];
% nsteps = 100;
% spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 qmax nsteps}))
% proj = [[1 0 0]' [0 1 0]' [0 0 1]'];
% dproj = [1, 1e-6, 1e-6];
% sw_spec2MDHisto(spec, proj, dproj, 'testmdh.nxs');
% Note that: 
% (1) In the call to `spinwave`, only one q-direction may be specified
%    e.g. the HKL specifier must be of the form {q0 q0+qdir nsteps} 
% (2) one column in the `proj` matrix must be the q-direction used in 
%    `spinwave` (e.g. `qdir`).


if nargin==0
    swhelp sw_spec2MDHisto
    return
end
[unit_cell,Bmat,proj_out,D,dat,proj,name] = read_struct(spectra,proj,dproj);
%check if hdf file exists and delete if it does.
if exist(filename,'file')
    delete(filename)
end

h5createnwrite(filename,'/MDHistoWorkspace/coordinate_system',3,0); %  None = 0, QLab = 1, QSample = 2, HKL = 3 
h5createnwrite(filename,'/MDHistoWorkspace/visual_normalization',0,0);
h5writeatt(filename,'/MDHistoWorkspace','NX_class','NXentry');
h5writeatt(filename,'/MDHistoWorkspace','Qconvention','Inelastic');
h5writeatt(filename,'/MDHistoWorkspace','SaveMDVersion', 2);
% write data
rtpth = NXScreategroup(filename,'/MDHistoWorkspace','data','NXdata');
% write D dimensions
Dszs=zeros(1,4);
Dstrcell={};
for idx=1:length(D)
    szd=size(D{idx});
    Dszs(idx)=szd(2)-1;
    Dstrcell{idx} = strcat('D',num2str(idx-1));
    Dpth = strcat(rtpth,'/',Dstrcell{idx});
    h5createnwritevec(filename,rtpth,Dstrcell{idx},D{idx});
    if idx<length(D)
        h5writeatt(filename,Dpth,'units','r.l.u');
        h5writeatt(filename,Dpth,'frame',mat2str(transpose(proj(:,idx))));
        h5writeatt(filename,Dpth,'long_name',mat2str(transpose(proj(:,idx))));
    else
        h5writeatt(filename,Dpth,'long_name','DeltaE');
        h5writeatt(filename,Dpth,'frame','General Frame');
        h5writeatt(filename,Dpth,'units','DeltaE');
    end
end
%write signal
fDszs = flip(Dszs);
signal = reshape(dat,fDszs); % change signal array dimensions to match the number of changing dimensions
h5createnwrite(filename,strcat(rtpth,'/signal'),signal,fDszs);
h5createnwrite(filename,strcat(rtpth,'/errors_squared'),zeros(fDszs),fDszs);
h5createnwrite(filename,strcat(rtpth,'/num_events'),zeros(fDszs),fDszs);
h5createnwrite(filename,strcat(rtpth,'/mask'),zeros(fDszs,'int8'),fDszs);
axesstr='';
for idx=1:length(Dstrcell)
    if idx<length(Dstrcell)
        axesstr=strcat(axesstr,Dstrcell{idx},':');
    else
        axesstr=strcat(axesstr,Dstrcell{idx});
    end
end
h5writeatt(filename,strcat(rtpth,'/signal'),'axes',axesstr);
h5writeatt(filename,strcat(rtpth,'/signal'),'signal',1);

%write experiment
exppth = NXScreategroup(filename,'/MDHistoWorkspace','experiment0','NXgroup');
h5writeatt(filename,exppth,'version', 1)
% write logs
log_pth = NXScreategroup(filename,exppth,'logs','NXgroup');
h5writeatt(filename,log_pth,'version',1)
writeNXlog(filename,log_pth,'W_MATRIX',proj_out,' ')
writeNXlog(filename,log_pth,'RUBW_MATRIX',proj_out,' ')
%write sample
smplpth = NXScreategroup(filename,exppth,'sample','NXsample');
h5writeatt(filename,smplpth,'version',1);
h5writeatt(filename,smplpth,'name',name);
h5writeatt(filename,smplpth,'shape_xml','<type name="userShape">  </type>');
h5createnwritevec(filename,smplpth,'num_oriented_lattice',int32(1))
h5createnwritevec(filename,smplpth,'num_other_samples',int32(0))
h5createnwritevec(filename,smplpth,'geom_height', 0)
h5createnwritevec(filename,smplpth,'geom_id', int32(0))
h5createnwritevec(filename,smplpth,'geom_thickness', 0)
h5createnwritevec(filename,smplpth,'geom_width' ,0)
%write material
mtl_pth = NXScreategroup(filename,smplpth,'material','NXdata');
h5writeatt(filename,mtl_pth,'formulaStyle','empty')
h5writeatt(filename,mtl_pth,'name',' ')
h5writeatt(filename,mtl_pth,'version',int32(2))
h5createnwritevec(filename,mtl_pth,'packing_fraction',1)
h5createnwritevec(filename,mtl_pth,'number_density',0)
h5createnwritevec(filename,mtl_pth,'pressure',0)
h5createnwritevec(filename,mtl_pth,'temperature',0)
%write crystal lattice
OL_pth = NXScreategroup(filename,smplpth,'oriented_lattice','NXcrystal');
%write lattice parameters
u_parm_names={'a','b','c','alpha','beta','gamma'};
for idx=1:length(u_parm_names)
    parm_name = strcat('unit_cell_',u_parm_names{idx});
    h5createnwritevec(filename,OL_pth,parm_name,unit_cell(idx))
end 
%write orientation matrix
om_path =strcat(OL_pth,'/orientation_matrix');
h5createnwrite(filename,om_path,Bmat,0);
%write instrument
instr_pth = NXScreategroup(filename,exppth,'instrument','NXinstrument');
h5writeatt(filename,instr_pth,'version',int32(1))
h5createnwritevec(filename,instr_pth,'name','SEQUOIA');
end

function h5createnwrite(filename,path,val,shp)
dtyp = class(val);
if shp ==0 
    shp = size(val)
end
h5create(filename,path,shp,'Datatype',dtyp);
h5write(filename,path,val);
end

function h5createnwritevec(filename,group,ds,val)
fh = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
sph = H5S.create_simple(1,length(val),length(val));
gsh = H5G.open(fh,group);
memtype = 'H5ML_DEFAULT';
switch class(val)
    case 'int32'
        H5type = 'H5T_NATIVE_INT';
        
    case 'int8'
        H5type = 'H5T_NATIVE_CHAR';
    case 'char'
        dims= size(val);
        sph = H5S.create_simple(1,fliplr(dims(1)),[]);
        H5type = H5T.copy('H5T_FORTRAN_S1');
        H5T.set_size(H5type,(dims(2)+1))
        memtype = H5T.copy ('H5T_C_S1');
        H5T.set_size(memtype,dims(2));
    otherwise
        H5type = 'H5T_NATIVE_DOUBLE';
end
dsh = H5D.create(gsh,ds,H5type,sph,'H5P_DEFAULT');
H5D.write(dsh,memtype,'H5S_ALL','H5S_ALL','H5P_DEFAULT',val)
H5S.close(sph)
H5D.close(dsh)
H5G.close(gsh)
H5F.close(fh)
end

function pthout = NXScreategroup(filename,pth,group,NX_class)
% ### Syntax

% pthout = NXScreategroup(filename,pth,group,NX_class)

% ### Description
%
% create a group with attributes of a nexus class
%
% ### Input Arguments
% filename, name of hdf5 file
% pth, path to where the group should be created
% group the group name
% NX_class a string containing a valid NX_class definition
%
% ### Output Arguments
% returns a path to the group

fh = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
pthh = H5G.open(fh,pth);
gh = H5G.create(pthh,group,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
pthout =strcat([pth,'/',group]);
H5G.close(gh)
H5G.close(pthh)
H5F.close(fh)
h5writeatt(filename,pthout,'NX_class',NX_class)
end

function writeNXlog(filename,log_pth,log_nm,value,units)
% ### Description
%
% create and write a Nexus log
%
% ### Input Arguments
%filename, name of hdf5 file
%log_pth, path to where the log should be created
%log_nm the name of the log
% value the vlaue of the log
% the units if any
pth = NXScreategroup(filename,log_pth,log_nm,'NXlog' );
h5createnwritevec(filename,pth,'value',value)
h5writeatt(filename,pth,'units',units)

end
function [latt_parms,Bmat,proj_out,D,signal,proj,name] = read_struct(dstruct,proj,dproj)

% ### Input Arguments
% dstruct: is the spinw structure
% proj: is the viewing projection matrix each column is a different axis.
%       One axis must be parallel to the drection of propogation in the
%       dstruct
% drproj: is the size of the bin in each direction (ignored for the direction
%         of the cut).
% 
% ### Output Arguments
% latt_parms: the lattice parameters from the spinw structure
% Bmat: a matrix for converting from inverse angstroms to rlu
% proj_out: 
% D: a cell array of the number of steps in each direction
% signal: the signal array from the spinw spec strcuture
% name : the chemical formula from the spinW file
    check_ortho(proj)
    fnames=fieldnames(dstruct);
    objnum = find(strcmp(fnames,'obj'));
    swobj = dstruct.(subsref(fnames,substruct('{}',{objnum})));
    latt_parms = abc(swobj);
    M = basisvector(swobj);
    Bmat = inv(M);
    name =formula(swobj).chemform;
    %proj_out = proj(:);
    hkls = dstruct.hkl;
    hkls_sz = size(hkls);
    % determine the direction where the hkl changes
    dir_vec = hkls(:,hkls_sz(2))-hkls(:,1);
    dir_vec = dir_vec/norm(dir_vec);
    %qout = hkls'/dir_vec';
    D={};
    % loop through each of the projection directions
    for qidx=1:3
        procjv = proj(:,qidx)/norm(proj(:,qidx));
       % if the projection direction is perpendicular to the propogation direction
       % then set the values to +/- dproj
       if abs(norm(cross(dir_vec,procjv)))> 1e-6
           dtmp = dot(hkls(:,2),procjv);
           D{qidx} = dtmp+dproj(qidx)/2.*[-1 1]; 
       else
           %assume it aslong the propogation direction
           hkl_proj = hkls'/proj(:,qidx)'; 
           dhkl = hkl_proj(2)-hkl_proj(1); % get the spacing along the q axis
           %hkl_proj = dot(hkls(:,1),procjv);
           D{qidx} = zeros([1,length(hkl_proj)+1]);
           D{qidx}(1:length(hkl_proj)) = hkl_proj-dhkl/2;
           D{qidx}(length(D{qidx})) = hkl_proj(length(hkl_proj))+dhkl/2;% change to bin boundaries
           %proj(:,qidx) = dir_vec; %set varying projection vector to spectra object 
       end  
    end
    D{4}=dstruct.Evect;
    signal = dstruct.swConv;
    proj_out = proj(:);
end

function check_ortho(mat)
% check if the three column vectors in proj are orthogonal to each other
sm = size(mat);
    for idx = 1:sm(2)
        for idx2 = 1:sm(2)
            if idx~=idx2
                if norm(dot(mat(:,idx),mat(:,idx2))) >1e-6
                    error("read_struct:nonorthogonal","the 3 vectors in proj must form an orthogonal basis set")
                end
            end
        end
    end
end