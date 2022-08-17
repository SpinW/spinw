%% Antiferromagnetic square lattice
% We define a square lattice in the ab plane, with Cu+ ions with S=1 spin.

AFsq = spinw;
AFsq.genlattice('lat_const',[3 3 6],'angled',[90 90 90],'sym',0)
AFsq.addatom('r',[0 0 0],'S', 1,'label','Cu1','color','b')
AFsq.table('atom')
plot(AFsq)
swplot.zoom(1.5)

%% Couplings
% We create first neighbor couplings in the ab plane and plot the bonds.
% You can click on the different bonds to get the value of the
% corresponding matrix.

AFsq.gencoupling('maxDistance',5)
AFsq.table('bond',[])

AFsq.addmatrix('label','J1','value',   1,'color','red')
AFsq.addmatrix('label','J2','value',-0.1,'color','green')
AFsq.addcoupling('mat','J1','bond',1)
AFsq.addcoupling('mat','J2','bond',2)
plot(AFsq,'range',[2 2 0.5])

%% Magnetic structure
% For weak second neighbor ferromagnetic interaction the magnetic structure
% is Neel type, with the following parameters:
%
% * ordering wave vector k = (1/2 1/2 0)
% * spin are in arbitrary plane, lets point along the S = (1 0 0) direction
% * normal to the spin vectors n = (0 0 1)
% * magnetic supercell is 2x2x1
%
% We use magnetic supercell, since the 2*k equal to a reciprocal lattice
% vector. In this case the spin wave code cannot use the incommensurate
% mode, thus we have to create a zero-k structure, that is a 2x2x1 magnetic
% supercell.Note that the spinw.genmagstr() automatically normalizes the spin
% vectors to the spin value of the magnetic atoms.

AFsq.genmagstr('mode','helical','k',[1/2 1/2 0],'n',[0 0 1], 'S',[1; 0; 0],'nExt',[1 1 1]);  
disp('Magnetic structure:')
AFsq.table('mag')

AFsq.energy
plot(AFsq,'range',[2 2 1])

%% Spin wave spectrum
% We calculate the spin wave spectrum and correlation function along
% several linear q-scans in reciprocal space defined by the Qcorner corner
% points. The last number in the cell defines the number of steps in each
% linear scan.

Qcorner = {[0 0 0] [1/2 0 0] [1/2 1/2 0] [0 0 0] 200};
sqSpec = AFsq.spinwave(Qcorner, 'hermit', false);
sqSpec = sw_neutron(sqSpec); 
sqSpec = sw_egrid(sqSpec,'Evect',linspace(0,6.5,500));
figure
sw_plotspec(sqSpec,'mode',3,'dashed',true,'dE',0.4,'qLabel',{'\Gamma' 'X' 'M' '\Gamma'})
caxis([0 4])

%%
%  Written by
%  Bjorn Fak & Sandor Toth
%  06-Jun-2014, 06-Feb-2017
