%% Description
% This tutorial fits a spinwave model to simulated powder data.

%% Define instrument resolution parameters for MARI
% This is a simplified parameterisation of the resolution of a
% Time-of-Flight spectrometer which will be used to broaden the simulated
% data. These parameters will depend on the instrument used - the instrument
% contact should be able to provide these.

Ei = 20; % meV

% Resolution for MARI with Gd chopper running at 200Hz and Ei=20meV (from PyChop)
eres = @(en) 2.1750e-04*sqrt((Ei-en).^3 .* ( (32.27217*(0.168+0.400*(Ei./(Ei-en)).^1.5)).^2 + (14.53577*(1.168+0.400*(Ei./(Ei-en)).^1.5)).^2) );
% Q-resolution (parameters for MARI)
e2k = @(en) sqrt( en .* (2*1.67492728e-27*1.60217653e-22) )./1.05457168e-34./1e10;
L1 = 11.8;  % Moderator to Sample distance in m
L2 = 4.0;   % Sample to detector distance in m
ws = 0.05;  % Width of sample in m
wm = 0.12;  % Width of moderator in m
wd = 0.025; % Width of detector in m
ki = e2k(Ei);
a1 = ws/L1; % Angular width of sample seen from moderator
a2 = wm/L1; % Angular width of moderator seen from sample
a3 = wd/L2; % Angular width of detector seen from sample
a4 = ws/L2; % Angular width of sample seen from detector
dQ = 2.35 * sqrt( (ki*a1)^2/12 + (ki*a2)^2/12 + (ki*a3)^2/12 + (ki*a4)^2/12 );

%% Generate data
% Simulate a powder spectrum of a triangular AFM with constant background
% and normally distributed noise. This can also be done using the powder
% fitting class if an empty data array is passed.


% intialise spinw model
tri = sw_model('triAF',1);
scale = 1000;

% calculate powder spectrum
E = linspace(0,4,80);
Q = linspace(0.1,3,60);
spec = tri.powspec(Q,'Evect',E,'nRand',2e2);
% scale data and add background
spec.swConv= scale*(spec.swConv + 0.1*mean(spec.swConv, "all", 'omitnan'));
% crop to kinematic range of instrument and apply some arbitrary resolution
spec = sw_instrument(spec, 'dQ', dQ, 'dE', eres, 'ThetaMin', 3.5, 'Ei', Ei);
% add noise
e = sqrt(spec.swConv);
spec.swConv = spec.swConv + e.*randn(size(spec.swConv));

% make xye struct containing data
data = struct('x', {{0.5*(E(1:end-1) + E(2:end)), Q}}, ...
              'y', spec.swConv', ...
              'e', e');

%% Initialise the sw_fitpowder class and plot data to inspect initial guess
% To initialise the powder fitting class (sw_fitpowder) object requires a 
% spinw model with a magnetic strucutre and couplings added. 
% Additional arguments are the data (can be xye structs for 2D data a cell
% array of xye structs for 1D data, or HORACE objects - see sw_fitpowder
% help for more details), a fit function which alters the matrices on the 
% spinw object and a vector of initial guesses for the fit parameters.
%
% Optinally a background strategy can be specified as a 5th positional 
% argument ("planar" background in energy transfer and |Q| or "independent"
% background for each 1D cut (1D only)).
%
% The class will call other spinw algorithms such as powspec and
% sw_instrument. The parameters used in these calls can be set by the user
% by modifying structs as shown below.

fit_func =  @(obj, p) matparser(obj, 'param', p, 'mat', {'J_1'}, 'init', true);
J1 = 0.85;

fitpow = sw_fitpowder(tri, data,  fit_func, [J1]);
fitpow.set_caching(true); % will store last calculated spectra

% set parameters passed to powspec
fitpow.powspec_args.dE = eres; % emergy resolution
fitpow.powspec_args.fastmode = true;
fitpow.powspec_args.neutron_output = true;
fitpow.powspec_args.nRand = 1e3; % low for speed (typically want > 1e3)
fitpow.powspec_args.hermit = true;
% set parameters passsed to sw_instrument
fitpow.sw_instrument_args = struct('dQ', dQ, 'ThetaMin', 3.5, 'Ei', Ei);


%% Estimate the background automatically
% Before fitting the data it is necessary to first get a good estimate for
% the background and the scale factor.
% There is a means or trying to automatically determine the background
% region - it is worth trying this first

fitpow.estimate_constant_background();
fitpow.estimate_scale_factor();

% view the background region on 2D colorfill plot of data
figure("color","white")
fitpow.plot_2d_data(fitpow.y);
fitpow.plot_background_region();

% view 1D cuts of data and the initial guess
qcens = 0.75:0.75:2.25;
dq = 0.1;
fitpow.plot_1d_cuts_of_2d_data(qcens-dq, qcens+dq, fitpow.params);

%% User supplied background region
% If the automatic background estimation doesn't look reasonable then it is
% possible to define regions of |Q| and energy transfer that are background
% and use these for fitting.
%
% General polynomial backgrounds in both |Q| and energy transfer are
% supported (governed by the attributes npoly_modQ and npoly_en).
% The default background is planar of the form
%
% Q1*|Q| + E1*en + E0
%
% The background parameters can be referenced by index or label 
% (e.g. "E0" for constant background).
% The labels are of the form "QN" or "EN" (for |Q| and energy
% polynomial coeficients repectively) where N represents the order of the 
% polynomial. The background parameters are indexed in the following order
%
% [QN, QN-1, ..., Q1, EM, EM-1, ..., E1, E0]
%
% where N = npoly_modQ and M=npoly_en. Note there is no Q0 parameter, this
% would correspond to a constant which would be perfectly correlated with 
% E0.


% set background regions
%-----------------------

fitpow.clear_background_region(); % in case previously set

% set background region with en > 3.75 meV
fitpow.set_bg_region(3.75, inf);
% set background region with
% 0.65 < |Q| < 0.9 Ang^-1
% 0.8 < en < 1.2 meV
fitpow.set_bg_region(0.8, 1.2, 0.65, 0.9);


% fit background (to voxels in background region)
%-----------------------------------------------------

% fit only a constant/flat background
fitpow.fix_bg_parameters(["Q1", "E1"]); % fix slopes to initial value (0)
fitpow.set_bg_parameter_bounds("E0", 0.0, []) % Set lb of constant bg = 0

% general polynomial backgrounds are supported the parameters are ordered

[pfit, ~, stat] = fitpow.fit_background();
fitpow.estimate_scale_factor();

% view 1D cuts of data and the initial guess
fitpow.plot_1d_cuts_of_2d_data(qcens-dq, qcens+dq, fitpow.params);

% print inital parameters
fitpow.disp_params();
% ---
% MODEL
% ---
% Index		Initial		
%        1	    0.85	
%    SCALE	 1316.58	
% ---
% BACKGROUND
% ---
% Label	Index	Initial	
%    Q1	       1	       0	
%    E1	       2	       0	
%    E0	       3	0.976825	

%% Cropping data
% Data can be cropped in |Q| and energy transfer - for example it is
% advisible not to include the low-energy transfer region close to the
% elastic line (relative to the resolution). 
%
% You may want to re-optimise the background and scale factor after
% cropping

% crop data
fitpow.crop_energy_range(1.5*eres(0), inf); % typically want to avoid elastic line
fitpow.crop_q_range(0.25, 3);

% plot the cropped data and the inital simulated spectrum
fig = figure("color","white"); hold all;
fig.Position(3) = 1.5*fig.Position(3);
subplot(1,2,1)
ax_dat = fitpow.plot_2d_data(fitpow.y);
subplot(1,2,2)
[ycalc, ~] = fitpow.calc_spinwave_spec(fitpow.params);
ax_calc = fitpow.plot_2d_data(ycalc);
ax_calc.CLim = ax_dat.CLim;
for ax = [ax_dat, ax_calc]
    ax.YLim(1) = 1.5*eres(0);
end
sgtitle(fig, "Initial guess")


%% Fit and inspect the result
% It can be desirable to set bounds on fit parameters and fix some
% parameters initially to get a reasonable fit with as few degrees of
% freedom as possible. For example one might assume the background is 
% constant in an initial fit.
%
% The cost function minimised in the fit can be chi-squared ("chisq") or 
% the unweighted sum of the squared residuals ("Rsq"). Also the minimiser
% can be set to any callable as long as it shares the API of the ndbase
% minimisers. Typically users would be recommended to use ndbase.lm4 or
% ndbase.simplex.

% set cost function and minimser
fitpow.cost_function = "Rsq";  % "Rsq" or "chisq"
fitpow.optimizer = @ndbase.simplex; % or e.g. ndbase.lm4

% Fix the constant background at the result of the previous fit
fitpow.fix_bg_parameters("E0"); % fix fitted constant bg

% set bounds on first (and only in this case) model parameter
fitpow.set_model_parameter_bounds(1, 0, []) % Force J_1 to be AFM

[pfit, cost_val, stat] = fitpow.fit();

% plot 2D colorfill
fitpow.plot_result(pfit, 10,  'EdgeAlpha', 0.9, 'LineWidth', 1 ,'EdgeColor', 'k')

% make 1D plots
fitpow.plot_1d_cuts_of_2d_data(qcens-dq, qcens+dq, pfit)

% print parameters
fitpow.disp_params(pfit);
% ---
% MODEL
% ---
% Index		Initial		Best Fit		
%        1	    0.85	 1.00185	
%    SCALE	 1316.58	 1024.67	
% ---
% BACKGROUND
% ---
% Label	Index	Initial	Best Fit	
%    Q1	       1	       0	       0	
%    E1	       2	       0	       0	
%    E0	       3	0.976825	0.976825	

%% Fit 1D cuts
% As discussed the sw_fitpowder class can be initialised with a cell array
% of xye structs correpsonding to 1D cuts at constant |Q|.
% Alternatively, having loaded in 2D data, the sw_fitpowder class has a 
% convenient method to take 1D cuts and replace the data. 
% 
% If 2D data are replaced with 1D cuts then the class evaluates the 
% cost function using the average of the model over the |Q| bins in the 
% 2D data within the integration limits. 
% If the class is initialised with 1D cuts then the class averages the model 
% over nQ points. 

% set inital parameters to those of 2D fit
fitpow.params = pfit(:);

% set data to 1D cuts
fitpow.replace_2D_data_with_1D_cuts(qcens-dq, qcens+dq);

% change background to quadratic in |Q|
fitpow.set_bg_npoly_modQ(2);  % zeros previously optimised background
% Fix linear terms of background to be 0
fitpow.fix_bg_parameters(["Q1", "E1"]);
% fit background and scale (fix model params) to all bins (not just bg
% region -0 only advisable if near good solution)
[pfit, ~, stat] = fitpow.fit_background_and_scale();
% fitpow.plot_result(fitpow.params)

% fit data with lm4 minimiser (which calculates parameter uncertainties)
% this requires a larger nRand as derivatives in LM algorithm sensitive to 
% noise introduced in powder averaging.
% if uncertainties come out imaginary try increasing diff_step
fitpow.optimizer = @ndbase.lm4;
fitpow.powspec_args.nRand = 1e4;
fit_kwargs = {'resid_handle', true, 'diff_step', 1e-3};

[pfit, cost_val, stat] = fitpow.fit(fit_kwargs{:});

fitpow.plot_result(pfit);

% print parameters
fitpow.disp_params(pfit, stat.sigP);
% ---
% MODEL
% ---
% Index		Initial		Best Fit		Error		
%        1	 1.00185	 1.00091	0.00357452	
%    SCALE	 1106.42	 1105.33	 19.8056	
% ---
% BACKGROUND
% ---
% Label	Index	Initial	Best Fit	Error	
%    Q2	       1	-0.618568	-0.613758	0.294827	
%    Q1	       2	       0	       0	       0	
%    E1	       3	       0	       0	       0	
%    E0	       4	0.444711	0.455551	0.497835	

%% Change background strategy
% The background for 1D cuts can be handled in two ways:
%
%  - A planar background as in the 2D case
%
%  - A linear background in energy transfer that's independent for each cut.
%
% The default is planar - but this can be changed after the class has been
% initialised. 
% The class will attempt to get reasonable guesses for the
% new background parameters, but note the background parameter bounds will
% be reset.

fitpow.params = pfit(:);
fitpow.set_background_strategy("independent");

[pfit, cost_val, stat] = fitpow.fit(fit_kwargs{:});

fitpow.plot_result(pfit)

% print parameters
fitpow.disp_params(pfit, stat.sigP);
% ---
% MODEL
% ---
% Index		Initial		Best Fit		Error		
%        1	 1.00091	 1.00181	0.00361181	
%    SCALE	 1105.33	 1213.84	 29.0031	
% ---
% BACKGROUND
% ---
% Label	Index	Cut Index	Initial	Best Fit	Error	
%    E1	1	       1	       0	-1.03769	0.359709	
%    E0	2	       1	0.0977161	 1.79735	0.857853	
%    E1	1	       2	       0	 1.85999	 0.41829	
%    E0	2	       2	-0.926966	-7.27229	 1.43701	
%    E1	1	       3	       0	 0.23947	0.337007	
%    E0	2	       3	-2.61892	-2.66401	0.904453	
