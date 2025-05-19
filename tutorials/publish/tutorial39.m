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
fitpow.estimate_scale_factor()

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
% default background is planar with parameters:
% [slope_energy, slope_modQ,  intercept]
fitpow.fix_bg_parameters(1:2); % fix slopes to initial value (0)
fitpow.set_bg_parameter_bounds(3, 0.0, []) % Set lb of constant bg = 0

[pfit, ~, stat] = fitpow.fit_background();
fitpow.estimate_scale_factor();

% view 1D cuts of data and the initial guess
fitpow.plot_1d_cuts_of_2d_data(qcens-dq, qcens+dq, fitpow.params);


%% Cropping data
% Data can be cropped in |Q| and energy transfer - for example it is
% advisible not to include the low-energy transfer region close to the
% elastic line (relative to the resolution). 

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
fitpow.optimizer = @ndbase.simplex; % or e.g. ndbase.simplex

% Fix the constant background at the result of the previous fit
fitpow.fix_bg_parameters(3); % fix fitted constant bg

% set bounds on first (and only in this case) model parameter
fitpow.set_model_parameter_bounds(1, 0, []) % Force J_1 to be AFM

[pfit, cost_val, stat] = fitpow.fit();

% plot 2D colorfill
fitpow.plot_result(pfit, 10,  'EdgeAlpha', 0.9, 'LineWidth', 1 ,'EdgeColor', 'k')

% make 1D plots
fitpow.plot_1d_cuts_of_2d_data(qcens-dq, qcens+dq, pfit)

%% Fit 1D cuts
% As discussed the sw_fitpowder class can be initialised with a cell array
% of xye structs correpsonding to 1D cuts at constant |Q|.
% Alternatively, having loaded in 2D data the sw_fitpowder class has a 
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
% Fix linear terms of background
% bg params ordered as en^n,...en^1, q^n,..q^1, const
fitpow.fix_bg_parameters(3);  % fix term ~|Q|
fitpow.fix_bg_parameters(1);  % fix term ~energy
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

% display result
par_names = ["J1", "Bg_En1", "Bg_Q2", "Bg_Q1", "Bg_Q0", "Scale"];
["Name", "Value", "Error"; par_names(:), pfit(:), stat.sigP(:)]

%     "Name"      "Value"           "Error"     
%     "J1"        "1.0027473"       "0.00360547"
%     "Bg_En1"    "0"               "0"         
%     "Bg_Q2"     "0.02066691"      "0.095079"  
%     "Bg_Q1"     "0"               "0"         
%     "Bg_Q0"     "-0.036883773"    "0.37194"   
%     "Scale"     "1082.117"        "19.405"    


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

% display result
bg_names = repmat(["Bg_En1_cut", "Bg_En0_cut"]', fitpow.ncuts, 1);
bg_names = bg_names + repelem([1:fitpow.ncuts]', fitpow.nparams_bg,1);
par_names = ["J1", bg_names', "Scale"];
["Name", "Value", "Error"; par_names(:), pfit(:), stat.sigP(:)]

%     "Name"           "Value"          "Error"     
%     "J1"             "1.0059739"      "0.00403765"
%     "Bg_En1_cut1"    "-0.85631989"    "0.376039"  
%     "Bg_En0_cut1"    "1.7546252"      "0.892913"  
%     "Bg_En1_cut2"    "0.52590868"     "0.433608"  
%     "Bg_En0_cut2"    "-2.3563309"     "1.49539"   
%     "Bg_En1_cut3"    "-0.41373845"    "0.350704"  
%     "Bg_En0_cut3"    "0.52287368"     "0.936986"  
%     "Scale"          "1131.9098"      "30.4495"   
