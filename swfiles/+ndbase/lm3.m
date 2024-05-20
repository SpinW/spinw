function [pOpt,yCalc,stat] = lm3(dat,func0,p0,varargin)
% least square refinement of parameters using Levenberg Marquardt method
%
% [pOpt,fVal,stat] = NDBASE.LM(dat,func,p0,'Option1','Value1',...)
%
% Levenberg Marquardt curve-fitting function, minimizes sum of weighted
% squared residuals.
%
% Input:
%
% dat       Data to be fitted stored in a structure with fields:
%               dat.x   vector of N independent variables,
%               dat.y   vector of N data values to be fitted,
%               dat.e   vector of N standard deviation (positive numbers)
%                       used to weight the fit. If zero or missing
%                       1./dat.y^2 will be assigned to each point, however
%                       in this case the refinement is reduced to an
%                       optimisation problem.
% func      Function handle with the following definition:
%               y = func(x,p)
%           where x is a vector of N independent variables, p are the
%           M parameters to be optimized and y is the simulated model.
% p0        Row vector of M initial parameters.
%
% Options:
%
% Options can be given using the modified output of optimset() or as option
% name string option value pairs.
%
% dp        Vector with N or 1 element, defines the fractional increment of
%           'p' when calculating the Jacobian matrix dFunc(x,p)/dp:
%               dp(j)>0     central differences calculated,
%               dp(j)<0     one sided 'backwards' differences calculated,
%               dp(j)=0     sets corresponding partials to zero, i.e. holds
%                           p(j) fixed.
%           Default value if 1e-3.
% vary      Vector with N elements, if an element is false, the
%           corresponding parameter will be fixed. Default value is
%           false(1,N).
% win       Limits for the independent variabel values where the function
%           is fitted. Default is [-inf inf].
% lb        Vector with N elements, lower boundary of the parameters.
%           Default value is -inf.
% ub        Vector with N elements, upper boundary of the parameters.
%           Default value is inf.
% MaxIter   Maximum number of iterations, default value is 10*M.
% MaxFunEvals Maximum number of function evaluations, default value is
%           100*M.
% TolX      Convergence tolerance for parameters, defines the maximum of
%           the relative chande of any parameter value. Default value is
%           1e-3.
% eps1      Convergence tolerance for gradient, default value is 1e-3.
% eps2      Convergence tolerance for reduced Chi-square, default value is
%           1e-2.
% eps3      Determines acceptance of a L-M step, default value is 0.1.
% lambda0   Initial value of L-M paramter, default value is 1e-2.
% nu0       Value that determines the speed of convergence. Default value
%           is 10. It should be larger than 1.
% lUp       Factor for increasing lambda, default value is 30.
% lDown     Factor for decreasing lambda, default value is 7.
% update    Type of parameter update:
%                   'lm'        Levenberg-Marquardt lambda update,
%                   'quadratic' Quadratic update,
%                   'nielsen'   Nielsen's lambda update equations (default).
% extraStat Calculates extra statistics: covariance matrix of parameters,
%           cofficient of multiple determination, asymptotic standard
%           error of the curve-fit and convergence history.
% confLev   Confidence level, where the error of the curve fit (sigY) is
%           calculated. Default is erf(1/sqrt(2))~0.6827 for standard
%           deviation (+/- 1 sigma).
%
% Output:
%
% pOpt      Value of the M optimal parameters.
% fVal      Value of the model function calculated with the optimal
%           parameters at the N independent values of x.
%
% stat      Structure, storing the detailed output of the calculation with
%           the following fields:
%               p       Least-squares optimal estimate of the parameter
%                       values.
%               redX2   Reduced Chi squared error criteria, its value
%                       should be close to 1. If the value is larger, the
%                       model is not a good description of the data. If the
%                       value is smaller, the model is overparameterized
%                       and fitting the statistical error of the data.
%               sigP    Asymptotic standard error of the parameters.
%               sigY    Asymptotic standard error of the curve-fit.
%               corrP   Correlation matrix of the parameters.
%               Rsq     R-squared cofficient of multiple determination.
%               cvgHst  Convergence history.
%               exitFlag The reason, why the code stopped:
%                           1       convergence in r.h.s. ("JtWdy"),
%                           2       convergence in parameters,
%                           3       convergence in reduced Chi-square,
%                           4       maximum Number of iterations reached
%                                   without convergence
%               message String, one of the above messages.
%               nIter   The number of iterations executed during the fit.
%               nFunEvals The number of function evaluations executed
%                       during the fit.
%
% See also NDBASE.PSO.
%

% T.G.Perring Jan 2009:
% ------------------------
% Generalise to arbitrary data objects which have a certain set of methods defined on them (see
% multifit.m for details)
%
% T.G.Perring 11-Jan-2007:
% ------------------------
% Core Levenberg-Marquardt minimisation method inspired by speclsqr.m from spec1d, but massively
% edited to make more memory efficient, remove redundant code, and especially rewrite the *AWFUL*
% estimation of errors on the parameter values (which needed a temporary
% array of size m^2, where m is the number of data values - 80GB RAM
% for m=100,000!). The error estimates were also a factor sqrt(ndat/(ndat-npfree))
% too large - as determined by comparing with analytical result for fitting to
% a straight line, from e.g. G.L.Squires, Practical Physics, (CUP ~1980). The
% current routine gives correct result.
%
% Previous history:
% -----------------
% Version 3.beta
% Levenberg-Marquardt nonlinear regression of f(x,p) to y(x)
% Richard I. Shrager (301)-496-1122
% Modified by A.Jutan (519)-679-2111
% Modified by Ray Muzic 14-Jul-1992

if nargin == 0
    swhelp ndbase.lm
    return
end

% number of parameters
Np    = numel(p0);
% vertical parameter vector
p0    = p0(:)';

inpForm.fname  = {'dp'   'lb'       'ub'      'MaxIter' 'eps1' 'TolX' 'MaxFunEvals' 'confLev'      'win'     };
inpForm.defval = {1e-3   -inf(1,Np) inf(1,Np) 10*Np     1e-3   1e-3   100*Np        erf(1/sqrt(2)) [-inf inf]};
inpForm.size   = {[1 -1] [1 Np]     [1 Np]    [1 1]     [1 1]  [1 1]  [1 1]         [1 1]          [1 2]     };

inpForm.fname  = [inpForm.fname  {'eps2' 'eps3' 'lambda0' 'lUp' 'lDown' 'update' 'extraStat' 'vary'      'nu0'}];
inpForm.defval = [inpForm.defval {1e-2    0.1    1e-2      30    7       'nielsen' true       true(1,Np) 10   }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]  [1 1]     [1 1] [1 1]   [1 -2]    [1 1]      [1 Np]      [1 1]}];

param = sw_readparam(inpForm, varargin{:});

% check input function
if ischar(func0)
    % convert to fuction handle
    func0 = str2func(func0);
end

if ~isa(func0,'function_handle')
    error('lm:WrongInput','The input function is neither a string, not function handle!');
end

if isempty(dat)
    % define fake data if is not given
    dat   = struct;
    dat.x = (1:Np)';
    dat.y = dat.x*0;
    dat.e = dat.x*0+1;
    
    func = @(x,p)func0(p)*ones(Np,1);
else
    % just use the given function
    func = func0;
    % make column vectors of all input
    dat.x = dat.x(:);
    dat.y = dat.y(:);
    
    %remove unnecessary data
    xKeep = dat.x>param.win(1) & dat.x<param.win(2);
    dat.x = dat.x(xKeep);
    dat.y = dat.y(xKeep);
    if isfield(dat,'e')
        dat.e = dat.e(:);
        dat.e = dat.e(xKeep);
    end
end

if numel(dat.x) ~= numel(dat.y)
    error('lm:WrongInput',['The number of independent vriables (x) is '...
        'not equal to the number of dependent variables (y)!']);
end

% define fit weights
if ~isfield(dat,'e') || isempty(dat.e) || ~any(dat.e(:))
    weight = 1./sqrt(abs(dat.y));
else
    if any(dat.e(:)<0)
        error('lm:WrongInput','Standard deviations have to be positive!')
    end
    weight = 1./dat.e(:);
end

if numel(param.dp) == 1
    param.dp = repmat(param.dp,[1 Np]);
end

if numel(param.vary) == 1
    param.vary = repmat(param.vary,[1 Np]);
end

% fixed the requested parameters
param.vary = logical(param.vary);
if any(~param.vary)
    param.dp(~param.vary) = 0;
end

% check that there are more data points than free parameters
Nval   = numel(dat.y);
NpFree = numel(p0);

% we allow for the case nval=npfree
nnorm  = max(Nval-NpFree,1);
if Nval<NpFree
    error('lm:WrongInput',['Number of data points must be greater than '...
        'or equal to the number of free parameters']);
end

% derivative step length
dp    = param.dp;
if dp>0
    nCallperJAC = Np;
else
    nCallperJAC = 2*Np;
end

% maximum number of iterations
Niter = param.MaxIter;
% convergence criterion
tol   = param.eps2;

if abs(dp)<1e-12
    error('lm:WrongParam','Derivative step length must be greater or equal to 10^-12')
end
if Niter<0
    error('lm:WrongParam','Number of iterations must be non-negative.')
end

% starting values of parameters and function values
yCalc = func(dat.x,p0);
nFunEval = 1;
resid = weight.*(dat.y-yCalc(:));

% un-normalised chi-squared
chi2Best = resid'*resid;

% best values for parameters at start
pBest = p0;
p     = pBest;

% function values at start
yBest = yCalc(:);
converged=false;
if Niter > 0
    % optimisation
    
    lambda = 1;
    lambda_table = [1 1 100 100 100 100];  %test this: [.1, 1, 1e2, 1e4, 1e6]
    converged = false;
    max_rescale_lambda = false;
    
    % iterate to find best solution
    for iter = 1:Niter
        % compute Jacobian matrix
        resid = weight.*(dat.y-yBest);
        jac   = dfdpf(dat,func,pBest,yCalc,dp);
        nFunEval = nFunEval + nCallperJAC;
        nrm   = zeros(NpFree,1);
        for j=1:NpFree
            jac(:,j) = weight.*jac(:,j);
            nrm(j)   = jac(:,j)'*jac(:,j);
            if nrm(j)>0
                nrm(j)=1/sqrt(nrm(j));
            end
            jac(:,j)=nrm(j)*jac(:,j);
        end
        
        [jac,s,v] = svd(jac,0);
        s = diag(s);
        g = jac'*resid;
        
        % Compute change in parameter values.
        % If the change does not improve chisqr  then increase the
        % Levenberg-Marquardt parameter until it does (up to a maximum
        % number of times gicven by the length of lambda_table).
        if tol>0
            chi2Goal = (1-tol)*chi2Best;  % Goal for improvement in chisqr
        else
            chi2Goal = chi2Best-abs(tol);
        end
        
        lambda = lambda/10;
        for itable = 1:numel(lambda_table)
            se    = sqrt((s.*s)+lambda);
            gse   = g./se;
            p_chg = ((v*gse).*nrm);   % compute change in parameter values
            if any(p_chg)  % there is a change in (at least one of) the parameters
                p     = pBest + p_chg(:)';
                yCalc = func(dat.x,p);
                nFunEval = nFunEval + 1;
                resid = weight.*(dat.y-yCalc(:));
                c     = resid'*resid;
                if c<chi2Best || c==0
                    pBest   = p;
                    yBest   = yCalc;
                    chi2Best = c;
                    break
                end
            end
            
            if itable == numel(lambda_table) % Gone to end of table without improving chisqr
                max_rescale_lambda = true;
                break
            end
            
            % Chisqr didn't improve - increase lambda and recompute step in parameters
            lambda = lambda*lambda_table(itable);
        end
        
        % If chisqr lowered, but not to goal, so converged; or chisqr==0 i.e. perfect fit; then exit loop
        if (chi2Best>chi2Goal) || (chi2Best==0)
            converged = true;
            break
        end
        
        % If multipled lambda to limit of the table, give up
        if max_rescale_lambda
            converged = false;
            break
        end
        
    end
end

% Wrap up for exit from fitting routine
if converged
    chisqr = chi2Best/nnorm;
    % Now get Jacobian matrix
    jac = dfdpf(dat,func,p,yCalc,dp);
    nFunEval = nFunEval + nCallperJAC;
    for j=1:NpFree
        jac(:,j)=weight.*jac(:,j);
    end
    
    %[jac,s,v]=svd(jac,0);
    [~,s,v]=svd(jac,0);
    s   = repmat((1./diag(s))',[NpFree,1]);
    v   = v.*s;
    cov = chisqr*(v*v');  % true covariance matrix;
    sigP = sqrt(diag(cov))';
    tmp = repmat(1./sqrt(diag(cov)),[1,NpFree]);
    cor = tmp.*cov.*tmp';
else
    sigP = [];
    iter=0;
    chisqr = chi2Best/nnorm;
    ok     = true;
    warning('WARNING: Convergence not achieved')
    cor = [];
end

% save output parameters
stat.p    = p;
stat.sigP = sigP;
% save counters
stat.nIter      = iter;
stat.nFunEvals  = nFunEval;
stat.algorithm  = 'Levenberg-Marquardt';
stat.func       = func;
stat.redX2      = chisqr;
% save the input parameters
stat.param = param;

% final parameter value
pOpt = p;

end

function jac = dfdpf(dat,func,p,f,dp)
% Calculate partial derivatives of function with respect to parameters
%
% jac = DFDPF(w,xye,func,bkdfunc,pin,bpin,p,p_info,f,dp)
%
% Input:
% ------
%   w       Cell array of data objects
%   xye     Logical array sye(i)==true if w{i} is x-y-e triple
%   func    Handle to global function
%   bkdfunc Cell array of handles to background functions
%   pin     Function arguments for global function
%   bpin    Cell array of function arguments for background functions
%   p       Parameter values of free parameters
%   p_info  Structure with information to convert free parameters to numerical
%           parameters needed for function evaluation
%   f       Function values at parameter values p sbove
%   dp      Fractional step change in p for calculation of partial derivatives
%                - if dp > 0    calculate as (f(p+h)-f(p))/h
%                - if dp < 0    calculate as (f(p+h)-f(p-h))/(2h)
%   listing Screen output control
%
% Output:
% -------
%   jac     Matrix of partial derivatives: m x n array where m=length(f) and
%           n = length(p)
%

% initialise Jacobian to zero
jac = zeros(numel(f),numel(p));

min_abs_del = 1e-12;

for j = 1:numel(p)
    % dp is fractional change in parameter
    del = dp(j)*p(j);
    % Ensure del non-zero
    if abs(del)<=min_abs_del
        if p(j)>=0
            del=min_abs_del;
        else
            del=-min_abs_del;
        end
    end
    if dp(j)>=0
        ppos=p;
        ppos(j)=p(j)+del;
        jac(:,j)=(func(dat.x,ppos)-f)/del;
    else
        ppos=p; ppos(j)=p(j)+del;
        pneg=p; pneg(j)=p(j)-del;
        jac(:,j)=(func(dat.x,ppos)-func(dat.x,pneg))/(2*del);
    end
end
end