function [CI_wvar, edof, Qeta, AJ] = modwt_wvar_ci(wvar, MJ, ci_method, ...
                                                   WJt, lbound, ubound, p)
% modwt_wvar_ci -- Calculate confidence interval of MODWT wavelet variance.
%
%****f* wmtsa.dwt/modwt_wvar_ci
%
% SYNOPSIS
%   modwt_wvar_ci -- Calculate confidence interval of MODWT wavelet variance.
%
% USAGE
%   [CI_wvar, edof, Qeta] = modwt_wvar_ci(wvar, MJ, [ci_method],
%                                         [WJt], [lbound], [ubound], [p])
%
% INPUTS
%   wvar         = wavelet variance (1xJ vector).
%   MJ          = number of coefficients used calculate the wavelet variance at
%                  each level (Jx1).
%   ci_method    = (optional) method for calculating confidence interval
%                  valid values:  'gaussian', 'chi2eta1', 'chi2eta3'
%                  default: 'chi2eta3'
%   WJt          = MODWT wavelet coefficients (NxJ array).
%                  where N = number of time intervals
%                        J = number of levels
%                  required for 'gaussian' and 'chi2eta1' methods.
%   lbound       = lower bound of range of WJt for calculating ACVS for each
%                  level (Jx1 vector).
%   ubound       = upper bound of range of WJt for calculating ACVS for each
%                  level (Jx1 vector).
%   p            = (optional) percentage point for chi2square distribution.
%                  default: 0.025 ==> 95% confidence interval
%
% OUTPUTS
%   CI_wvar      = confidence interval of wavelet variance  (Jx2 array).
%                  lower bound (column 1) and upper bound (column 2).
%   edof         = equivalent degrees of freedom (Jx1 vector).
%   Qeta         = p x 100% percentage point of chi-square(eta) distribution (Jx2 array).
%                  lower bound (column 1) and upper bound (column 2).
%   AJ           = integral of squared SDF for WJt (Jx1 vector).
%
% SIDE EFFECTS
%
%
% DESCRIPTION
%   MJ is vector containing the number of coefficients used to calculate the 
%   wavelet variance at each level. 
%   For the unbiased estimator, MJ = MJ for j=1:J0, where MJ is the number 
%   of nonboundary MODWT coefficients at each level.
%   For the biased estimator, MJ = N for all levels.
%   For the weaklybiased estimator, MJ = MJ(Haar), for j=1:J0, where MJ(Haar) 
%   is the number of nonboundary MODWT coefficients for Haar filter at each level.
%
% EXAMPLE
%
%
% ERRORS  
%   WMTSA:InvalidNumArguments
%   WMTSA:WVAR:InvalidCIMethod
%
% NOTES
%   The output argument edof (equivalent degrees of freedom) is returned for
%   the chi2 confidence interval methods.  For the gaussian method, a null
%   value is returned for edof.
%  
%
% ALGORITHM
%   See section 8.4 of WMTSA on pages 311-315.
%
% REFERENCES
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%     Time Series Analysis. Cambridge: Cambridge University Press.
%   
% SEE ALSO
%   wmtsa_acvs
%
% TOOLBOX
%   wmtsa/wmtsa
%
% CATEGORY
%   ANOVA:WVAR:MODWT
%


% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2003-04-23
%
% CREDITS
% 
%
% COPYRIGHT
%
%
% CREDITS
%
%
% REVISION
%   $Revision: 612 $
%
%***

% $Id: modwt_wvar_ci.m 612 2005-10-28 21:42:24Z ccornish $

valid_ci_methods = {'gaussian', 'chi2eta1', 'chi2eta3', 'none'};

default_ci_method = 'chi2eta3';
default_p = 0.025;
  
usage_str = ['[CI_wvar, edof, Qeta, AJ] = ', mfilename, ...
               '(wvar, MJ, [ci_method], [WJt], [lbound], [ubound], [p]'];
  
%%  Check input arguments and set defaults.
error(nargerr(mfilename, nargin, [2:7], nargout, [0:4], 1, usage_str, 'struct'));


if (~exist('ci_method', 'var') || isempty(ci_method))
  ci_method = default_ci_method;
else
  if (isempty(strmatch( ci_method, valid_ci_methods, 'exact')))
    error('WMTSA:WVAR:InvalidCIMethod', ...
          [ci_method, ' is not a valid confidence interval method.']);
  end
end


if (strmatch(ci_method, {'gaussian', 'chi2eta1'}, 'exact'))
  if (~exist('WJt', 'var') || isempty(WJt) || ...
      ~exist('lbound', 'var') || isempty(lbound) || ...
      ~exist('ubound', 'var') || isempty(ubound))
    error('WMTSA:MissingRequiredArgument', ...
          wmtsa_encode_errmsg('WMTSA:MissingRequiredArgument', ...
                        'WJt, lbound, ubound arguments required', ...
                       [' when using ci method (', ci_method,').']));
  end
end


if (~exist('p', 'var') || isempty(p))
  p = default_p;
end

sz = size(wvar);
J = sz(1);
nsets = sz(2);

% Initialize output
CI_wvar = [];
edof = [];
AJ = [];
Qeta = [];

j_range = (1:J)';

% Calculate ACVS and AJ used by gaussian and chi2eta1 ci methods.
switch lower(ci_method)
 case {'gaussian', 'chi2eta1'}

  % For 'gaussian', 'chi2eta1' ci methods, 
  % compute AJ (integral of squared SDF for WJt) 
  % via calculating ACVS of WJt. (see page 312 of WMTSA).

%  AJ = zeros(J,1,nsets);

  for (i = 1:nsets)
    for j = 1:J
      if (MJ(j) > 0)
        ACVS = wmtsa_acvs(WJt(lbound(j):ubound(j),j,i), 'biased', 0, '', 1);
        % Get lags for taus >= 0.
        ACVS = ACVS(MJ(j):2*MJ(j)-1);
        AJ(j,i) = ACVS(1).^2 / 2 + sum(ACVS(2:MJ(j)).^2);
      else
        AJ(j,i) = NaN;
      end
    end
  end
end  % end switch ci_method



% Calculate confidence intervals

switch lower(ci_method)
 case 'gaussian'
  % Gaussian method, Eq. 311 of WMSTA
  VARwvar = 2 .* AJ ./ MJ;
  CI_wvar(:,1) = wvar - norminv(1-p) * sqrt(VARwvar);
  CI_wvar(:,2) = wvar + norminv(1-p) * sqrt(VARwvar);
 
 case 'chi2eta1'
  % Chi-square equivalent degree of freedom method 1
  % Eqns. 313c & 313d of WMSTA
  eta1 = (repmat(MJ, [1,size(wvar,2)]) .* wvar.^2);
  eta1 = eta1 ./ AJ;
  [Qeta, CI_wvar] = calc_ci_via_edof(wvar, eta1, p);
  edof = eta1;

 case 'chi2eta2'
  error('WMTSA:WVAR:InvalidCIMethod', [ci_method, ' is not currently implemented']);
 
 case 'chi2eta3'
  % Chi-square wquivalent degree of freedom method 3
  % Eqns. 314c & 313c of WMSTA
  
  eta3 = max(MJ ./ (2.^j_range), 1);

  [Qeta, CI_wvar] = calc_ci_via_edof(wvar, eta3, p);
  edof = eta3;
 
 otherwise
  error('WMTSA:WVAR:InvalidCIMethod', ['Unknown confidence interval method', ci_method]);
end

return


function [Qeta, CI_wvar] = calc_ci_via_edof(wvar, eta, p)
  Qeta(:,1,:) = chi2inv(1-p, eta);
  Qeta(:,2,:) = chi2inv(p, eta);
  CI_wvar(:,1,:) = eta .* wvar ./ squeeze(Qeta(:,1,:));
  CI_wvar(:,2,:) = eta .* wvar ./ squeeze(Qeta(:,2,:));
return
  
