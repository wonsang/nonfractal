function [C, L, U] = bfn_wave_var_dwt(X, I)
% BFN_WAVE_VAR_DWT  Compute wavelet variance with approximate 95% confidence interval.
%
% Syntax:
%   [C, L, U] = bfn_wave_var_dwt(X, I)
%
% Description:
%   It computes the cross-wavelet product sum.
%
% Input Arguments:
%   X - Matrix containing wavelet coefficients with appropriate 
%       boundary condition
%   L - lower boundary
%   U - upper boundary
%
% Output Arguments:
%   C - Matrix containing the wavelet variance (column 1), lower 
%       95% quantile for confidence interval, upper 95% quantile 
%       for confidence interval
%
% References:
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%   Time Series Analysis. Cambridge: Cambridge University Press.
%
%__________________________________________________________________________
% Original from WMTSA toolbox, Modified by Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_wave_var_dwt.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

J = length(X);
Y = []; NX = zeros(1,J);
for j = 1:J    
    Xj = X{j};
    NX(j) = sum(~isnan(Xj));
    XNaN = Xj(~isnan(Xj))';
    Y = [Y sum(XNaN.^2)/sum(~isnan(Xj))];   % unbiased
    
    %NJ = length(Xj);
    %Y = [Y sum(XNaN.^2)/NJ];   % biased
end

%%% Computes the variance of the wavelet variance via Chapter 8.
%{
VARnu = []; 
for j = 1:J    
    Xj = X{j};
    XjNaN = Xj(~isnan(Xj));
    ACFj = wmtsa_acvs(XjNaN);
    ACFNaN = ACFj(~isnan(ACFj))';
    if numel(ACFNaN)
        A = sum(ACFNaN.^2) - ACFNaN(1).^2 ./ 2;
        A = 2 * A / I;
    else
        A = NaN;
    end
    VARnu = [VARnu A];
end
%}
%C = [Y; Y - norminv(0.975) * sqrt(VARnu); Y + norminv(0.975) * sqrt(VARnu)]';

%%% Computes eta_3 from Chapter 8.

format short e
ETA3 = max(NX ./ 2.^(1:J), 1);

C = Y;
L = ETA3 .* Y ./ chi2inv(0.975, ETA3);
U = ETA3 .* Y ./ chi2inv(0.025, ETA3);

