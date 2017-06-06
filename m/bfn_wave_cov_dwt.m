function [C, L, U] = bfn_wave_cov_dwt(X, Y)
% BFN_WAVE_COV_DWT  Compute wavelet covariance with approximate 95% confidence interval
%
% Syntax:
%   [C, L, U] = bfn_wave_cov_dwt(X, Y)
%
% Description:
%   It computes wavelet covariance with approximate 95% confidence
%   interval.
%
% Input Arguments:
%   X,Y - Matrix containing wavelet coefficients with appropriate 
%         boundary condition
%
% Output Arguments:
%   C - Matrix containing the wavelet covariance (column 1), lower 
%       95% quantile for confidence interval, upper 95% quantile 
%       for confidence interval
%
% References:
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%   Time Series Analysis. Cambridge: Cambridge University Press.
%
%__________________________________________________________________________
% Original from WMTSA toolbox, Modified by Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_wave_cov_dwt.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

J = length(X);

XY = {};
Z = []; NX = zeros(1,J);
for j = 1:J
    NX(j) = sum(~isnan(X{j}));
    XY{j} = X{j} .* Y{j};
    
    XYj = XY{j};
    XYNaN = XYj(~isnan(XYj))';
    Z = [Z sum(XYNaN)/sum(~isnan(XYj))];   % unbiased
    
    %NJ = length(XYj);
    %Z = [Z sum(XYNaN)/NJ];   % biased
end

VARgamma = [];
for j = 1:J
    Xj = X{j};
    XjNaN = Xj(~isnan(Xj));
    XACFj = wmtsa_acvs(XjNaN);
    XACFNaN = XACFj(~isnan(XACFj))';

    Yj = Y{j};
    YjNaN = Yj(~isnan(Yj));
    YACFj = wmtsa_acvs(YjNaN);
    YACFNaN = YACFj(~isnan(YACFj))';

    SUMXYCCFj = wmtsa_ccvs(XjNaN,YjNaN);
    
    CCFj = wmtsa_ccvs(XjNaN,YjNaN); 
    CCFNaN = CCFj(~isnan(CCFj))';
    if numel(XACFNaN) == numel(YACFNaN)
        A = sum(XACFNaN .* YACFNaN) + CCFNaN(1).^2 ./ 2;
        A = A/(2 * NX(j));
    else
        A = NaN;
    end
    
    VARgamma = [VARgamma A];
end

C = Z;
L = Z - norminv(0.975) .* sqrt(VARgamma);
U = Z + norminv(0.975) .* sqrt(VARgamma);
