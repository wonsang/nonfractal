function [C, L, U] = bfn_wave_cor_dwt(X, Y, I)
% bfn_WAVE_COR_DWT  Compute wavelet correlation with approximate 95% confidence interval
%
% Syntax:
%   [C, L, U] = bfn_wave_cor_dwt(X, Y, I)
%
% Description:
%   It computes wavelet correlation with approximate 95% confidence
%   interval.
%
% Input Arguments:
%   X,Y - Matrix containing wavelet coefficients with appropriate 
%         boundary condition
%
% Output Arguments:
%   C - Matrix containing the wavelet correlation (column 1), lower 
%       95% quantile for confidence interval, upper 95% quantile 
%       for confidence interval
%
% References:
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%   Time Series Analysis. Cambridge: Cambridge University Press.
%
%__________________________________________________________________________
% Original from WMTSA toolbox, Modified by Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_wave_cor_dwt.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

J = length(X);

XY = {};
SSX = []; SSY = []; SSXY = [];
for j = 1:J
    XY{j} = X{j} .* Y{j};
    
    Xj = X{j};
    XNaN = Xj(~isnan(Xj))';
    SSX = [SSX sum(XNaN.^2)/sum(~isnan(Xj))];
    
    %NJ = length(Xj);
    %SSX = [SSX sum(XNaN.^2)/NJ];   % biased
    
    Yj = Y{j};
    YNaN = Yj(~isnan(Yj))';
    SSY = [SSY sum(YNaN.^2)/sum(~isnan(Yj))];
    
    %NJ = length(Yj);
    %SSY = [SSY sum(YNaN.^2)/NJ];   % biased
    
    XYj = XY{j};
    XYNaN = XYj(~isnan(XYj))';
    SSXY = [SSXY sum(XYNaN)/sum(~isnan(XYj))];
    
    %NJ = length(XYj);
    %SSXY = [SSXY sum(XYNaN)/NJ];   % biased
end
COR = (SSXY) ./ (sqrt(SSX) .* sqrt(SSY));

NDWT = floor(I ./ 2.^(1:J));

C = COR;
L = tanh(atanh(COR) - norminv(0.975) ./ sqrt(NDWT-3));
U = tanh(atanh(COR) + norminv(0.975) ./ sqrt(NDWT-3));

