function [Z, NJ] = bfn_wave_sum_dwt(X, Y)
% BFN_WAVE_SUM_DWT  Compute the cross-wavelet product sum.
%
% Syntax:
%   [Z, NJ] = bfn_wave_sum_dwt(X, Y)
%
% Description:
%   It computes the cross-wavelet product sum.
%
% Input Arguments:
%   X,Y - Matrix containing wavelet coefficients with appropriate 
%         boundary condition
%
% Output Arguments:
%   Z - Matrix containing the cross-wavelet product sum.
%   NJ - Vector containing the number of coefficients in each scale.
%
% References:
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
%   Time Series Analysis. Cambridge: Cambridge University Press.
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_wave_sum_dwt.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

J  = length(X);

XY = {};
Z  = []; NJ = []; NX = zeros(1,J);
for j = 1:J
    NX(j) = sum(~isnan(X{j}));
    XY{j} = X{j} .* Y{j};
    
    XYj = XY{j};
    XYNaN = XYj(~isnan(XYj))';
    Z = [Z sum(XYNaN)];
    NJ = [NJ sum(~isnan(XYj))];
end
