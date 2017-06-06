function C = bfn_dwt_brick_wall(X, wavelet, N)
% BFN_DWT_BRICK_WALL  Sets the first N_j DWT coefficients to NaN.
%
% Syntax:
%   cor = bfn_CovToCor(cov)
%
% Description:
%   It creates the correlation matrix from a given covariance matrix.
%
% Input Arguments:
%   X - Matrix containing wavelet coefficients with appropriate 
%     boundary condition
%   wavelet - Character string; 'haar', 'd4', 'la8', 'la16'
%   N - Length of original vector of observations
%
% Output Arguments:
%   C - Matrix containing wavelet coefficients where ones affected
%       by boundary conditions are replaced with NaNs
%
% References:
%   Lindsay et al. (1996).  The Discrete Wavelet Transform
%   and the Scale Anlaysis of the Surface Properties of Sea
%   Ice.  IEEE Trans. on Geo. and Rem. Sen., 34(3), pp. 771-787.
%
%__________________________________________________________________________
% Original from WMTSA toolbox, Modified by Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_dwt_brick_wall.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

wtf = dwt_filter(wavelet);
L = wtf.L;
J = size(X,1);

C = X;
for j = 1:(J-1)
    n = ceil((L-2)*(1 - 1/2^j)); % for DWT
    n = min(n,length(X{j}));
    C{j}(1:n) = NaN;
end
C{j+1}(1:n) = NaN;
