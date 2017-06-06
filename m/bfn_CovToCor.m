function cor = bfn_CovToCor(cov, varargin)
% BFN_COVTOCOR  Creating correlation matrix from covariance matrix.
%
% Syntax:
%   cor = bfn_CovToCor(cov)
%   cor = bfn_CovToCor(cov, zeroDiag)
%   cor = bfn_CovToCor(cov, <Property name 1>,<Property value 1>,...)
%
% Description:
%   It creates the correlation matrix from a given covariance matrix.
%
% Input Arguments:
%   cov - a covariance matrix
%
% Output Arguments:
%   cor - the output correlation matrix
%
% Examples:
%   cov = [1, 0.7; 0.3, 1];
%   cor = bfn_CovToCor(cov,1);
%
% See also modwt, modwt_wcov
%
% References:
%   Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time
%   Series Analysis. Cambridge: Cambridge University Press
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_CovToCor.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

params = struct('zeroDiag'      ,1); 
[params, args] = bfn_parseArgs(varargin,params);

% set range
if nargin < 2
    zeroDiag = 1;
else
    if numel(args) >= 1
        zeroDiag = args{1};
    else
        zeroDiag = params.zeroDiag;
    end
end

N = size(cov,1);

cor = zeros(N);
for i=1:N
    for j=1:N
        cor(i,j) = cov(i,j)/sqrt(abs(cov(i,i)*cov(j,j)));
    end
end

if zeroDiag
    cor = cor - eye(N);
end
