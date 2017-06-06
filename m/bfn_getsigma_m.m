function Sigma = bfn_getsigma_m(d, C)
% BFN_GETSIGMA_M Getting the covariance matrix of a multivariate long memory process
%
% Syntax:
%   Sigma = bfn_getsigma_m(d, C)
%
% Description:
%   It computes the covariance matrix of fractal parameters of a
%   multivariate long memory process. The wavelet covariance of a long
%   memory process has a linearship  in log scale. The fractal parameter c
%   is identical to the y-axis  intercept in the plot of log-scale wavelet
%   variances with log wavelet scales.
%
% Input Arguments:
%   d - a vector of Hurst exponents
%   C - a matrix of fractal parameter c
%
% Output Arguments:
%   Sigma - the covariance matrix
%
% Examples:
%   c = bfn_getc(.8, 1);
%
% See also getsigma
%
% References:
%   Achard, S., Bassett, D. S., Meyer-Lindenberg, A., & Bullmore, E.
%   (2008). Fractal connectivity of long-memory networks. 
%   Physical Review E, 77(3), 1-12.
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_getsigma_m.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.
    
N = length(d);  
if N==1
    
    B1    = (1 - 2^(2*d(1)-1))/(1-2*d(1));
    Sigma = 2^(C(1,1)-1)/B1/(2*pi)^(1-2*d(1));
    
elseif N>1
    
    Sigma = zeros(N,N);
    for q1 = 1:N
        for q2 = q1:N
            d1 = d(q1); d2 = d(q2);
            B1 = (1 - 2^(d1+d2-1))/(1-d1-d2);
            Sigma(q1,q2) = 2^(C(q1,q2)-1)/B1/(2*pi)^(-d1-d2)/cos((d1-d2)*pi/2);
            Sigma(q2,q1) = Sigma(q1,q2);
        end
    end
    
end
