function sdf = bfn_finsdf_b(d1, d2, range, SNR)
% BFN_FINSDF_B Computation of scale-dependent variances of a bivariate fractionally integrated noise.
%
% Syntax:
%   sdf = bfn_finsdf_b(d1, d2, range, SNR)
%
% Description:
%   It computes the scale-dependent covariances of a bivariate
%   fractionally integrated noise (FIN) given its long memory
%   parameter. Also, it supports computation for a bFIN process
%   corrupted by white noise. 
%
% Input Arguments:
%   d1,d2 : a scalar of Hurst exponent
%   range : the range of wavelet scales
%   SNR : the signal-to-noise ratio
%
% Output Arguments:
%   sdf : a vector of scale-dependent covariances
%
% References:
%   Achard, S., Bassett, D. S., Meyer-Lindenberg, A., & Bullmore, E.
%    (2008). Fractal connectivity of long-memory networks. 
%    Physical Review E, 77(3), 1-12.
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_finsdf_b.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

dsum = d1+d2;    
ddif = d1-d2;
Js   = (range(1):range(2))';

B1   = (1 - 2^(dsum-1))/(1-dsum);
B2   = (1 - 2^(dsum-2))/(2-dsum);
B3   = (1 - 2^(dsum-3))/(3-dsum);
A    = -ddif^2/8 + dsum/24;  
a0   = cos(ddif*pi/2);
a1   = 2*pi*sin(ddif*pi/2)*ddif/2*B2/B1;
a2   = (2*pi)^2*A*B3/B1*cos(ddif*pi/2);
K    = 2*B1*(2*pi)^(1-dsum);

sdf  = K*2.^(Js*dsum)*a0;

if ~isnan(SNR)
    varRatio = 10^(-SNR/10);
    alpha    = 2*B1*(2*pi)^(-dsum)*(a0*sum(2.^((dsum-1)*Js)) + ...
                a1*sum(2.^((dsum-2)*Js)) + ...
                a2*sum(2.^((dsum-3)*Js)));
    sdf      = sdf + varRatio/2*alpha;
end
