function sdf = bfn_finsdf(d, range, SNR)
% BFN_FINSDF Computation of scale-dependent variances of a fractionally integrated noise.
%
% Syntax:
%   sdf = bfn_finsdf(d, range, SNR)
%
% Description:
%   It computes the scale-dependent variances of a fractionally
%   integrated noise (FIN) given its long memory parameter. Also, it
%   supports computation for a FIN process corrupted by white noise.
%
% Input Arguments:
%   d : a scalar or vector of memory parameters
%   range : the range of wavelet scales
%   SNR : the signal-to-noise ratio
%
% Output Arguments:
%   sdf : a vector of scale-dependent variances
%
% References:
%   Achard, S., Bassett, D. S., Meyer-Lindenberg, A., & Bullmore, E.
%   (2008). Fractal connectivity of long-memory networks. 
%   Physical Review E, 77(3), 1-12.
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_finsdf.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

Js   = (range(1):range(2))';

dd = 2*d;    

B1   = (1 - 2^(dd-1))/(1-dd);
B3   = (1 - 2^(dd-3))/(3-dd);
A    = d/12;   
a2   = (2*pi)^2*A*B3/B1;
K    = 2*B1*(2*pi)^(-dd);

sdf  = K*2.^(Js*dd).*(1 + a2./(2.^(2*Js)));

if ~isnan(SNR)   
    varRatio = 10^(-SNR/10);
    sdf      = sdf + varRatio;
end
