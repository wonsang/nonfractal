function A = bfn_acorr(H, Omega)
% BFN_ACORR Computing asymptotic correlation for a bivariate FIN process.
%
% Syntax:
%   A = bfn_acorr(H, Omega)
%
% Description:
%   It computes the asymptotic correlation of a bivariate fractionally
%   integrated noise from its long memory parameters.
%
% Input Arguments:
%   H - a vector of Hurst exponents.
%   Omega - a covariance matrix.
%
% Output Arguments:
%   A - a matrix of asymptotic correlation.
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_acorr.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

nSeries = length(H);
A = zeros(nSeries,nSeries);
for s1 = 1:nSeries
    for s2 = s1:nSeries
        d = [H(s1),H(s2)] - .5;   
        omega = [Omega(s1,s1),Omega(s1,s2);Omega(s2,s1),Omega(s2,s2)];
        A(s1,s2) = acorrb(d,omega);            
        A(s2,s1) = A(s1,s2);
    end
end

end

function R = acorrb(d, Omega)

    % Asymtotic correlation based on FIN model
    
    d1 = d(1);
    d2 = d(2);

    dif1 = 1-2*d1;
    dif2 = 1-2*d2;
    dif12 = 1-d1-d2;

    B12 = (1 - 1/(2^dif12)) / dif12;
    B11 = (1 - 1/(2^dif1)) / dif1;
    B22 = (1 - 1/(2^dif2)) / dif2;

    b0 = cos((d1-d2)*pi/2);

    K = (Omega(1,2) / ((Omega(1,1)*Omega(2,2))^0.5)) * (B12 / ((B11*B22)^0.5));
    R = K*b0;
end
