function [H, nfcor, fcor, att] = bfn_mfin_ml(X, varargin)
% BFN_MFIN_ML Univariate maximum likelihood estimation for a multivariate fractionally integrated noise
%
% Syntax:
%   [H, nfcor, fcor, att] = bfn_mfin_ml(X, 'Property Name 1','Property value 1',...)
%
% Description:
%   It estimates the Hurst exponent and covariance matrix of a
%   multivariate fractionally integrated noise via the univariate maximum
%   likelihood method.
%
% Input Arguments:
%   X : a matrix of multiple time series
%   
% Options:
%   method - the estimation method. Default is 'ML'.
%     'ML' - the standard ML method
%     'MR' - the reduced ML method
%   range - a vector of scale range. Default is [1 100].
%   wavelet - the type of wavelets. Default is 'dwt'.
%     'dwt' - discrete wavelet transform
%     'modwt' - maximal overlap discrete wavelet transform
%   boundary - the mode of boundary coefficients. Default is 'periodic'.
%   filter - The mode of wavelet filter. The list of filter
%            modes are shown below. Default is 'la8'.
%
%     Extremal phase filters
%       - Haar, D4, D6, D8, D10, D12, D14, D16, D18, D20
%     Least asymmetric filters
%       - LA8, LA10, LA12, LA14, LA16, LA18, LA20
%     Best Localized filters
%       - BL14, BL18, BL20
%     Coiflet filters
%       - C6, C12, C18, C24, C30
%
%   abstol - the absolute tolerance.
%   maxit - the maximum iterations.
%   init - the initial vector for maximum likelihood estimation.
%   lb - the two-element vector for the lower bound
%   ub - the two-element vector for the upper bound
%   omegamode - the covariance computation mode. Default is 'cov'.
%     'lin' - based on the linearity of wavelet covariances
%     'sdf' - based on the spectral density.
%     'cov' - based on the covariance.
%   d_only - only estimate the Hurst exponents (H).
%   verbose - whether to display the text messages during
%                   processing. Default is 1.
%
% Output Arguments:
%   H     - the Hurst exponent
%   nfcor - the correlation matrix of short memory
%   fcor  - the asymptotic wavelet correlation
%   att   - the list of estimator attributes
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_mfin_ml.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

warning('off','all')
 
params = struct('method'    ,'ML',...
                'wavelet'   ,'dwt',...
                'filter'    ,'la8',...
                'boundary'  ,'periodic',...
                'range'     ,[1 100],...
                'abstol'    ,1e-10,...
                'maxit'     ,100,...
                'init'      ,[.1 5],...
                'lb'        ,[-0.5 0],...
                'ub'        ,[0.5 10],...
                'options'   ,[],...
                'omegamode' ,'cov',...
                'd_only'    ,0,...
                'verbose'   ,1);
if nargin > 1
    if iscell(varargin{1})
        params = bfn_parseArgs(varargin{1},params);
    else
        params = bfn_parseArgs(varargin,params);
    end
end

Q = size(X,2);

%% wavelet coefficients
switch params.wavelet
    case 'dwt'
        Jmax = floor(log2(size(X,1)));
        N    = 2^Jmax;
        X    = X(1:N,:);

        x_dwt = dwt(X(:,1), params.filter, 'conservative', params.boundary);
        x_bw  = bfn_dwt_brick_wall(x_dwt, params.filter, N);
        J     = length(x_bw);
        W     = cell(1,J);
        for j=1:J
            NJ(j) = length(x_bw{j});
            W{j}(1:NJ(j),1) = x_bw{j};
        end

        for q=2:Q
            x_dwt = dwt(X(:,q), params.filter, 'conservative', params.boundary);
            x_bw  = bfn_dwt_brick_wall(x_dwt, params.filter, N);
            for j=1:J
                W{j}(1:NJ(j),q) = x_bw{j};
            end
        end 
        
    case 'modwt'
        
        x_dwt = modwt(X(:,1), params.filter, 'conservative', 'circular');
        
        N     = size(x_dwt,1);
        J     = size(x_dwt,2);
        
        Jmax  = floor(log2(size(X,1)));
        NJ    = 2.^(Jmax-1:-1:Jmax-J);
        
        x_bw  = bfn_modwt_brick_wall(x_dwt, params.filter, N);

        W = cell(1,J);
        for j=1:J
            W{j} = zeros(N,Q);
            W{j}(:,1) = squeeze(x_bw(:,j));  
        end
        
        for q=2:Q
            x_dwt = modwt(X(:,q), params.filter, 'conservative', 'circular');
            x_bw  = bfn_modwt_brick_wall(x_dwt, params.filter, N);
            for j=1:J
                W{j}(:,q) = squeeze(x_bw(:,j));
            end
        end 

    otherwise
end

J1    = max(1,params.range(1));
J2    = min(J,params.range(2));  
range = [J1, J2];
Js    = J1:J2;

%% wavelet variances
switch params.wavelet
    case {'dwt'}
        switch params.omegamode
            case 'sdf'
                wsum = zeros(J,Q,Q);
                for q1 = 1:Q
                    for q2 = q1:Q
                        W1 = getWq(W,q1);
                        W2 = getWq(W,q2);

                        wsum(:,q1,q2) = abs(bfn_wave_sum_dwt(W1, W2));
                        wsum(:,q2,q1) = wsum(:,q1,q2);
                    end
                end
            case 'lin'
                
                if ~params.d_only
                    wcov = zeros(J,Q,Q);
                    for q1 = 1:Q
                        for q2 = q1:Q
                            W1 = getWq(W,q1);
                            W2 = getWq(W,q2);
                            wcov(:,q1,q2) = bfn_wave_cov_dwt(W1, W2);
                            wcov(:,q2,q1) = wcov(:,q1,q2);
                        end
                    end  
                end

                wsum = zeros(J,Q);
                for q = 1:Q
                    W1 = getWq(W,q);
                    wsum(:,q) = bfn_wave_sum_dwt(W1, W1);
                end

            case 'cov'
                wsum = zeros(J,Q);
                for q = 1:Q
                    W1 = getWq(W,q);
                    wsum(:,q) = abs(bfn_wave_sum_dwt(W1, W1));
                end
                
                if ~params.d_only
                    wsum2 = zeros(J,Q,Q);
                    for q1 = 1:Q
                        for q2 = q1:Q
                            W1 = getWq(W,q1);
                            W2 = getWq(W,q2);

                            wsum2(:,q1,q2) = abs(bfn_wave_sum_dwt(W1, W2));
                            wsum2(:,q2,q1) = wsum2(:,q1,q2);
                        end
                    end
                end
        end
        
    case {'modwt'}
        
        switch params.omegamode
            case 'sdf'
                wsum = zeros(J,Q);
                wsum2 = zeros(J,Q,Q);
                for q1 = 1:Q
                    W1 = getWq(W,q1);
                    W1 = cell2mat(W1);
                    wsum2(:,q1,q1) = modwt_wvar(W1,'gaussian','unbiased',params.filter);
                    wsum(:,q1) = wsum2(:,q1,q1);
                    
                    if ~params.d_only
                        for q2 = q1:Q
                            if q1 ~= q2
                                W2 = getWq(W,q2);
                                W2 = cell2mat(W2);
                                wsum2(:,q1,q2) = modwt_wcov(W1, W2,'gaussian','unbiased',params.filter);
                                wsum2(:,q2,q1) = wsum2(:,q1,q2);
                            end
                        end
                    end
                end
            case 'lin'
                wsum = zeros(J,Q);
                wcov = zeros(J,Q,Q);
                for q1 = 1:Q
                    W1 = getWq(W,q1);
                    W1 = cell2mat(W1);
                    wcov(:,q1,q1) = modwt_wvar(W1,'gaussian','unbiased',params.filter);
                    wsum(:,q1) = wcov(:,q1,q1);
                    
                    if ~params.d_only
                        for q2 = q1:Q
                            if q1 ~= q2
                                W2 = getWq(W,q2);
                                W2 = cell2mat(W2);
                                wcov(:,q1,q2) = modwt_wcov(W1, W2,'gaussian','unbiased',params.filter);
                                wcov(:,q2,q1) = wcov(:,q1,q2);
                            end
                        end
                    end
                end
            case 'cov'
                wsum = zeros(J,Q);
                wsum2 = zeros(J,Q,Q);
                for q1 = 1:Q
                    W1 = getWq(W,q1);
                    W1 = cell2mat(W1);
                    wsum2(:,q1,q1) = modwt_wvar(W1,'gaussian','unbiased',params.filter);
                    wsum(:,q1) = wsum2(:,q1,q1);
                    
                    if ~params.d_only
                        for q2 = q1:Q
                            if q1 ~= q2
                                W2 = getWq(W,q2);
                                W2 = cell2mat(W2);
                                wsum2(:,q1,q2) = modwt_wcov(W1, W2,'gaussian','unbiased',params.filter);
                                wsum2(:,q2,q1) = wsum2(:,q1,q2);
                            end
                        end
                    end
                end
        end

    otherwise
end

%% compute memory parameters
d = zeros(1,Q);
switch params.method
    case 'ML'            
        for q = 1:Q
            ws = squeeze(wsum(:,q));

            va = var(X(:,q));
            d(q) = fin_ml(ws, va, N, NJ, range, params);

            if params.verbose   
                disp(strcat('q= ',num2str(q),': H= ',num2str(d(q)+.5)));
            end
        end

    case 'MR'
        for q = 1:Q
            ws = squeeze(wsum(:,q));

            va = var(X(:,q));
            d(q) = fin_mr(ws, va, N, NJ, range, params);

            if params.verbose   
                disp(strcat('q= ',num2str(q),': H= ',num2str(d(q)+.5)));
            end
        end

    otherwise
end

H = d+.5;

%% compute short memory covariance
if ~params.d_only
    switch params.omegamode
        case 'lin'
            C = zeros(Q,Q);             
            for i = 1:Q
                for j=i:Q
                    lwc = log2(abs(squeeze(wcov(Js,i,j))));
                    C(i,j) = mean(lwc' - (d(i)+d(j))*Js);
                    C(j,i) = C(i,j);
                end
            end
            Omega = bfn_getsigma_m(d,C);

        case 'cov'
            %Xcov = cov(X);             
            Omega = zeros(Q);
            for i = 1:Q
                for j=i:Q
                    B1 = (1 - 2.^(d(i)+d(j)-1))./(1-d(i)-d(j)); 
                    ws = squeeze(wsum2(:,i,j))';
                    xcov = sum(ws(Js)./2.^Js)/N;    
                    Omega(i,j) = xcov/(2*B1*cos(pi/2*(d(i)-d(j)))*(2*pi)^(1-d(i)-d(j))*sum(2.^((d(i)+d(j)-1)*Js)));
                    %Omega(i,j) = Xcov(i,j)/(2*B1*(2*pi)^(1-d(i)-d(j))*sum(2.^((d(i)+d(j)-1)*(1:10))));
                    Omega(j,i) = Omega(i,j);
                end
            end
        case 'sdf'
            wsteps = J;
            wsize  = N - N/2^wsteps;

            Omega = zeros(Q);
            for i = 1:Q
                for j=i:Q
                    sdf = bfn_finsdf_b(d(i),d(j),[1,J],NaN);
                    ws = squeeze(wsum2(:,i,j));
                    Omega(i,j) = sum(ws(1:J)./sdf)/wsize; 
                    Omega(j,i) = Omega(i,j);
                end
            end     
        otherwise
    end
    
    nfcor = bfn_CovToCor(Omega);
    fcor  = bfn_acorr(H, Omega);
else
    nfcor = [];
    fcor  = [];
end

att = struct('method'   ,params.method,...
            'wavelet'   ,params.wavelet,...
            'range'     ,params.range,...
            'filter'    ,params.filter,...
            'boundary'  ,params.boundary,...
            'abstol'    ,params.abstol,...
            'maxit'     ,params.maxit,...
            'init'      ,params.init,...
            'lb'        ,params.lb,...
            'ub'        ,params.ub,...
            'options'   ,params.options,...
            'omegamode' ,params.omegamode);

warning('on','all')
        
end

%% function: fin_ml
function [d, sigma] = fin_ml(wsum, va, N, NJ, range, params)
        
    J1 = range(1);
    J2 = range(2);     
    Js = J1:J2;

    wsize = N;
    
    %wsum = wsum/max(wsum(:));

    % initial parameters
    lb = params.lb(1);
    ub = params.ub(1);

    dnew = params.init(1);
    switch params.omegamode
        case {'cov','lin'}
            B1 = (1 - 2.^(2*dnew-1))./(1-2*dnew);  
            signew2 = va/(2*B1*(2*pi)^(1-2*dnew)*sum(2.^((2*dnew-1)*(Js))));
        case 'sdf'
            sdf = bfn_finsdf(dnew,[J1,J2],NaN);
            signew2 = sum(wsum(Js)./sdf)/wsize;
    end

    % iteration
    it = 1;
    sigold2 = signew2/(1-params.abstol)+10;                    
    while ((abs(1-signew2/sigold2)>params.abstol)) && (it<=params.maxit) 
        sigold2 = signew2;
        if isempty(params.options)
            options = optimset('MaxFunEvals'    ,1000,...
                               'MaxIter'        ,500,...
                               'TolFun'         ,1e-10,...
                               'TolX'           ,1e-10,...
                               'Diagnostics'    ,'off',...
                               'Display'        ,'off');
        else
            options = params.options;
        end
        f = @(x)likelihood(x,wsum,NJ,[J1,J2],signew2);
        [est,fval] = fminbnd(f,lb,ub,options);
        dnew = est(1);

        switch params.omegamode
            case {'cov','lin'}
                B1 = (1 - 2.^(2*dnew-1))./(1-2*dnew);    
                signew2 = va/(2*B1*(2*pi)^(1-2*dnew)*sum(2.^((2*dnew-1)*(Js))));
            case 'sdf'
                sdf = bfn_finsdf(dnew,[J1,J2],NaN);
                signew2 = sum(wsum(Js)./sdf)/wsize; 
        end

        if params.verbose                            
            disp(strcat('Rep. ',num2str(it),...
                        ': Signew2= ',num2str(signew2),...
                        ': fval= ',num2str(fval),...
                        ', d= ',num2str(dnew)));
        end
        it = it+1;
    end

    if it>1
        d = dnew;
        sigma = signew2; 
    else
        d = NaN;
        sigma = NaN;
    end
    
end

%% function : fin_mr
function [d, sigma] = fin_mr(wsum, va, N, NJ, range, params)

    J1 = range(1);
    J2 = range(2); 
    Js = J1:J2;

    wsize = N; 

    % initial parameters
    lb = params.lb(1);
    ub = params.ub(1);

    % estimation
    if isempty(params.options)
        options = optimset('MaxFunEvals'    ,300,...
                           'MaxIter'        ,50,...
                           'TolFun'         ,1e-10,...
                           'TolX'           ,1e-10,...
                           'Diagnostics'    ,'off',...
                           'Display'        ,'off');                      
    else
        options = params.options;
    end            
    f = @(x)likelihood_reduced(x,wsum,NJ,N,[J1,J2]);
    est = fminbnd(f,lb,ub,options);
    d = est(1);

    switch params.omegamode
        case 'sdf'
            sdf   = bfn_finsdf(d,[J1,J2],NaN);
            sigma = sum(wsum(Js)./sdf)/wsize; 
        case 'cov'            
            B1    = (1 - 2.^(2*d-1))./(1-2*d);    
            sigma = va/(2*B1*(2*pi)^(1-2*d)*sum(2.^((2*d-1)*(Js))));
    end                      
end

%% function : likelihood
function L = likelihood(x,wsum,Nj,range,sig)

    d   = x(1);  
    
    sdf = bfn_finsdf(d,range,NaN); 
    Js  = (range(1):range(2))';
    S   = sum(wsum(Js)./sdf);
    L   = sum(Nj(Js)'.*log(sdf)) + S/sig;
    
end

%% function : likelihood_reduced
function L = likelihood_reduced(x,wsum,Nj,N,range)

    d   = x(1);  
    
    sdf = bfn_finsdf(d,range,NaN); 
    Js  = (range(1):range(2))';
    S   = sum(wsum(Js)./sdf);    
    L   = N*log(S/N) + sum(Nj(Js)'.*log(sdf));
    
end

function Wq = getWq(W,q)

    J = length(W);
    
    Wq = {};
    for j=1:J
        Wj    = W{j};
        Wq{j} = squeeze(Wj(:,q));
    end

end
