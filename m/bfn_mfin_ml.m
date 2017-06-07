function [H, cor, att] = bfn_mfin_ml(X, varargin)
% BFN_MFIN_ML Univariate maximum likelihood estimation for a multivariate fractionally integrated noise
%
% Syntax:
%   [H, cor, att] = bfn_mfin_ml(X, 'Property Name 1','Property value 1',...)
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
%   verbose - whether to display the text messages during
%                   processing. Default is 1.
%
% Output Arguments:
%   H - the Hurst exponent
%   cor - the correlation matrix of short memory
%   att - the list of estimator attributes
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_mfin_ml.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

 
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
                'verbose'   ,1);
if nargin > 0
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

    otherwise
end

J1    = params.range(1);
J2    = min(J,params.range(2));  
range = [J1, J2];
Js    = J1:J2;

%% wavelet variances
switch params.wavelet
    case {'dwt','dwt-ext'}
        switch params.omegamode
            case 'sdf'
                wsum = zeros(J,Q,Q);
                for q1 = 1:Q
                    for q2 = q1:Q
                        W1 = getWq(W,q1);
                        W2 = getWq(W,q2);

                        wsum(:,q1,q2) = bfn_wave_sum_dwt(W1, W2);
                        wsum(:,q2,q1) = wsum(:,q1,q2);
                    end
                end
            case 'lin'
                wcov = zeros(J,Q,Q);
                for q1 = 1:Q
                    for q2 = q1:Q
                        W1 = getWq(W,q1);
                        W2 = getWq(W,q2);
                        wcov(:,q1,q2) = bfn_wave_cov_dwt(W1, W2);
                        wcov(:,q2,q1) = wcov(:,q1,q2);
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
                    wsum(:,q) = bfn_wave_sum_dwt(W1, W1);
                end
                
                wsum2 = zeros(J,Q,Q);
                for q1 = 1:Q
                    for q2 = q1:Q
                        W1 = getWq(W,q1);
                        W2 = getWq(W,q2);

                        wsum2(:,q1,q2) = bfn_wave_sum_dwt(W1, W2);
                        wsum2(:,q2,q1) = wsum2(:,q1,q2);
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
            switch params.omegamode
                case 'sdf'
                    ws = squeeze(wsum(:,q,q));
                otherwise
                    ws = squeeze(wsum(:,q));
            end

            va = var(X(:,q));
            d(q) = fin_ml(ws, va, N, NJ, range, params);

            if params.verbose   
                disp(strcat('q= ',num2str(q),': H= ',num2str(d(q)+.5)));
            end
        end

    case 'MR'
        for q = 1:Q
            switch params.omegamode
                case 'sdf'
                    ws = squeeze(wsum(:,q,q));
                otherwise
                    ws = squeeze(wsum(:,q));
            end

            va = var(X(:,q));
            d(q) = fin_mr(ws, va, N, NJ, range, params);

            if params.verbose   
                disp(strcat('q= ',num2str(q),': H= ',num2str(d(q)+.5)));
            end
        end

    otherwise
end

%% compute short memory covariance
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
                ws = squeeze(wsum(:,i,j));
                Omega(i,j) = sum(ws(1:J)./sdf)/wsize; 
                Omega(j,i) = Omega(i,j);
            end
        end     
    otherwise
end

H = d+.5;
cor = bfn_CovToCor(Omega);
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

end

%% function: fin_ml
function [d, sigma] = fin_ml(wsum, va, N, NJ, range, params)
        
    J1 = range(1);
    J2 = range(2);     
    Js = J1:J2;

    wsize = N;

    % initial parameters
    lb = params.lb(1);
    ub = params.ub(1);

    sdf = bfn_finsdf(params.init(1),[J1,J2],NaN);
    signew2 = sum(wsum(Js)./sdf)/wsize;     

    % iteration
    it = 1;
    sigold2 = signew2+params.abstol+10;                    
    while ((abs(1-signew2/sigold2)>params.abstol)) && (it<=params.maxit) 
        sigold2 = signew2;
        if isempty(params.options)
            options = optimset('Algorithm'      ,'trust-region-reflective',...
                               'MaxFunEvals'    ,300,...
                               'MaxIter'        ,50,...
                               'TolFun'         ,1e-5,...
                               'TolX'           ,1e-7,...
                               'Diagnostics'    ,'off',...
                               'Display'        ,'off');
        else
            options = params.options;
        end
        f = @(x)likelihood(x,wsum,NJ,[J1,J2],signew2);
        [est,fval] = fminbnd(f,lb,ub,options);
        dnew = est(1);

        switch params.omegamode
            case 'lin'
            case 'cov'
                B1 = (1 - 2.^(2*dnew-1))./(1-2*dnew);    
                signew2 = va/(2*B1*(2*pi)^(1-2*dnew)*sum(2.^((2*dnew-1)*(1:10))));
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

    d = dnew;
    sigma = signew2; 
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
                           'TolFun'         ,1e-5,...
                           'TolX'           ,1e-7,...
                           'Algorithm'      ,'trust-region-reflective',...
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
            sigma = va/(2*B1*(2*pi)^(1-2*d)*sum(2.^((2*d-1)*(1:10))));
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
