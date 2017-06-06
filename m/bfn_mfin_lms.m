function [H, Omega, att] = bfn_mfin_lms(X, varargin)
% BFN_MFIN_LMS  Univariate LMS estimation of a multivariate fractionally integrated noise
%
% Syntax:
%   [H, Omega, att] = bfn_mfin_lms(X, 'Property Name 1','Property value 1',...)
%
% Description:
%   It estimates the Hurst exponent and covariance matrix of a
%   multivariate fractionally integrated noise via the univariate
%   least-mean squares (LMS) method.
%
% Input Arguments:
%   X : a matrix of multiple time series
%   
% Options:
%   method - the estimation method. Default is 'US'.
%     'US' - the univariate LMS method
%   range - a vector of scale range. Default is [1 100].
%   wavelet - the type of wavelets. Default is 'modwt'.
%     'dwt' - discrete wavelet transform
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
%   fixMaxScale - whether to fix the maximum wavelet scale. 
%                   Default is 1.
%   minScaleRange - the minimum interval between maximum and minimum
%                   scales. Default is 2.
%   maxerr - the maximum allowance of LMS error. Default is 1e-2.
%   discardCor - whether to discard the scale correlation from LMS
%                   factors. Default is 0.
%   omegamode - the covariance computation mode.
%     'lin' - based on the linearity of wavelet covariances
%     'sdf' - based on the spectral density.
%     'cov' - based on the covariance.
%   verbose - whether to display the text messages during
%                   processing. Default is 1.
%   isplot - whether display the plot. The plot is not displayed 
%                   if isplot equal to 0. Default is 1.
%
% Output Arguments:
%   H - the Hurst exponent
%   Omega - the covariance matrix of short memory
%   att - the list of estimator attributes
%
% Examples:
%   [H, Omega, att] = bfn_mfin_lms(x,'method','LN','range',[1 100],...
%                    'wavelet','modwt','filter','la8');
%
% References:
%   Achard, S., Bassett, D. S., Meyer-Lindenberg, A., & Bullmore, E.
%    (2008). Fractal connectivity of long-memory networks. 
%    Physical Review E, 77(3), 1-12.
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_mfin_lms.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

params = struct('method'        ,'US',...
                'wavelet'       ,'dwt',...
                'filter'        ,'la8',...
                'boundary'      ,'periodic',...
                'range'         ,[1 100],...
                'fixMaxScale'   ,0,...
                'minScaleRange' ,2,...  
                'maxerr'        ,1e-1,...
                'omegamode'     ,'sdf',...
                'verbose'       ,1,...
                'isplot'        ,0);
params = bfn_parseArgs(varargin,params);

N  = size(X,1);
Q  = size(X,2);
Np = 2*Q;       % the number of parameters

% wavelet coefficients
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
            W{j}(1:NJ(j),1) = x_bw{j};  %x_dwt for biased
        end

        for q=2:Q
            x_dwt = dwt(X(:,q), params.filter, 'conservative', params.boundary);
            x_bw  = bfn_dwt_brick_wall(x_dwt, params.filter, N);
            for j=1:J
                W{j}(1:NJ(j),q) = x_bw{j}; %x_dwt for biased
            end
        end 
    otherwise
end

J1 = params.range(1);
J2 = min(J,params.range(2));
Js = J1:J2;

% wavelet variances
switch params.wavelet
    case {'dwt','dwt-ext'}
        wvar = zeros(J,Q);
        for q = 1:Q
            Wq = getWq(W,q);
            wvar(:,q) = bfn_wave_var_dwt(Wq, N);
        end

        wcov = zeros(J,Q,Q);
        wcor = zeros(J,Q,Q);
        for q1 = 1:Q-1
            for q2 = q1+1:Q
                W1 = getWq(W,q1);
                W2 = getWq(W,q2);

                wcov(:,q1,q2) = bfn_wave_cov_dwt(W1, W2);
                wcor(:,q1,q2) = bfn_wave_cor_dwt(W1, W2, N);

                wcov(:,q2,q1) = wcov(:,q1,q2);
                wcor(:,q2,q1) = wcor(:,q1,q2);
            end
        end

        for q = 1:Q
            wcov(:,q,q) = wvar(:,q);
            wcor(:,q,q) = 1;
        end 

        % wsum
        switch params.omegamode
            case {'cov','sdf'}
                wsum = zeros(J,Q,Q);
                for q1 = 1:Q
                    for q2 = q1:Q
                        W1 = getWq(W,q1);
                        W2 = getWq(W,q2);

                        wsum(:,q1,q2) = bfn_wave_sum_dwt(W1, W2);
                        wsum(:,q2,q1) = wsum(:,q1,q2);
                    end
                end
            otherwise
        end
    otherwise
end

% main
switch params.method
    case {'US'}     
        
        % compute LMS over wavelet scale ranges.
        lms     = zeros(J2,J2);
        par     = zeros(Np,J2,J2);
        lmsMax  = 0;        
        for m = J1:J2+1-params.minScaleRange
            
            if params.fixMaxScale
                nStart = J2;
            else
                nStart = m+params.minScaleRange-1;
            end
            
            for n = nStart:J2
                [lmp, lms(m,n)]  = lmsuni(wcov,[m,n]);                        
                par(1:Np,m,n)    = lmp;

                if lms(m,n) > lmsMax
                    lmsMax = lms(m,n);
                end
            end
        end
        
        lms(isnan(lms)) = 0;
        lms(~lms) = lmsMax;

        % find optimal wavelet scale range.
        if lmsMax == 0
            optRange = [J1 J2];
            optLMS = 0;
            lmsNorm = 0;
        else
            lmsNorm  = lms/lmsMax;
            maxInter = -10;
            optRange = [0,0];
            lmsMin   = 1e10;
            cnt      = 0;
            for m = J1:J2-1
                
                if params.fixMaxScale
                    nStart = J2;
                else
                    nStart = m+params.minScaleRange-1;
                end
                
                for n = nStart:J2
                    if (lmsNorm(m,n)>0)&&(lmsNorm(m,n)<params.maxerr)
                        if (n-m+1>maxInter) || ((n-m+1==maxInter)&&(m>optRange(1)))
                            maxInter    = n-m+1;
                            optRange    = [m,n];
                            optLMS      = lmsNorm(m,n);
                            cnt         = cnt+1;
                        end
                    end

                    if lmsNorm(m,n) < lmsMin
                        lmsMin    = lmsNorm(m,n);
                        optRange2 = [m,n];
                        optLMS2   = lmsNorm(m,n);
                    end

                    if params.verbose
                        disp(strcat( ...               
                                'Scale Range = (',num2str(m),',',num2str(n), ...            
                                '),   LMS = ',num2str(lmsNorm(m,n)) ...
                                ));
                    end
                end
            end 

            if cnt<1
                optRange    = optRange2;
                optLMS      = optLMS2;
                if params.verbose
                    disp('Warning: There is no pair which satisfies the original criteria.');
                end
            end
        end

        % obtain d
        optJ1 = optRange(1); optJ2 = optRange(2);
        minlsq = lms(optJ1,optJ2);
        optPar = par(:,optJ1,optJ2);
        d = optPar(1:Q);        

        % obtain covariance
        switch params.omegamode                
            case 'lin'                   
                C = zeros(Q,Q);                 
                for i = 1:Q
                    for j=i:Q
                        lwc = log2(squeeze(abs(wcov(Js,i,j))));
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
                        ws = squeeze(wsum(:,i,j))';
                        xcov = sum(ws(Js)./2.^Js)/N;    
                        Omega(i,j) = xcov/(2*B1*cos(pi/2*(d(i)-d(j)))*(2*pi)^(1-d(i)-d(j))*sum(2.^((d(i)+d(j)-1)*Js)));
                        %Omega(i,j) = Xcov(i,j)/(2*B1*(2*pi)^(1-d(i)-d(j))*sum(2.^((d(i)+d(j)-1)*(1:10))));
                        Omega(j,i) = Omega(i,j);
                    end
                end
            case 'sdf'
                wsteps = J;
                wsize = N - N/2^wsteps;

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

        % messaging
        if params.verbose
            disp('== Selected Scale Range ==');
            disp(strcat( ...               
                        'Scale Range = (',num2str(optJ1),...
                        ',',num2str(optJ2), ...  
                        '),   Normalized LMS = ',num2str(optLMS) ...   
                        ));
        end

        % plot
        if params.isplot
            figure;
            Z = lms(J1:J2-1,J1+1:J2);
            imagesc(Z);
            xlabel('wavelet scale, end'); ylabel('wavelet scale, start');
            set(gca, 'XTick',[1:J2-J1],'YTick',[1:J2-J1],...
                'XTickLabel',[J1+1:J2],'YTickLabel',[J1:J2-1],...
                'XAxisLocation','bottom','YAxisLocation','right','YDir','normal');
            colorbar;
        end

    otherwise
end

H = d+.5;
att = struct('method'       ,params.method,...
            'wavelet'       ,params.wavelet,...
            'range'         ,params.range,...
            'fixMaxScale'   ,params.fixMaxScale,...
            'minScaleRange' ,params.minScaleRange,...
            'maxerr'        ,params.maxerr,...
            'omegamode'     ,params.omegamode,...
            'filter'        ,params.filter,...
            'boundary'      ,params.boundary,...
            'optpar'        ,optPar,...
            'minlsq'        ,minlsq,...
            'par'           ,par,...
            'lsq'           ,lmsNorm);

end

function Wq = getWq(W,q)

    J = length(W);
    
    Wq = {};
    for j=1:J
        Wj = W{j};
        Wq{j} = squeeze(Wj(:,q));
    end

end

function [lmp, lms] = lmsuni(wcov, range)

J   = size(wcov,1);
J1  = range(1);
J2  = min(J,range(2));
Jp  = J1:J2;
dif = J2-J1+1;

Q = size(wcov,2);   % the number of elements
d = zeros(Q,1);
c = zeros(Q,1);

for q = 1:Q
    lwc  = log2(squeeze(wcov(Jp,q,q)));
    p    = polyfit(Jp',lwc,1);                    
    d(q) = p(1)/2;
    c(q) = p(2); 
end   
lmp = [d;c];

% compute the least mean square (LMS)
L = 0; 
for q = 1:Q
    va = squeeze(wcov(:,q,q));
    lv = log2(abs(va));
    lv(isnan(lv))=0;
    L = L + sum((lv(Jp) - 2*d(q)*Jp' - c(q)).^2);
end
lms = L/dif;

end

