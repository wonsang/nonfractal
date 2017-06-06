%% Estimation of nonfractal connectivity by simulation studies
% This demo shows an example of estimating nonfractal connectivity of
% an fMRI BOLD signals of the brain.
% We applied the ML-LIN estimator to a resting state fMRI data of the rat
% brain so as to estimate the nonfractal connectivity. The following figure
% shows the comparison of Pearson correlation and estimated nonfractal
% connectivity.
%
% References
% Wonsang You, Sophie Achard, Joerg Stadler, Bernd Bruekner, and Udo
% Seiffert, "Fractal analysis of resting state functional connectivity of
% the brain," in 2012 International Joint Conference on Neural Networks,
% 2012. 

warning off

fprintf('Processing an fMRI data...\n');
Y = bfn_readtable('bfn_dat_fmri_roi.txt');

[H_FMRI, Omega_FMRI] = bfn_mfin_ml(Y', ...         
                    'range'     ,[2 100],...
                    'abstol'    ,1e-10,...
                    'maxit'     ,100,...
                    'omegamode' ,'lin',...
                    'verbose'   ,0);

C_FMRI      = corr(Y') - eye(15);
R_FMRI      = bfn_CovToCor(Omega_FMRI);

fprintf('Creating plots...\n');
figure('position',[100 100 450 450]);
colormap Jet

subplot(1,2,1); imagesc(R_FMRI,[0 1]);
title('Nonfractal connectivity (fMRI)'  ,'FontSize',10); set(gca,'PlotBoxAspectRatio',[1.2 1 1])

subplot(1,2,2); imagesc(C_FMRI,[0 1]);
title('Pearson correlation (fMRI)'      ,'FontSize',10); set(gca,'PlotBoxAspectRatio',[1.2 1 1])

drawnow
set(gcf, 'Name', 'Nonfractal connectivity in fMRI data');

fprintf(' DONE.\n');

warning on