
clear all; close all; clc;

%% Prep

% Load data
load('matfiles/SCA_Processed.mat');

% Load template
load('fsaverage_7.1.0.mat');
load('glasser_fsaverage_7.1.0.mat');
load('Vasa_G360.mat');

% Colormaps
loadColorMaps;

% Significance threshold
thresh = 0.05;
nPerm = 10000;

% Save subfolder
figDir = 'sca';

%% Compute statistics for X and Y dosage effects

% Compute z-scores for XXY cohort
XXY_CV = zscore(XXY.CV);
XXY_SA = zscore(XXY.SA);
XXY_CT = zscore(XXY.CT);

% Compute z-scores for XYY cohort
XYY_CV = zscore(XYY.CV);
XYY_SA = zscore(XYY.SA);
XYY_CT = zscore(XYY.CT);

% Compute global phenotypes
XXY.CV_GS = sum(XXY.CV,2);
XXY.SA_GS = sum(XXY.SA,2);
XXY.CT_GS = sum(XXY.CT .* (XXY.SA ./ sum(XXY.SA,2)),2);

XYY.CV_GS = sum(XYY.CV,2);
XYY.SA_GS = sum(XYY.SA,2);
XYY.CT_GS = sum(XYY.CT .* (XYY.SA ./ sum(XYY.SA,2)),2);

% Prepare terms
XXY_SEX = term(XXY.SEX);
XXY_AGE = term(XXY.AGE - mean(XXY.AGE));
XXY_EUL = term(XXY.EUL - mean(XXY.EUL));
XXY_G_CV = term(XXY.CV_GS - mean(XXY.CV_GS));
XXY_G_SA = term(XXY.SA_GS - mean(XXY.SA_GS));
XXY_G_CT = term(XXY.CT_GS - mean(XXY.CT_GS));

XYY_SEX = term(XYY.SEX);
XYY_AGE = term(XYY.AGE - mean(XYY.AGE));
XYY_EUL = term(XYY.EUL - mean(XYY.EUL));
XYY_G_CV = term(XYY.CV_GS - mean(XYY.CV_GS));
XYY_G_SA = term(XYY.SA_GS - mean(XYY.SA_GS));
XYY_G_CT = term(XYY.CT_GS - mean(XYY.CT_GS));

% Define GLM
XXY_M_CV = 1 + XXY_SEX + XXY_AGE + XXY_EUL + XXY_G_CV;
XXY_M_SA = 1 + XXY_SEX + XXY_AGE + XXY_EUL + XXY_G_SA;
XXY_M_CT = 1 + XXY_SEX + XXY_AGE + XXY_EUL + XXY_G_CT;

XYY_M_CV = 1 + XYY_SEX + XYY_AGE + XYY_EUL + XYY_G_CV;
XYY_M_SA = 1 + XYY_SEX + XYY_AGE + XYY_EUL + XYY_G_SA;
XYY_M_CT = 1 + XYY_SEX + XYY_AGE + XYY_EUL + XYY_G_CT;

% Compute stats for XXY
XXY_CV_SEX = SurfStatT(SurfStatLinMod(XXY_CV, XXY_M_CV), XXY_SEX);
XXY_CV_SEX_P = 2*(1 - tcdf(abs(XXY_CV_SEX.t), XXY_CV_SEX.df));
[~,~,~,XXY_CV_SEX_Q] = fdr_bh(XXY_CV_SEX_P);

XXY_SA_SEX = SurfStatT(SurfStatLinMod(XXY_SA, XXY_M_SA), XXY_SEX);
XXY_SA_SEX_P = 2*(1 - tcdf(abs(XXY_SA_SEX.t), XXY_SA_SEX.df));
[~,~,~,XXY_SA_SEX_Q] = fdr_bh(XXY_SA_SEX_P);

XXY_CT_SEX = SurfStatT(SurfStatLinMod(XXY_CT, XXY_M_CT), XXY_SEX);
XXY_CT_SEX_P = 2*(1 - tcdf(abs(XXY_CT_SEX.t), XXY_CT_SEX.df));
[~,~,~,XXY_CT_SEX_Q] = fdr_bh(XXY_CT_SEX_P);

% Compute stats for XYY
XYY_CV_SEX = SurfStatT(SurfStatLinMod(XYY_CV, XYY_M_CV), XYY_SEX);
XYY_CV_SEX_P = 2*(1 - tcdf(abs(XYY_CV_SEX.t), XYY_CV_SEX.df));
[~,~,~,XYY_CV_SEX_Q] = fdr_bh(XYY_CV_SEX_P);

XYY_SA_SEX = SurfStatT(SurfStatLinMod(XYY_SA, XYY_M_SA), XYY_SEX);
XYY_SA_SEX_P = 2*(1 - tcdf(abs(XYY_SA_SEX.t), XYY_SA_SEX.df));
[~,~,~,XYY_SA_SEX_Q] = fdr_bh(XYY_SA_SEX_P);

XYY_CT_SEX = SurfStatT(SurfStatLinMod(XYY_CT, XYY_M_CT), XYY_SEX);
XYY_CT_SEX_P = 2*(1 - tcdf(abs(XYY_CT_SEX.t), XYY_CT_SEX.df));
[~,~,~,XYY_CT_SEX_Q] = fdr_bh(XYY_CT_SEX_P);

% Save group maps
save('matfiles/sca_results.mat', 'XXY_CV_SEX', 'XXY_CV_SEX_Q', 'XXY_CV_SEX_P', ...
                                 'XXY_SA_SEX', 'XXY_SA_SEX_Q', 'XXY_SA_SEX_P', ...
                                 'XXY_CT_SEX', 'XXY_CT_SEX_Q', 'XXY_CT_SEX_P', ...
                                 'XYY_CV_SEX', 'XYY_CV_SEX_Q', 'XYY_CV_SEX_P', ...
                                 'XYY_SA_SEX', 'XYY_SA_SEX_Q', 'XYY_SA_SEX_P', ...
                                 'XYY_CT_SEX', 'XYY_CT_SEX_Q', 'XYY_CT_SEX_P', ...
                                 '-v7.3');

%% Visualize X and Y dosage effect maps

% XXY
h = HyoSurfacePlotFunction(XXY_CV_SEX.coef(2,:) .* (XXY_CV_SEX_Q < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xxy_cv_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(XXY_SA_SEX.coef(2,:) .* (XXY_SA_SEX_Q < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xxy_sa_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(XXY_CT_SEX.coef(2,:) .* (XXY_CT_SEX_Q < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xxy_ct_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(XXY_CV_SEX.coef(2,:), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xxy_cv.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(XXY_SA_SEX.coef(2,:), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xxy_sa.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(XXY_CT_SEX.coef(2,:), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xxy_ct.png'), 'Resolution', 300); close(h);

% XYY
h = HyoSurfacePlotFunction(XYY_CV_SEX.coef(2,:) .* (XYY_CV_SEX_Q < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xyy_cv_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(XYY_SA_SEX.coef(2,:) .* (XYY_SA_SEX_Q < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xyy_sa_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(XYY_CT_SEX.coef(2,:) .* (XYY_CT_SEX_Q < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xyy_ct_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(XYY_CV_SEX.coef(2,:), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xyy_cv.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(XYY_SA_SEX.coef(2,:), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xyy_sa.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(XYY_CT_SEX.coef(2,:), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'xyy_ct.png'), 'Resolution', 300); close(h);

%% Regional congruence between HCP and XXY cohorts

% Conjunction map for CV
XXY_CV_AD = zeros(1,360);
HCP_CV_AD = zeros(1,360);
CV_AD = zeros(1,360);

% Annote regions based on direction and significance of X-dosage effects
XXY_CV_AD(XXY_CV_SEX.coef(2,:) > 0) =  1;
XXY_CV_AD(XXY_CV_SEX.coef(2,:) < 0) = -1;
XXY_CV_AD((XXY_CV_SEX_Q < thresh) & (XXY_CV_SEX.coef(2,:) > 0))  =  2;
XXY_CV_AD((XXY_CV_SEX_Q < thresh) & (XXY_CV_SEX.coef(2,:) < 0))  = -2;

% Annote regions based on direction and significance of sex effects
HCP_CV_AD(HCP_CV_SEX.coef(2,:) > 0) =  1;
HCP_CV_AD(HCP_CV_SEX.coef(2,:) < 0) = -1;
HCP_CV_AD((HCP_CV_SEX_Q < thresh) & (HCP_CV_SEX.coef(2,:) > 0))  =  2;
HCP_CV_AD((HCP_CV_SEX_Q < thresh) & (HCP_CV_SEX.coef(2,:) < 0))  = -2;

% Compute regional conjunction
CV_AD((HCP_CV_AD < 2) & (HCP_CV_AD > -2)) = 0;

CV_AD((HCP_CV_AD == 2) & (XXY_CV_AD == 2)) = -1.9; % Congruent
CV_AD((HCP_CV_AD == 2) & (XXY_CV_AD == 1)) = -1.6; % Congruent

CV_AD((HCP_CV_AD == -2) & (XXY_CV_AD == 2)) = -0.9; % Incongruent
CV_AD((HCP_CV_AD == -2) & (XXY_CV_AD == 1)) = -0.6; % Incongruent

CV_AD((HCP_CV_AD == 2) & (XXY_CV_AD == -2)) = 0.3; % Incongruent
CV_AD((HCP_CV_AD == 2) & (XXY_CV_AD == -1)) = 0.6; % Incongruent

CV_AD((HCP_CV_AD == -2) & (XXY_CV_AD == -2)) = 1.3; % Congruent
CV_AD((HCP_CV_AD == -2) & (XXY_CV_AD == -1)) = 1.6; % Congruent

% Visualize
h = HyoSurfacePlotFunction(CV_AD, IS, G360, cmap.quad_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_xxy_CV.png'), 'Resolution', 300); close(h);

% Now for SA
XXY_SA_AD = zeros(1,360);
HCP_SA_AD = zeros(1,360);
SA_AD = zeros(1,360);

% Annote regions based on direction and significance of X-dosage effects
XXY_SA_AD(XXY_SA_SEX.coef(2,:) > 0) =  1;
XXY_SA_AD(XXY_SA_SEX.coef(2,:) < 0) = -1;
XXY_SA_AD((XXY_SA_SEX_Q < thresh) & (XXY_SA_SEX.coef(2,:) > 0))  =  2;
XXY_SA_AD((XXY_SA_SEX_Q < thresh) & (XXY_SA_SEX.coef(2,:) < 0))  = -2;

% Annote regions based on direction and significance of sex effects
HCP_SA_AD(HCP_SA_SEX.coef(2,:) > 0) =  1;
HCP_SA_AD(HCP_SA_SEX.coef(2,:) < 0) = -1;
HCP_SA_AD((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) > 0))  =  2;
HCP_SA_AD((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) < 0))  = -2;

% Compute regional conjunction
SA_AD((HCP_SA_AD < 2) & (HCP_SA_AD > -2)) = 0;

SA_AD((HCP_SA_AD == 2) & (XXY_SA_AD == 2)) = -1.9;
SA_AD((HCP_SA_AD == 2) & (XXY_SA_AD == 1)) = -1.6;

SA_AD((HCP_SA_AD == -2) & (XXY_SA_AD == 2)) = -0.9;
SA_AD((HCP_SA_AD == -2) & (XXY_SA_AD == 1)) = -0.6;

SA_AD((HCP_SA_AD == 2) & (XXY_SA_AD == -2)) = 0.3;
SA_AD((HCP_SA_AD == 2) & (XXY_SA_AD == -1)) = 0.6;

SA_AD((HCP_SA_AD == -2) & (XXY_SA_AD == -2)) = 1.3;
SA_AD((HCP_SA_AD == -2) & (XXY_SA_AD == -1)) = 1.6;

% Visualize
h = HyoSurfacePlotFunction(SA_AD, IS, G360, cmap.quad_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_xxy_SA.png'), 'Resolution', 300); close(h);

% Finally for CT
XXY_CT_AD = zeros(1,360);
HCP_CT_AD = zeros(1,360);
CT_AD = zeros(1,360);

% Annote regions based on direction and significance of X-dosage effects
XXY_CT_AD(XXY_CT_SEX.coef(2,:) > 0) =  1;
XXY_CT_AD(XXY_CT_SEX.coef(2,:) < 0) = -1;
XXY_CT_AD((XXY_CT_SEX_Q < thresh) & (XXY_CT_SEX.coef(2,:) > 0))  =  2;
XXY_CT_AD((XXY_CT_SEX_Q < thresh) & (XXY_CT_SEX.coef(2,:) < 0))  = -2;

% Annote regions based on direction and significance of sex effects
HCP_CT_AD(HCP_CT_SEX.coef(2,:) > 0) =  1;
HCP_CT_AD(HCP_CT_SEX.coef(2,:) < 0) = -1;
HCP_CT_AD((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) > 0))  =  2;
HCP_CT_AD((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) < 0))  = -2;

% Compute regional conjunction
CT_AD((HCP_CT_AD < 2) & (HCP_CT_AD > -2)) = 0;

CT_AD((HCP_CT_AD == 2) & (XXY_CT_AD == 2)) = -1.9;
CT_AD((HCP_CT_AD == 2) & (XXY_CT_AD == 1)) = -1.6;

CT_AD((HCP_CT_AD == -2) & (XXY_CT_AD == 2)) = -0.9;
CT_AD((HCP_CT_AD == -2) & (XXY_CT_AD == 1)) = -0.6;

CT_AD((HCP_CT_AD == 2) & (XXY_CT_AD == -2)) = 0.3;
CT_AD((HCP_CT_AD == 2) & (XXY_CT_AD == -1)) = 0.6;

CT_AD((HCP_CT_AD == -2) & (XXY_CT_AD == -2)) = 1.3;
CT_AD((HCP_CT_AD == -2) & (XXY_CT_AD == -1)) = 1.6;

% Visualize
h = HyoSurfacePlotFunction(CT_AD, IS, G360, cmap.quad_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_xxy_CT.png'), 'Resolution', 300); close(h);

% Bar graphs
CV_AD_Count = [sum(CV_AD == 1.3) sum(CV_AD == 1.6) 0 ...
               sum(CV_AD == -1.9) sum(CV_AD == -1.6) 0 ...
               sum(CV_AD == -0.9) sum(CV_AD == -0.6) 0 ...
               sum(CV_AD == 0.3) sum(CV_AD == 0.6) 0];
h = figure;
b = bar(1, CV_AD_Count, 'stacked');
for i = 1:length(CV_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_xxy_CV_bar.png'), 'Resolution', 300); close(h);

SA_AD_Count = [sum(SA_AD == 1.3) sum(SA_AD == 1.6) 0 ...
               sum(SA_AD == -1.9) sum(SA_AD == -1.6) 0 ...
               sum(SA_AD == -0.9) sum(SA_AD == -0.6) 0 ...
               sum(SA_AD == 0.3) sum(SA_AD == 0.6) 0];
h = figure;
b = bar(1, SA_AD_Count, 'stacked');
for i = 1:length(SA_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_xxy_SA_bar.png'), 'Resolution', 300); close(h);

CT_AD_Count = [sum(CT_AD == 1.3) sum(CT_AD == 1.6) 0 ...
               sum(CT_AD == -1.9) sum(CT_AD == -1.6) 0 ...
               sum(CT_AD == -0.9) sum(CT_AD == -0.6) 0 ...
               sum(CT_AD == 0.3) sum(CT_AD == 0.6) 0];
h = figure;
b = bar(1, CT_AD_Count, 'stacked');
for i = 1:length(CT_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_xxy_CT_bar.png'), 'Resolution', 300); close(h);

% Calculate number of significantly sex-biased regions for each feature
numSigCV = sum((HCP_CV_AD == 2) | (HCP_CV_AD == -2));
numSigSA = sum((HCP_SA_AD == 2) | (HCP_SA_AD == -2));
numSigCT = sum((HCP_CT_AD == 2) | (HCP_CT_AD == -2));

disp('Proportion of congruent regions')
[(sum(abs(CV_AD) > 1) ./ numSigCV) (sum(abs(SA_AD) > 1) ./ numSigSA) (sum(abs(CT_AD) > 1) ./ numSigCT)]

% --- Spin permutations --- (start)
for i = 1:nPerm
    SPIN_HCP_CV_AD(:,i) = HCP_CV_AD(Vasa(:,i));
    SPIN_HCP_SA_AD(:,i) = HCP_SA_AD(Vasa(:,i));
    SPIN_HCP_CT_AD(:,i) = HCP_CT_AD(Vasa(:,i));
end

% Initialize
SPIN_CV_AD = zeros(1,360*nPerm);
SPIN_SA_AD = zeros(1,360*nPerm);
SPIN_CT_AD = zeros(1,360*nPerm);

% Compute agreement maps
SPIN_CV_AD((repmat(XXY_CV_AD, nPerm, 1)' == 2) & (SPIN_HCP_CV_AD == 2)) = -1.9;
SPIN_CV_AD((repmat(XXY_CV_AD, nPerm, 1)' == 1) & (SPIN_HCP_CV_AD == 2)) = -1.6;
SPIN_CV_AD((repmat(XXY_CV_AD, nPerm, 1)' == 2) & (SPIN_HCP_CV_AD == -2)) = -0.9;
SPIN_CV_AD((repmat(XXY_CV_AD, nPerm, 1)' == 1) & (SPIN_HCP_CV_AD == -2)) = -0.6;
SPIN_CV_AD((repmat(XXY_CV_AD, nPerm, 1)' == -2) & (SPIN_HCP_CV_AD == 2)) = 0.3;
SPIN_CV_AD((repmat(XXY_CV_AD, nPerm, 1)' == -1) & (SPIN_HCP_CV_AD == 2)) = 0.6;
SPIN_CV_AD((repmat(XXY_CV_AD, nPerm, 1)' == -2) & (SPIN_HCP_CV_AD == -2)) = 1.3;
SPIN_CV_AD((repmat(XXY_CV_AD, nPerm, 1)' == -1) & (SPIN_HCP_CV_AD == -2)) = 1.6;

SPIN_SA_AD((repmat(XXY_SA_AD, nPerm, 1)' == 2) & (SPIN_HCP_SA_AD == 2)) = -1.9;
SPIN_SA_AD((repmat(XXY_SA_AD, nPerm, 1)' == 1) & (SPIN_HCP_SA_AD == 2)) = -1.6;
SPIN_SA_AD((repmat(XXY_SA_AD, nPerm, 1)' == 2) & (SPIN_HCP_SA_AD == -2)) = -0.9;
SPIN_SA_AD((repmat(XXY_SA_AD, nPerm, 1)' == 1) & (SPIN_HCP_SA_AD == -2)) = -0.6;
SPIN_SA_AD((repmat(XXY_SA_AD, nPerm, 1)' == -2) & (SPIN_HCP_SA_AD == 2)) = 0.3;
SPIN_SA_AD((repmat(XXY_SA_AD, nPerm, 1)' == -1) & (SPIN_HCP_SA_AD == 2)) = 0.6;
SPIN_SA_AD((repmat(XXY_SA_AD, nPerm, 1)' == -2) & (SPIN_HCP_SA_AD == -2)) = 1.3;
SPIN_SA_AD((repmat(XXY_SA_AD, nPerm, 1)' == -1) & (SPIN_HCP_SA_AD == -2)) = 1.6;

SPIN_CT_AD((repmat(XXY_CT_AD, nPerm, 1)' == 2) & (SPIN_HCP_CT_AD == 2)) = -1.9;
SPIN_CT_AD((repmat(XXY_CT_AD, nPerm, 1)' == 1) & (SPIN_HCP_CT_AD == 2)) = -1.6;
SPIN_CT_AD((repmat(XXY_CT_AD, nPerm, 1)' == 2) & (SPIN_HCP_CT_AD == -2)) = -0.9;
SPIN_CT_AD((repmat(XXY_CT_AD, nPerm, 1)' == 1) & (SPIN_HCP_CT_AD == -2)) = -0.6;
SPIN_CT_AD((repmat(XXY_CT_AD, nPerm, 1)' == -2) & (SPIN_HCP_CT_AD == 2)) = 0.3;
SPIN_CT_AD((repmat(XXY_CT_AD, nPerm, 1)' == -1) & (SPIN_HCP_CT_AD == 2)) = 0.6;
SPIN_CT_AD((repmat(XXY_CT_AD, nPerm, 1)' == -2) & (SPIN_HCP_CT_AD == -2)) = 1.3;
SPIN_CT_AD((repmat(XXY_CT_AD, nPerm, 1)' == -1) & (SPIN_HCP_CT_AD == -2)) = 1.6;

% Reshape
SPIN_CV_AD = reshape(SPIN_CV_AD, [360 nPerm]);
SPIN_SA_AD = reshape(SPIN_SA_AD, [360 nPerm]);
SPIN_CT_AD = reshape(SPIN_CT_AD, [360 nPerm]);

% Calculate proportions
SPIN_CV = sum(abs(SPIN_CV_AD) > 1) ./ numSigCV;
SPIN_SA = sum(abs(SPIN_SA_AD) > 1) ./ numSigSA;
SPIN_CT = sum(abs(SPIN_CT_AD) > 1) ./ numSigCT;

% Calculate p value
CV_SPIN_P = sum((sum(abs(CV_AD) > 1) ./ numSigCV) < SPIN_CV) ./ nPerm;
SA_SPIN_P = sum((sum(abs(SA_AD) > 1) ./ numSigSA) < SPIN_SA) ./ nPerm;
CT_SPIN_P = sum((sum(abs(CT_AD) > 1) ./ numSigCT) < SPIN_CT) ./ nPerm;

% Bonferroni correction
CV_SPIN_Q = 3*CV_SPIN_P;
SA_SPIN_Q = 3*SA_SPIN_P;
CT_SPIN_Q = 3*CT_SPIN_P;

% --- Spin permutations --- (end)

%% Regional congruence between HCP and XYY cohorts

% Conjunction map for CV
XYY_CV_AD = zeros(1,360);
HCP_CV_AD = zeros(1,360);
CV_AD = zeros(1,360);

% Annote regions based on direction and significance of Y-dosage effects
XYY_CV_AD(XYY_CV_SEX.coef(2,:) > 0) =  1;
XYY_CV_AD(XYY_CV_SEX.coef(2,:) < 0) = -1;
XYY_CV_AD((XYY_CV_SEX_Q < thresh) & (XYY_CV_SEX.coef(2,:) > 0))  =  2;
XYY_CV_AD((XYY_CV_SEX_Q < thresh) & (XYY_CV_SEX.coef(2,:) < 0))  = -2;

% Annote regions based on direction and significance of sex effects
HCP_CV_AD(HCP_CV_SEX.coef(2,:) > 0) =  1;
HCP_CV_AD(HCP_CV_SEX.coef(2,:) < 0) = -1;
HCP_CV_AD((HCP_CV_SEX_Q < thresh) & (HCP_CV_SEX.coef(2,:) > 0))  =  2;
HCP_CV_AD((HCP_CV_SEX_Q < thresh) & (HCP_CV_SEX.coef(2,:) < 0))  = -2;

% Compute regional conjunction
CV_AD((HCP_CV_AD < 2) & (HCP_CV_AD > -2)) = 0;

CV_AD((HCP_CV_AD == 2) & (XYY_CV_AD == -2)) = -1.9; % Congruent
CV_AD((HCP_CV_AD == 2) & (XYY_CV_AD == -1)) = -1.6; % Congruent

CV_AD((HCP_CV_AD == -2) & (XYY_CV_AD == -2)) = -0.9; % Incongruent
CV_AD((HCP_CV_AD == -2) & (XYY_CV_AD == -1)) = -0.6; % Incongruent

CV_AD((HCP_CV_AD == 2) & (XYY_CV_AD == 2)) = 0.3; % Incongruent
CV_AD((HCP_CV_AD == 2) & (XYY_CV_AD == 1)) = 0.6; % Incongruent

CV_AD((HCP_CV_AD == -2) & (XYY_CV_AD == 2)) = 1.3; % Congruent
CV_AD((HCP_CV_AD == -2) & (XYY_CV_AD == 1)) = 1.6; % Congruent

% Visualize
h = HyoSurfacePlotFunction(CV_AD, IS, G360, cmap.quad_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_xyy_CV.png'), 'Resolution', 300); close(h);

% Now for SA
XYY_SA_AD = zeros(1,360);
HCP_SA_AD = zeros(1,360);
SA_AD = zeros(1,360);

% Annote regions based on direction and significance of Y-dosage effects
XYY_SA_AD(XYY_SA_SEX.coef(2,:) > 0) =  1;
XYY_SA_AD(XYY_SA_SEX.coef(2,:) < 0) = -1;
XYY_SA_AD((XYY_SA_SEX_Q < thresh) & (XYY_SA_SEX.coef(2,:) > 0))  =  2;
XYY_SA_AD((XYY_SA_SEX_Q < thresh) & (XYY_SA_SEX.coef(2,:) < 0))  = -2;

% Annote regions based on direction and significance of sex effects
HCP_SA_AD(HCP_SA_SEX.coef(2,:) > 0) =  1;
HCP_SA_AD(HCP_SA_SEX.coef(2,:) < 0) = -1;
HCP_SA_AD((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) > 0))  =  2;
HCP_SA_AD((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) < 0))  = -2;

% Compute regional conjunction
SA_AD((HCP_SA_AD < 2) & (HCP_SA_AD > -2)) = 0;

SA_AD((HCP_SA_AD == 2) & (XYY_SA_AD == -2)) = -1.9;
SA_AD((HCP_SA_AD == 2) & (XYY_SA_AD == -1)) = -1.6;

SA_AD((HCP_SA_AD == -2) & (XYY_SA_AD == -2)) = -0.9;
SA_AD((HCP_SA_AD == -2) & (XYY_SA_AD == -1)) = -0.6;

SA_AD((HCP_SA_AD == 2) & (XYY_SA_AD == 2)) = 0.3;
SA_AD((HCP_SA_AD == 2) & (XYY_SA_AD == 1)) = 0.6;

SA_AD((HCP_SA_AD == -2) & (XYY_SA_AD == 2)) = 1.3;
SA_AD((HCP_SA_AD == -2) & (XYY_SA_AD == 1)) = 1.6;

% Visualize
h = HyoSurfacePlotFunction(SA_AD, IS, G360, cmap.quad_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_xyy_SA.png'), 'Resolution', 300); close(h);

% Finally for CT
XYY_CT_AD = zeros(1,360);
HCP_CT_AD = zeros(1,360);
CT_AD = zeros(1,360);

% Annote regions based on direction and significance of Y-dosage effects
XYY_CT_AD(XYY_CT_SEX.coef(2,:) > 0) =  1;
XYY_CT_AD(XYY_CT_SEX.coef(2,:) < 0) = -1;
XYY_CT_AD((XYY_CT_SEX_Q < thresh) & (XYY_CT_SEX.coef(2,:) > 0))  =  2;
XYY_CT_AD((XYY_CT_SEX_Q < thresh) & (XYY_CT_SEX.coef(2,:) < 0))  = -2;

% Annote regions based on direction and significance of sex effects
HCP_CT_AD(HCP_CT_SEX.coef(2,:) > 0) =  1;
HCP_CT_AD(HCP_CT_SEX.coef(2,:) < 0) = -1;
HCP_CT_AD((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) > 0))  =  2;
HCP_CT_AD((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) < 0))  = -2;

% Generate agreement maps
CT_AD((HCP_CT_AD < 2) & (HCP_CT_AD > -2)) = 0;

CT_AD((HCP_CT_AD == 2) & (XYY_CT_AD == -2)) = -1.9;
CT_AD((HCP_CT_AD == 2) & (XYY_CT_AD == -1)) = -1.6;

CT_AD((HCP_CT_AD == -2) & (XYY_CT_AD == -2)) = -0.9;
CT_AD((HCP_CT_AD == -2) & (XYY_CT_AD == -1)) = -0.6;

CT_AD((HCP_CT_AD == 2) & (XYY_CT_AD == 2)) = 0.3;
CT_AD((HCP_CT_AD == 2) & (XYY_CT_AD == 1)) = 0.6;

CT_AD((HCP_CT_AD == -2) & (XYY_CT_AD == 2)) = 1.3;
CT_AD((HCP_CT_AD == -2) & (XYY_CT_AD == 1)) = 1.6;

% Visualize
h = HyoSurfacePlotFunction(CT_AD, IS, G360, cmap.quad_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_xyy_CT.png'), 'Resolution', 300); close(h);

% Bar graphs
CV_AD_Count = [sum(CV_AD == 1.3) sum(CV_AD == 1.6) 0 ...
               sum(CV_AD == -1.9) sum(CV_AD == -1.6) 0 ...
               sum(CV_AD == -0.9) sum(CV_AD == -0.6) 0 ...
               sum(CV_AD == 0.3) sum(CV_AD == 0.6) 0];
h = figure;
b = bar(1, CV_AD_Count, 'stacked');
for i = 1:length(CV_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_xyy_CV_bar.png'), 'Resolution', 300); close(h);

SA_AD_Count = [sum(SA_AD == 1.3) sum(SA_AD == 1.6) 0 ...
               sum(SA_AD == -1.9) sum(SA_AD == -1.6) 0 ...
               sum(SA_AD == -0.9) sum(SA_AD == -0.6) 0 ...
               sum(SA_AD == 0.3) sum(SA_AD == 0.6) 0];
h = figure;
b = bar(1, SA_AD_Count, 'stacked');
for i = 1:length(SA_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_xyy_SA_bar.png'), 'Resolution', 300); close(h);

CT_AD_Count = [sum(CT_AD == 1.3) sum(CT_AD == 1.6) 0 ...
               sum(CT_AD == -1.9) sum(CT_AD == -1.6) 0 ...
               sum(CT_AD == -0.9) sum(CT_AD == -0.6) 0 ...
               sum(CT_AD == 0.3) sum(CT_AD == 0.6) 0];
h = figure;
b = bar(1, CT_AD_Count, 'stacked');
for i = 1:length(CT_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_xyy_CT_bar.png'), 'Resolution', 300); close(h);

% Number of significant ROIs
numSigCV = sum((HCP_CV_AD == 2) | (HCP_CV_AD == -2));
numSigSA = sum((HCP_SA_AD == 2) | (HCP_SA_AD == -2));
numSigCT = sum((HCP_CT_AD == 2) | (HCP_CT_AD == -2));

disp('Proportion of congruent regions')
[(sum(abs(CV_AD) > 1) ./ numSigCV) (sum(abs(SA_AD) > 1) ./ numSigSA) (sum(abs(CT_AD) > 1) ./ numSigCT)]

% --- Spin permutations --- (start)
for i = 1:nPerm
    SPIN_HCP_CV_AD(:,i) = HCP_CV_AD(Vasa(:,i));
    SPIN_HCP_SA_AD(:,i) = HCP_SA_AD(Vasa(:,i));
    SPIN_HCP_CT_AD(:,i) = HCP_CT_AD(Vasa(:,i));
end

SPIN_CV_AD = zeros(1,360*nPerm);
SPIN_SA_AD = zeros(1,360*nPerm);
SPIN_CT_AD = zeros(1,360*nPerm);

SPIN_CV_AD((repmat(XYY_CV_AD, nPerm, 1)' == -2) & (SPIN_HCP_CV_AD == 2)) = -1.9;
SPIN_CV_AD((repmat(XYY_CV_AD, nPerm, 1)' == -1) & (SPIN_HCP_CV_AD == 2)) = -1.6;
SPIN_CV_AD((repmat(XYY_CV_AD, nPerm, 1)' == -2) & (SPIN_HCP_CV_AD == -2)) = -0.9;
SPIN_CV_AD((repmat(XYY_CV_AD, nPerm, 1)' == -1) & (SPIN_HCP_CV_AD == -2)) = -0.6;
SPIN_CV_AD((repmat(XYY_CV_AD, nPerm, 1)' == 2) & (SPIN_HCP_CV_AD == 2)) = 0.3;
SPIN_CV_AD((repmat(XYY_CV_AD, nPerm, 1)' == 1) & (SPIN_HCP_CV_AD == 2)) = 0.6;
SPIN_CV_AD((repmat(XYY_CV_AD, nPerm, 1)' == 2) & (SPIN_HCP_CV_AD == -2)) = 1.3;
SPIN_CV_AD((repmat(XYY_CV_AD, nPerm, 1)' == 1) & (SPIN_HCP_CV_AD == -2)) = 1.6;

SPIN_SA_AD((repmat(XYY_SA_AD, nPerm, 1)' == -2) & (SPIN_HCP_SA_AD == 2)) = -1.9;
SPIN_SA_AD((repmat(XYY_SA_AD, nPerm, 1)' == -1) & (SPIN_HCP_SA_AD == 2)) = -1.6;
SPIN_SA_AD((repmat(XYY_SA_AD, nPerm, 1)' == -2) & (SPIN_HCP_SA_AD == -2)) = -0.9;
SPIN_SA_AD((repmat(XYY_SA_AD, nPerm, 1)' == -1) & (SPIN_HCP_SA_AD == -2)) = -0.6;
SPIN_SA_AD((repmat(XYY_SA_AD, nPerm, 1)' == 2) & (SPIN_HCP_SA_AD == 2)) = 0.3;
SPIN_SA_AD((repmat(XYY_SA_AD, nPerm, 1)' == 1) & (SPIN_HCP_SA_AD == 2)) = 0.6;
SPIN_SA_AD((repmat(XYY_SA_AD, nPerm, 1)' == 2) & (SPIN_HCP_SA_AD == -2)) = 1.3;
SPIN_SA_AD((repmat(XYY_SA_AD, nPerm, 1)' == 1) & (SPIN_HCP_SA_AD == -2)) = 1.6;

SPIN_CT_AD((repmat(XYY_CT_AD, nPerm, 1)' == -2) & (SPIN_HCP_CT_AD == 2)) = -1.9;
SPIN_CT_AD((repmat(XYY_CT_AD, nPerm, 1)' == -1) & (SPIN_HCP_CT_AD == 2)) = -1.6;
SPIN_CT_AD((repmat(XYY_CT_AD, nPerm, 1)' == -2) & (SPIN_HCP_CT_AD == -2)) = -0.9;
SPIN_CT_AD((repmat(XYY_CT_AD, nPerm, 1)' == -1) & (SPIN_HCP_CT_AD == -2)) = -0.6;
SPIN_CT_AD((repmat(XYY_CT_AD, nPerm, 1)' == 2) & (SPIN_HCP_CT_AD == 2)) = 0.3;
SPIN_CT_AD((repmat(XYY_CT_AD, nPerm, 1)' == 1) & (SPIN_HCP_CT_AD == 2)) = 0.6;
SPIN_CT_AD((repmat(XYY_CT_AD, nPerm, 1)' == 2) & (SPIN_HCP_CT_AD == -2)) = 1.3;
SPIN_CT_AD((repmat(XYY_CT_AD, nPerm, 1)' == 1) & (SPIN_HCP_CT_AD == -2)) = 1.6;

% Reshape
SPIN_CV_AD = reshape(SPIN_CV_AD, [360 nPerm]);
SPIN_SA_AD = reshape(SPIN_SA_AD, [360 nPerm]);
SPIN_CT_AD = reshape(SPIN_CT_AD, [360 nPerm]);

% Calculate proportions
SPIN_CV = sum(abs(SPIN_CV_AD) > 1) ./ numSigCV;
SPIN_SA = sum(abs(SPIN_SA_AD) > 1) ./ numSigSA;
SPIN_CT = sum(abs(SPIN_CT_AD) > 1) ./ numSigCT;

% Calculate p value
CV_SPIN_P = sum((sum(abs(CV_AD) > 1) ./ numSigCV) < SPIN_CV) ./ nPerm;
SA_SPIN_P = sum((sum(abs(SA_AD) > 1) ./ numSigSA) < SPIN_SA) ./ nPerm;
CT_SPIN_P = sum((sum(abs(CT_AD) > 1) ./ numSigCT) < SPIN_CT) ./ nPerm;

% Bonferroni correction
CV_SPIN_Q = 3*CV_SPIN_P;
SA_SPIN_Q = 3*SA_SPIN_P;
CT_SPIN_Q = 3*CT_SPIN_P;

% --- Spin permutations --- (end)
