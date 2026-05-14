
clear all; close all; clc;

%% Prepare

% Load data
load('matfiles/IGD_Processed.mat');

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
figDir = 'igd';

%% Statistics for group maps

% Compute z-score
IGD_CV = zscore(IGD.CV);
IGD_SA = zscore(IGD.SA);
IGD_CT = zscore(IGD.CT);

% Compute global phenotyeps
IGD.CV_GS = sum(IGD.CV,2);
IGD.SA_GS = sum(IGD.SA,2);
IGD.CT_GS = sum(IGD.CT .* (IGD.SA ./ sum(IGD.SA,2)),2);

% Prepare terms (mean-centered)
IGD_GRP = term(IGD.GRP); % The contrast is healthy males vs IGD males. Positive direction indicates effect of intact testicular hormone production.
IGD_SEX = term(IGD.SEX);
IGD_AGE = term(IGD.AGE - mean(IGD.AGE));
IGD_EUL = term(IGD.EUL' - mean(IGD.EUL));
IGD_G_CV = term(IGD.CV_GS - mean(IGD.CV_GS));
IGD_G_SA = term(IGD.SA_GS - mean(IGD.SA_GS));
IGD_G_CT = term(IGD.CT_GS - mean(IGD.CT_GS));

% Define GLM
IGD_M_CV = 1 + IGD_GRP + IGD_AGE + IGD_EUL + IGD_G_CV;
IGD_M_SA = 1 + IGD_GRP + IGD_AGE + IGD_EUL + IGD_G_SA;
IGD_M_CT = 1 + IGD_GRP + IGD_AGE + IGD_EUL + IGD_G_CT;

% Compute stats
IGD_CV_GRP = SurfStatT(SurfStatLinMod(IGD_CV, IGD_M_CV), IGD_GRP);
IGD_CV_GRP_P = 2*(1 - tcdf(abs(IGD_CV_GRP.t), IGD_CV_GRP.df));
[~,~,~,IGD_CV_GRP_Q] = fdr_bh(IGD_CV_GRP_P);

IGD_SA_GRP = SurfStatT(SurfStatLinMod(IGD_SA, IGD_M_SA), IGD_GRP);
IGD_SA_GRP_P = 2*(1 - tcdf(abs(IGD_SA_GRP.t), IGD_SA_GRP.df));
[~,~,~,IGD_SA_GRP_Q] = fdr_bh(IGD_SA_GRP_P);

IGD_CT_GRP = SurfStatT(SurfStatLinMod(IGD_CT, IGD_M_CT), IGD_GRP);
IGD_CT_GRP_P = 2*(1 - tcdf(abs(IGD_CT_GRP.t), IGD_CT_GRP.df));
[~,~,~,IGD_CT_GRP_Q] = fdr_bh(IGD_CT_GRP_P);

% Save group maps
save('igd_results.mat', 'IGD_CV_GRP', 'IGD_CV_GRP_Q', 'IGD_CV_GRP_P', ...
                        'IGD_SA_GRP', 'IGD_SA_GRP_Q', 'IGD_SA_GRP_P', ...
                        'IGD_CT_GRP', 'IGD_CT_GRP_Q', 'IGD_CT_GRP_P', ...
                        'IGD', '-v7.3');

%% Group maps

% Map of purberty effect (positive: control-biased, negative: IGD-biased)
h = HyoSurfacePlotFunction(IGD_CV_GRP.coef(2,:) .* (IGD_CV_GRP_Q < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'igd_cv_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(IGD_SA_GRP.coef(2,:) .* (IGD_SA_GRP_Q < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'igd_sa_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(IGD_CT_GRP.coef(2,:) .* (IGD_CT_GRP_Q < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'igd_ct_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(IGD_CV_GRP.coef(2,:) .* (IGD_CV_GRP_P < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'igd_cv_nominal.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(IGD_SA_GRP.coef(2,:) .* (IGD_SA_GRP_P < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'igd_sa_nominal.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(IGD_CT_GRP.coef(2,:) .* (IGD_CT_GRP_P < thresh), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'igd_ct_nominal.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(IGD_CV_GRP.coef(2,:), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'igd_cv.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(IGD_SA_GRP.coef(2,:), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'igd_sa.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(IGD_CT_GRP.coef(2,:), IS, G360, cmap.PuOr, '', [-1 1], 'square');
exportgraphics(h, fullfile(figDir, 'igd_ct.png'), 'Resolution', 300); close(h);

%% Regional congruence between HCP and IGD cohorts

% Compute regional conjunction for CV
IGD_CV_AD = zeros(1,360);
HCP_CV_AD = zeros(1,360);
CV_AD = zeros(1,360);

% Annote regions based on direction and significance of testicular hormone production effects
IGD_CV_AD(IGD_CV_GRP.coef(2,:) > 0) =  1;
IGD_CV_AD(IGD_CV_GRP.coef(2,:) < 0) = -1;
IGD_CV_AD((IGD_CV_GRP_P < thresh) & (IGD_CV_GRP.coef(2,:) > 0))  =  2;
IGD_CV_AD((IGD_CV_GRP_P < thresh) & (IGD_CV_GRP.coef(2,:) < 0))  = -2;

% Annote regions based on direction and significance of sex effects
HCP_CV_AD(HCP_CV_SEX.coef(2,:) > 0) =  1;
HCP_CV_AD(HCP_CV_SEX.coef(2,:) < 0) = -1;
HCP_CV_AD((HCP_CV_SEX_Q < thresh) & (HCP_CV_SEX.coef(2,:) > 0))  =  2;
HCP_CV_AD((HCP_CV_SEX_Q < thresh) & (HCP_CV_SEX.coef(2,:) < 0))  = -2;

% Compute regional conjunction
CV_AD((HCP_CV_AD < 2) & (HCP_CV_AD > -2)) = 0;

CV_AD((HCP_CV_AD == 2) & (IGD_CV_AD == -2)) = -1.9; % Congruent
CV_AD((HCP_CV_AD == 2) & (IGD_CV_AD == -1)) = -1.6; % Congruent

CV_AD((HCP_CV_AD == -2) & (IGD_CV_AD == -2)) = -0.9; % Incongruent
CV_AD((HCP_CV_AD == -2) & (IGD_CV_AD == -1)) = -0.6; % Incongruent

CV_AD((HCP_CV_AD == 2) & (IGD_CV_AD == 2)) = 0.3; % Incongruent
CV_AD((HCP_CV_AD == 2) & (IGD_CV_AD == 1)) = 0.6; % Incongruent

CV_AD((HCP_CV_AD == -2) & (IGD_CV_AD == 2)) = 1.3; % Congruent
CV_AD((HCP_CV_AD == -2) & (IGD_CV_AD == 1)) = 1.6; % Congruent

% Visualize
h = HyoSurfacePlotFunction(CV_AD, IS, G360, cmap.quad_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_IGD_CV.png'), 'Resolution', 300); close(h);

% Now for SA
IGD_SA_AD = zeros(1,360);
HCP_SA_AD = zeros(1,360);
SA_AD = zeros(1,360);

% Annote regions based on direction and significance of testicular hormone production effects
IGD_SA_AD(IGD_SA_GRP.coef(2,:) > 0) =  1;
IGD_SA_AD(IGD_SA_GRP.coef(2,:) < 0) = -1;
IGD_SA_AD((IGD_SA_GRP_P < thresh) & (IGD_SA_GRP.coef(2,:) > 0))  =  2;
IGD_SA_AD((IGD_SA_GRP_P < thresh) & (IGD_SA_GRP.coef(2,:) < 0))  = -2;

% Annote regions based on direction and significance of sex effects
HCP_SA_AD(HCP_SA_SEX.coef(2,:) > 0) =  1;
HCP_SA_AD(HCP_SA_SEX.coef(2,:) < 0) = -1;
HCP_SA_AD((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) > 0))  =  2;
HCP_SA_AD((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) < 0))  = -2;

% Compute regional conjunction
SA_AD((HCP_SA_AD < 2) & (HCP_SA_AD > -2)) = 0;

SA_AD((HCP_SA_AD == 2) & (IGD_SA_AD == -2)) = -1.9; % Congruent
SA_AD((HCP_SA_AD == 2) & (IGD_SA_AD == -1)) = -1.6; % Congruent

SA_AD((HCP_SA_AD == -2) & (IGD_SA_AD == -2)) = -0.9; % Incongruent
SA_AD((HCP_SA_AD == -2) & (IGD_SA_AD == -1)) = -0.6; % Incongruent

SA_AD((HCP_SA_AD == 2) & (IGD_SA_AD == 2)) = 0.3; % Incongruent
SA_AD((HCP_SA_AD == 2) & (IGD_SA_AD == 1)) = 0.6; % Incongruent

SA_AD((HCP_SA_AD == -2) & (IGD_SA_AD == 2)) = 1.3; % Congruent
SA_AD((HCP_SA_AD == -2) & (IGD_SA_AD == 1)) = 1.6; % Congruent

% Visualize
h = HyoSurfacePlotFunction(SA_AD, IS, G360, cmap.quad_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_IGD_SA.png'), 'Resolution', 300); close(h);

% Finally for CT
IGD_CT_AD = zeros(1,360);
HCP_CT_AD = zeros(1,360);
CT_AD = zeros(1,360);

% Annote regions based on direction and significance of testicular hormone production effects
IGD_CT_AD(IGD_CT_GRP.coef(2,:) > 0) =  1;
IGD_CT_AD(IGD_CT_GRP.coef(2,:) < 0) = -1;
IGD_CT_AD((IGD_CT_GRP_P < thresh) & (IGD_CT_GRP.coef(2,:) > 0))  =  2;
IGD_CT_AD((IGD_CT_GRP_P < thresh) & (IGD_CT_GRP.coef(2,:) < 0))  = -2;

% Annote regions based on direction and significance of sex effects
HCP_CT_AD(HCP_CT_SEX.coef(2,:) > 0) =  1;
HCP_CT_AD(HCP_CT_SEX.coef(2,:) < 0) = -1;
HCP_CT_AD((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) > 0))  =  2;
HCP_CT_AD((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) < 0))  = -2;

% Generate agreement maps
CT_AD((HCP_CT_AD < 2) & (HCP_CT_AD > -2)) = 0;

CT_AD((HCP_CT_AD == 2) & (IGD_CT_AD == -2)) = -1.9; % Congruent
CT_AD((HCP_CT_AD == 2) & (IGD_CT_AD == -1)) = -1.6; % Congruent

CT_AD((HCP_CT_AD == -2) & (IGD_CT_AD == -2)) = -0.9; % Incongruent
CT_AD((HCP_CT_AD == -2) & (IGD_CT_AD == -1)) = -0.6; % Incongruent

CT_AD((HCP_CT_AD == 2) & (IGD_CT_AD == 2)) = 0.3; % Incongruent
CT_AD((HCP_CT_AD == 2) & (IGD_CT_AD == 1)) = 0.6; % Incongruent

CT_AD((HCP_CT_AD == -2) & (IGD_CT_AD == 2)) = 1.3; % Congruent
CT_AD((HCP_CT_AD == -2) & (IGD_CT_AD == 1)) = 1.6; % Congruent

% Visualize
h = HyoSurfacePlotFunction(CT_AD, IS, G360, cmap.quad_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_IGD_CT.png'), 'Resolution', 300); close(h);

% Number of significant ROIs
numSigCV = sum((HCP_CV_AD == 2) | (HCP_CV_AD == -2));
numSigSA = sum((HCP_SA_AD == 2) | (HCP_SA_AD == -2));
numSigCT = sum((HCP_CT_AD == 2) | (HCP_CT_AD == -2));

% CV
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
exportgraphics(h, fullfile(figDir, 'conjunction_IGD_CV_bar.png'), 'Resolution', 300); close(h);

% SA
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
exportgraphics(h, fullfile(figDir, 'conjunction_IGD_SA_bar.png'), 'Resolution', 300); close(h);

% CT
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
exportgraphics(h, fullfile(figDir, 'conjunction_IGD_CT_bar.png'), 'Resolution', 300); close(h);

% Display proportion of expected directions
[(sum(abs(CV_AD) > 1) ./ numSigCV) (sum(abs(SA_AD) > 1) ./ numSigSA) (sum(abs(CT_AD) > 1) ./ numSigCT)]

% --- Spin permutations --- (start)
for i = 1:nPerm
    SPIN_HCP_CV_AD(:,i) = HCP_CV_AD(Vasa(:,i));
    SPIN_HCP_SA_AD(:,i) = HCP_SA_AD(Vasa(:,i));
    SPIN_HCP_CT_AD(:,i) = HCP_CT_AD(Vasa(:,i));
end

% Number of significant ROIs
numSigCV = sum((HCP_CV_AD == 2) | (HCP_CV_AD == -2));
numSigSA = sum((HCP_SA_AD == 2) | (HCP_SA_AD == -2));
numSigCT = sum((HCP_CT_AD == 2) | (HCP_CT_AD == -2));

% Initialize
SPIN_CV_AD = zeros(1,360*nPerm);
SPIN_SA_AD = zeros(1,360*nPerm);
SPIN_CT_AD = zeros(1,360*nPerm);

% Compute agreement maps
SPIN_CV_AD((repmat(IGD_CV_AD, nPerm, 1)' == -2) & (SPIN_HCP_CV_AD == 2)) = -1.9; % Congruent
SPIN_CV_AD((repmat(IGD_CV_AD, nPerm, 1)' == -1) & (SPIN_HCP_CV_AD == 2)) = -1.6; % Congruent
SPIN_CV_AD((repmat(IGD_CV_AD, nPerm, 1)' == -2) & (SPIN_HCP_CV_AD == -2)) = -0.9; % Incongruent
SPIN_CV_AD((repmat(IGD_CV_AD, nPerm, 1)' == -1) & (SPIN_HCP_CV_AD == -2)) = -0.6; % Incongruent
SPIN_CV_AD((repmat(IGD_CV_AD, nPerm, 1)' == 2) & (SPIN_HCP_CV_AD == 2)) = 0.3; % Incongruent
SPIN_CV_AD((repmat(IGD_CV_AD, nPerm, 1)' == 1) & (SPIN_HCP_CV_AD == 2)) = 0.6; % Incongruent
SPIN_CV_AD((repmat(IGD_CV_AD, nPerm, 1)' == 2) & (SPIN_HCP_CV_AD == -2)) = 1.3; % Congruent
SPIN_CV_AD((repmat(IGD_CV_AD, nPerm, 1)' == 1) & (SPIN_HCP_CV_AD == -2)) = 1.6; % Congruent

SPIN_SA_AD((repmat(IGD_SA_AD, nPerm, 1)' == -2) & (SPIN_HCP_SA_AD == 2)) = -1.9; % Congruent
SPIN_SA_AD((repmat(IGD_SA_AD, nPerm, 1)' == -1) & (SPIN_HCP_SA_AD == 2)) = -1.6; % Congruent
SPIN_SA_AD((repmat(IGD_SA_AD, nPerm, 1)' == -2) & (SPIN_HCP_SA_AD == -2)) = -0.9; % Incongruent
SPIN_SA_AD((repmat(IGD_SA_AD, nPerm, 1)' == -1) & (SPIN_HCP_SA_AD == -2)) = -0.6; % Incongruent
SPIN_SA_AD((repmat(IGD_SA_AD, nPerm, 1)' == 2) & (SPIN_HCP_SA_AD == 2)) = 0.3; % Incongruent
SPIN_SA_AD((repmat(IGD_SA_AD, nPerm, 1)' == 1) & (SPIN_HCP_SA_AD == 2)) = 0.6; % Incongruent
SPIN_SA_AD((repmat(IGD_SA_AD, nPerm, 1)' == 2) & (SPIN_HCP_SA_AD == -2)) = 1.3; % Congruent
SPIN_SA_AD((repmat(IGD_SA_AD, nPerm, 1)' == 1) & (SPIN_HCP_SA_AD == -2)) = 1.6; % Congruent

SPIN_CT_AD((repmat(IGD_CT_AD, nPerm, 1)' == -2) & (SPIN_HCP_CT_AD == 2)) = -1.9; % Congruent
SPIN_CT_AD((repmat(IGD_CT_AD, nPerm, 1)' == -1) & (SPIN_HCP_CT_AD == 2)) = -1.6; % Congruent
SPIN_CT_AD((repmat(IGD_CT_AD, nPerm, 1)' == -2) & (SPIN_HCP_CT_AD == -2)) = -0.9; % Incongruent
SPIN_CT_AD((repmat(IGD_CT_AD, nPerm, 1)' == -1) & (SPIN_HCP_CT_AD == -2)) = -0.6; % Incongruent
SPIN_CT_AD((repmat(IGD_CT_AD, nPerm, 1)' == 2) & (SPIN_HCP_CT_AD == 2)) = 0.3; % Incongruent
SPIN_CT_AD((repmat(IGD_CT_AD, nPerm, 1)' == 1) & (SPIN_HCP_CT_AD == 2)) = 0.6; % Incongruent
SPIN_CT_AD((repmat(IGD_CT_AD, nPerm, 1)' == 2) & (SPIN_HCP_CT_AD == -2)) = 1.3; % Congruent
SPIN_CT_AD((repmat(IGD_CT_AD, nPerm, 1)' == 1) & (SPIN_HCP_CT_AD == -2)) = 1.6; % Congruent

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