
clear all; close all; clc;

%% Prep

% Load data
load('matfiles/UKB_Processed.mat');
load('matfiles/HCP_Processed.mat');

% Load template
load('fsaverage_7.1.0.mat');
load('glasser_fsaverage_7.1.0.mat');
load('Vasa_G360.mat');

% Age selection (for 45 - 50 years)
indAGE = UKB.AGE < 50;

% Misc
loadColorMaps;
nPerm = 10000;
thresh = 0.05;
figDir = 'ukb';

%% Statistical testing for sex differences in cortical anatomy

% Compute global phenotypes
UKB.CV_GS = sum(UKB.CV(indAGE,:),2);
UKB.SA_GS = sum(UKB.SA(indAGE,:),2);
UKB.CT_GS = sum(UKB.CT(indAGE,:) .* (UKB.SA(indAGE,:) ./ sum(UKB.SA(indAGE,:),2)),2);

% z-score
CV = zscore(UKB.CV(indAGE,:));
SA = zscore(UKB.SA(indAGE,:));
CT = zscore(UKB.CT(indAGE,:));

% Prepare terms
SEX = term(UKB.SEX(indAGE));
AGE = term(UKB.AGE(indAGE) - mean(UKB.AGE(indAGE)));
EUL = term(UKB.EUL(indAGE) - mean(UKB.EUL(indAGE)));
SITE_1 = term(UKB.SITE_3_CAT_1(indAGE));
SITE_2 = term(UKB.SITE_3_CAT_2(indAGE));
SITE_3 = term(UKB.SITE_3_CAT_3(indAGE));
G_CV = term(UKB.CV_GS - mean(UKB.CV_GS));
G_SA = term(UKB.SA_GS - mean(UKB.SA_GS));
G_CT = term(UKB.CT_GS - mean(UKB.CT_GS));

% Define GLM
M_CV = 1 + SEX + AGE + EUL + SITE_1 + SITE_2 + SITE_3 + G_CV;
M_SA = 1 + SEX + AGE + EUL + SITE_1 + SITE_2 + SITE_3 + G_SA;
M_CT = 1 + SEX + AGE + EUL + SITE_1 + SITE_2 + SITE_3 + G_CT;

% Compute stats
UKB_CV = SurfStatLinMod(CV, M_CV);
UKB_CV_SEX = SurfStatT(UKB_CV, SEX);
UKB_CV_SEX_P = 2*(1 - tcdf(abs(UKB_CV_SEX.t), UKB_CV_SEX.df));
[~,~,~,UKB_CV_SEX_Q] = fdr_bh(UKB_CV_SEX_P);

UKB_SA = SurfStatLinMod(SA, M_SA);
UKB_SA_SEX = SurfStatT(UKB_SA, SEX);
UKB_SA_SEX_P = 2*(1 - tcdf(abs(UKB_SA_SEX.t), UKB_SA_SEX.df));
[~,~,~,UKB_SA_SEX_Q] = fdr_bh(UKB_SA_SEX_P);

UKB_CT = SurfStatLinMod(CT, M_CT);
UKB_CT_SEX = SurfStatT(UKB_CT, SEX);
UKB_CT_SEX_P = 2*(1 - tcdf(abs(UKB_CT_SEX.t), UKB_CT_SEX.df));
[~,~,~,UKB_CT_SEX_Q] = fdr_bh(UKB_CT_SEX_P);

% Save group maps
save('matfiles/ukb_results.mat', 'UKB_CV_SEX', 'UKB_CV_SEX_Q', 'UKB_SA_SEX', 'UKB_SA_SEX_Q', 'UKB_CT_SEX', 'UKB_CT_SEX_Q', '-v7.3');

%% Visualize regional sex difference maps

% Standardized effect size maps
h = HyoSurfacePlotFunction(UKB_CV_SEX.coef(2,:) .* (UKB_CV_SEX_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'ukb_cv_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(UKB_SA_SEX.coef(2,:) .* (UKB_SA_SEX_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'ukb_sa_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(UKB_CT_SEX.coef(2,:) .* (UKB_CT_SEX_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'ukb_ct_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(UKB_CV_SEX.coef(2,:), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'ukb_cv.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(UKB_SA_SEX.coef(2,:), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'ukb_sa.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(UKB_CT_SEX.coef(2,:), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'ukb_ct.png'), 'Resolution', 300); close(h);

%% Scatter plots showing correlations of sex differences between HCP and UKB

% Calculate correlations
[r_cv ~] = corr(HCP_CV_SEX.coef(2,:)', UKB_CV_SEX.coef(2,:)');
[r_sa ~] = corr(HCP_SA_SEX.coef(2,:)', UKB_SA_SEX.coef(2,:)');
[r_ct ~] = corr(HCP_CT_SEX.coef(2,:)', UKB_CT_SEX.coef(2,:)');

% Vectors
CMP_CV = [HCP_CV_SEX.coef(2,:); UKB_CV_SEX.coef(2,:)];
CMP_SA = [HCP_SA_SEX.coef(2,:); UKB_SA_SEX.coef(2,:)];
CMP_CT = [HCP_CT_SEX.coef(2,:); UKB_CT_SEX.coef(2,:)];
IDL = [-1:1e-5:1; -1:1e-5:1];

% Calculate residuals from the identity line
for i = 1:size(CMP_CV,2)
    ID_CMP_CV(i) = min(sqrt(sum((CMP_CV(:,i) - IDL).^2)));
    ID_CMP_SA(i) = min(sqrt(sum((CMP_SA(:,i) - IDL).^2)));
    ID_CMP_CT(i) = min(sqrt(sum((CMP_CT(:,i) - IDL).^2)));
end

% Generate scatter plots
h = figure('WindowState', 'maximized'); hold on;
scatter(HCP_CV_SEX.coef(2,:), UKB_CV_SEX.coef(2,:), 200, [0 0.5 0.5], 'o', 'filled', 'MarkerEdgeColor', [0 0 0]); xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out'); b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--'; l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xticklabels({}); yticklabels({});
axis square;
exportgraphics(h, fullfile(figDir, 'scatter_cv.eps'), 'ContentType', 'vector');

h = figure('WindowState', 'maximized'); hold on;
scatter(HCP_SA_SEX.coef(2,:), UKB_SA_SEX.coef(2,:), 200, [0 0.5 0.5], 'o', 'filled', 'MarkerEdgeColor', [0 0 0]); xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out'); b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--'; l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xticklabels({}); yticklabels({});
axis square;
exportgraphics(h, fullfile(figDir, 'scatter_sa.eps'), 'ContentType', 'vector');

h = figure('WindowState', 'maximized'); hold on;
scatter(HCP_CT_SEX.coef(2,:), UKB_CT_SEX.coef(2,:), 200, [0 0.5 0.5], 'o', 'filled', 'MarkerEdgeColor', [0 0 0]); xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out'); b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--'; l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xticklabels({}); yticklabels({});
axis square;
exportgraphics(h, fullfile(figDir, 'scatter_ct.eps'), 'ContentType', 'vector');

close all;

% Perform spin permutation tests
p_spin_cv = doPerms(HCP_CV_SEX_Spin, UKB_CV_SEX.coef(2,:), r_cv, 'top');
p_spin_sa = doPerms(HCP_SA_SEX_Spin, UKB_SA_SEX.coef(2,:), r_sa, 'top');
p_spin_ct = doPerms(HCP_CT_SEX_Spin, UKB_CT_SEX.coef(2,:), r_ct, 'top');

% Perform sex permutation tests
p_perm_cv = doPerms(HCP_CV_SEX_Perm, UKB_CV_SEX.coef(2,:), r_cv, 'top');
p_perm_sa = doPerms(HCP_SA_SEX_Perm, UKB_SA_SEX.coef(2,:), r_sa, 'top');
p_perm_ct = doPerms(HCP_CT_SEX_Perm, UKB_CT_SEX.coef(2,:), r_ct, 'top');

%% Agreement/Disagreement between HCP and UKB sex differences

% Compute a conjunction map for HCP vs. UKB sex differences in CV
UKB_CV_AD = zeros(1,360);
HCP_CV_AD = zeros(1,360);
CV_AD = zeros(1,360);

% Annotate cortical region for direction and significance of sex effcts
UKB_CV_AD(UKB_CV_SEX.coef(2,:) > 0) =  1;
UKB_CV_AD(UKB_CV_SEX.coef(2,:) < 0) = -1;
UKB_CV_AD((UKB_CV_SEX_Q < thresh) & (UKB_CV_SEX.coef(2,:) > 0))  =  2;
UKB_CV_AD((UKB_CV_SEX_Q < thresh) & (UKB_CV_SEX.coef(2,:) < 0))  = -2;

HCP_CV_AD(HCP_CV_SEX.coef(2,:) > 0) =  1;
HCP_CV_AD(HCP_CV_SEX.coef(2,:) < 0) = -1;
HCP_CV_AD((HCP_CV_SEX_Q < thresh) & (HCP_CV_SEX.coef(2,:) > 0))  =  2;
HCP_CV_AD((HCP_CV_SEX_Q < thresh) & (HCP_CV_SEX.coef(2,:) < 0))  = -2;

% Compute regional conjucntion
CV_AD((UKB_CV_AD == 2) & (HCP_CV_AD == 2)) = -1.9;
CV_AD((UKB_CV_AD == 2) & (HCP_CV_AD == 1)) = -1.5;
CV_AD((UKB_CV_AD == 1) & (HCP_CV_AD == 2)) = -1.5;
CV_AD((UKB_CV_AD == 1) & (HCP_CV_AD == 1)) = -1.1;

CV_AD((UKB_CV_AD == 2) & (HCP_CV_AD == -2)) = -0.9;
CV_AD((UKB_CV_AD == 2) & (HCP_CV_AD == -1)) = -0.5;
CV_AD((UKB_CV_AD == 1) & (HCP_CV_AD == -2)) = -0.5;
CV_AD((UKB_CV_AD == 1) & (HCP_CV_AD == -1)) = -0.1;

CV_AD((UKB_CV_AD == -2) & (HCP_CV_AD == 2)) = 0.1;
CV_AD((UKB_CV_AD == -2) & (HCP_CV_AD == 1)) = 0.5;
CV_AD((UKB_CV_AD == -1) & (HCP_CV_AD == 2)) = 0.5;
CV_AD((UKB_CV_AD == -1) & (HCP_CV_AD == 1)) = 0.9;

CV_AD((UKB_CV_AD == -2) & (HCP_CV_AD == -2)) = 1.1;
CV_AD((UKB_CV_AD == -2) & (HCP_CV_AD == -1)) = 1.5;
CV_AD((UKB_CV_AD == -1) & (HCP_CV_AD == -2)) = 1.5;
CV_AD((UKB_CV_AD == -1) & (HCP_CV_AD == -1)) = 1.9;

% Visualize
h = HyoSurfacePlotFunction(CV_AD, IS, G360, cmap.tri_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_cv.png'), 'Resolution', 300); close(h);

% Now for SA
UKB_SA_AD = zeros(1,360);
HCP_SA_AD = zeros(1,360);
SA_AD = zeros(1,360);

UKB_SA_AD(UKB_SA_SEX.coef(2,:) > 0) =  1;
UKB_SA_AD(UKB_SA_SEX.coef(2,:) < 0) = -1;
UKB_SA_AD((UKB_SA_SEX_Q < thresh) & (UKB_SA_SEX.coef(2,:) > 0))  =  2;
UKB_SA_AD((UKB_SA_SEX_Q < thresh) & (UKB_SA_SEX.coef(2,:) < 0))  = -2;

HCP_SA_AD(HCP_SA_SEX.coef(2,:) > 0) =  1;
HCP_SA_AD(HCP_SA_SEX.coef(2,:) < 0) = -1;
HCP_SA_AD((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) > 0))  =  2;
HCP_SA_AD((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) < 0))  = -2;

% Compute regional conjucntion
SA_AD((UKB_SA_AD == 2) & (HCP_SA_AD == 2)) = -1.9;
SA_AD((UKB_SA_AD == 2) & (HCP_SA_AD == 1)) = -1.5;
SA_AD((UKB_SA_AD == 1) & (HCP_SA_AD == 2)) = -1.5;
SA_AD((UKB_SA_AD == 1) & (HCP_SA_AD == 1)) = -1.1;

SA_AD((UKB_SA_AD == 2) & (HCP_SA_AD == -2)) = -0.9;
SA_AD((UKB_SA_AD == 2) & (HCP_SA_AD == -1)) = -0.5;
SA_AD((UKB_SA_AD == 1) & (HCP_SA_AD == -2)) = -0.5;
SA_AD((UKB_SA_AD == 1) & (HCP_SA_AD == -1)) = -0.1;

SA_AD((UKB_SA_AD == -2) & (HCP_SA_AD == 2)) = 0.1;
SA_AD((UKB_SA_AD == -2) & (HCP_SA_AD == 1)) = 0.5;
SA_AD((UKB_SA_AD == -1) & (HCP_SA_AD == 2)) = 0.5;
SA_AD((UKB_SA_AD == -1) & (HCP_SA_AD == 1)) = 0.9;

SA_AD((UKB_SA_AD == -2) & (HCP_SA_AD == -2)) = 1.1;
SA_AD((UKB_SA_AD == -2) & (HCP_SA_AD == -1)) = 1.5;
SA_AD((UKB_SA_AD == -1) & (HCP_SA_AD == -2)) = 1.5;
SA_AD((UKB_SA_AD == -1) & (HCP_SA_AD == -1)) = 1.9;

% Visualize
h = HyoSurfacePlotFunction(SA_AD, IS, G360, cmap.tri_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_sa.png'), 'Resolution', 300); close(h);

% Finally for CT
UKB_CT_AD = zeros(1,360);
HCP_CT_AD = zeros(1,360);
CT_AD = zeros(1,360);

UKB_CT_AD(UKB_CT_SEX.coef(2,:) > 0) =  1;
UKB_CT_AD(UKB_CT_SEX.coef(2,:) < 0) = -1;
UKB_CT_AD((UKB_CT_SEX_Q < thresh) & (UKB_CT_SEX.coef(2,:) > 0))  =  2;
UKB_CT_AD((UKB_CT_SEX_Q < thresh) & (UKB_CT_SEX.coef(2,:) < 0))  = -2;

HCP_CT_AD(HCP_CT_SEX.coef(2,:) > 0) =  1;
HCP_CT_AD(HCP_CT_SEX.coef(2,:) < 0) = -1;
HCP_CT_AD((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) > 0))  =  2;
HCP_CT_AD((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) < 0))  = -2;

% Compute regional conjucntion
CT_AD((UKB_CT_AD == 2) & (HCP_CT_AD == 2)) = -1.9;
CT_AD((UKB_CT_AD == 2) & (HCP_CT_AD == 1)) = -1.5;
CT_AD((UKB_CT_AD == 1) & (HCP_CT_AD == 2)) = -1.5;
CT_AD((UKB_CT_AD == 1) & (HCP_CT_AD == 1)) = -1.1;

CT_AD((UKB_CT_AD == 2) & (HCP_CT_AD == -2)) = -0.9;
CT_AD((UKB_CT_AD == 2) & (HCP_CT_AD == -1)) = -0.5;
CT_AD((UKB_CT_AD == 1) & (HCP_CT_AD == -2)) = -0.5;
CT_AD((UKB_CT_AD == 1) & (HCP_CT_AD == -1)) = -0.1;

CT_AD((UKB_CT_AD == -2) & (HCP_CT_AD == 2)) = 0.1;
CT_AD((UKB_CT_AD == -2) & (HCP_CT_AD == 1)) = 0.5;
CT_AD((UKB_CT_AD == -1) & (HCP_CT_AD == 2)) = 0.5;
CT_AD((UKB_CT_AD == -1) & (HCP_CT_AD == 1)) = 0.9;

CT_AD((UKB_CT_AD == -2) & (HCP_CT_AD == -2)) = 1.1;
CT_AD((UKB_CT_AD == -2) & (HCP_CT_AD == -1)) = 1.5;
CT_AD((UKB_CT_AD == -1) & (HCP_CT_AD == -2)) = 1.5;
CT_AD((UKB_CT_AD == -1) & (HCP_CT_AD == -1)) = 1.9;

% Visualize
h = HyoSurfacePlotFunction(CT_AD, IS, G360, cmap.tri_div, '', [-2 2], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction_ct.png'), 'Resolution', 300); close(h);

%% Bar graphs for agreements/disagreements

CV_AD_Count = [sum(CV_AD == 1.1) sum(CV_AD == 1.5) sum(CV_AD == 1.9) ...
               sum(CV_AD == -1.9) sum(CV_AD == -1.5) sum(CV_AD == -1.1) ...
               sum((CV_AD == -0.9) | (CV_AD == 0.1)) sum((CV_AD == -0.5) | (CV_AD == 0.5)) sum((CV_AD == -0.1) | (CV_AD == 0.9))];
h = figure;
b = bar(1, CV_AD_Count, 'stacked');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_bar_cv.png'), 'Resolution', 300); close(h);

SA_AD_Count = [sum(SA_AD == 1.1) sum(SA_AD == 1.5) sum(SA_AD == 1.9) ...
               sum(SA_AD == -1.9) sum(SA_AD == -1.5) sum(SA_AD == -1.1) ...
               sum((SA_AD == -0.9) | (SA_AD == 0.1)) sum((SA_AD == -0.5) | (SA_AD == 0.5)) sum((SA_AD == -0.1) | (SA_AD == 0.9))];
h = figure;
b = bar(1, SA_AD_Count, 'stacked');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_bar_sa.png'), 'Resolution', 300); close(h);

CT_AD_Count = [sum(CT_AD == 1.1) sum(CT_AD == 1.5) sum(CT_AD == 1.9) ...
               sum(CT_AD == -1.9) sum(CT_AD == -1.5) sum(CT_AD == -1.1) ...
               sum((CT_AD == -0.9) | (CT_AD == 0.1)) sum((CT_AD == -0.5) | (CT_AD == 0.5)) sum((CT_AD == -0.1) | (CT_AD == 0.9))];
h = figure;
b = bar(1, CT_AD_Count, 'stacked');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_bar_ct.png'), 'Resolution', 300); close(h);

[sum(CV_AD_Count(1:6)) sum(SA_AD_Count(1:6)) sum(CT_AD_Count(1:6))] ./ 360

%% Testing if proportion of agreement is significantly elvated against
%  the null model based on spin permutations of HCP sex difference maps

for i = 1:nPerm
    SPIN_HCP_CV_AD(:,i) = HCP_CV_AD(Vasa(:,i));
    SPIN_HCP_SA_AD(:,i) = HCP_SA_AD(Vasa(:,i));
    SPIN_HCP_CT_AD(:,i) = HCP_CT_AD(Vasa(:,i));
end

SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == 2) & (SPIN_HCP_CV_AD == 2)) = -1.9;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == 2) & (SPIN_HCP_CV_AD == 1)) = -1.5;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == 1) & (SPIN_HCP_CV_AD == 2)) = -1.5;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == 1) & (SPIN_HCP_CV_AD == 1)) = -1.1;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == 2) & (SPIN_HCP_CV_AD == -2)) = -0.9;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == 2) & (SPIN_HCP_CV_AD == -1)) = -0.5;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == 1) & (SPIN_HCP_CV_AD == -2)) = -0.5;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == 1) & (SPIN_HCP_CV_AD == -1)) = -0.1;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == -2) & (SPIN_HCP_CV_AD == 2)) = 0.1;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == -2) & (SPIN_HCP_CV_AD == 1)) = 0.5;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == -1) & (SPIN_HCP_CV_AD == 2)) = 0.5;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == -1) & (SPIN_HCP_CV_AD == 1)) = 0.9;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == -2) & (SPIN_HCP_CV_AD == -2)) = 1.1;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == -2) & (SPIN_HCP_CV_AD == -1)) = 1.5;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == -1) & (SPIN_HCP_CV_AD == -2)) = 1.5;
SPIN_CV_AD((repmat(UKB_CV_AD, nPerm, 1)' == -1) & (SPIN_HCP_CV_AD == -1)) = 1.9;

SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == 2) & (SPIN_HCP_SA_AD == 2)) = -1.9;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == 2) & (SPIN_HCP_SA_AD == 1)) = -1.5;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == 1) & (SPIN_HCP_SA_AD == 2)) = -1.5;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == 1) & (SPIN_HCP_SA_AD == 1)) = -1.1;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == 2) & (SPIN_HCP_SA_AD == -2)) = -0.9;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == 2) & (SPIN_HCP_SA_AD == -1)) = -0.5;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == 1) & (SPIN_HCP_SA_AD == -2)) = -0.5;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == 1) & (SPIN_HCP_SA_AD == -1)) = -0.1;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == -2) & (SPIN_HCP_SA_AD == 2)) = 0.1;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == -2) & (SPIN_HCP_SA_AD == 1)) = 0.5;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == -1) & (SPIN_HCP_SA_AD == 2)) = 0.5;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == -1) & (SPIN_HCP_SA_AD == 1)) = 0.9;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == -2) & (SPIN_HCP_SA_AD == -2)) = 1.1;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == -2) & (SPIN_HCP_SA_AD == -1)) = 1.5;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == -1) & (SPIN_HCP_SA_AD == -2)) = 1.5;
SPIN_SA_AD((repmat(UKB_SA_AD, nPerm, 1)' == -1) & (SPIN_HCP_SA_AD == -1)) = 1.9;

SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == 2) & (SPIN_HCP_CT_AD == 2)) = -1.9;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == 2) & (SPIN_HCP_CT_AD == 1)) = -1.5;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == 1) & (SPIN_HCP_CT_AD == 2)) = -1.5;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == 1) & (SPIN_HCP_CT_AD == 1)) = -1.1;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == 2) & (SPIN_HCP_CT_AD == -2)) = -0.9;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == 2) & (SPIN_HCP_CT_AD == -1)) = -0.5;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == 1) & (SPIN_HCP_CT_AD == -2)) = -0.5;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == 1) & (SPIN_HCP_CT_AD == -1)) = -0.1;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == -2) & (SPIN_HCP_CT_AD == 2)) = 0.1;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == -2) & (SPIN_HCP_CT_AD == 1)) = 0.5;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == -1) & (SPIN_HCP_CT_AD == 2)) = 0.5;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == -1) & (SPIN_HCP_CT_AD == 1)) = 0.9;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == -2) & (SPIN_HCP_CT_AD == -2)) = 1.1;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == -2) & (SPIN_HCP_CT_AD == -1)) = 1.5;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == -1) & (SPIN_HCP_CT_AD == -2)) = 1.5;
SPIN_CT_AD((repmat(UKB_CT_AD, nPerm, 1)' == -1) & (SPIN_HCP_CT_AD == -1)) = 1.9;

% Reshape
SPIN_CV_AD = reshape(SPIN_CV_AD, [360 nPerm]);
SPIN_SA_AD = reshape(SPIN_SA_AD, [360 nPerm]);
SPIN_CT_AD = reshape(SPIN_CT_AD, [360 nPerm]);

% Calculate proportions
SPIN_CV = sum(abs(SPIN_CV_AD) > 1) ./ 360;
SPIN_SA = sum(abs(SPIN_SA_AD) > 1) ./ 360;
SPIN_CT = sum(abs(SPIN_CT_AD) > 1) ./ 360;

% Calculate p value
CV_SPIN_P = sum((sum(abs(CV_AD) > 1) ./ 360) < SPIN_CV) ./ nPerm;
SA_SPIN_P = sum((sum(abs(SA_AD) > 1) ./ 360) < SPIN_SA) ./ nPerm;
CT_SPIN_P = sum((sum(abs(CT_AD) > 1) ./ 360) < SPIN_CT) ./ nPerm;

% --- Spin permutations --- (end)
