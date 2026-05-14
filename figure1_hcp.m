
clear all; close all; clc;

%% Prep

% Load processed data
load('matfiles/HCP_Processed.mat');

% Load surface templates
load('fsaverage_7.1.0.mat');
load('glasser_fsaverage_7.1.0.mat');
load('Vasa_G360.mat');

% Colormaps for figures
loadColorMaps;

% Misc
nPerm = 10000;
thresh = 0.05;
figDir = 'hcp';

%% Statistical testing for sex differences in cortical anatomy

% Compute Z-score across individuals
CV = zscore(MRI.CV);
SA = zscore(MRI.SA);
CT = zscore(MRI.CT);

% Prepare terms for statistical models
SEX = term(MRI.SEX);
AGE = term(MRI.AGE - mean(MRI.AGE)); % mean-centered age
EUL = term(MRI.EUL - mean(MRI.EUL)); % mean-centered Euler
G_CV = term(MRI.CV_GS - mean(MRI.CV_GS)); % mean-centered total CV
G_SA = term(MRI.SA_GS - mean(MRI.SA_GS)); % mean-centered total SA
G_CT = term(MRI.CT_GS - mean(MRI.CT_GS)); % mean-centered mean CT weighted by SA

% Define design matrices
M_CV = 1 + SEX + AGE + EUL + G_CV;
M_SA = 1 + SEX + AGE + EUL + G_SA;
M_CT = 1 + SEX + AGE + EUL + G_CT;

% Compute statiscal parameters
HCP_CV = SurfStatLinMod(CV, M_CV);
HCP_CV_SEX = SurfStatT(HCP_CV, SEX); % t score
HCP_CV_SEX_P = 2*(1 - tcdf(abs(HCP_CV_SEX.t), HCP_CV_SEX.df)); % two-tailed t test
[~,~,~,HCP_CV_SEX_Q] = fdr_bh(HCP_CV_SEX_P); % FDR across regions

HCP_SA = SurfStatLinMod(SA, M_SA);
HCP_SA_SEX = SurfStatT(HCP_SA, SEX); % t score
HCP_SA_SEX_P = 2*(1 - tcdf(abs(HCP_SA_SEX.t), HCP_SA_SEX.df)); % two-tailed t test
[~,~,~,HCP_SA_SEX_Q] = fdr_bh(HCP_SA_SEX_P); % FDR across regions

HCP_CT = SurfStatLinMod(CT, M_CT);
HCP_CT_SEX = SurfStatT(HCP_CT, SEX); % t score
HCP_CT_SEX_P = 2*(1 - tcdf(abs(HCP_CT_SEX.t), HCP_CT_SEX.df)); % two-tailed t test
[~,~,~,HCP_CT_SEX_Q] = fdr_bh(HCP_CT_SEX_P); % FDR across regions

% Permutation method 1: Spin permutations using Vasa's approach for parcellated cortical data
for i = 1:nPerm
    % Effect sizes for sex difference
    HCP_CV_SEX_Spin(:,i) = HCP_CV_SEX.coef(2,Vasa(:,i));
    HCP_SA_SEX_Spin(:,i) = HCP_SA_SEX.coef(2,Vasa(:,i));
    HCP_CT_SEX_Spin(:,i) = HCP_CT_SEX.coef(2,Vasa(:,i));

    % FDR corrected p values
    HCP_CV_SEX_Spin_Q(:,i) = HCP_CV_SEX_Q(Vasa(:,i));
    HCP_SA_SEX_Spin_Q(:,i) = HCP_SA_SEX_Q(Vasa(:,i));
    HCP_CT_SEX_Spin_Q(:,i) = HCP_CT_SEX_Q(Vasa(:,i));
end

% Permutation method 2: Sex permutations
for i = 1:nPerm
    disp(i)
    % Permutation index
    [~, ind] = datasample(CV, size(CV,1), 'Replace', false);

    % Define design matrices with sex permutation
    SEX_Perm = term(MRI.SEX(ind));
    M_CV_Perm = 1 + SEX_Perm + AGE + EUL + G_CV;
    M_SA_Perm = 1 + SEX_Perm + AGE + EUL + G_SA;
    M_CT_Perm = 1 + SEX_Perm + AGE + EUL + G_CT;

    % Compute statistical parameters
    HCP_CV_Perm = SurfStatT(SurfStatLinMod(CV, M_CV_Perm), SEX_Perm);
    HCP_SA_Perm = SurfStatT(SurfStatLinMod(SA, M_SA_Perm), SEX_Perm);
    HCP_CT_Perm = SurfStatT(SurfStatLinMod(CT, M_CT_Perm), SEX_Perm);

    HCP_CV_SEX_Perm_P(:,i) = 2*(1 - tcdf(abs(HCP_CV_Perm.t), HCP_CV_Perm.df));
    [~,~,~,HCP_CV_SEX_Perm_Q(:,i)] = fdr_bh(HCP_CV_SEX_Perm_P(:,i));

    HCP_SA_SEX_Perm_P(:,i) = 2*(1 - tcdf(abs(HCP_SA_Perm.t), HCP_SA_Perm.df));
    [~,~,~,HCP_SA_SEX_Perm_Q(:,i)] = fdr_bh(HCP_SA_SEX_Perm_P(:,i));

    HCP_CT_SEX_Perm_P(:,i) = 2*(1 - tcdf(abs(HCP_CT_Perm.t), HCP_CT_Perm.df));
    [~,~,~,HCP_CT_SEX_Perm_Q(:,i)] = fdr_bh(HCP_CT_SEX_Perm_P(:,i));

    % Extract effect sizes
    HCP_CV_SEX_Perm(:,i) = HCP_CV_Perm.coef(2,:);
    HCP_SA_SEX_Perm(:,i) = HCP_SA_Perm.coef(2,:);
    HCP_CT_SEX_Perm(:,i) = HCP_CT_Perm.coef(2,:);
end

% Save the results
save('matfiles/hcp_results.mat', 'HCP_CV_SEX', 'HCP_CV_SEX_Q', 'HCP_CV_SEX_P', ...
                                 'HCP_CV_SEX_Spin', 'HCP_CV_SEX_Spin_Q', ...
                                 'HCP_CV_SEX_Perm', 'HCP_CV_SEX_Perm_Q', ...
                                 'HCP_SA_SEX', 'HCP_SA_SEX_Q', 'HCP_SA_SEX_P', ...
                                 'HCP_SA_SEX_Spin', 'HCP_SA_SEX_Spin_Q', ...
                                 'HCP_SA_SEX_Perm', 'HCP_SA_SEX_Perm_Q', ...
                                 'HCP_CT_SEX', 'HCP_CT_SEX_Q', 'HCP_CT_SEX_P', ...
                                 'HCP_CT_SEX_Spin', 'HCP_CT_SEX_Spin_Q', ...
                                 'HCP_CT_SEX_Perm', 'HCP_CT_SEX_Perm_Q', ...
                                 '-v7.3');

%% Visualize regional sex difference maps

% Standardized sex effects after FDR thresholding
h = HyoSurfacePlotFunction(HCP_CV_SEX.coef(2,:) .* (HCP_CV_SEX_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_cv_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(HCP_SA_SEX.coef(2,:) .* (HCP_SA_SEX_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_sa_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(HCP_CT_SEX.coef(2,:) .* (HCP_CT_SEX_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_ct_fdr.png'), 'Resolution', 300); close(h);

% Standardized sex effects
h = HyoSurfacePlotFunction(HCP_CV_SEX.coef(2,:), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_cv.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(HCP_SA_SEX.coef(2,:), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_sa.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(HCP_CT_SEX.coef(2,:), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_ct.png'), 'Resolution', 300); close(h);

%% Scatter plots showing inter-feature relationships

% Calculate cortex-wide correlations between cortical features
[r_cv_sa ~] = corr(HCP_CV_SEX.coef(2,:)', HCP_SA_SEX.coef(2,:)');
[r_cv_ct ~] = corr(HCP_CV_SEX.coef(2,:)', HCP_CT_SEX.coef(2,:)');
[r_sa_ct ~] = corr(HCP_SA_SEX.coef(2,:)', HCP_CT_SEX.coef(2,:)');

% Generate scatter plots
h = figure('WindowState', 'maximized'); hold on;
scatter(HCP_SA_SEX.coef(2,:), HCP_CV_SEX.coef(2,:), 200, 'o', 'filled', 'MarkerFaceColor',[0 .5 .5], 'MarkerEdgeColor', [0 0 0]);
xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
xticklabels({}); yticklabels({});
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2; axis square;
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
axis square;
exportgraphics(h, fullfile(figDir, 'scatter_cv_sa.eps'), 'ContentType', 'vector');

h = figure('WindowState', 'maximized'); hold on;
scatter(HCP_CT_SEX.coef(2,:), HCP_CV_SEX.coef(2,:), 200, 'o', 'filled', 'MarkerFaceColor',[0 .5 .5], 'MarkerEdgeColor', [0 0 0]);
xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
xticklabels({}); yticklabels({});
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2; axis square;
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
axis square;
exportgraphics(h, fullfile(figDir, 'scatter_cv_ct.eps'), 'ContentType', 'vector');

h = figure('WindowState', 'maximized'); hold on;
scatter(HCP_CT_SEX.coef(2,:), HCP_SA_SEX.coef(2,:), 200, 'o', 'filled', 'MarkerFaceColor',[0 .5 .5], 'MarkerEdgeColor', [0 0 0]);
xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
xticklabels({}); yticklabels({});
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2; axis square;
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
axis square;
exportgraphics(h, fullfile(figDir, 'scatter_sa_ct.eps'), 'ContentType', 'vector');

% Test if correlations are significantly elevated against spin permutations
% (right-tailed)
p_spin_cv_sa = doPerms(HCP_CV_SEX_Spin, HCP_SA_SEX.coef(2,:), r_cv_sa, 'top');
p_spin_cv_ct = doPerms(HCP_CV_SEX_Spin, HCP_CT_SEX.coef(2,:), r_cv_ct, 'top');
p_spin_sa_ct = doPerms(HCP_SA_SEX_Spin, HCP_CT_SEX.coef(2,:), r_sa_ct, 'top');

% Test if correlations are significantly elevated against sex permutations
% (right-tailed)
p_perm_cv_sa = doPerms(HCP_CV_SEX_Perm, HCP_SA_SEX.coef(2,:), r_cv_sa, 'top');
p_perm_cv_ct = doPerms(HCP_CV_SEX_Perm, HCP_CT_SEX.coef(2,:), r_cv_ct, 'top');
p_perm_sa_ct = doPerms(HCP_SA_SEX_Perm, HCP_CT_SEX.coef(2,:), r_sa_ct, 'top');

%% Composite log-ratio map for regional contributions of SA and CT to sex-biased CV

% Compute log composite map
composite = log10(abs(HCP_CV_SEX.coef(2,:) - HCP_CT_SEX.coef(2,:)) ./ abs(HCP_CV_SEX.coef(2,:) - HCP_SA_SEX.coef(2,:)));

% Composite map for female- and male-biased regions (log ratios are clipped
% at 1.45 and off-centered at 1.5 so that they can be displayed on two
% different colormaps ranging from -3 to 0 for male-biases and from 0 to 3
% for female-biases)
composite_female = composite(HCP_CV_SEX.coef(2,:) > 0);
composite_female(composite_female < -1.45) = -1.45;
composite_female(composite_female > 1.45) = 1.45;
composite_female = composite_female + 1.5;

composite_male = composite(HCP_CV_SEX.coef(2,:) < 0);
composite_male(composite_male < -1.45) = -1.45;
composite_male(composite_male > 1.45) = 1.45;
composite_male = composite_male - 1.5;

% Compute the female and male composite maps into one
composite_overlaid = zeros(1,360);
composite_overlaid(HCP_CV_SEX.coef(2,:) > 0) = composite_female;
composite_overlaid(HCP_CV_SEX.coef(2,:) < 0) = composite_male;

% Mask the brain regions that are not significantly sex biased
composite_overlaid_thr = composite_overlaid;
composite_overlaid_thr((HCP_CV_SEX_Q > thresh)) = -1.5;

% Display on two distinct colormaps
h = HyoSurfacePlotFunction(composite_overlaid_thr, IS, G360, cmap.comp, '', [-3 3], 'square');
exportgraphics(h, fullfile(figDir, 'log_ratio.png'), 'Resolution', 300); close(h);

%% Anotate cortical regions based on significances in CV, SA and CT

% Regional conjunction map
conjunction = zeros(1,360);

% Significant CV regions (gray)
conjunction(HCP_CV_SEX_Q < thresh) = 0.5;

% Insignificant in CV (white)
conjunction(HCP_CV_SEX_Q > thresh) = 1.5;

% SA significant only
conjunction((HCP_CV_SEX_Q > thresh) &  (HCP_SA_SEX_Q < thresh) & (HCP_CT_SEX_Q > thresh)) = 2.5;

% CT significant only
conjunction((HCP_CV_SEX_Q > thresh) &  (HCP_SA_SEX_Q > thresh) & (HCP_CT_SEX_Q < thresh)) = 3.5;

% SA and CT significant in opposite direction
conjunction((HCP_CV_SEX_Q > thresh) &  ((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) > 0)) & ((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) < 0))) = 4.5;
conjunction((HCP_CV_SEX_Q > thresh) &  ((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) < 0)) & ((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) > 0))) = 4.5;

% SA and CT significant in same direction
conjunction((HCP_CV_SEX_Q > thresh) &  ((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) < 0)) & ((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) < 0))) = 5.5;
conjunction((HCP_CV_SEX_Q > thresh) &  ((HCP_SA_SEX_Q < thresh) & (HCP_SA_SEX.coef(2,:) > 0)) & ((HCP_CT_SEX_Q < thresh) & (HCP_CT_SEX.coef(2,:) > 0))) = 5.5;

% Visualize map
h = HyoSurfacePlotFunction(conjunction, IS, G360, cmap.qual, '', [0 6], 'square');
exportgraphics(h, fullfile(figDir, 'conjunction.png'), 'Resolution', 300); close(h);

% Display percentages
insigCV = [100 * sum(conjunction == 0.5) / 360
100 * sum(conjunction == 1.5) / 360
100 * sum(conjunction == 2.5) / 360
100 * sum(conjunction == 3.5) / 360
100 * sum(conjunction == 4.5) / 360
100 * sum(conjunction == 5.5) / 360];

% Show bar plot
h = figure;
b = bar(1, insigCV([1 3 4 5 6 2]), 'stacked');
c = 1;
for i = [1 3 4 5 6 2]
    b(c).BarWidth = 0.1;
    b(c).FaceColor = cmap.qual(i,:);
    c = c + 1;
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_bar.png'), 'Resolution', 300); close(h);
