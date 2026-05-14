
clear all; close all; clc;

%% Prep

% Load HCP results
load('matfiles/hcp_results.mat');

% Load template
load('fsaverage_7.1.0.mat');
load('glasser_fsaverage_7.1.0.mat');
load('Vasa_G360.mat');
 
% Parcellations
YeoKrienen = vertexToParcel(fetch_parcellation('fsaverage', 'yeo', 17), G360, 'mode');
[YeoKrienenKrienen_annot YeoKrienen_cmap] = fetch_yeo_networks_metadata(17);
YeoKrienen_cmap = [0 0 0; YeoKrienen_cmap];

% Hippocampus and piriform cortex are assigned to limbic network
YeoKrienen([110 120 290 300]) = 9;

% Colormaps for this figure
loadColorMaps;

% Misc
nPerm = 10000;
nTests = 3; % across 3 features
thresh = 0.05;

% Figure path
figDir = 'function';

%% Spin-based enrichment test

% Spin functional network maps
indPerm = randperm(nPerm);
for i = 1:nPerm
    YeoKrienen_Spin(:,i) = YeoKrienen(Vasa(:,indPerm(i)));
end

% Compute mean effect size in HCP sex differences for each network
for i = 1:nNetworks
    CV_YeoKrienen_ES(i) = mean(HCP_CV_SEX.coef(2,YeoKrienen==i));
    SA_YeoKrienen_ES(i) = mean(HCP_SA_SEX.coef(2,YeoKrienen==i));
    CT_YeoKrienen_ES(i) = mean(HCP_CT_SEX.coef(2,YeoKrienen==i));
    
    % Compute mean effect size in spin-permuted YeoKrienen maps using
    % indendently spun HCP sex differences
    for j = 1:nPerm
        CV_YeoKrienen_ES_Null(i,j) = mean(HCP_CV_SEX_Spin(YeoKrienen_Spin(:,j)==i,j));
        SA_YeoKrienen_ES_Null(i,j) = mean(HCP_SA_SEX_Spin(YeoKrienen_Spin(:,j)==i,j));
        CT_YeoKrienen_ES_Null(i,j) = mean(HCP_CT_SEX_Spin(YeoKrienen_Spin(:,j)==i,j));
    end
end

% Calculate empirical p value (two-tailed)
CV_YeoKrienen_P = 2 * min([sum(CV_YeoKrienen_ES_Null > CV_YeoKrienen_ES', 2) sum(CV_YeoKrienen_ES_Null < CV_YeoKrienen_ES', 2)]') / nPerm;
SA_YeoKrienen_P = 2 * min([sum(SA_YeoKrienen_ES_Null > SA_YeoKrienen_ES', 2) sum(SA_YeoKrienen_ES_Null < SA_YeoKrienen_ES', 2)]') / nPerm;
CT_YeoKrienen_P = 2 * min([sum(CT_YeoKrienen_ES_Null > CT_YeoKrienen_ES', 2) sum(CT_YeoKrienen_ES_Null < CT_YeoKrienen_ES', 2)]') / nPerm;

% Bonferroni correction across cortical features
CV_YeoKrienen_Q = CV_YeoKrienen_P * nTests;
SA_YeoKrienen_Q = SA_YeoKrienen_P * nTests;
CT_YeoKrienen_Q = CT_YeoKrienen_P * nTests;

%% Visualize

% Bar plots
[~, ind] = sort(CV_YeoKrienen_ES, 'descend', 'MissingPlacement', 'last');
h = figure; hold on;
for i = 1:nNetworks
    b(i) = bar(i, CV_YeoKrienen_ES(ind(i)));
    b(i).FaceColor = YeoKrienen_cmap(ind(i)+1,:);
    b(i).LineStyle = 'none';
    if CV_YeoKrienen_Q(ind(i)) <= 0.05
        b(i).LineStyle = '-';
        b(i).LineWidth = 6;
    elseif CV_YeoKrienen_P(ind(i)) <= 0.05
        b(i).LineStyle = '--';
        b(i).LineWidth = 4;
    else
        b(i).FaceAlpha = 0.2;
    end
end
xticks(1:nNetworks);
xticklabels(YeoKrienen_annot(ind));
yticks(-0.2:0.1:0.2);
xticklabels({}); yticklabels({});
ylim([-0.2 0.2]);
ax = gca;
ax.YGrid = 'on';
ax.TickDir = 'out';
ax.Box = 'off';
h.Position = [0 0 2000 1000];
exportgraphics(h, fullfile(figDir, 'YeoKrienen_CV_Spin.eps'), 'ContentType', 'vector');


[~, ind] = sort(SA_YeoKrienen_ES, 'descend', 'MissingPlacement', 'last');
h = figure; hold on;
for i = 1:nNetworks
    b(i) = bar(i, SA_YeoKrienen_ES(ind(i)));
    b(i).FaceColor = YeoKrienen_cmap(ind(i)+1,:);
    b(i).LineStyle = 'none';
    if SA_YeoKrienen_Q(ind(i)) <= 0.05
        b(i).LineStyle = '-';
        b(i).LineWidth = 6;
    elseif SA_YeoKrienen_P(ind(i)) <= 0.05
        b(i).LineStyle = '--';
        b(i).LineWidth = 4;
    else
        b(i).FaceAlpha = 0.2;
    end
end
xticks(1:nNetworks);
xticklabels(YeoKrienen_annot(ind));
yticks(-0.2:0.1:0.2);
ylim([-0.2 0.2]);
xticklabels({}); yticklabels({});
ax = gca;
ax.YGrid = 'on';
ax.TickDir = 'out';
ax.Box = 'off';
h.Position = [0 0 2000 1000];
exportgraphics(h, fullfile(figDir, 'YeoKrienen_SA_Spin.eps'), 'ContentType', 'vector');


[~, ind] = sort(CT_YeoKrienen_ES, 'descend', 'MissingPlacement', 'last');
h = figure; hold on;
for i = 1:nNetworks
    b(i) = bar(i, CT_YeoKrienen_ES(ind(i)));
    b(i).FaceColor = YeoKrienen_cmap(ind(i)+1,:);
    b(i).LineStyle = 'none';
    if CT_YeoKrienen_Q(ind(i)) <= 0.05
        b(i).LineStyle = '-';
        b(i).LineWidth = 6;
    elseif CT_YeoKrienen_P(ind(i)) <= 0.05
        b(i).LineStyle = '--';
        b(i).LineWidth = 4;
    else
        b(i).FaceAlpha = 0.2;
    end
end
xticks(1:nNetworks);
xticklabels(YeoKrienen_annot(ind));
yticks(-0.2:0.1:0.2);
ylim([-0.2 0.2]);
xticklabels({}); yticklabels({});
ax = gca;
ax.YGrid = 'on';
ax.TickDir = 'out';
ax.Box = 'off';
h.Position = [0 0 2000 1000];
exportgraphics(h, fullfile(figDir, 'YeoKrienen_CT_Spin.eps'), 'ContentType', 'vector');