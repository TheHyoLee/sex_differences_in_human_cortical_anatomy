
clear all; close all; clc;

%% Prep

% Load results
load('matfiles/hcp_results.mat');

% Load template
load('fsaverage_7.1.0.mat');
load('glasser_fsaverage_7.1.0.mat');
load('Vasa_G360.mat');

% Misc
loadColorMaps;
thresh = 0.05;
nPerm = 10000;
nTests = 6; % = 3 cortical features x 2 directions for sex bias

% Set up paths
figDir = 'gsea';

%% Prep MAGICC

% Load MAGICC-AHBA data mapped onto Glasser parcellation (Wagstyl et
% al.,2023, eLife)
load('matfiles/ahba_glasser.mat', 'GO_GE', 'GO_GE_ID', 'GE_Glasser', 'GE_ID');
GE_Glasser = repmat(GE_Glasser, 2, 1);

% Select cortex-coding genes
indCortex = GO_GE(:,53) == 1;
GE_Glasser = GE_Glasser(:,indCortex);
GE_ID = GE_ID(indCortex);

% Select relevent columns
indCol = [2:52 54:57];
GO_GE = GO_GE(indCortex,indCol) == 1;
GO_GE_ID = GO_GE_ID(indCol); 

%% Compute correlation between sex difference and spatial transcripts (MAGICC)

% Precompute indices for female and male biased regions
indCV_F = HCP_CV_SEX.coef(2,:) > 0;
indSA_F = HCP_SA_SEX.coef(2,:) > 0;
indCT_F = HCP_CT_SEX.coef(2,:) > 0;
indCV_M = HCP_CV_SEX.coef(2,:) < 0;
indSA_M = HCP_SA_SEX.coef(2,:) < 0;
indCT_M = HCP_CT_SEX.coef(2,:) < 0;

% Compute imaging-transriptomic associations
for i = 1:size(GE_Glasser,2) % For each gene
    % Show progress
    disp(['Iteration ' num2str(i)]);

    % Calculate spatial correlation
    [MAGICC_R_CV_F(i) ~] = corr(GE_Glasser(indCV_F,i), HCP_CV_SEX.coef(2,indCV_F)');
    [MAGICC_R_SA_F(i) ~] = corr(GE_Glasser(indSA_F,i), HCP_SA_SEX.coef(2,indSA_F)');
    [MAGICC_R_CT_F(i) ~] = corr(GE_Glasser(indCT_F,i), HCP_CT_SEX.coef(2,indCT_F)');
    [MAGICC_R_CV_M(i) ~] = corr(GE_Glasser(indCV_M,i), -HCP_CV_SEX.coef(2,indCV_M)');
    [MAGICC_R_SA_M(i) ~] = corr(GE_Glasser(indSA_M,i), -HCP_SA_SEX.coef(2,indSA_M)');
    [MAGICC_R_CT_M(i) ~] = corr(GE_Glasser(indCT_M,i), -HCP_CT_SEX.coef(2,indCT_M)');
end

% Indices for spin permutations
indCV_F_Spin = HCP_CV_SEX_Spin > 0;
indSA_F_Spin = HCP_SA_SEX_Spin > 0;
indCT_F_Spin = HCP_CT_SEX_Spin > 0;
indCV_M_Spin = HCP_CV_SEX_Spin < 0;
indSA_M_Spin = HCP_SA_SEX_Spin < 0;
indCT_M_Spin = HCP_CT_SEX_Spin < 0;

parfor i = 1:size(GE_Glasser,2) % For each gene
    % Show progress
    disp(['Iteration ' num2str(i)]);

    for j = 1:nPerm
        % Calculate spatial correlation between gene expression and sex
        % difference
        [MAGICC_R_CV_F_Spin(i,j) ~] = corr(GE_Glasser(indCV_F,i), HCP_CV_SEX_Spin(indCV_F_Spin(:,j),j));
        [MAGICC_R_SA_F_Spin(i,j) ~] = corr(GE_Glasser(indSA_F,i), HCP_SA_SEX_Spin(indSA_F_Spin(:,j),j));
        [MAGICC_R_CT_F_Spin(i,j) ~] = corr(GE_Glasser(indCT_F,i), HCP_CT_SEX_Spin(indCT_F_Spin(:,j),j));
        [MAGICC_R_CV_M_Spin(i,j) ~] = corr(GE_Glasser(indCV_M,i), -HCP_CV_SEX_Spin(indCV_M_Spin(:,j),j));
        [MAGICC_R_SA_M_Spin(i,j) ~] = corr(GE_Glasser(indSA_M,i), -HCP_SA_SEX_Spin(indSA_M_Spin(:,j),j));
        [MAGICC_R_CT_M_Spin(i,j) ~] = corr(GE_Glasser(indCT_M,i), -HCP_CT_SEX_Spin(indCT_M_Spin(:,j),j));
    end
end

% Write out the spin permuted results
for i = 1:nPerm
    % Show progress
    disp(['Iteration ' num2str(i)]);

    writematrix(MAGICC_R_CV_F_Spin(:,i), fullfile(figDir, ['/GSEA_R/input/MAGICC_R_CV_F_Spin_BS' num2str(i) '.txt']));
    writematrix(MAGICC_R_SA_F_Spin(:,i), fullfile(figDir, ['/GSEA_R/input/MAGICC_R_SA_F_Spin_BS' num2str(i) '.txt']));
    writematrix(MAGICC_R_CT_F_Spin(:,i), fullfile(figDir, ['/GSEA_R/input/MAGICC_R_CT_F_Spin_BS' num2str(i) '.txt']));
    writematrix(MAGICC_R_CV_M_Spin(:,i), fullfile(figDir, ['/GSEA_R/input/MAGICC_R_CV_M_Spin_BS' num2str(i) '.txt']));
    writematrix(MAGICC_R_SA_M_Spin(:,i), fullfile(figDir, ['/GSEA_R/input/MAGICC_R_SA_M_Spin_BS' num2str(i) '.txt']));
    writematrix(MAGICC_R_CT_M_Spin(:,i), fullfile(figDir, ['/GSEA_R/input/MAGICC_R_CT_M_Spin_BS' num2str(i) '.txt']));
end

%% Perform GSEA in R (see the accompanying R scripts)


%% Read the GSEA results done in R

% Pathways
MAGICC_terms = sort(GO_GE_ID);

% Read GSEA results from R
MAGICC_FCV_pos = xlsread(['gsea/output_main/MAGICC_GSEA_F_CV_pos.xlsx']);
MAGICC_FSA_pos = xlsread(['gsea/output_main/MAGICC_GSEA_F_SA_pos.xlsx']);
MAGICC_FCT_pos = xlsread(['gsea/output_main/MAGICC_GSEA_F_CT_pos.xlsx']);
MAGICC_MCV_pos = xlsread(['gsea/output_main/MAGICC_GSEA_M_CV_pos.xlsx']);
MAGICC_MSA_pos = xlsread(['gsea/output_main/MAGICC_GSEA_M_SA_pos.xlsx']);
MAGICC_MCT_pos = xlsread(['gsea/output_main/MAGICC_GSEA_M_CT_pos.xlsx']);

% Read the results based on spin-permuted HCP sex difference maps from R
parfor iBS = 1:nPerm
    MAGICC_FCV_pos_bs(:,:,iBS) = xlsread(['gsea/output_bs/MAGICC_GSEA_F_CV_pos_BS' num2str(iBS) '.xlsx']);
    MAGICC_FSA_pos_bs(:,:,iBS) = xlsread(['gsea/output_bs/MAGICC_GSEA_F_SA_pos_BS' num2str(iBS) '.xlsx']);
    MAGICC_FCT_pos_bs(:,:,iBS) = xlsread(['gsea/output_bs/MAGICC_GSEA_F_CT_pos_BS' num2str(iBS) '.xlsx']);
    MAGICC_MCV_pos_bs(:,:,iBS) = xlsread(['gsea/output_bs/MAGICC_GSEA_M_CV_pos_BS' num2str(iBS) '.xlsx']);
    MAGICC_MSA_pos_bs(:,:,iBS) = xlsread(['gsea/output_bs/MAGICC_GSEA_M_SA_pos_BS' num2str(iBS) '.xlsx']);
    MAGICC_MCT_pos_bs(:,:,iBS) = xlsread(['gsea/output_bs/MAGICC_GSEA_M_CT_pos_BS' num2str(iBS) '.xlsx']);
end
save('matfiles/MAGICC_BS.mat', 'MAGICC_FCV_pos_bs', 'MAGICC_FSA_pos_bs', 'MAGICC_FCT_pos_bs', ...
                               'MAGICC_MCV_pos_bs', 'MAGICC_MSA_pos_bs', 'MAGICC_MCT_pos_bs', '-v7.3');

% Extract the normalize enrichment scores and p-values
MAGICC_NES = [MAGICC_FCV_pos(:,5)'; ...
              MAGICC_MCV_pos(:,5)'; ...
              zeros(1,55);
              MAGICC_FSA_pos(:,5)'; ...
              MAGICC_MSA_pos(:,5)'; ...
              zeros(1,55);
              MAGICC_FCT_pos(:,5)'; ...
              MAGICC_MCT_pos(:,5)'];

% Compute the p-values for normalized enrichment scores against the ones
% from the null distribution (derived from spin permutations of HCP sex
% differences)
MAGICC_P_PERM = [sum(MAGICC_FCV_pos(:,1) > squeeze(MAGICC_FCV_pos_bs(:,1,:)),2)'; ...
                 sum(MAGICC_MCV_pos(:,1) > squeeze(MAGICC_MCV_pos_bs(:,1,:)),2)'; ...
                 zeros(1,55);
                 sum(MAGICC_FSA_pos(:,1) > squeeze(MAGICC_FSA_pos_bs(:,1,:)),2)'; ...
                 sum(MAGICC_MSA_pos(:,1) > squeeze(MAGICC_MSA_pos_bs(:,1,:)),2)'; ...
                 zeros(1,55);
                 sum(MAGICC_FCT_pos(:,1) > squeeze(MAGICC_FCT_pos_bs(:,1,:)),2)'; ...
                 sum(MAGICC_MCT_pos(:,1) > squeeze(MAGICC_MCT_pos_bs(:,1,:)),2)'] ./ nPerm;

% Bonf. correction
MAGICC_Q_PERM = MAGICC_P_PERM .* nTests;

%% Visualize the GSEA results

% Misc
MAGICC_NES(:,end+1) = 0; % Dummy column
MAGICC_P_PERM(:,end+1) = 0; % Dummy column
MAGICC_Q_PERM(:,end+1) = 0; % Dummy column

% Orders of ontological terms
% 1:3 Epochs
% 4:10 Cells
% 11:21 Cell compartment
% 22:28 Fetal layers
% 29:44 Fetal cells
% 45:50 Layers
% 51:54 Rare diseases
% 55: SynGO
MAGICC_terms{56} = '';

% Plot order
ordind = [1:3 56 22:28 56 29:44 56 45:50 56 4:10 56 11:21 55 56];

% Visualize
h = figure;
imagesc(abs(MAGICC_NES(:,ordind)), [0 3.5]); hold on;
colormap(cmap.Purples);
xticks(1:length(ordind));
xticklabels(strrep(MAGICC_terms(ordind), '_', ' '));
xticklabels({}); yticklabels({});
xticks([]); yticks([]);
ax = gca;
ax.TickDir = 'out';
ax.Box = 'off';
h.Position = [0 0 3000 1000];
for i = ordind
    for j = [1:2 4:5 7:8] % only pos enrichments
        if i < 56
            if MAGICC_Q_PERM(j,i) < thresh % Dashed box
                vertices = 0.5 + [find(ordind==i)-1 j-1; find(ordind==i) j-1; find(ordind==i) j; find(ordind==i)-1 j];
                patch('Faces', 1:4, 'Vertices', vertices, 'FaceColor', 'none', 'LineWidth', 6, 'EdgeColor', [0.7 0.7 0]);
            elseif MAGICC_P_PERM(j,i) < thresh % Solid box
                vertices = 0.5 + [find(ordind==i)-1 j-1; find(ordind==i) j-1; find(ordind==i) j; find(ordind==i)-1 j];
                patch('Faces', 1:4, 'Vertices', vertices, 'FaceColor', 'none', 'LineWidth', 6, 'EdgeColor', [0.7 0.7 0], 'LineStyle', ':');
            end
        end
    end
end
exportgraphics(h, fullfile(figDir, 'MAGICC_GSEA_unipos_BS_pPerm.eps'), 'ContentType', 'vector');

close all;
