
clear all; close all; clc;

%% Prep

% Load data
load('matfiles/HCP_Processed.mat');

% Load template
load('fsaverage_7.1.0.mat');
load('glasser_fsaverage_7.1.0.mat');

% Sampling parameters
nRep = 10000;
sRep = 50:50:450;

% Misc
figDir = 'replication_within_hcp';
thresh = 0.05;

%% With TBV and EUL as a linear covariate

% Bootstrap sampling
for iS = 1:length(sRep)

    disp(['Current sample size: ' num2str(sRep(iS))]);

    for iBS = 1:nRep % Iteration

        % Sex stratified bootstrap sampling
        [~, indBS_F] = datasample(MRI.CV(MRI.SEX==1,:), sRep(iS), 'Replace', false);
        [~, indBS_M] = datasample(MRI.CV(MRI.SEX==0,:), sRep(iS), 'Replace', false);
        tempF = find(MRI.SEX==1);
        tempM = find(MRI.SEX==0);
        indBS1 = [tempF(indBS_F(1:end/2)); tempM(indBS_M(1:end/2))];
        indBS2 = [tempF(indBS_F(end/2+1:end)); tempM(indBS_M(end/2+1:end))];

        % Prepare features
        iCV1 = CV(indBS1,:);
        iSA1 = SA(indBS1,:);
        iCT1 = CT(indBS1,:);
        iCV2 = CV(indBS2,:);
        iSA2 = SA(indBS2,:);
        iCT2 = CT(indBS2,:);

        % Prepare terms
        SEX1 = term(MRI.SEX(indBS1));
        AGE1 = term(MRI.AGE(indBS1) - mean(MRI.AGE(indBS1)));
        GS_CV1 = term(MRI.CV_GS(indBS1) - mean(MRI.CV_GS(indBS1)));
        GS_SA1 = term(MRI.SA_GS(indBS1) - mean(MRI.SA_GS(indBS1)));
        GS_CT1 = term(MRI.CT_GS(indBS1) - mean(MRI.CT_GS(indBS1)));
        EUL1 = term(MRI.EUL(indBS1) - mean(MRI.EUL(indBS1)));

        SEX2 = term(MRI.SEX(indBS2));
        AGE2 = term(MRI.AGE(indBS2) - mean(MRI.AGE(indBS2)));
        GS_CV2 = term(MRI.CV_GS(indBS2) - mean(MRI.CV_GS(indBS2)));
        GS_SA2 = term(MRI.SA_GS(indBS2) - mean(MRI.SA_GS(indBS2)));
        GS_CT2 = term(MRI.CT_GS(indBS2) - mean(MRI.CT_GS(indBS2)));
        EUL2 = term(MRI.EUL(indBS2) - mean(MRI.EUL(indBS2)));

        % Define GLM
        M1_CV = 1 + SEX1 + AGE1 + GS_CV1 + EUL1;
        M1_SA = 1 + SEX1 + AGE1 + GS_SA1 + EUL1;
        M1_CT = 1 + SEX1 + AGE1 + GS_CT1 + EUL1;
        M2_CV = 1 + SEX2 + AGE2 + GS_CV2 + EUL2;
        M2_SA = 1 + SEX2 + AGE2 + GS_SA2 + EUL2;
        M2_CT = 1 + SEX2 + AGE2 + GS_CT2 + EUL2;

        % Compute stats (first half)
        SLM1_CV = SurfStatLinMod(iCV1, M1_CV);
        SLM1_CV_SEX = SurfStatT(SLM1_CV, SEX1);

        SLM1_SA = SurfStatLinMod(iSA1, M1_SA);
        SLM1_SA_SEX = SurfStatT(SLM1_SA, SEX1);

        SLM1_CT = SurfStatLinMod(iCT1, M1_CT);
        SLM1_CT_SEX = SurfStatT(SLM1_CT, SEX1);

        % Compute stats (second half)
        SLM2_CV = SurfStatLinMod(iCV2, M2_CV);
        SLM2_CV_SEX = SurfStatT(SLM2_CV, SEX2);

        SLM2_SA = SurfStatLinMod(iSA2, M2_SA);
        SLM2_SA_SEX = SurfStatT(SLM2_SA, SEX2);

        SLM2_CT = SurfStatLinMod(iCT2, M2_CT);
        SLM2_CT_SEX = SurfStatT(SLM2_CT, SEX2);

        % Extract effect sizes for sex difference (first half)
        CV_B(:,iBS,iS,1) = SLM1_CV_SEX.coef(2,:);
        SA_B(:,iBS,iS,1) = SLM1_SA_SEX.coef(2,:);
        CT_B(:,iBS,iS,1) = SLM1_CT_SEX.coef(2,:);

        % Extract effect sizes for sex difference (second half)
        CV_B(:,iBS,iS,2) = SLM2_CV_SEX.coef(2,:);
        SA_B(:,iBS,iS,2) = SLM2_SA_SEX.coef(2,:);
        CT_B(:,iBS,iS,2) = SLM2_CT_SEX.coef(2,:);

        % Cortex-wide spatial correlation between the two halves
        CV_R(iBS,iS) = corr(CV_B(:,iBS,iS,1), CV_B(:,iBS,iS,2));
        SA_R(iBS,iS) = corr(SA_B(:,iBS,iS,1), SA_B(:,iBS,iS,2));
        CT_R(iBS,iS) = corr(CT_B(:,iBS,iS,1), CT_B(:,iBS,iS,2));
    end
end

% Save results
save('matfiles/hcp_replication.mat', 'CV_B', 'SA_B', 'CT_B', 'CV_R', 'SA_R', 'CT_R', 'nRep', 'sRep', '-v7.3');

%% Raincloud plot for visualization

data{1} = CV_R;
data{2} = SA_R;
data{3} = CT_R;
group_inx = [ones(1,nRep), 2.*ones(1,nRep) 3.*ones(1,nRep) 4.*ones(1,nRep) 5.*ones(1,nRep) 6.*ones(1,nRep) 7.*ones(1,nRep) 8.*ones(1,nRep) 9.*ones(1,nRep)];
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30]; 
group_names = {'CV', 'SA', 'CT'};
condition_names = {'50', '100', '150', '200', '250', '300', '350', '400', '450'};

h = figure; d = daviolinplot(data, 'groups', group_inx(1:90), 'colors', c, 'box', 3,...
    'boxcolor', 'w', 'scatter', 2, 'jitter', 1, 'scattercolor', 'same',...
    'scattersize', 10, 'scatteralpha', 0.7, 'linkline', 1, 'outliers', 0, ...
    'xtlabels', condition_names,...
    'legend', group_names(1:2)); 
% ylabel('Performance');
set(d.sc,'MarkerEdgeColor','none');
set(gca,'FontSize',10);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]);
ylim([-0.3 0.9]); yticks(-0.3:0.1:0.9);
xticklabels({}); yticklabels({});
ax = gca;
ax.YGrid = 'on';
ax.TickDir = 'out';
h.Position = [0 0 2000 1000];
legend('off')
exportgraphics(h, fullfile(figDir, ['raincloud_R_' num2str(nRep) '.eps']), 'ContentType', 'vector');

close all;
