%% Permutation DFA of bat ID for 7 acoustic groups

LocalDataDir = '/Users/elie/Documents/DeafBats/Data';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
SorKOnly = 2; % Set to 1 to only run calculations with saline Bats, 2 to only run on Kanamycin,  0 otherwise (run on both)
%%
% Loading previous data

load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls2.mat'), 'CallType', 'BatID')
load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls.mat'),'MicAudioGood','TmicAll7');
MicAudioGood01 = MicAudioGood;
MicAudioGood01(isnan(MicAudioGood01)) = 0;
MicAudioGood01 = logical(MicAudioGood01);
BatIDMic = BatID(MicAudioGood01);
%% %% Load MPS for this part and Restrict the MPS spectral frequency to 3 cycles/kHz and normalized all
% MPS

load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls2.mat'),'MPS_mic', 'MPS_mic_wf', 'MPS_mic_wt')
IndWf = (MPS_mic_wf*1000)<=3; % Max Spectral frequency taken into account in cycles/kHz
NVoc = sum(MicAudioGood01);
MPS_mic_norm = cell(1,NVoc);
MPS_mic_mean = reshape(mean(MPS_mic(MicAudioGood01,:), 'omitnan')',length(MPS_mic_wf), length(MPS_mic_wt));
IndMicAudioGood = find(MicAudioGood01);
for nn=1:NVoc
    vv=IndMicAudioGood(nn);
    MPS_local = reshape(MPS_mic(vv,:)',length(MPS_mic_wf), length(MPS_mic_wt));
    MPS_local = MPS_local(IndWf,:);
    MPS_mic_norm{nn} = reshape(MPS_local./MPS_mic_mean(IndWf,:), numel(MPS_local),1);
end
MPS_mic_norm = [MPS_mic_norm{:}]';
%%
% Get color vectors ready
ColorCode = [get(groot, 'DefaultAxesColorOrder'); 0 1 1; 0.5 0.5 0.5; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0 1; 0 0 0];

% Get the color vector ready for BatName, Sex and Deafness
GGPath = dir('/Users/elie/Google Drive*');
Path2RecordingTable = fullfile(GGPath.folder, GGPath.name,'My Drive/JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx');
[~,~,RecTableData]=xlsread(Path2RecordingTable,2,'A1:k3','basic');
BatName = cell2mat(RecTableData(3,2:11));
BatSexDeaf = cell(size(BatName));
CSexDeaf =  nan(length(BatID),3);
SexDeaf =  cell(length(BatID),1);
Sex =  cell(length(BatID),1);
Deaf = cell(length(BatID),1);
CBat = nan(length(BatID),3);
UCBatName = nan(length(BatName),3);
for bat = 1:length(BatName)
    ColID = 1+ bat;
    CBat(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(bat,:),sum(contains(BatID,num2str(BatName(bat)))),1);
    UCBatName(bat,:) = ColorCode(bat,:);
    if strcmp(RecTableData(1,ColID), 'K') && strcmp(RecTableData(2,ColID), 'M')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(1,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'DM'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'M'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'D'};
        BatSexDeaf{bat} = 'Deaf Male';
    elseif strcmp(RecTableData(1,ColID), 'S') && strcmp(RecTableData(2,ColID), 'M')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(6,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'HM'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'M'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'H'};
        BatSexDeaf{bat} = 'Hearing Male';
    elseif strcmp(RecTableData(1,ColID), 'S') && strcmp(RecTableData(2,ColID), 'F')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(3,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'HF'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'F'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'H'};
        BatSexDeaf{bat} = 'Hearing Female';
    elseif strcmp(RecTableData(1,ColID), 'K') && strcmp(RecTableData(2,ColID), 'F')
        CSexDeaf(contains(BatID,num2str(BatName(bat))),:) = repmat(ColorCode(2,:),sum(contains(BatID,num2str(BatName(bat)))),1);
        SexDeaf(contains(BatID,num2str(BatName(bat))),:) = {'DF'};
        Sex(contains(BatID,num2str(BatName(bat))),:) = {'F'};
        Deaf(contains(BatID,num2str(BatName(bat))),:) = {'D'};
        BatSexDeaf{bat} = 'Deaf Female';
    end

end
USexDeaf = {'Deaf Male' 'Hearing Male' 'Deaf Female' 'Hearing Female'};
UCSexDeaf = ColorCode([1 6 2 3], :);
BatName = BatName([1:7 10 8:9]);
BatSexDeaf = BatSexDeaf([1:7 10 8:9]);
UCBatName = UCBatName([1:7 10 8:9],:);
UCBatNameSexDeaf = ColorCode([2 2 2 3 3 3 1 1 6 6], :);
SexDeafMic = SexDeaf(MicAudioGood01);
%% Parameters for random permutations
NRandPerm=100;
NumPCsup = 10;

%% For All Male Calls run a PCA  and then a permutation DFA for males
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------All Male Calls----------------------------</strong>\n')
if SorKOnly==1
    MaleLogical = contains(SexDeafMic, 'HM');
elseif SorKOnly == 2
    MaleLogical = contains(SexDeafMic, 'DM');
else
    MaleLogical = contains(SexDeafMic, 'M');
end
BatIDMic_M=BatIDMic(MaleLogical);
NumInd = length(unique(BatIDMic_M));
[~,Score_M,~,~,VarExpl_M] = pca(MPS_mic_norm(MaleLogical,:));
CSVarExpl_M = cumsum(VarExpl_M);
NPC90var_M = find(cumsum(VarExpl_M)>90,1);
%% DFA for male calls
% Find the optimal regularization parameters that gives the best value
% of discrimination in cross-validation
fprintf(1, 'Male calls Finding the optimal regularization parameter (= # PC) according to cross-validated discrimination performance\n')
if (round(NPC90var_M/100)*100)>=500
    PC_val_M = [10:10:90 100:20:490 500:100:round(NPC90var_M/100)*100];
elseif (round(NPC90var_M/100)*100)>=100
    PC_val_M = [10:10:90 100:20:round(NPC90var_M/100)*100];
else
    PC_val_M = 10:10:100;
end

L_males = cell(length(PC_val_M),1);
PCC_males_mean_std = nan(length(PC_val_M),2);
for npc = 1:length(PC_val_M)
    fprintf(1, 'Male calls #PC = %d (%d/%d)\n', PC_val_M(npc), npc, length(PC_val_M))
    Mdl = fitcdiscr(Score_M(:, 1:PC_val_M(npc)),BatIDMic_M, 'CrossVal','on','KFold',10, 'Prior', 1/NumInd.*ones(1,NumInd), 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_males{npc} = kfoldLoss(Mdl, 'Mode', 'Individual');
    PCC_males_mean_std(npc,1) = mean(100*(1-L_males{npc}));
    PCC_males_mean_std(npc,2) = std(100*(1-L_males{npc}));
end
%% Determine the optimal number of PC as the first value where the...
% increase in performance  weighted by the increase of % variance
% explained by the change of number of PC is smaller than half the average standard
% deviation of the performance
Objective = diff(PCC_males_mean_std(:,1)) .* diff(CSVarExpl_M(PC_val_M));
NPC_opt_ind = find(Objective-mean(PCC_males_mean_std(:,2))/2<0, 1, 'first')+1;
NPC_opt_M = PC_val_M(NPC_opt_ind);
% Plot the figure
FIG=figure();
FIG.Position(3) = 2*FIG.Position(3);
FIG.PaperPosition(3) = 2*FIG.PaperPosition(3);
subplot(1,2,1)
shadedErrorBar(PC_val_M, PCC_males_mean_std(:,1),PCC_males_mean_std(:,2),{'-','color',ColorCode(5,:), 'LineWidth',2})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [0 100])
hold on
H=hline(1/NumInd*100, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_M, 'b--', sprintf('Optimal #PC = %d', NPC_opt_M));
V.LineWidth = 2;
%% Permutation tests for male calls
% Permutation of ID irrespective of Sex
fprintf(1, 'Male calls Permutation test')
PC_val_rand = PC_val_M(1:min(length(PC_val_M),(NPC_opt_ind+NumPCsup)));
% PC_val_rand = PC_val_M;
Lrand_males = cell(length(PC_val_rand),1);
for bb=1:NRandPerm
    fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
    RandInd = randperm(length(BatIDMic_M));
    for npc=1:length(PC_val_rand)
        fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
        if bb==1
            Lrand_males{npc} = nan(NRandPerm,10);
        end
        Mdlrand = fitcdiscr(Score_M(:,1:PC_val_rand(npc)),BatIDMic_M(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', 1/NumInd .* ones(1,NumInd), 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lrand_males{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end
end
PCCrand_males_mean_std = nan(length(PC_val_rand),2);
for npc=1:length(PC_val_rand)
    PCCrand_males_mean_std(npc,1) = mean(reshape(100*(1-Lrand_males{npc}),numel(Lrand_males{npc}),1));
    PCCrand_males_mean_std(npc,2) = std(reshape(100*(1-Lrand_males{npc}),numel(Lrand_males{npc}),1));
end

% Add the permutation values to the figure
FIG;
subplot(1,2,2)
shadedErrorBar(PC_val_rand, PCC_males_mean_std(1:length(PC_val_rand),1),PCC_males_mean_std(1:length(PC_val_rand),2),{'-','color',ColorCode(5,:), 'LineWidth',2})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [0 100])
hold on
H=hline(1/NumInd*100, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_M, 'b--', sprintf('Optimal #PC = %d', NPC_opt_M));
V.LineWidth = 2;
hold on
shadedErrorBar(PC_val_rand, PCCrand_males_mean_std(:,1),PCCrand_males_mean_std(:,2),{'-', 'Color', 'k', 'LineWidth',2})
hold off
XLim = get(gca,'XLim');
text(XLim(2)/2,97, 'Observed data', 'Color', ColorCode(5,:))
text(XLim(2)/2,93, 'Full permutation', 'Color', 'k')


% %% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
% fprintf(1, 'Male calls Getting the DF Axis\n')
% Mdl_M = fitcdiscr(Score_M(:, 1:NPC_opt_M),SexDeafMic_M, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
% % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
% % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
% %  (default) — Matrix of size p-by-p, where p is the number of predictors
% % The axis of the DFA is the largest eigenvector of inv(Sigma) *
% % BetweenSigma which in the case of 2 classes...
% % is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
% % each category
% D_mean = mean(Score_M(contains(SexDeafMic_M, 'D'),1:NPC_opt_M));
% H_mean = mean(Score_M(contains(SexDeafMic_M, 'H'),1:NPC_opt_M));
% EigenVec_DFA_M = Mdl_M.Sigma \ (D_mean - H_mean)';
% PC_DF_males = PC_M(:,1:NPC_opt_M) * EigenVec_DFA_M(:,1);
% PC_DF_males = reshape(PC_DF_males, sum(IndWf), length(MPS_mic_wt));

% fprintf(1, 'Classification performance between male K and male S calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\n#PC:%d \n', PCC_males_mean_std(NPC_opt_ind,1), PCC_males_mean_std(NPC_opt_ind,2),PCCperm_males_mean_std(NPC_opt_ind,1), PCCperm_males_mean_std(NPC_opt_ind,2),PCCrand_males_mean_std(NPC_opt_ind,1), PCCrand_males_mean_std(NPC_opt_ind,2), NPC_opt_M)
% save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFA_MaleCalls.mat'),  'PC_val_M', 'NPC_opt_M','L_males', 'PCC_males_mean_std', 'Lperm_males', 'PCCperm_males_mean_std', 'Lrand_males', 'PCCrand_males_mean_std', 'PC_DF_males', 'EigenVec_DFA_M', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append' )
fprintf(1, 'Classification performance between male ID calls: %.1f%% +/-%.1f%%\nPermutation value (irrespective of Sex): %.1f%% +/-%.1f%%\n#PC:%d \n', PCC_males_mean_std(NPC_opt_ind,1), PCC_males_mean_std(NPC_opt_ind,2),PCCrand_males_mean_std(NPC_opt_ind,1), PCCrand_males_mean_std(NPC_opt_ind,2), NPC_opt_M)
if SorKOnly==1
    save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAIDS_MaleCalls.mat'),  'PC_val_M', 'NPC_opt_M','L_males', 'PCC_males_mean_std', 'Lrand_males', 'PCCrand_males_mean_std')
elseif SorKOnly==2
    save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAIDK_MaleCalls.mat'),  'PC_val_M', 'NPC_opt_M','L_males', 'PCC_males_mean_std', 'Lrand_males', 'PCCrand_males_mean_std')
else
    save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAID_MaleCalls.mat'),  'PC_val_M', 'NPC_opt_M','L_males', 'PCC_males_mean_std', 'Lrand_males', 'PCCrand_males_mean_std')
end

% % plot the positive direction of the DF1 axis in the MPS space
% FIG;
% subplot(1,3,3)
% plot_mps(PC_DF_males, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
% title('DF1 K vs S axis for males calls')
% Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
%     flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
% colormap(Cmap)
% Axis = caxis();
% Lim = max(abs(Axis));
% caxis([-Lim Lim])
suplabel('Male Calls DFA performance ID', 't');
for cc=1:length(FIG.Children)
    FIG.Children(cc).FontSize=12;
end
if SorKOnly==1
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAIDS_MaleCalls.png') , '-dpng')
elseif SorKOnly==2
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAIDK_MaleCalls.png') , '-dpng')
else
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAID_MaleCalls.png') , '-dpng')
end

clear Score_M PCC_M Mdl_M Mdlrand
%% For All female Calls run a PCA  and then a permutation DFA for females
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------All Female Calls----------------------------</strong>\n')
if SorKOnly==1
    FemaleLogical = contains(SexDeafMic, 'HF');
elseif SorKOnly==2
    FemaleLogical = contains(SexDeafMic, 'DF');
else
    FemaleLogical = contains(SexDeafMic, 'F');
end
BatIDMic_F=BatIDMic(FemaleLogical);
NumInd = length(unique(BatIDMic_F));
[~,Score_F,~,~,VarExpl_F] = pca(MPS_mic_norm(FemaleLogical,:));
CSVarExpl_F = cumsum(VarExpl_F);
NPC90var_F = find(cumsum(VarExpl_F)>90,1);
%% DFA for female calls
% Find the optimal regularization parameters that gives the best value
% of discrimination in cross-validation
fprintf(1, 'Female calls Finding the optimal regularization parameter (= # PC) according to cross-validated discrimination performance\n')
if (round(NPC90var_F/100)*100)>=500
    PC_val_F = [10:10:90 100:20:490 500:100:round(NPC90var_F/100)*100];
elseif (round(NPC90var_F/100)*100)>=100
    PC_val_F = [10:10:90 100:20:round(NPC90var_F/100)*100];
else
    PC_val_F = 10:10:100;
end

L_females = cell(length(PC_val_F),1);
PCC_females_mean_std = nan(length(PC_val_F),2);
for npc = 1:length(PC_val_F)
    fprintf(1, 'Female calls #PC = %d (%d/%d)\n', PC_val_F(npc), npc, length(PC_val_F))
    Mdl = fitcdiscr(Score_F(:, 1:PC_val_F(npc)),BatIDMic_F, 'CrossVal','on','KFold',10, 'Prior', 1/NumInd*ones(1,NumInd), 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_females{npc} = kfoldLoss(Mdl, 'Mode', 'Individual');
    PCC_females_mean_std(npc,1) = mean(100*(1-L_females{npc}));
    PCC_females_mean_std(npc,2) = std(100*(1-L_females{npc}));
end
%% Determine the optimal number of PC as the first value where the...
% increase in performance  weighted by the increase of % variance
% explained by the change of number of PC is smaller than half the average standard
% deviation of the performance
Objective = diff(PCC_females_mean_std(:,1)) .* diff(CSVarExpl_F(PC_val_F));
NPC_opt_ind = find(Objective-mean(PCC_females_mean_std(:,2))/2<0, 1, 'first')+1;
NPC_opt_F = PC_val_F(NPC_opt_ind);
% Plot the figure
FIG=figure();
FIG.Position(3) = 2*FIG.Position(3);
FIG.PaperPosition(3) = 2*FIG.PaperPosition(3);
subplot(1,2,1)
shadedErrorBar(PC_val_F, PCC_females_mean_std(:,1),PCC_females_mean_std(:,2),{'-','color',ColorCode(5,:), 'LineWidth',2})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [0 100])
hold on
H=hline(1/NumInd*100, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_F, 'b--', sprintf('Optimal #PC = %d', NPC_opt_F));
V.LineWidth = 2;
%% Permutation tests for female calls
% Permutation of ID irrespective of Sex
fprintf(1, 'Female calls Permutation test irrespective of Sex')
PC_val_rand = PC_val_F(1:min(length(PC_val_F),(NPC_opt_ind+NumPCsup)));
% PC_val_rand = PC_val_F;
Lrand_females = cell(length(PC_val_rand),1);
for bb=1:NRandPerm
    fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
    RandInd = randperm(length(BatIDMic_F));
    for npc=1:length(PC_val_rand)
        fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
        if bb==1
            Lrand_females{npc} = nan(NRandPerm,10);
        end
        Mdlrand = fitcdiscr(Score_F(:,1:PC_val_rand(npc)),BatIDMic_F(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', 1/NumInd*ones(1,NumInd), 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lrand_females{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end
end
PCCrand_females_mean_std = nan(length(PC_val_rand),2);
for npc=1:length(PC_val_rand)
    PCCrand_females_mean_std(npc,1) = mean(reshape(100*(1-Lrand_females{npc}),numel(Lrand_females{npc}),1));
    PCCrand_females_mean_std(npc,2) = std(reshape(100*(1-Lrand_females{npc}),numel(Lrand_females{npc}),1));
end

% Add the permutation values to the figure
FIG;
subplot(1,2,2)
shadedErrorBar(PC_val_rand, PCC_females_mean_std(1:length(PC_val_rand),1),PCC_females_mean_std(1:length(PC_val_rand),2),{'-','color',ColorCode(5,:), 'LineWidth',2})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [0 100])
hold on
H=hline(1/NumInd*100, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_F, 'b--', sprintf('Optimal #PC = %d', NPC_opt_F));
V.LineWidth = 2;
hold on
shadedErrorBar(PC_val_rand, PCCrand_females_mean_std(:,1),PCCrand_females_mean_std(:,2),{'-','Color','k', 'LineWidth',2})
hold off
XLim = get(gca,'XLim');
text(XLim(2)/2,97, 'Observed data', 'Color', ColorCode(5,:))
text(XLim(2)/2,93, 'Full permutation', 'Color', 'k')


% %% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
% fprintf(1, 'Female calls Getting the DF Axis\n')
% Mdl_F = fitcdiscr(Score_F(:, 1:NPC_opt_F),SexDeafMic_F, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
% % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
% % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
% %  (default) — Matrix of size p-by-p, where p is the number of predictors
% % The axis of the DFA is the largest eigenvector of inv(Sigma) *
% % BetweenSigma which in the case of 2 classes...
% % is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
% % each category
% D_mean = mean(Score_F(contains(SexDeafMic_F, 'D'),1:NPC_opt_F));
% H_mean = mean(Score_F(contains(SexDeafMic_F, 'H'),1:NPC_opt_F));
% EigenVec_DFA_F = Mdl_F.Sigma \ (D_mean - H_mean)';
% PC_DF_females = PC_M(:,1:NPC_opt_F) * EigenVec_DFA_F(:,1);
% PC_DF_females = reshape(PC_DF_females, sum(IndWf), length(MPS_mic_wt));

% save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFA_FemaleCalls.mat'),  'PC_val_F', 'NPC_opt_F','L_females', 'PCC_females_mean_std', 'Lperm_females', 'PCCperm_females_mean_std', 'Lrand_females', 'PCCrand_females_mean_std', 'PC_DF_females','EigenVec_DFA_F', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append' )
fprintf(1, 'Classification performance between female ID calls: %.1f%% +/-%.1f%%\nPermutation value (irrespective of Sex): %.1f%% +/-%.1f%%\n#PC:%d \n', PCC_females_mean_std(NPC_opt_ind,1), PCC_females_mean_std(NPC_opt_ind,2),PCCrand_females_mean_std(NPC_opt_ind,1), PCCrand_females_mean_std(NPC_opt_ind,2), NPC_opt_F)
if SorKOnly==1
    save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAIDS_FemaleCalls.mat'),  'PC_val_F', 'NPC_opt_F','L_females', 'PCC_females_mean_std', 'Lrand_females', 'PCCrand_females_mean_std')
elseif SorKOnly==2
    save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAIDK_FemaleCalls.mat'),  'PC_val_F', 'NPC_opt_F','L_females', 'PCC_females_mean_std', 'Lrand_females', 'PCCrand_females_mean_std')
else
    save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAID_FemaleCalls.mat'),  'PC_val_F', 'NPC_opt_F','L_females', 'PCC_females_mean_std', 'Lrand_females', 'PCCrand_females_mean_std')
end

% plot the positive direction of the DF1 axis in the MPS space
% FIG;
% subplot(1,3,3)
% plot_mps(PC_DF_females, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
% title('DF1 K vs S axis for females calls')
% Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
%     flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
% colormap(Cmap)
% Axis = caxis();
% Lim = max(abs(Axis));
% caxis([-Lim Lim])
suplabel('Female Calls DFA performance ID', 't');
for cc=1:length(FIG.Children)
    FIG.Children(cc).FontSize=12;
end
if SorKOnly==1
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAIDS_FemaleCalls.png') , '-dpng')
elseif SorKOnly==2
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAIDK_FemaleCalls.png') , '-dpng')
else
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAID_FemaleCalls.png') , '-dpng')
end

clear Score_F PCC_F Mdl_F Mdlrand

%% Now same analysis for each acoustic group
UGroup = unique(TmicAll7);
for gg=1:length(UGroup)
    GR = UGroup(gg);

    %% For All Male Calls run a PCA  and then a permutation DFA for males
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-----------------------Acoustic Group %d Male Calls----------------------------</strong>\n', GR)
    if SorKOnly==1
        MaleLogical = logical(contains(SexDeafMic, 'HM').*(TmicAll7==GR));
    elseif SorKOnly==2
        MaleLogical = logical(contains(SexDeafMic, 'DM').*(TmicAll7==GR));
    else
        MaleLogical = logical(contains(SexDeafMic, 'M').*(TmicAll7==GR));
    end
    BatIDMic_M = BatIDMic(MaleLogical);
    NumInd = length(unique(BatIDMic_M));
    [~,Score_M,~, ~, VarExpl_M,~] = pca(MPS_mic_norm(MaleLogical,:));
    NPC90var_M = find(cumsum(VarExpl_M)>90,1);
    CSVarExpl_M = cumsum(VarExpl_M);

    %% DFA for male calls
    % Find the optimal regularization parameters that gives the best value
    % of discrimination in cross-validation
    fprintf(1, 'Male calls Ac Grp %d Finding the optimal regularization parameter (= # PC) according to cross-validated discrimination performance\n', GR)
    if (round(NPC90var_M/100)*100)>=500
        PC_val_M = [10:10:90 100:20:490 500:100:round(NPC90var_M/100)*100];
    elseif (round(NPC90var_M/100)*100)>=100
        PC_val_M = [10:10:90 100:20:round(NPC90var_M/100)*100];
    else
        PC_val_M = 10:10:100;
    end

    L_males = cell(length(PC_val_M),1);
    PCC_males_mean_std = nan(length(PC_val_M),2);
    for npc = 1:length(PC_val_M)
        fprintf(1, 'Male calls #PC = %d (%d/%d)\n', PC_val_M(npc), npc, length(PC_val_M))
        Mdl = fitcdiscr(Score_M(:, 1:PC_val_M(npc)),BatIDMic_M, 'CrossVal','on','KFold',10, 'Prior', ones(1,NumInd)./NumInd, 'SaveMemory', 'on', 'FillCoeffs', 'off');
        L_males{npc} = kfoldLoss(Mdl, 'Mode', 'Individual');
        PCC_males_mean_std(npc,1) = mean(100*(1-L_males{npc}));
        PCC_males_mean_std(npc,2) = std(100*(1-L_males{npc}));
    end
    %% Determine the optimal number of PC as the first value where the...
    % increase in performance  weighted by the increase of % variance
    % explained by the change of number of PC is smaller than half the average standard
    % deviation of the performance
    Objective = diff(PCC_males_mean_std(:,1)) .* diff(CSVarExpl_M(PC_val_M));
    NPC_opt_ind = find(Objective-mean(PCC_males_mean_std(:,2))/2<0, 1, 'first')+1;
    NPC_opt_M = PC_val_M(NPC_opt_ind);
    % Plot the figure
    FIG=figure();
    FIG.Position(3) = 2*FIG.Position(3);
    FIG.PaperPosition(3) = 2*FIG.PaperPosition(3);
    subplot(1,2,1)
    shadedErrorBar(PC_val_M, PCC_males_mean_std(:,1),PCC_males_mean_std(:,2),{'-','color',ColorCode(5,:), 'LineWidth',2})
    xlabel('# PC')
    ylabel('% Classification correct')
    set(gca, 'YLim', [0 100])
    hold on
    H=hline(1/NumInd*100, 'k--','Chance level');
    H.LineWidth = 2;
    hold on
    V = vline(NPC_opt_M, 'b--', sprintf('Optimal #PC = %d', NPC_opt_M));
    V.LineWidth = 2;
    %% Permutation tests for male calls
    % Permutation of ID irrespective of Sex
    fprintf(1, 'Male calls Ac Grp%d Permutation test', GR)
    PC_val_rand = PC_val_M(1:min(length(PC_val_M),(NPC_opt_ind+NumPCsup)));
    %     PC_val_rand = PC_val_M;
    Lrand_males = cell(length(PC_val_rand),1);
    for bb=1:NRandPerm
        fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
        RandInd = randperm(length(BatIDMic_M));
        for npc=1:length(PC_val_rand)
            fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
            if bb==1
                Lrand_males{npc} = nan(NRandPerm,10);
            end
            Mdlrand = fitcdiscr(Score_M(:,1:PC_val_rand(npc)),BatIDMic_M(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', ones(1,NumInd)./NumInd, 'SaveMemory', 'on', 'FillCoeffs', 'off');
            Lrand_males{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        end
    end
    PCCrand_males_mean_std = nan(length(PC_val_rand),2);
    for npc=1:length(PC_val_rand)
        PCCrand_males_mean_std(npc,1) = mean(reshape(100*(1-Lrand_males{npc}),numel(Lrand_males{npc}),1));
        PCCrand_males_mean_std(npc,2) = std(reshape(100*(1-Lrand_males{npc}),numel(Lrand_males{npc}),1));
    end

    % Add the permutation values to the figure
    FIG;
    subplot(1,2,2)
    shadedErrorBar(PC_val_rand, PCC_males_mean_std(1:length(PC_val_rand),1),PCC_males_mean_std(1:length(PC_val_rand),2),{'-','color',ColorCode(5,:), 'LineWidth',2})
    xlabel('# PC')
    ylabel('% Classification correct')
    set(gca, 'YLim', [0 100])
    hold on
    H=hline(1/NumInd*100, 'k--','Chance level');
    H.LineWidth = 2;
    hold on
    V = vline(NPC_opt_M, 'b--', sprintf('Optimal #PC = %d', NPC_opt_M));
    V.LineWidth = 2;
    hold on
    shadedErrorBar(PC_val_rand, PCCrand_males_mean_std(:,1),PCCrand_males_mean_std(:,2), {'-','Color','k', 'LineWidth',2})
    hold off
    XLim = get(gca,'XLim');
    text(XLim(2)/2,97, 'Observed data', 'Color', ColorCode(5,:))
    text(XLim(2)/2,93, 'Full permutation', 'Color', 'k')


    %     %% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
    %     fprintf(1, 'Male calls Ac Grp %d Getting the DF Axis\n', GR)
    %     Mdl_M = fitcdiscr(Score_M(:, 1:NPC_opt_M),SexDeafMic_M, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
    %     % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
    %     % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
    %     %  (default) — Matrix of size p-by-p, where p is the number of predictors
    %     % The axis of the DFA is the largest eigenvector of inv(Sigma) *
    %     % BetweenSigma which in the case of 2 classes...
    %     % is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
    %     % each category
    %     D_mean = mean(Score_M(contains(SexDeafMic_M, 'D'),1:NPC_opt_M));
    %     H_mean = mean(Score_M(contains(SexDeafMic_M, 'H'),1:NPC_opt_M));
    %     EigenVec_DFA_M = Mdl_M.Sigma \ (D_mean - H_mean)';
    %     PC_DF_males = PC_M(:,1:NPC_opt_M) * EigenVec_DFA_M(:,1);
    %     PC_DF_males = reshape(PC_DF_males, sum(IndWf), length(MPS_mic_wt));

    %     save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFA_MaleAcGroup%d.mat', GR)),  'PC_val_M', 'NPC_opt_M','L_males', 'PCC_males_mean_std', 'Lperm_males', 'PCCperm_males_mean_std', 'Lrand_males', 'PCCrand_males_mean_std', 'PC_DF_males', 'EigenVec_DFA_M', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append' )
    fprintf(1, 'Classification performance between male ID calls in Ac Grp %d: %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\n#PC:%d \n', GR, PCC_males_mean_std(NPC_opt_ind,1), PCC_males_mean_std(NPC_opt_ind,2),PCCrand_males_mean_std(NPC_opt_ind,1), PCCrand_males_mean_std(NPC_opt_ind,2), NPC_opt_M)
    if SorKOnly==1
        save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAIDS_MaleAcGroup%d.mat', GR)),  'PC_val_M', 'NPC_opt_M','L_males', 'PCC_males_mean_std', 'Lrand_males', 'PCCrand_males_mean_std')
    elseif SorKOnly==2
        save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAIDK_MaleAcGroup%d.mat', GR)),  'PC_val_M', 'NPC_opt_M','L_males', 'PCC_males_mean_std', 'Lrand_males', 'PCCrand_males_mean_std')
    else
        save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAID_MaleAcGroup%d.mat', GR)),  'PC_val_M', 'NPC_opt_M','L_males', 'PCC_males_mean_std', 'Lrand_males', 'PCCrand_males_mean_std')
    end

    %     % plot the positive direction of the DF1 axis in the MPS space
    %     FIG;
    %     subplot(1,,3)
    %     plot_mps(PC_DF_males, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
    %     title('DF1 K vs S axis for males calls')
    %     Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
    %         flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
    %     colormap(Cmap)
    %     Axis = caxis();
    %     Lim = max(abs(Axis));
    %     caxis([-Lim Lim])
    suplabel(sprintf('Acoustid Group %d Male Calls DFA performance ID', GR), 't');
    for cc=1:length(FIG.Children)
        FIG.Children(cc).FontSize=12;
    end
    if SorKOnly==1
        print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFAIDS_MaleCallsAcGrp%d.png',GR)) , '-dpng')
    elseif SorKOnly==2
        print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFAIDK_MaleCallsAcGrp%d.png',GR)) , '-dpng')
    else
        print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFAID_MaleCallsAcGrp%d.png',GR)) , '-dpng')
    end

    clear Score_M PCC_M Mdl_M Mdlrand

    %% For All female Calls run a PCA  and then a permutation DFA for females
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>---------------------Female Calls Acoustic Group %d----------------------------</strong>\n', GR)
    if SorKOnly==1
        FemaleLogical = logical(contains(SexDeafMic, 'HF').*(TmicAll7==GR));
    elseif SorKOnly==2
        FemaleLogical = logical(contains(SexDeafMic, 'DF').*(TmicAll7==GR));
    else
        FemaleLogical = logical(contains(SexDeafMic, 'F').*(TmicAll7==GR));
    end
    BatIDMic_F = BatIDMic(FemaleLogical);
    [~,Score_F,~, ~, VarExpl_F,~] = pca(MPS_mic_norm(FemaleLogical,:));
    NPC90var_F = find(cumsum(VarExpl_F)>90,1);
    CSVarExpl_F = cumsum(VarExpl_F);
    NumInd = length(unique(BatIDMic_F));
    %% DFA for female calls
    % Find the optimal regularization parameters that gives the best value
    % of discrimination in cross-validation
    fprintf(1, 'Female calls Ac Grp %d Finding the optimal regularization parameter (= # PC) according to cross-validated discrimination performance\n', GR)
    if (round(NPC90var_F/100)*100)>=500
        PC_val_F = [10:10:90 100:20:490 500:100:round(NPC90var_F/100)*100];
    elseif (round(NPC90var_F/100)*100)>=100
        PC_val_F = [10:10:90 100:20:round(NPC90var_F/100)*100];
    else
        PC_val_F = 10:10:100;
    end

    L_females = cell(length(PC_val_F),1);
    PCC_females_mean_std = nan(length(PC_val_F),2);
    for npc = 1:length(PC_val_F)
        fprintf(1, 'Female calls #PC = %d (%d/%d)\n', PC_val_F(npc), npc, length(PC_val_F))
        Mdl = fitcdiscr(Score_F(:, 1:PC_val_F(npc)),BatIDMic_F, 'CrossVal','on','KFold',10, 'Prior', ones(1,NumInd)./NumInd, 'SaveMemory', 'on', 'FillCoeffs', 'off');
        L_females{npc} = kfoldLoss(Mdl, 'Mode', 'Individual');
        PCC_females_mean_std(npc,1) = mean(100*(1-L_females{npc}));
        PCC_females_mean_std(npc,2) = std(100*(1-L_females{npc}));
    end
    %% Determine the optimal number of PC as the first value where the...
    % increase in performance  weighted by the increase of % variance
    % explained by the change of number of PC is smaller than half the average standard
    % deviation of the performance
    Objective = diff(PCC_females_mean_std(:,1)) .* diff(CSVarExpl_F(PC_val_F));
    NPC_opt_ind = find(Objective-mean(PCC_females_mean_std(:,2))/2<0, 1, 'first')+1;
    NPC_opt_F = PC_val_F(NPC_opt_ind);
    % Plot the figure
    FIG=figure();
    FIG.Position(3) = 2*FIG.Position(3);
    FIG.PaperPosition(3) = 2*FIG.PaperPosition(3);
    subplot(1,2,1)
    shadedErrorBar(PC_val_F, PCC_females_mean_std(:,1),PCC_females_mean_std(:,2),{'-','color',ColorCode(5,:), 'LineWidth',2})
    xlabel('# PC')
    ylabel('% Classification correct')
    set(gca, 'YLim', [0 100])
    hold on
    H=hline(1/NumInd*100, 'k--','Chance level');
    H.LineWidth = 2;
    hold on
    V = vline(NPC_opt_F, 'b--', sprintf('Optimal #PC = %d', NPC_opt_F));
    V.LineWidth = 2;
    %% Permutation tests for female calls
    % Permutation of ID irrespective of Sex
    fprintf(1, 'Female calls Ac Grp%d Permutation test', GR)
    PC_val_rand = PC_val_F(1:min(length(PC_val_F),(NPC_opt_ind+NumPCsup)));
    %     PC_val_rand = PC_val_F;
    Lrand_females = cell(length(PC_val_rand),1);
    for bb=1:NRandPerm
        fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
        RandInd = randperm(length(BatIDMic_F));
        for npc=1:length(PC_val_rand)
            fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
            if bb==1
                Lrand_females{npc} = nan(NRandPerm,10);
            end
            Mdlrand = fitcdiscr(Score_F(:,1:PC_val_rand(npc)),BatIDMic_F(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', ones(1,NumInd)./NumInd, 'SaveMemory', 'on', 'FillCoeffs', 'off');
            Lrand_females{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        end
    end
    PCCrand_females_mean_std = nan(length(PC_val_rand),2);
    for npc=1:length(PC_val_rand)
        PCCrand_females_mean_std(npc,1) = mean(reshape(100*(1-Lrand_females{npc}),numel(Lrand_females{npc}),1));
        PCCrand_females_mean_std(npc,2) = std(reshape(100*(1-Lrand_females{npc}),numel(Lrand_females{npc}),1));
    end

    % Add the permutation values to the figure
    FIG;
    subplot(1,2,2)
    shadedErrorBar(PC_val_rand, PCC_females_mean_std(1:length(PC_val_rand),1),PCC_females_mean_std(1:length(PC_val_rand),2),{'-','color',ColorCode(5,:), 'LineWidth',2})
    xlabel('# PC')
    ylabel('% Classification correct')
    set(gca, 'YLim', [0 100])
    hold on
    H=hline(1/NumInd*100, 'k--','Chance level');
    H.LineWidth = 2;
    hold on
    V = vline(NPC_opt_F, 'b--', sprintf('Optimal #PC = %d', NPC_opt_F));
    V.LineWidth = 2;
    hold on
    shadedErrorBar(PC_val_rand, PCCrand_females_mean_std(:,1),PCCrand_females_mean_std(:,2), {'-', 'Color', 'k', 'LineWidth',2})
    hold off
    XLim = get(gca,'XLim');
    text(XLim(2)/2,97, 'Observed data', 'Color', ColorCode(5,:))
    text(XLim(2)/2,93, 'Full permutation', 'Color', 'k')


    %     %% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
    %     fprintf(1, 'Female calls Ac Grp%d Getting the DF Axis\n', GR)
    %     Mdl_F = fitcdiscr(Score_F(:, 1:NPC_opt_F),SexDeafMic_F, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
    %     % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
    %     % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
    %     %  (default) — Matrix of size p-by-p, where p is the number of predictors
    %     % The axis of the DFA is the largest eigenvector of inv(Sigma) *
    %     % BetweenSigma which in the case of 2 classes...
    %     % is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
    %     % each category
    %     D_mean = mean(Score_F(contains(SexDeafMic_F, 'D'),1:NPC_opt_F));
    %     H_mean = mean(Score_F(contains(SexDeafMic_F, 'H'),1:NPC_opt_F));
    %     EigenVec_DFA_F = Mdl_F.Sigma \ (D_mean - H_mean)';
    %     PC_DF_females = PC_M(:,1:NPC_opt_F) * EigenVec_DFA_F(:,1);
    %     PC_DF_females = reshape(PC_DF_females, sum(IndWf), length(MPS_mic_wt));

    fprintf(1, 'Classification performance between female ID calls in Acoustic Group %d: %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\n#PC:%d \n', GR, PCC_females_mean_std(NPC_opt_ind,1), PCC_females_mean_std(NPC_opt_ind,2),PCCrand_females_mean_std(NPC_opt_ind,1), PCCrand_females_mean_std(NPC_opt_ind,2), NPC_opt_F)
    if SorKOnly==1
        save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAIDS_FemaleAcGroup%d.mat', GR)),  'PC_val_F', 'NPC_opt_F','L_females', 'PCC_females_mean_std', 'Lrand_females', 'PCCrand_females_mean_std')
    elseif SorKOnly==2
        save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAIDK_FemaleAcGroup%d.mat', GR)),  'PC_val_F', 'NPC_opt_F','L_females', 'PCC_females_mean_std', 'Lrand_females', 'PCCrand_females_mean_std')
    else
        save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAID_FemaleAcGroup%d.mat', GR)),  'PC_val_F', 'NPC_opt_F','L_females', 'PCC_females_mean_std', 'Lrand_females', 'PCCrand_females_mean_std')
    end

    % % plot the positive direction of the DF1 axis in the MPS space
    % FIG;
    % subplot(1,3,3)
    % plot_mps(PC_DF_females, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
    % title('DF1 K vs S axis for females calls')
    % Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
    %     flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
    % colormap(Cmap)
    % Axis = caxis();
    % Lim = max(abs(Axis));
    % caxis([-Lim Lim])
    suplabel(sprintf('Female Calls DFA performance ID in Ac Group%d', GR), 't');
    for cc=1:length(FIG.Children)
        FIG.Children(cc).FontSize=12;
    end
    if SorKOnly==1
        print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFAIDS_FemaleAcGroup%d.png', GR)) , '-dpng')
    elseif SorKOnly==2
        print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFAIDK_FemaleAcGroup%d.png', GR)) , '-dpng')
    else
        print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFAID_FemaleAcGroup%d.png', GR)) , '-dpng')
    end

    clear Score_F PCC_F Mdl_F Mdlrand
end
%% Obtain the confusion matrix for male calls in Ac Group 5
% rerrun de PCA on all male calls from group 5
GR=5;
MaleLogical = logical(contains(SexDeafMic, 'M').*(TmicAll7==GR));
BatIDMic_M = BatIDMic(MaleLogical);
NumInd = length(unique(BatIDMic_M));
[~,Score_M,~, ~, VarExpl_M,~] = pca(MPS_mic_norm(MaleLogical,:));
NPC90var_M = find(cumsum(VarExpl_M)>90,1);
CSVarExpl_M = cumsum(VarExpl_M);

%% DFA for male calls
% Find the optimal regularization parameters that gives the best value
% of discrimination in cross-validation
fprintf(1, 'Male calls Ac Grp %d Finding the optimal regularization parameter (= # PC) according to cross-validated discrimination performance\n', GR)
if (round(NPC90var_M/100)*100)>=500
    PC_val_M = [10:10:90 100:20:490 500:100:round(NPC90var_M/100)*100];
elseif (round(NPC90var_M/100)*100)>=100
    PC_val_M = [10:10:90 100:20:round(NPC90var_M/100)*100];
else
    PC_val_M = 10:10:100;
end

L_males = cell(length(PC_val_M),1);
PCC_males_mean_std = nan(length(PC_val_M),2);
for npc = 1:length(PC_val_M)
    fprintf(1, 'Male calls #PC = %d (%d/%d)\n', PC_val_M(npc), npc, length(PC_val_M))
    Mdl = fitcdiscr(Score_M(:, 1:PC_val_M(npc)),BatIDMic_M, 'CrossVal','on','KFold',10, 'Prior', ones(1,NumInd)./NumInd, 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_males{npc} = kfoldLoss(Mdl, 'Mode', 'Individual');
    PCC_males_mean_std(npc,1) = mean(100*(1-L_males{npc}));
    PCC_males_mean_std(npc,2) = std(100*(1-L_males{npc}));
end
%% Determine the optimal number of PC as the first value where the...
% increase in performance  weighted by the increase of % variance
% explained by the change of number of PC is smaller than half the average standard
% deviation of the performance
Objective = diff(PCC_males_mean_std(:,1)) .* diff(CSVarExpl_M(PC_val_M));
NPC_opt_ind = find(Objective-mean(PCC_males_mean_std(:,2))/2<0, 1, 'first')+1;
NPC_opt_M = PC_val_M(NPC_opt_ind);
% Plot the figure
FIG=figure();
FIG.Position(3) = 2*FIG.Position(3);
FIG.PaperPosition(3) = 2*FIG.PaperPosition(3);
subplot(1,2,1)
shadedErrorBar(PC_val_M, PCC_males_mean_std(:,1),PCC_males_mean_std(:,2),{'-','color',ColorCode(5,:), 'LineWidth',2})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [0 100])
hold on
H=hline(1/NumInd*100, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_M, 'b--', sprintf('Optimal #PC = %d', NPC_opt_M));
V.LineWidth = 2;
%% Calculate the DFA over all calls and get the confusion matrix
Mdl = fitcdiscr(Score_M(:, 1:NPC_opt_M),BatIDMic_M, 'CrossVal','off', 'Prior', ones(1,NumInd)./NumInd);
[Label,Score] = predict(Mdl,Score_M(:, 1:NPC_opt_M));
UID = unique(BatIDMic_M);
ConfMat = nan(length(UID)*ones(1,2));
for aa=1:length(UID)
    Ind = strcmp(BatIDMic_M, UID{aa});
    for pp=1:length(UID)
        ConfMat(aa,pp) = sum(strcmp(Label(Ind), UID{pp}));
    end
end
ConfpMat = ConfMat./repmat(sum(ConfMat,2),1,4);
subplot(1,2,2)
imagesc(ConfpMat)
colorbar()
colormap('gray')

%%
function [CB,Im,MPS4plot, Wt_local, Wf_local]=plot_mps(MPS, Wf, Wt, DBNOISE,CLim,Log, YLim, XLim, TwoDfilter)
if nargin<4
    DBNOISE=60;
end

if nargin<6
    Log=1;
end

if nargin<7
    YLim = [0 max(Wf*10^3)];
end
if nargin<8
    XLim = [min(Wt) max(Wt)];
end
if nargin<5
    CLim = nan;
end
if nargin<9
    TwoDfilter=0;
end

Wf_i = logical((Wf*10^3>=YLim(1)).* (Wf*10^3<=YLim(2)));
Wt_i = logical((Wt>=XLim(1)).* (Wt<=XLim(2)));
MPS = MPS(Wf_i, Wt_i);
Wf_local = Wf(Wf_i);
Wt_local = Wt(Wt_i);
if Log
    MPS4plot = 10*log10(MPS);
    MaxMPS = max(max(MPS4plot));
    MinMPS = MaxMPS-DBNOISE;
    MPS4plot(MPS4plot < MinMPS) = MinMPS;
else
    MPS4plot = MPS;
end

if TwoDfilter
    MPS4plot = imgaussfilt(MPS4plot,TwoDfilter);
end

%             imagesc(LogMPS, interpolation='nearest', aspect='auto', origin='lower', cmap=cmap, extent=ex)
Im = imagesc(MPS4plot);
axis xy
colormap(Im.Parent,'jet');
CB=colorbar();
%     CB.YTick=0;
%     CB.YTickLabel=0;
CB.Units = 'inches';
CB.Location='north';
CB.Position(4) = 0.1;
CB.Position(2) = CB.Position(2) + 0.2;
xlabel('Temporal Frequency (Hz)')
ylabel('Spectral Frequency (Cycles/kHz)')
% get nice X and y tick labels
MaxWf = max(floor(Wf_local*10^3));
YTickLabel=0:MaxWf;
YTick = nan(length(YTickLabel),1);
for yy=1:length(YTick)
    YTick(yy) = find(floor(Wf_local*10^3)==YTickLabel(yy),1);
end
set(gca, 'YTickLabel', YTickLabel, 'YTick', YTick)
MaxWt = max(ceil(Wt_local*10^-1))*10^1;
XTickLabel = [-MaxWt 0 MaxWt];
XTick = nan(length(XTickLabel),1);
for xx=1:length(XTick)
    [~,XTick(xx)] = min(abs(round(Wt_local)-XTickLabel(xx)));
end
set(gca, 'XTickLabel',XTickLabel, 'XTick', XTick)
if ~isnan(CLim)
    PlotMPS = gca;
    PlotMPS.CLim = CLim;
end
end