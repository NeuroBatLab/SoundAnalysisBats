%% Permutation DFA of K vs S for 6 acoustic groups
% here we are using now the clustering obtained with the 10 dimensions UMAP
% calculated in Notebbok 8ter instead of the clustering obtained with the 3
% dimensions UMAP calculated in Notebook 8. The optimization of the number
% of diemnsions for the UMAP clustering performance was done by looking at
% median values of Silhouette (notebook 8ter)

LocalDataDir = '/Users/elie/Documents/DeafBats/Data';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
GGPath = dir('/Users/elie/Google Drive*');
Path2Paper = fullfile(GGPath.folder, GGPath.name, 'My Drive', 'BatmanData', 'Deaf Paper');
%%
% Loading previous data

load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls2.mat'), 'CallType', 'BatID','MicAudioGood','MPS_mic', 'MPS_mic_wf', 'MPS_mic_wt')
%load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls.mat'),'MicAudioGood','TmicAll7');
load(fullfile(BaseDataDir, 'Data4_DeafBats_CatCalls_UMAPMic.mat'), 'TmicAll6');
MicAudioGood01 = MicAudioGood;
MicAudioGood01(isnan(MicAudioGood01)) = 0;
MicAudioGood01 = logical(MicAudioGood01);
%%
% Get color vectors ready

% Get the color vector ready for call type
ColorCode = [get(groot, 'DefaultAxesColorOrder'); 0 1 1; 0.5 0.5 0.5; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0 1; 0 0 0];
UCT = unique(CallType(~cellfun(@isempty, CallType)));
UCT = UCT([1:4 9 5:8 10:14]);
UCTFull = UCT;
UCTFull{strcmp(UCT, 'Tr')} = 'Trill';
UCTFull{strcmp(UCT, 'Ba')} = 'Bark';
UCTFull{strcmp(UCT, 'Pi')} = 'Pitchy Call';
UCTFull{strcmp(UCT, 'LPi')} = 'Loud Pitchy Call';
UCTFull{strcmp(UCT, 'Bu')} = 'Low Buzz';
UCTFull{strcmp(UCT, 'Pa')} = 'Panting';
UCTFull{strcmp(UCT, 'LT')} = 'Low Tuck';
UCTFull{strcmp(UCT, 'Sq')} = 'Squeal';
UCTFull{strcmp(UCT, 'Ra')} = 'Rattle';
UCTFull{strcmp(UCT, 'Ch')} = 'Chuckles';
UCTFull{strcmp(UCT, 'BB')} = 'Bark Buzz';
UCTFull{strcmp(UCT, 'PB')} = 'Pitchy Call Buzz';
UCTFull{strcmp(UCT, 'SB')} = 'Squeal Buzz';
UCTFull{strcmp(UCT, 'Un')} = 'Unknown';
CCT = nan(length(CallType),3);
IndNoCall = cellfun(@isempty, CallType);
CallType(IndNoCall) = {'Un'};
for ct=1:length(UCT)
    CCT(contains(CallType,UCT(ct)),:) = repmat(ColorCode(ct,:),sum(contains(CallType,UCT(ct))),1);
end
%%
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

%% Restrict the MPS spectral frequency to 3 cycles/kHz and normalized all
% MPS

addpath(genpath('/Users/elie/Documents/CODE/GitHub/umapDistribution'))
IndWf = (MPS_mic_wf*1000)<=3; % Max Spectral frequency taken into account in cycles/kHz
NVoc = sum(MicAudioGood01);
MPS_mic_norm = cell(1,NVoc);
MPS_mic_focus = cell(1,NVoc);
MPS_mic_mean = reshape(mean(MPS_mic(MicAudioGood01,:), 'omitnan')',length(MPS_mic_wf), length(MPS_mic_wt));
IndMicAudioGood = find(MicAudioGood01);
for nn=1:NVoc
    vv=IndMicAudioGood(nn);
    MPS_local = reshape(MPS_mic(vv,:)',length(MPS_mic_wf), length(MPS_mic_wt));
    MPS_local = MPS_local(IndWf,:);
    MPS_mic_norm{nn} = reshape(MPS_local./MPS_mic_mean(IndWf,:), numel(MPS_local),1);
    MPS_mic_focus{nn} = reshape(MPS_local, numel(MPS_local),1);
end
MPS_mic_norm = [MPS_mic_norm{:}]';
MPS_mic_focus = [MPS_mic_focus{:}]';
CCT_mic = CCT(MicAudioGood01,:);

%% For AllCalls run a permutation DFA for both male and
% females
NRandPerm=100;
PLim = 0.01;
NumPCsup = 10;

fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------All Calls----------------------------</strong>\n')

%% DFA for both male and female calls
DeafMic = Deaf(IndMicAudioGood);

% Find the optimal regularization parameters that gives the best value
% of discrimination in cross-validation
fprintf(1, 'All calls calculating the cross-validated discrimination performance\n')
PCC_all_mean_std = nan(1,2);
fprintf(1, 'All calls MPS_Norm as input\n')
Mdl = fitcdiscr(MPS_mic_norm,DeafMic, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
L_all = kfoldLoss(Mdl, 'Mode', 'Individual');
PCC_all_mean_std(1) = mean(100*(1-L_all));
PCC_all_mean_std(2) = std(100*(1-L_all));

%% Permutation tests for both male and female calls
% permutation of HD respecting ID
fprintf(1, '\nAll calls Permutation test respecting ID')
C=nchoosek(1:10,5);
Perm_all = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
PFisher = nan(size(Perm_all,1),1); % p-value of the right tail Fisher's exact test to determine if the
            % odds of correct classification is higher in the observed case
            % than the permutation case when we use the optimal number of
            % PC
BatSexDeaf4DFA = BatSexDeaf([1:3, 7:8, 4:6, 9:10]); % reorder to place Deaf bats first in the list, then hearing bats
BatName4DFA = BatName([1:3, 7:8, 4:6, 9:10]);
BatIDMicGR = BatID(IndMicAudioGood);
Lperm_all = nan(size(Perm_all,1),10);
parfor bb=1:size(Perm_all,1)
    fprintf(1, '\n Permutation %d/%d', bb, size(Perm_all,1))
    DeafMic_temp = DeafMic;
    BatSexDeaf_local = BatSexDeaf4DFA(Perm_all(bb,:));
    for bat=1:length(BatName4DFA)
        if sum(contains(BatIDMicGR,num2str(BatName4DFA(bat))))
            Ind = find(contains(BatIDMicGR,num2str(BatName4DFA(bat))));
            for ii=1:length(Ind)
                DeafMic_temp{Ind(ii)} = strtok(BatSexDeaf_local{bat});
            end
        end
    end
    Mdlrand = fitcdiscr(MPS_mic_norm, DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    Lperm_all(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    % perform a right tail Fisher's exact test to determine if the
    % odds of correct classification is higher in the observed case
    % than the permutation case
    CorrectRand = round((1-mean(Lperm_all(bb,:))).*size(MPS_mic_norm,1));
    ErrorRand = round(mean(Lperm_all(bb,:)).*size(MPS_mic_norm,1));
    Tbl = table([CorrectObs ; CorrectRand],[ErrorObs ; ErrorRand], 'VariableNames', {'Correct', 'Error'}, 'RowNames', {'Observed', 'Permutation'});
    [~,PFisher(bb),~] = fishertest(Tbl, 'Tail','right');
end

PCCperm_all_mean_std = nan(length(PC_val_rand),2);
PCCperm_all_mean_std(1) = mean(reshape(100*(1-Lperm_all),numel(Lperm_all),1));
PCCperm_all_mean_std(2) = std(reshape(100*(1-Lperm_all),numel(Lperm_all),1));
fprintf('# Significant exact Fisher test (obs vs perm) at optimal #PC: %d/%d', sum(PFisher<PLim), length(PFisher))
fprintf(1, 'Classification performance between K and S calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\n', PCC_all_mean_std(1), PCC_all_mean_std(2),PCCperm_all_mean_std(1), PCCperm_all_mean_std(2),PCCrand_all_mean_std(1), PCCrand_all_mean_std(2))
save(fullfile(LocalDataDir, 'DeafBats_NonRegularizedPermDFA_AllCalls.mat'),  'L_all', 'PCC_all_mean_std', 'Lperm_all', 'PCCperm_all_mean_std', 'PFisher')

% %% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
% fprintf(1, 'All calls Getting the DF Axis\n')
% Mdl_all = fitcdiscr(MPS_mic_norm,DeafMic, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
% % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
% % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
% %  (default) — Matrix of size p-by-p, where p is the number of predictors
% % The axis of the DFA is the largest eigenvector of inv(Sigma) *
% % BetweenSigma which in the case of 2 classes...
% % is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
% % each category
% D_mean = mean(MPS_mic_norm(contains(DeafMic, 'D'),:));
% H_mean = mean(MPS_mic_norm(contains(DeafMic, 'H'),:));
% EigenVec_DFA_all = Mdl_all.Sigma \ (D_mean - H_mean)';
% PC_DF_all = reshape(EigenVec_DFA_all, sum(IndWf), length(MPS_mic_wt));
% 
% 
% 
% % plot the positive direction of the DF1 axis in the MPS space
% FIG;
% subplot(1,3,3)
% plot_mps(PC_DF_all, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
% title('DF1 K vs S axis for all calls')
% Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
%     flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
% colormap(Cmap)
% Axis = caxis();
% Lim = max(abs(Axis));
% caxis([-Lim Lim])
% suplabel('All Calls DFA performance Saline vs Kanamycin irrespective of sex', 't');
% for cc=1:length(FIG.Children)
%     FIG.Children(cc).FontSize=12;
% end
% print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFA_AllCalls.png') , '-dpng')
% clear PC Mdl_all Mdlrand

%% For All Male Calls run a permutation DFA for males
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------All Male Calls----------------------------</strong>\n')
SexDeafMic = SexDeaf(IndMicAudioGood);
MaleLogical = contains(SexDeafMic, 'M');
SexDeafMic_M=SexDeafMic(MaleLogical);
MPS_mic_norm_local = MPS_mic_norm(MaleLogical,:);
%% DFA for male calls
% Find the optimal regularization parameters that gives the best value
% of discrimination in cross-validation
fprintf(1, 'Male calls calculating the cross-validated discrimination performance\n')
PCC_males_mean_std = nan(1,2);
Mdl = fitcdiscr(MPS_mic_norm_local,SexDeafMic_M, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
L_males = kfoldLoss(Mdl, 'Mode', 'Individual');
PCC_males_mean_std(1) = mean(100*(1-L_males{npc}));
PCC_males_mean_std(2) = std(100*(1-L_males{npc}));
%% Permutation tests for male calls
% permutation of HD respecting ID
fprintf(1, '\nMale calls Permutation test respecting ID')
C=nchoosek(1:4,2);
Perm_males = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
Perm_males = Perm_males(2:end,:); % The first row is the actual true order
Lperm_males = nan(size(Perm_males,1),10);
PFisher_males = nan(size(Perm_males,1),1);
BatSexDeaf4DFA = BatSexDeaf(7:10); % only take male bats
BatName4DFA = BatName(7:10); % only take male bats
BatIDMicGR = BatID(IndMicAudioGood(MaleLogical));

parfor bb=1:size(Perm_males,1)
    fprintf(1, '\n Permutation %d/%d', bb, size(Perm_males,1))
    DeafMic_temp = SexDeafMic_M;
    BatSexDeaf_local = BatSexDeaf4DFA(Perm_males(bb,:));
    for bat=1:length(BatName4DFA)
        if sum(contains(BatIDMicGR,num2str(BatName4DFA(bat))))
            Ind = find(contains(BatIDMicGR,num2str(BatName4DFA(bat))));
            for ii=1:length(Ind)
                DeafMic_temp{Ind(ii)} = strtok(BatSexDeaf_local{bat});
            end
        end
    end
    Mdlrand = fitcdiscr(MPS_mic_norm_local, DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    Lperm_males(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    % perform a right tail Fisher's exact test to determine if the
    % odds of correct classification is higher in the observed case
    % than the permutation case
    CorrectRand = round((1-mean(Lperm_males(bb,:))).*sum(MaleLogical));
    ErrorRand = round(mean(Lperm_males(bb,:)).*sum(MaleLogical));
    Tbl = table([CorrectObs ; CorrectRand],[ErrorObs ; ErrorRand], 'VariableNames', {'Correct', 'Error'}, 'RowNames', {'Observed', 'Permutation'});
    [~,PFisher_males(bb),~] = fishertest(Tbl, 'Tail','right');
end


PCCperm_males_mean_std(1) = mean(reshape(100*(1-Lperm_males),numel(Lperm_males),1));
PCCperm_males_mean_std(2) = std(reshape(100*(1-Lperm_males),numel(Lperm_males),1));

fprintf(1,'# Significant exact Fisher test (obs vs perm) at optimal #PC: %d/%d', sum(PFisher_males<PLim), length(PFisher_males))
fprintf(1, 'Classification performance between male K and male S calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\n', PCC_males_mean_std(NPC_opt_ind,1), PCC_males_mean_std(NPC_opt_ind,2),PCCperm_males_mean_std(NPC_opt_ind,1), PCCperm_males_mean_std(NPC_opt_ind,2))
save(fullfile(LocalDataDir, 'DeafBats_NonRegularizedPermDFA_MaleCalls.mat'), 'L_males', 'PCC_males_mean_std', 'Lperm_males','PFisher_males', 'PCCperm_males_mean_std')

clear MPS_mic_norm_local
%% For All female Calls run a permutation DFA for females
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------All Female Calls----------------------------</strong>\n')
SexDeafMic = SexDeaf(IndMicAudioGood);
FemaleLogical = contains(SexDeafMic, 'F');
SexDeafMic_F=SexDeafMic(FemaleLogical);
MPS_mic_norm_local= MPS_mic_norm(FemaleLogical,:);

%% DFA for female calls
fprintf(1, 'Female calls Calculating the cross-validated discrimination performance\n')

PCC_females_mean_std = nan(1,2);
Mdl = fitcdiscr(MPS_mic_norm_local,SexDeafMic_F, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
L_females = kfoldLoss(Mdl, 'Mode', 'Individual');
PCC_females_mean_std(1) = mean(100*(1-L_females));
PCC_females_mean_std(2) = std(100*(1-L_females));

%% Permutation tests for female calls
% permutation of HD respecting ID
fprintf(1, '\nFemale calls Permutation test respecting ID')
C=nchoosek(1:6,3);
Perm_females = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
Perm_females = Perm_females(2:end,:); % first row is in True order
Lperm_females = nan(size(Perm_females,1),10);
PFisher_females = nan(size(Perm_females,1),1);
BatSexDeaf4DFA = BatSexDeaf(1:6); % only take female bats
BatName4DFA = BatName(1:6); % only take female bats
BatIDMicGR = BatID(IndMicAudioGood(FemaleLogical));
for bb=1:size(Perm_females,1)
    fprintf(1, '\n Permutation %d/%d', bb, size(Perm_females,1))
    DeafMic_temp = SexDeafMic_F;
    BatSexDeaf_local = BatSexDeaf4DFA(Perm_females(bb,:));
    for bat=1:length(BatName4DFA)
        if sum(contains(BatIDMicGR,num2str(BatName4DFA(bat))))
            Ind = find(contains(BatIDMicGR,num2str(BatName4DFA(bat))));
            for ii=1:length(Ind)
                DeafMic_temp{Ind(ii)} = strtok(BatSexDeaf_local{bat});
            end
        end
    end
    Mdlrand = fitcdiscr(MPS_mic_norm_local, DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    Lperm_females(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    % perform a right tail Fisher's exact test to determine if the
    % odds of correct classification is higher in the observed case
    % than the permutation case
    CorrectRand = round((1-mean(Lperm_females(bb,:))).*sum(FemaleLogical));
    ErrorRand = round(mean(Lperm_females(bb,:)).*sum(FemaleLogical));
    Tbl = table([CorrectObs ; CorrectRand],[ErrorObs ; ErrorRand], 'VariableNames', {'Correct', 'Error'}, 'RowNames', {'Observed', 'Permutation'});
    [~,PFisher_females(bb),~] = fishertest(Tbl, 'Tail','right');
end

PCCperm_females_mean_std(1) = mean(reshape(100*(1-Lperm_females),numel(Lperm_females),1));
PCCperm_females_mean_std(2) = std(reshape(100*(1-Lperm_females),numel(Lperm_females),1));
fprintf(1,'# Significant exact Fisher test (obs vs perm) at optimal #PC: %d/%d', sum(PFisher_females<PLim), length(PFisher_females))
fprintf(1, 'Classification performance between female K and female S calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\n', PCC_females_mean_std(NPC_opt_ind,1), PCC_females_mean_std(NPC_opt_ind,2),PCCperm_females_mean_std(NPC_opt_ind,1), PCCperm_females_mean_std(NPC_opt_ind,2))
save(fullfile(LocalDataDir, 'DeafBats_NonRegularizedPermDFA_FemaleCalls.mat'),'L_females','PFisher_females', 'PCC_females_mean_std', 'Lperm_females', 'PCCperm_females_mean_std')

clear MPS_mic_norm_local


%% Now same analysis for each acoustic group

UGroup = unique(TmicAll6);
for gg=1:length(UGroup)
    GR = UGroup(gg);
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>------------------------Acoustic Group %d All calls----------------------------</strong>\n', GR)

    %% For All calls in that acoustic group
    MPS_mic_norm_local = MPS_mic_norm(TmicAll6==GR,:);
    %% DFA for both male and female calls
    DeafMic = Deaf(IndMicAudioGood(TmicAll6==GR));
    PCC_all_mean_std = nan(length(PC_val),2);
    fprintf(1, 'All calls Ac Grp %d \n', GR)
    Mdl = fitcdiscr(MPS_mic_norm_local,DeafMic, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_all = kfoldLoss(Mdl, 'Mode', 'Individual');
    PCC_all_mean_std(1) = mean(100*(1-L_all));
    PCC_all_mean_std(2) = std(100*(1-L_all));
    %% Permutation tests for both male and female calls
    % permutation of HD respecting ID
    fprintf(1, '\nAll calls Ac Grp %d Permutation test respecting ID', GR)
    C=nchoosek(1:10,5);
    Perm_all = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
    Lperm_all = nan(size(Perm_all,1),10);
    PFisher = nan(size(Perm_all,1),1);
    BatSexDeaf4DFA = BatSexDeaf([1:3, 7:8, 4:6, 9:10]); % reorder to place Deaf bats first in the list, then hearing bats
    BatName4DFA = BatName([1:3, 7:8, 4:6, 9:10]);
    BatIDMicGR = BatID(IndMicAudioGood(TmicAll6==GR));
    for bb=1:size(Perm_all,1)
        fprintf(1, '\n Permutation %d/%d', bb, size(Perm_all,1))
        DeafMic_temp = DeafMic;
        BatSexDeaf_local = BatSexDeaf4DFA(Perm_all(bb,:));
        for bat=1:length(BatName4DFA)
            if sum(contains(BatIDMicGR,num2str(BatName4DFA(bat))))
                Ind = find(contains(BatIDMicGR,num2str(BatName4DFA(bat))));
                for ii=1:length(Ind)
                    DeafMic_temp{Ind(ii)} = strtok(BatSexDeaf_local{bat});
                end
            end
        end
        Mdlrand = fitcdiscr(MPS_mic_norm_local, DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lperm_all(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        CorrectRand = round((1-mean(Lperm_all(bb,:))).*sum(TmicAll6==GR));
        ErrorRand = round(mean(Lperm_all(bb,:)).*sum(TmicAll6==GR));
        Tbl = table([CorrectObs ; CorrectRand],[ErrorObs ; ErrorRand], 'VariableNames', {'Correct', 'Error'}, 'RowNames', {'Observed', 'Permutation'});
        [~,PFisher(bb),~] = fishertest(Tbl, 'Tail','right');
    end

    PCCperm_all_mean_std(1) = mean(reshape(100*(1-Lperm_all),numel(Lperm_all),1));
    PCCperm_all_mean_std(2) = std(reshape(100*(1-Lperm_all),numel(Lperm_all),1));
    fprintf(1,'# Significant exact Fisher test (obs vs perm) at optimal #PC: %d/%d', sum(PFisher<PLim), length(PFisher))
    fprintf(1, 'Classification performance between K and S calls in Ac Grp %d: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\n', GR, PCC_all_mean_std(NPC_opt_ind,1), PCC_all_mean_std(NPC_opt_ind,2),PCCperm_all_mean_std(NPC_opt_ind,1), PCCperm_all_mean_std(NPC_opt_ind,2))
    save(fullfile(LocalDataDir, sprintf('DeafBats_NonRegularizedPermDFA_AcGroup%d.mat', GR)), 'L_all', 'PFisher', 'PCC_all_mean_std', 'Lperm_all', 'PCCperm_all_mean_std')

    clear MPS_mic_norm_local

    %% For All Male Calls run a PCA  and then a permutation DFA for males
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-----------------------Acoustic Group %d Male Calls----------------------------</strong>\n', GR)
    SexDeafMic = SexDeaf(IndMicAudioGood);
    MaleLogical = logical(contains(SexDeafMic, 'M').*(TmicAll6==GR));
    SexDeafMic_M=SexDeafMic(MaleLogical);
    MPS_mic_norm_local=MPS_mic_norm(MaleLogical,:);

    %% DFA for male calls
    fprintf(1, 'Male calls Ac Grp %d Calculating cross-validated discrimination performance\n', GR)
    Mdl = fitcdiscr(MPS_mic_norm_local,SexDeafMic_M, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_males = kfoldLoss(Mdl, 'Mode', 'Individual');
    PCC_males_mean_std(1) = mean(100*(1-L_males));
    PCC_males_mean_std(2) = std(100*(1-L_males));

    %% Permutation tests for male calls
    % permutation of HD respecting ID
    fprintf(1, '\nMale calls Ag Grp%d Permutation test respecting ID', GR)
    C=nchoosek(1:4,2);
    Perm_males = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
    Perm_males = Perm_males(2:end,:); % The first row is the actual true order
    Lperm_males = nan(size(Perm_males,1),10);
    PFisher_males = nan(size(Perm_males,1),1);
    BatSexDeaf4DFA = BatSexDeaf(7:10); % only take male bats
    BatName4DFA = BatName(7:10); % only take male bats
    BatIDMicGR = BatID(IndMicAudioGood(MaleLogical));
    for bb=1:size(Perm_males,1)
        fprintf(1, '\n Permutation %d/%d', bb, size(Perm_males,1))
        DeafMic_temp = SexDeafMic_M;
        BatSexDeaf_local = BatSexDeaf4DFA(Perm_males(bb,:));
        for bat=1:length(BatName4DFA)
            if sum(contains(BatIDMicGR,num2str(BatName4DFA(bat))))
                Ind = find(contains(BatIDMicGR,num2str(BatName4DFA(bat))));
                for ii=1:length(Ind)
                    DeafMic_temp{Ind(ii)} = strtok(BatSexDeaf_local{bat});
                end
            end
        end
        Mdlrand = fitcdiscr(MPC_mic_norm_local, DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lperm_males(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');

        % perform a right tail Fisher's exact test to determine if the
        % odds of correct classification is higher in the observed case
        % than the permutation case
        CorrectRand = round((1-mean(Lperm_males(bb,:))).*sum(MaleLogical));
        ErrorRand = round(mean(Lperm_males(bb,:)).*sum(MaleLogical));
        Tbl = table([CorrectObs ; CorrectRand],[ErrorObs ; ErrorRand], 'VariableNames', {'Correct', 'Error'}, 'RowNames', {'Observed', 'Permutation'});
        [~,PFisher_males(bb),~] = fishertest(Tbl, 'Tail','right');
    end

    PCCperm_males_mean_std(1) = mean(reshape(100*(1-Lperm_males),numel(Lperm_males),1));
    PCCperm_males_mean_std(2) = std(reshape(100*(1-Lperm_males),numel(Lperm_males),1));

    fprintf(1,'# Significant exact Fisher test (obs vs perm) at optimal #PC: %d/%d', sum(PFisher_males<PLim), length(PFisher_males))
    fprintf(1, 'Classification performance between male K and male S calls in Ac Grp %d: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\n', GR, PCC_males_mean_std(NPC_opt_ind,1), PCC_males_mean_std(NPC_opt_ind,2),PCCperm_males_mean_std(NPC_opt_ind,1), PCCperm_males_mean_std(NPC_opt_ind,2))
    save(fullfile(LocalDataDir, sprintf('DeafBats_NonRegularizedPermDFA_MaleAcGroup%d.mat', GR)), 'L_males', 'PFisher_males', 'PCC_males_mean_std', 'Lperm_males', 'PCCperm_males_mean_std')

    clear MPS_mic_norm_local

    %% For All female Calls run a permutation DFA for females
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>---------------------Female Calls Acoustic Group %d----------------------------</strong>\n', GR)
    SexDeafMic = SexDeaf(IndMicAudioGood);
    FemaleLogical = logical(contains(SexDeafMic, 'F').*(TmicAll6 == GR));
    SexDeafMic_F=SexDeafMic(FemaleLogical);
    MPS_mic_norm_local = MPS_mic_norm(FemaleLogical,:);


    %% DFA for female calls

    Mdl = fitcdiscr(MPS_mic_norm_local,SexDeafMic_F, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_females = kfoldLoss(Mdl, 'Mode', 'Individual');
    PCC_females_mean_std(1) = mean(100*(1-L_females));
    PCC_females_mean_std(2) = std(100*(1-L_females));

    %% Permutation tests for female calls
    % permutation of HD respecting ID
    fprintf(1, '\nFemale calls Ac Grp %d Permutation test respecting ID', GR)
    C=nchoosek(1:6,3);
    Perm_females = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
    Perm_females = Perm_females(2:end,:); % first row is in True order
    Lperm_females = nan(size(Perm_females,1),10);
    PFisher_females = nan(size(Perm_females,1),1);
    BatSexDeaf4DFA = BatSexDeaf(1:6); % only take female bats
    BatName4DFA = BatName(1:6); % only take female bats
    BatIDMicGR = BatID(IndMicAudioGood(FemaleLogical));
    for bb=1:size(Perm_females,1)
        fprintf(1, '\n Permutation %d/%d', bb, size(Perm_females,1))
        DeafMic_temp = SexDeafMic_F;
        BatSexDeaf_local = BatSexDeaf4DFA(Perm_females(bb,:));
        for bat=1:length(BatName4DFA)
            if sum(contains(BatIDMicGR,num2str(BatName4DFA(bat))))
                Ind = find(contains(BatIDMicGR,num2str(BatName4DFA(bat))));
                for ii=1:length(Ind)
                    DeafMic_temp{Ind(ii)} = strtok(BatSexDeaf_local{bat});
                end
            end
        end
        Mdlrand = fitcdiscr(MPS_mic_norm_local, DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lperm_females(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        % perform a right tail Fisher's exact test to determine if the
        % odds of correct classification is higher in the observed case
        % than the permutation case
        CorrectRand = round((1-mean(Lperm_females(bb,:))).*sum(FemaleLogical));
        ErrorRand = round(mean(Lperm_females(bb,:)).*sum(FemaleLogical));
        Tbl = table([CorrectObs ; CorrectRand],[ErrorObs ; ErrorRand], 'VariableNames', {'Correct', 'Error'}, 'RowNames', {'Observed', 'Permutation'});
        [~,PFisher_females(bb),~] = fishertest(Tbl, 'Tail','right');

    end

    PCCperm_females_mean_std(1) = mean(reshape(100*(1-Lperm_females),numel(Lperm_females),1));
    PCCperm_females_mean_std(2) = std(reshape(100*(1-Lperm_females),numel(Lperm_females),1));


    fprintf(1,'# Significant exact Fisher test (obs vs perm) at optimal #PC: %d/%d', sum(PFisher_females<PLim), length(PFisher_females))
    fprintf(1, 'Classification performance between female K and female S calls in Acoustic Group %d: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\n', GR, PCC_females_mean_std(NPC_opt_ind,1), PCC_females_mean_std(NPC_opt_ind,2),PCCperm_females_mean_std(NPC_opt_ind,1), PCCperm_females_mean_std(NPC_opt_ind,2))
    save(fullfile(LocalDataDir, sprintf('DeafBats_NonRegularizedPermDFA_FemaleAcGroup%d.mat', GR)),  'L_females', 'PFisher_females', 'PCC_females_mean_std', 'Lperm_females', 'PCCperm_females_mean_std')

    clear MPS_mic_norm_local
end
