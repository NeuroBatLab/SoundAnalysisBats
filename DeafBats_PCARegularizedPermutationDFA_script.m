%% Permutation DFA of K vs S for 7 acoustic groups

LocalDataDir = '/Users/elie/Documents/DeafBats/Data';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
%%
% Loading previous data

load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls2.mat'), 'CallType', 'BatID','MPS_mic', 'MPS_mic_wf', 'MPS_mic_wt')
load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls.mat'),'MicAudioGood','TmicAll7');
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

%% For AllCalls run a PCA  and then a permutation DFA for both male and
% females
NRandPerm=10;
NumPCsup = 5;

fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------All Calls----------------------------</strong>\n')
[PC,Score,~, ~, VarExpl,~] = pca(MPS_mic_norm);

% Plot the % variance explained by the PC

figure()
CSVarExpl = cumsum(VarExpl);
plot(CSVarExpl, 'Linewidth',2)
xlabel('# PC')
ylabel('%variance explained for All Calls')
NPC90var = find(cumsum(VarExpl)>90,1);
text(100,95, sprintf('90%% variance explained with %d PC',NPC90var ))
text(100,90, sprintf('%.1f%% variance explained with 100 PC', CSVarExpl(100)))
save(fullfile(LocalDataDir, 'DeafBats_PCA_AllCalls.mat'), 'PC', 'Score', 'VarExpl', 'NPC90var')
save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFA_AllCalls.mat'),  'NPC90var')

%% DFA for both male and female calls
DeafMic = Deaf(IndMicAudioGood);

% Find the optimal regularization parameters that gives the best value
% of discrimination in cross-validation
fprintf(1, 'All calls Finding the optimal regularization parameter (= # PC) according to cross-validated discrimination performance\n')
if (round(NPC90var/100)*100)>=500
    PC_val = [10:10:90 100:20:490 500:100:round(NPC90var/100)*100];
elseif (round(NPC90var/100)*100)>=100
    PC_val = [10:10:90 100:20:round(NPC90var/100)*100];
else
    PC_val = 10:10:100;
end

L_all = cell(length(PC_val),1);
PCC_all_mean_std = nan(length(PC_val),2);
for npc = 1:length(PC_val)
    fprintf(1, 'All calls #PC = %d (%d/%d)\n', PC_val(npc), npc, length(PC_val))
    Mdl = fitcdiscr(Score(:, 1:PC_val(npc)),DeafMic, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_all{npc} = kfoldLoss(Mdl, 'Mode', 'Individual');
    PCC_all_mean_std(npc,1) = mean(100*(1-L_all{npc}));
    PCC_all_mean_std(npc,2) = std(100*(1-L_all{npc}));
end
%% Determine the optimal number of PC as the first value where the...
% increase in performance  weighted by the increase of % variance
% explained by the change of number of PC is smaller than half the average standard
% deviation of the performance
Objective = diff(PCC_all_mean_std(:,1)) .* diff(CSVarExpl(PC_val));
NPC_opt_ind = find(Objective-mean(PCC_all_mean_std(:,2))/2<0, 1, 'first')+1;
NPC_opt_all = PC_val(NPC_opt_ind);
% Plot the figure
FIG=figure();
FIG.Position(3) = 2*FIG.Position(3);
FIG.PaperPosition(3) = 2*FIG.PaperPosition(3);
subplot(1,3,1)
shadedErrorBar(PC_val, PCC_all_mean_std(:,1),PCC_all_mean_std(:,2),{'-','color',ColorCode(5,:)})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [40 100])
hold on
H=hline(50, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_all, 'b--', sprintf('Optimal #PC = %d', NPC_opt_all));
V.LineWidth = 2;
%% Permutation tests for both male and female calls
% Permutation of HD irrespective of ID
fprintf(1, 'All calls Permutation test irrespective of ID')
%PC_val_rand = PC_val(1:(NPC_opt_ind+NumPCsup));
PC_val_rand = PC_val;
Lrand_all = cell(length(PC_val_rand),1);
for bb=1:NRandPerm
    fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
    RandInd = randperm(length(DeafMic));
    for npc=1:length(PC_val_rand)
        fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
        if bb==1
            Lrand_all{npc} = nan(NRandPerm,10);
        end
        Mdlrand = fitcdiscr(Score(:,1:PC_val_rand(npc)),DeafMic(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lrand_all{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end
end
PCCrand_all_mean_std = nan(length(PC_val_rand),2);
for npc=1:length(PC_val_rand)
    PCCrand_all_mean_std(npc,1) = mean(reshape(100*(1-Lrand_all{npc}),numel(Lrand_all{npc}),1));
    PCCrand_all_mean_std(npc,2) = std(reshape(100*(1-Lrand_all{npc}),numel(Lrand_all{npc}),1));
end


% permutation of HD respecting ID
fprintf(1, '\nAll calls Permutation test respecting ID')
C=nchoosek(1:10,5);
Perm_all = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
Lperm_all = cell(length(PC_val_rand),1);

BatSexDeaf4DFA = BatSexDeaf([1:3, 7:8, 4:6, 9:10]); % reorder to place Deaf bats first in the list, then hearing bats
BatName4DFA = BatName([1:3, 7:8, 4:6, 9:10]);
BatIDMicGR = BatID(IndMicAudioGood);
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
    for npc=1:length(PC_val_rand)
        fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
        if bb==1
            Lperm_all{npc} = nan(size(Perm_all,1),10);
        end
        Mdlrand = fitcdiscr(Score(:,1:PC_val_rand(npc)), DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lperm_all{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end
end

PCCperm_all_mean_std = nan(length(PC_val_rand),2);
for npc=1:length(PC_val_rand)
    PCCperm_all_mean_std(npc,1) = mean(reshape(100*(1-Lperm_all{npc}),numel(Lperm_all{npc}),1));
    PCCperm_all_mean_std(npc,2) = std(reshape(100*(1-Lperm_all{npc}),numel(Lperm_all{npc}),1));
end

% Add the permutation values to the figure
FIG;
subplot(1,3,2)
shadedErrorBar(PC_val_rand, PCC_all_mean_std(1:(NPC_opt_ind+1),1),PCC_all_mean_std(1:(NPC_opt_ind+1),2),{'-','color',ColorCode(5,:)})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [40 100])
hold on
H=hline(50, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_all, 'b--', sprintf('Optimal #PC = %d', NPC_opt_all));
V.LineWidth = 2;
hold on
shadedErrorBar(PC_val_rand, PCCrand_all_mean_std(:,1),PCCrand_all_mean_std(:,2))
hold on
shadedErrorBar(PC_val_rand, PCCperm_all_mean_std(:,1),PCCperm_all_mean_std(:,2))
hold off

text(10,55, 'Full permutation')
text(10,77, 'Permutation respecting ID')


%% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
fprintf(1, 'All calls Getting the DF Axis\n')
Mdl_all = fitcdiscr(Score(:, 1:NPC_opt_all),DeafMic, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
% Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
% Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
%  (default) — Matrix of size p-by-p, where p is the number of predictors
% The axis of the DFA is the largest eigenvector of inv(Sigma) *
% BetweenSigma which in the case of 2 classes...
% is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
% each category
D_mean = mean(Score(contains(DeafMic, 'D'),1:NPC_opt_all));
H_mean = mean(Score(contains(DeafMic, 'H'),1:NPC_opt_all));
EigenVec_DFA_all = Mdl_all.Sigma \ (D_mean - H_mean)';
PC_DF_all = PC(:,1:NPC_opt_all) * EigenVec_DFA_all(:,1);
PC_DF_all = reshape(PC_DF_all, sum(IndWf), length(MPS_mic_wt));

fprintf(1, 'Classification performance between K and S calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\n#PC:%d \n', PCC_all_mean_std(NPC_opt_ind,1), PCC_all_mean_std(NPC_opt_ind,2),PCCperm_all_mean_std(NPC_opt_ind,1), PCCperm_all_mean_std(NPC_opt_ind,2),PCCrand_all_mean_std(NPC_opt_ind,1), PCCrand_all_mean_std(NPC_opt_ind,2), NPC_opt_all)
save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFA_AllCalls.mat'),  'PC_val', 'NPC_opt_all','L_all', 'PCC_all_mean_std', 'Lperm_all', 'PCCperm_all_mean_std', 'Lrand_all', 'PCCrand_all_mean_std', 'PC_DF_all', 'EigenVec_DFA_all', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append' )

% plot the positive direction of the DF1 axis in the MPS space
FIG;
subplot(1,3,3)
plot_mps(PC_DF_all, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
title('DF1 K vs S axis for all calls')
Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
    flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
colormap(Cmap)
Axis = caxis();
Lim = max(abs(Axis));
caxis([-Lim Lim])
suplabel('All Calls DFA performance Saline vs Kanamycin irrespective of sex', 't')
print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFA_AllCalls.png') , '-dpng')
clear Score PC Mdl_all Mdlrand
%% For All Male Calls run a PCA  and then a permutation DFA for males
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------All Male Calls----------------------------</strong>\n')
SexDeafMic = SexDeaf(IndMicAudioGood);
MaleLogical = contains(SexDeafMic, 'M');
SexDeafMic_M=SexDeafMic(MaleLogical);
[PC_M,Score_M,~, ~, VarExpl_M,~] = pca(MPS_mic_norm(MaleLogical,:));

% Plot the % variance explained by the PC

figure()
CSVarExpl_M = cumsum(VarExpl_M);
plot(CSVarExpl_M, 'Linewidth',2)
xlabel('# PC')
ylabel('%variance explained for All Male Calls')
NPC90var_M = find(cumsum(VarExpl_M)>90,1);
text(100,95, sprintf('90%% variance explained with %d PC',NPC90var_M ))
text(100,90, sprintf('%.1f%% variance explained with 100 PC', CSVarExpl_M(100)))
save(fullfile(LocalDataDir, 'DeafBats_PCA_MaleCalls.mat'), 'PC_M', 'VarExpl_M', 'NPC90var_M')
save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFA_MaleCalls.mat'),  'NPC90var_M')

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
    Mdl = fitcdiscr(Score_M(:, 1:PC_val_M(npc)),SexDeafMic_M, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
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
subplot(1,3,1)
shadedErrorBar(PC_val_M, PCC_males_mean_std(:,1),PCC_males_mean_std(:,2),{'-','color',ColorCode(5,:)})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [40 100])
hold on
H=hline(50, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_M, 'b--', sprintf('Optimal #PC = %d', NPC_opt_M));
V.LineWidth = 2;
%% Permutation tests for male calls
% Permutation of HD irrespective of ID
fprintf(1, 'Male calls Permutation test irrespective of ID')
% PC_val_rand = PC_val_M(1:(NPC_opt_ind+NumPCsup));
PC_val_rand = PC_val_M;
Lrand_males = cell(length(PC_val_rand),1);
for bb=1:NRandPerm
    fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
    RandInd = randperm(length(SexDeafMic_M));
    for npc=1:length(PC_val_rand)
        fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
        if bb==1
            Lrand_males{npc} = nan(NRandPerm,10);
        end
        Mdlrand = fitcdiscr(Score_M(:,1:PC_val_rand(npc)),SexDeafMic_M(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lrand_males{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end
end
PCCrand_males_mean_std = nan(length(PC_val_rand),2);
for npc=1:length(PC_val_rand)
    PCCrand_males_mean_std(npc,1) = mean(reshape(100*(1-Lrand_males{npc}),numel(Lrand_males{npc}),1));
    PCCrand_males_mean_std(npc,2) = std(reshape(100*(1-Lrand_males{npc}),numel(Lrand_males{npc}),1));
end


% permutation of HD respecting ID
fprintf(1, '\nMale calls Permutation test respecting ID')
C=nchoosek(1:4,2);
Perm_males = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
Perm_males = Perm_males(2:end,:); % The first row is the actual true order
Lperm_males = cell(length(PC_val_rand),1);

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
    for npc=1:length(PC_val_rand)
        fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
        if bb==1
            Lperm_males{npc} = nan(size(Perm_males,1),10);
        end
        Mdlrand = fitcdiscr(Score_M(:,1:PC_val_rand(npc)), DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lperm_males{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end
end

PCCperm_males_mean_std = nan(length(PC_val_rand),2);
for npc=1:length(PC_val_rand)
    PCCperm_males_mean_std(npc,1) = mean(reshape(100*(1-Lperm_males{npc}),numel(Lperm_males{npc}),1));
    PCCperm_males_mean_std(npc,2) = std(reshape(100*(1-Lperm_males{npc}),numel(Lperm_males{npc}),1));
end

% Add the permutation values to the figure
FIG;
subplot(1,3,2)
shadedErrorBar(PC_val_rand, PCC_males_mean_std(1:(NPC_opt_ind+1),1),PCC_males_mean_std(1:(NPC_opt_ind+1),2),{'-','color',ColorCode(5,:)})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [40 100])
hold on
H=hline(50, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_M, 'b--', sprintf('Optimal #PC = %d', NPC_opt_M));
V.LineWidth = 2;
hold on
shadedErrorBar(PC_val_rand, PCCrand_males_mean_std(:,1),PCCrand_males_mean_std(:,2))
hold on
shadedErrorBar(PC_val_rand, PCCperm_males_mean_std(:,1),PCCperm_males_mean_std(:,2))
hold off

text(10,55, 'Full permutation')
text(10,77, 'Permutation respecting ID')


%% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
fprintf(1, 'Male calls Getting the DF Axis\n')
Mdl_M = fitcdiscr(Score_M(:, 1:NPC_opt_M),SexDeafMic_M, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
% Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
% Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
%  (default) — Matrix of size p-by-p, where p is the number of predictors
% The axis of the DFA is the largest eigenvector of inv(Sigma) *
% BetweenSigma which in the case of 2 classes...
% is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
% each category
D_mean = mean(Score_M(contains(SexDeafMic_M, 'D'),1:NPC_opt_M));
H_mean = mean(Score_M(contains(SexDeafMic_M, 'H'),1:NPC_opt_M));
EigenVec_DFA_M = Mdl_M.Sigma \ (D_mean - H_mean)';
PC_DF_males = PC_M(:,1:NPC_opt_M) * EigenVec_DFA_M(:,1);
PC_DF_males = reshape(PC_DF_males, sum(IndWf), length(MPS_mic_wt));

fprintf(1, 'Classification performance between male K and male S calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\n#PC:%d \n', PCC_males_mean_std(NPC_opt_ind,1), PCC_males_mean_std(NPC_opt_ind,2),PCCperm_males_mean_std(NPC_opt_ind,1), PCCperm_males_mean_std(NPC_opt_ind,2),PCCrand_males_mean_std(NPC_opt_ind,1), PCCrand_males_mean_std(NPC_opt_ind,2), NPC_opt_M)
save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFA_MaleCalls.mat'),  'PC_val_M', 'NPC_opt_M','L_males', 'PCC_males_mean_std', 'Lperm_males', 'PCCperm_males_mean_std', 'Lrand_males', 'PCCrand_males_mean_std', 'PC_DF_males', 'EigenVec_DFA_M', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append' )

% plot the positive direction of the DF1 axis in the MPS space
FIG;
subplot(1,3,3)
plot_mps(PC_DF_males, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
title('DF1 K vs S axis for males calls')
Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
    flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
colormap(Cmap)
Axis = caxis();
Lim = max(abs(Axis));
caxis([-Lim Lim])
suplabel('Male Calls DFA performance Saline vs Kanamycin', 't')
print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFA_MaleCalls.png') , '-dpng')

clear Score_M PCC_M Mdl_M Mdlrand
%% For All female Calls run a PCA  and then a permutation DFA for females
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------All Female Calls----------------------------</strong>\n')
SexDeafMic = SexDeaf(IndMicAudioGood);
FemaleLogical = contains(SexDeafMic, 'F');
SexDeafMic_F=SexDeafMic(FemaleLogical);
[PC_F,Score_F,~, ~, VarExpl_F,~] = pca(MPS_mic_norm(FemaleLogical,:));

% Plot the % variance explained by the PC

figure()
CSVarExpl_F = cumsum(VarExpl_F);
plot(CSVarExpl_F, 'Linewidth',2)
xlabel('# PC')
ylabel('%variance explained for All Female Calls')
NPC90var_F = find(cumsum(VarExpl_F)>90,1);
text(100,95, sprintf('90%% variance explained with %d PC',NPC90var_F ))
text(100,90, sprintf('%.1f%% variance explained with 100 PC', CSVarExpl_F(100)))
save(fullfile(LocalDataDir, 'DeafBats_PCA_FemaleCalls.mat'), 'PC_F', 'VarExpl_F', 'NPC90var_F')
save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFA_FemaleCalls.mat'),  'NPC90var_F')

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
    Mdl = fitcdiscr(Score_F(:, 1:PC_val_F(npc)),SexDeafMic_F, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
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
subplot(1,3,1)
shadedErrorBar(PC_val_F, PCC_females_mean_std(:,1),PCC_females_mean_std(:,2),{'-','color',ColorCode(5,:)})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [40 100])
hold on
H=hline(50, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_F, 'b--', sprintf('Optimal #PC = %d', NPC_opt_F));
V.LineWidth = 2;
%% Permutation tests for female calls
% Permutation of HD irrespective of ID
fprintf(1, 'Female calls Permutation test irrespective of ID')
%PC_val_rand = PC_val_F(1:(NPC_opt_ind+NumPCsup));
PC_val_rand = PC_val_F;
Lrand_females = cell(length(PC_val_rand),1);
for bb=1:NRandPerm
    fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
    RandInd = randperm(length(SexDeafMic_F));
    for npc=1:length(PC_val_rand)
        fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
        if bb==1
            Lrand_females{npc} = nan(NRandPerm,10);
        end
        Mdlrand = fitcdiscr(Score_F(:,1:PC_val_rand(npc)),SexDeafMic_F(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lrand_females{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end
end
PCCrand_females_mean_std = nan(length(PC_val_rand),2);
for npc=1:length(PC_val_rand)
    PCCrand_females_mean_std(npc,1) = mean(reshape(100*(1-Lrand_females{npc}),numel(Lrand_females{npc}),1));
    PCCrand_females_mean_std(npc,2) = std(reshape(100*(1-Lrand_females{npc}),numel(Lrand_females{npc}),1));
end


% permutation of HD respecting ID
fprintf(1, '\nFemale calls Permutation test respecting ID')
C=nchoosek(1:6,3);
Perm_females = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
Perm_females = Perm_females(2:end,:); % first row is in True order
Lperm_females = cell(length(PC_val_rand),1);

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
    for npc=1:length(PC_val_rand)
        fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
        if bb==1
            Lperm_females{npc} = nan(size(Perm_females,1),10);
        end
        Mdlrand = fitcdiscr(Score_F(:,1:PC_val_rand(npc)), DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lperm_females{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end
end

PCCperm_females_mean_std = nan(length(PC_val_rand),2);
for npc=1:length(PC_val_rand)
    PCCperm_females_mean_std(npc,1) = mean(reshape(100*(1-Lperm_females{npc}),numel(Lperm_females{npc}),1));
    PCCperm_females_mean_std(npc,2) = std(reshape(100*(1-Lperm_females{npc}),numel(Lperm_females{npc}),1));
end

% Add the permutation values to the figure
FIG;
subplot(1,3,2)
shadedErrorBar(PC_val_rand, PCC_females_mean_std(1:(NPC_opt_ind+1),1),PCC_females_mean_std(1:(NPC_opt_ind+1),2),{'-','color',ColorCode(5,:)})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [40 100])
hold on
H=hline(50, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_F, 'b--', sprintf('Optimal #PC = %d', NPC_opt_F));
V.LineWidth = 2;
hold on
shadedErrorBar(PC_val_rand, PCCrand_females_mean_std(:,1),PCCrand_females_mean_std(:,2))
hold on
shadedErrorBar(PC_val_rand, PCCperm_females_mean_std(:,1),PCCperm_females_mean_std(:,2))
hold off

text(10,55, 'Full permutation')
text(10,77, 'Permutation respecting ID')


%% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
fprintf(1, 'Female calls Getting the DF Axis\n')
Mdl_F = fitcdiscr(Score_F(:, 1:NPC_opt_F),SexDeafMic_F, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
% Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
% Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
%  (default) — Matrix of size p-by-p, where p is the number of predictors
% The axis of the DFA is the largest eigenvector of inv(Sigma) *
% BetweenSigma which in the case of 2 classes...
% is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
% each category
D_mean = mean(Score_F(contains(SexDeafMic_F, 'D'),1:NPC_opt_F));
H_mean = mean(Score_F(contains(SexDeafMic_F, 'H'),1:NPC_opt_F));
EigenVec_DFA_F = Mdl_F.Sigma \ (D_mean - H_mean)';
PC_DF_females = PC_M(:,1:NPC_opt_F) * EigenVec_DFA_F(:,1);
PC_DF_females = reshape(PC_DF_females, sum(IndWf), length(MPS_mic_wt));

fprintf(1, 'Classification performance between female K and female S calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\n#PC:%d \n', PCC_females_mean_std(NPC_opt_ind,1), PCC_females_mean_std(NPC_opt_ind,2),PCCperm_females_mean_std(NPC_opt_ind,1), PCCperm_females_mean_std(NPC_opt_ind,2),PCCrand_females_mean_std(NPC_opt_ind,1), PCCrand_females_mean_std(NPC_opt_ind,2), NPC_opt_F)
save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFA_FemaleCalls.mat'),  'PC_val_F', 'NPC_opt_F','L_females', 'PCC_females_mean_std', 'Lperm_females', 'PCCperm_females_mean_std', 'Lrand_females', 'PCCrand_females_mean_std', 'PC_DF_females','EigenVec_DFA_F', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append' )

% plot the positive direction of the DF1 axis in the MPS space
FIG;
subplot(1,3,3)
plot_mps(PC_DF_females, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
title('DF1 K vs S axis for females calls')
Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
    flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
colormap(Cmap)
Axis = caxis();
Lim = max(abs(Axis));
caxis([-Lim Lim])
suplabel('Female Calls DFA performance Saline vs Kanamycin', 't')
print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFA_FemaleCalls.png') , '-dpng')

clear Score_F PCC_F Mdl_F Mdlrand


%% Now same analysis for each acoustic group

NumPCsup = 5;
UGroup = unique(TmicAll7);
for gg=1:length(UGroup)
    GR = UGroup(gg);
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>------------------------Acoustic Group %d All calls----------------------------</strong>\n', GR)
    
    %% For All calls PCA for all calls in that acoustic group
    [PC,Score,~, ~, VarExpl,~] = pca(MPS_mic_norm(TmicAll7==GR,:));

    % Plot the % variance explained by the PC
    figure()
    CSVarExpl = cumsum(VarExpl);
    plot(CSVarExpl, 'Linewidth',2)
    xlabel('# PC')
    ylabel(sprintf('%%variance explained for All Calls Acoustic group %d', GR))
    NPC90var = find(cumsum(VarExpl)>90,1);
    text(100,95, sprintf('90%% variance explained with %d PC',NPC90var ))
    text(100,90, sprintf('%.1f%% variance explained with 100 PC', CSVarExpl(100)))
    save(fullfile(LocalDataDir, sprintf('DeafBats_PCA_AcGroup%d.mat', GR)), 'PC', 'Score', 'VarExpl', 'NPC90var')
    save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFA_AcGroup%d.mat', GR)),  'NPC90var')

    %% DFA for both male and female calls
    DeafMic = Deaf(IndMicAudioGood(TmicAll7==GR));

    % Find the optimal regularization parameters that gives the best value
    % of discrimination in cross-validation
    fprintf(1, 'All calls Acsoutic Group %d Finding the optimal regularization parameter (= # PC) according to cross-validated discrimination performance\n', GR)
    if (round(NPC90var/100)*100)>=500
        PC_val = [10:10:90 100:20:490 500:100:round(NPC90var/100)*100];
    elseif (round(NPC90var/100)*100)>=100
        PC_val = [10:10:90 100:20:round(NPC90var/100)*100];
    else
        PC_val = 10:10:100;
    end

    L_all = cell(length(PC_val),1);
    PCC_all_mean_std = nan(length(PC_val),2);
    for npc = 1:length(PC_val)
        fprintf(1, 'All calls Ac Grp %d #PC = %d (%d/%d)\n', GR, PC_val(npc), npc, length(PC_val))
        Mdl = fitcdiscr(Score(:, 1:PC_val(npc)),DeafMic, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        L_all{npc} = kfoldLoss(Mdl, 'Mode', 'Individual');
        PCC_all_mean_std(npc,1) = mean(100*(1-L_all{npc}));
        PCC_all_mean_std(npc,2) = std(100*(1-L_all{npc}));
    end
    %% Determine the optimal number of PC as the first value where the...
    % increase in performance  weighted by the increase of % variance
    % explained by the change of number of PC is smaller than half the average standard
    % deviation of the performance
    Objective = diff(PCC_all_mean_std(:,1)) .* diff(CSVarExpl(PC_val));
    NPC_opt_ind = find(Objective-mean(PCC_all_mean_std(:,2))/2<0, 1, 'first')+1;
    NPC_opt_all = PC_val(NPC_opt_ind);
    % Plot the figure
    FIG=figure();
    FIG.Position(3) = 2*FIG.Position(3);
    FIG.PaperPosition(3) = 2*FIG.PaperPosition(3);
    subplot(1,3,1)
    shadedErrorBar(PC_val, PCC_all_mean_std(:,1),PCC_all_mean_std(:,2),{'-','color',ColorCode(5,:)})
    xlabel('# PC')
    ylabel('% Classification correct')
    set(gca, 'YLim', [40 100])
    hold on
    H=hline(50, 'k--','Chance level');
    H.LineWidth = 2;
    hold on
    V = vline(NPC_opt_all, 'b--', sprintf('Optimal #PC = %d', NPC_opt_all));
    V.LineWidth = 2;
    %% Permutation tests for both male and female calls
    % Permutation of HD irrespective of ID
    fprintf(1, 'All calls Ac Grp %d Permutation test irrespective of ID', GR)
    %PC_val_rand = PC_val(1:(NPC_opt_ind+NumPCsup));
    PC_val_rand = PC_val;
    Lrand_all = cell(length(PC_val_rand),1);
    for bb=1:NRandPerm
        fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
        RandInd = randperm(length(DeafMic));
        for npc=1:length(PC_val_rand)
            fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
            if bb==1
                Lrand_all{npc} = nan(NRandPerm,10);
            end
            Mdlrand = fitcdiscr(Score(:,1:PC_val_rand(npc)),DeafMic(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
            Lrand_all{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        end
    end
    PCCrand_all_mean_std = nan(length(PC_val_rand),2);
    for npc=1:length(PC_val_rand)
        PCCrand_all_mean_std(npc,1) = mean(reshape(100*(1-Lrand_all{npc}),numel(Lrand_all{npc}),1));
        PCCrand_all_mean_std(npc,2) = std(reshape(100*(1-Lrand_all{npc}),numel(Lrand_all{npc}),1));
    end


    % permutation of HD respecting ID
    fprintf(1, '\nAll calls Ac Grp %d Permutation test respecting ID', GR)
    C=nchoosek(1:10,5);
    Perm_all = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
    Lperm_all = cell(length(PC_val_rand),1);

    BatSexDeaf4DFA = BatSexDeaf([1:3, 7:8, 4:6, 9:10]); % reorder to place Deaf bats first in the list, then hearing bats
    BatName4DFA = BatName([1:3, 7:8, 4:6, 9:10]);
    BatIDMicGR = BatID(IndMicAudioGood(TmicAll7==GR));
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
        for npc=1:length(PC_val_rand)
            fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
            if bb==1
                Lperm_all{npc} = nan(size(Perm_all,1),10);
            end
            Mdlrand = fitcdiscr(Score(:,1:PC_val_rand(npc)), DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
            Lperm_all{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        end
    end

    PCCperm_all_mean_std = nan(length(PC_val_rand),2);
    for npc=1:length(PC_val_rand)
        PCCperm_all_mean_std(npc,1) = mean(reshape(100*(1-Lperm_all{npc}),numel(Lperm_all{npc}),1));
        PCCperm_all_mean_std(npc,2) = std(reshape(100*(1-Lperm_all{npc}),numel(Lperm_all{npc}),1));
    end

    % Add the permutation values to the figure
    FIG;
    subplot(1,3,2)
    shadedErrorBar(PC_val_rand, PCC_all_mean_std(1:(NPC_opt_ind+NumPCsup),1),PCC_all_mean_std(1:(NPC_opt_ind+NumPCsup),2),{'-','color',ColorCode(5,:)})
    xlabel('# PC')
    ylabel('% Classification correct')
    set(gca, 'YLim', [40 100])
    hold on
    H=hline(50, 'k--','Chance level');
    H.LineWidth = 2;
    hold on
    V = vline(NPC_opt_all, 'b--', sprintf('Optimal #PC = %d', NPC_opt_all));
    V.LineWidth = 2;
    hold on
    shadedErrorBar(PC_val_rand, PCCrand_all_mean_std(:,1),PCCrand_all_mean_std(:,2))
    hold on
    shadedErrorBar(PC_val_rand, PCCperm_all_mean_std(:,1),PCCperm_all_mean_std(:,2))
    hold off
    
    text(10,55, 'Full permutation')
    text(10,77, 'Permutation respecting ID')


    %% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
    fprintf(1, 'All calls Getting the DF Axis\n')
    Mdl_all = fitcdiscr(Score(:, 1:NPC_opt_all),DeafMic, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
    % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
    % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
    %  (default) — Matrix of size p-by-p, where p is the number of predictors
    % The axis of the DFA is the largest eigenvector of inv(Sigma) *
    % BetweenSigma which in the case of 2 classes...
    % is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
    % each category
    D_mean = mean(Score(contains(DeafMic, 'D'),1:NPC_opt_all));
    H_mean = mean(Score(contains(DeafMic, 'H'),1:NPC_opt_all));
    EigenVec_DFA_all = Mdl_all.Sigma \ (D_mean - H_mean)';
    PC_DF_all = PC(:,1:NPC_opt_all) * EigenVec_DFA_all(:,1);
    PC_DF_all = reshape(PC_DF_all, sum(IndWf), length(MPS_mic_wt));

    fprintf(1, 'Classification performance between K and S calls in Ac Grp %d: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\n#PC:%d \n', GR, PCC_all_mean_std(NPC_opt_ind,1), PCC_all_mean_std(NPC_opt_ind,2),PCCperm_all_mean_std(NPC_opt_ind,1), PCCperm_all_mean_std(NPC_opt_ind,2),PCCrand_all_mean_std(NPC_opt_ind,1), PCCrand_all_mean_std(NPC_opt_ind,2), NPC_opt_all)
    save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFA_AcGroup%d.mat', GR)),  'PC_val', 'NPC_opt_all','L_all', 'PCC_all_mean_std', 'Lperm_all', 'PCCperm_all_mean_std', 'Lrand_all', 'PCCrand_all_mean_std', 'PC_DF_all', 'EigenVec_DFA_all', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append' )

    % plot the positive direction of the DF1 axis in the MPS space
    FIG;
    subplot(1,3,3)
    plot_mps(PC_DF_all, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
    title('DF1 K vs S axis for all calls')
    Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
    flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
    colormap(Cmap)
    Axis = caxis();
    Lim = max(abs(Axis));
    caxis([-Lim Lim])
    suplabel(sprintf('Acoustid Group %d DFA performance Saline vs Kanamycin irrespective of sex', GR), 't')
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFA_AcGrp%d.png', GR)) , '-dpng')
    clear Score PC Mdl_all Mdlrand

    %% For All Male Calls run a PCA  and then a permutation DFA for males
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-----------------------Acoustic Group %d Male Calls----------------------------</strong>\n', GR)
    SexDeafMic = SexDeaf(IndMicAudioGood);
    MaleLogical = logical(contains(SexDeafMic, 'M').*(TmicAll7==GR));
    SexDeafMic_M=SexDeafMic(MaleLogical);
    [PC_M,Score_M,~, ~, VarExpl_M,~] = pca(MPS_mic_norm(MaleLogical,:));
    
    % Plot the % variance explained by the PC
    figure()
    CSVarExpl_M = cumsum(VarExpl_M);
    plot(CSVarExpl_M, 'Linewidth',2)
    xlabel('# PC')
    ylabel(sprintf('%%variance explained for Male Calls in Acoutsic Group %d', GR))
    NPC90var_M = find(cumsum(VarExpl_M)>90,1);
    text(100,95, sprintf('90%% variance explained with %d PC',NPC90var_M ))
    text(100,90, sprintf('%.1f%% variance explained with 100 PC', CSVarExpl_M(100)))
    if GR==5
        save(fullfile(LocalDataDir, sprintf('DeafBats_PCA_MaleAcGroup%d.mat', GR)), 'PC_M', 'VarExpl_M', 'NPC90var_M','Score_M')
    else
        save(fullfile(LocalDataDir, sprintf('DeafBats_PCA_MaleAcGroup%d.mat', GR)), 'PC_M', 'VarExpl_M', 'NPC90var_M')
    end

    save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFA_MaleAcGroup%d.mat',GR)),  'NPC90var_M')

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
        Mdl = fitcdiscr(Score_M(:, 1:PC_val_M(npc)),SexDeafMic_M, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
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
    subplot(1,3,1)
    shadedErrorBar(PC_val_M, PCC_males_mean_std(:,1),PCC_males_mean_std(:,2),{'-','color',ColorCode(5,:)})
    xlabel('# PC')
    ylabel('% Classification correct')
    set(gca, 'YLim', [40 100])
    hold on
    H=hline(50, 'k--','Chance level');
    H.LineWidth = 2;
    hold on
    V = vline(NPC_opt_M, 'b--', sprintf('Optimal #PC = %d', NPC_opt_M));
    V.LineWidth = 2;
    %% Permutation tests for male calls
    % Permutation of HD irrespective of ID
    fprintf(1, 'Male calls Ac Grp%d Permutation test irrespective of ID', GR)
    %PC_val_rand = PC_val_M(1:(NPC_opt_ind+NumPCsup));
    PC_val_rand = PC_val_M;
    Lrand_males = cell(length(PC_val_rand),1);
    for bb=1:NRandPerm
        fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
        RandInd = randperm(length(SexDeafMic_M));
        for npc=1:length(PC_val_rand)
            fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
            if bb==1
                Lrand_males{npc} = nan(NRandPerm,10);
            end
            Mdlrand = fitcdiscr(Score_M(:,1:PC_val_rand(npc)),SexDeafMic_M(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
            Lrand_males{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        end
    end
    PCCrand_males_mean_std = nan(length(PC_val_rand),2);
    for npc=1:length(PC_val_rand)
        PCCrand_males_mean_std(npc,1) = mean(reshape(100*(1-Lrand_males{npc}),numel(Lrand_males{npc}),1));
        PCCrand_males_mean_std(npc,2) = std(reshape(100*(1-Lrand_males{npc}),numel(Lrand_males{npc}),1));
    end


    % permutation of HD respecting ID
    fprintf(1, '\nMale calls Ag Grp%d Permutation test respecting ID', GR)
    C=nchoosek(1:4,2);
    Perm_males = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
    Perm_males = Perm_males(2:end,:); % The first row is the actual true order
    Lperm_males = cell(length(PC_val_rand),1);

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
        for npc=1:length(PC_val_rand)
            fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
            if bb==1
                Lperm_males{npc} = nan(size(Perm_males,1),10);
            end
            Mdlrand = fitcdiscr(Score_M(:,1:PC_val_rand(npc)), DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
            Lperm_males{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        end
    end

    PCCperm_males_mean_std = nan(length(PC_val_rand),2);
    for npc=1:length(PC_val_rand)
        PCCperm_males_mean_std(npc,1) = mean(reshape(100*(1-Lperm_males{npc}),numel(Lperm_males{npc}),1));
        PCCperm_males_mean_std(npc,2) = std(reshape(100*(1-Lperm_males{npc}),numel(Lperm_males{npc}),1));
    end

    % Add the permutation values to the figure
    FIG;
    subplot(1,3,2)
    shadedErrorBar(PC_val_rand, PCC_males_mean_std(1:(NPC_opt_ind+NumPCsup),1),PCC_males_mean_std(1:(NPC_opt_ind+NumPCsup),2),{'-','color',ColorCode(5,:)})
    xlabel('# PC')
    ylabel('% Classification correct')
    set(gca, 'YLim', [40 100])
    hold on
    H=hline(50, 'k--','Chance level');
    H.LineWidth = 2;
    hold on
    V = vline(NPC_opt_M, 'b--', sprintf('Optimal #PC = %d', NPC_opt_M));
    V.LineWidth = 2;
    hold on
    shadedErrorBar(PC_val_rand, PCCrand_males_mean_std(:,1),PCCrand_males_mean_std(:,2))
    hold on
    shadedErrorBar(PC_val_rand, PCCperm_males_mean_std(:,1),PCCperm_males_mean_std(:,2))
    hold off
    
    text(10,55, 'Full permutation')
    text(10,77, 'Permutation respecting ID')


    %% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
    fprintf(1, 'Male calls Ac Grp %d Getting the DF Axis\n', GR)
    Mdl_M = fitcdiscr(Score_M(:, 1:NPC_opt_M),SexDeafMic_M, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
    % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
    % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
    %  (default) — Matrix of size p-by-p, where p is the number of predictors
    % The axis of the DFA is the largest eigenvector of inv(Sigma) *
    % BetweenSigma which in the case of 2 classes...
    % is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
    % each category
    D_mean = mean(Score_M(contains(SexDeafMic_M, 'D'),1:NPC_opt_M));
    H_mean = mean(Score_M(contains(SexDeafMic_M, 'H'),1:NPC_opt_M));
    EigenVec_DFA_M = Mdl_M.Sigma \ (D_mean - H_mean)';
    PC_DF_males = PC_M(:,1:NPC_opt_M) * EigenVec_DFA_M(:,1);
    PC_DF_males = reshape(PC_DF_males, sum(IndWf), length(MPS_mic_wt));

    fprintf(1, 'Classification performance between male K and male S calls in Ac Grp %d: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\n#PC:%d \n', GR, PCC_males_mean_std(NPC_opt_ind,1), PCC_males_mean_std(NPC_opt_ind,2),PCCperm_males_mean_std(NPC_opt_ind,1), PCCperm_males_mean_std(NPC_opt_ind,2),PCCrand_males_mean_std(NPC_opt_ind,1), PCCrand_males_mean_std(NPC_opt_ind,2), NPC_opt_M)
    save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFA_MaleAcGroup%d.mat', GR)),  'PC_val_M', 'NPC_opt_M','L_males', 'PCC_males_mean_std', 'Lperm_males', 'PCCperm_males_mean_std', 'Lrand_males', 'PCCrand_males_mean_std', 'PC_DF_males', 'EigenVec_DFA_M', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append' )

    % plot the positive direction of the DF1 axis in the MPS space
    FIG;
    subplot(1,3,3)
    plot_mps(PC_DF_males, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
    title('DF1 K vs S axis for males calls')
    Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
        flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
    colormap(Cmap)
    Axis = caxis();
    Lim = max(abs(Axis));
    caxis([-Lim Lim])
    suplabel(sprintf('Acoustid Group %d Male Calls DFA performance Saline vs Kanamycin', GR), 't')
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFA_MaleCallsAcGrp%d.png',GR)) , '-dpng')
    
    clear Score_M PCC_M Mdl_M Mdlrand

    %% For All female Calls run a PCA  and then a permutation DFA for females
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>---------------------Female Calls Acoustic Group %d----------------------------</strong>\n', GR)
    SexDeafMic = SexDeaf(IndMicAudioGood);
    FemaleLogical = logical(contains(SexDeafMic, 'F').*(TmicAll7 == GR));
    SexDeafMic_F=SexDeafMic(FemaleLogical);
    [PC_F,Score_F,~, ~, VarExpl_F,~] = pca(MPS_mic_norm(FemaleLogical,:));

    % Plot the % variance explained by the PC
    figure()
    CSVarExpl_F = cumsum(VarExpl_F);
    plot(CSVarExpl_F, 'Linewidth',2)
    xlabel('# PC')
    ylabel(sprintf('%%variance explained for Female Calls in Ac Grp %d',GR))
    NPC90var_F = find(cumsum(VarExpl_F)>90,1);
    text(100,95, sprintf('90%% variance explained with %d PC',NPC90var_F ))
    text(100,90, sprintf('%.1f%% variance explained with 100 PC', CSVarExpl_F(100)))
    if GR==1
        save(fullfile(LocalDataDir, sprintf('DeafBats_PCA_FemaleAcGroup%d.mat', GR)), 'PC_F', 'VarExpl_F', 'NPC90var_F', 'Score_F')
    else
        save(fullfile(LocalDataDir, sprintf('DeafBats_PCA_FemaleAcGroup%d.mat', GR)), 'PC_F', 'VarExpl_F', 'NPC90var_F')
    end
    save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFA_FemaleAcGroup%d.mat',GR)),  'NPC90var_F')

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
        Mdl = fitcdiscr(Score_F(:, 1:PC_val_F(npc)),SexDeafMic_F, 'CrossVal','on','KFold',10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
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
    subplot(1,3,1)
    shadedErrorBar(PC_val_F, PCC_females_mean_std(:,1),PCC_females_mean_std(:,2),{'-','color',ColorCode(5,:)})
    xlabel('# PC')
    ylabel('% Classification correct')
    set(gca, 'YLim', [40 100])
    hold on
    H=hline(50, 'k--','Chance level');
    H.LineWidth = 2;
    hold on
    V = vline(NPC_opt_F, 'b--', sprintf('Optimal #PC = %d', NPC_opt_F));
    V.LineWidth = 2;
    %% Permutation tests for female calls
    % Permutation of HD irrespective of ID
    fprintf(1, 'Female calls Ac Grp%d Permutation test irrespective of ID', GR)
    %PC_val_rand = PC_val_F(1:(NPC_opt_ind+NumPCsup));
    PC_val_rand = PC_val_F;
    Lrand_females = cell(length(PC_val_rand),1);
    for bb=1:NRandPerm
        fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
        RandInd = randperm(length(SexDeafMic_F));
        for npc=1:length(PC_val_rand)
            fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
            if bb==1
                Lrand_females{npc} = nan(NRandPerm,10);
            end
            Mdlrand = fitcdiscr(Score_F(:,1:PC_val_rand(npc)),SexDeafMic_F(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
            Lrand_females{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        end
    end
    PCCrand_females_mean_std = nan(length(PC_val_rand),2);
    for npc=1:length(PC_val_rand)
        PCCrand_females_mean_std(npc,1) = mean(reshape(100*(1-Lrand_females{npc}),numel(Lrand_females{npc}),1));
        PCCrand_females_mean_std(npc,2) = std(reshape(100*(1-Lrand_females{npc}),numel(Lrand_females{npc}),1));
    end


    % permutation of HD respecting ID
    fprintf(1, '\nFemale calls Ac Grp %d Permutation test respecting ID', GR)
    C=nchoosek(1:6,3);
    Perm_females = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:))];
    Perm_females = Perm_females(2:end,:); % first row is in True order
    Lperm_females = cell(length(PC_val_rand),1);
    
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
        for npc=1:length(PC_val_rand)
            fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
            if bb==1
                Lperm_females{npc} = nan(size(Perm_females,1),10);
            end
            Mdlrand = fitcdiscr(Score_F(:,1:PC_val_rand(npc)), DeafMic_temp, 'CrossVal','on', 'KFold', 10,'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
            Lperm_females{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        end
    end

    PCCperm_females_mean_std = nan(length(PC_val_rand),2);
    for npc=1:length(PC_val_rand)
        PCCperm_females_mean_std(npc,1) = mean(reshape(100*(1-Lperm_females{npc}),numel(Lperm_females{npc}),1));
        PCCperm_females_mean_std(npc,2) = std(reshape(100*(1-Lperm_females{npc}),numel(Lperm_females{npc}),1));
    end

    % Add the permutation values to the figure
    FIG;
    subplot(1,3,2)
    shadedErrorBar(PC_val_rand, PCC_females_mean_std(1:(NPC_opt_ind+NumPCsup),1),PCC_females_mean_std(1:(NPC_opt_ind+NumPCsup),2),{'-','color',ColorCode(5,:)})
    xlabel('# PC')
    ylabel('% Classification correct')
    set(gca, 'YLim', [40 100])
    hold on
    H=hline(50, 'k--','Chance level');
    H.LineWidth = 2;
    hold on
    V = vline(NPC_opt_F, 'b--', sprintf('Optimal #PC = %d', NPC_opt_F));
    V.LineWidth = 2;
    hold on
    shadedErrorBar(PC_val_rand, PCCrand_females_mean_std(:,1),PCCrand_females_mean_std(:,2))
    hold on
    shadedErrorBar(PC_val_rand, PCCperm_females_mean_std(:,1),PCCperm_females_mean_std(:,2))
    hold off
    
    text(10,55, 'Full permutation')
    text(10,77, 'Permutation respecting ID')


    %% Rerun the DFA with all the data and the optimal number of PC to obtain the filter
    fprintf(1, 'Female calls Ac Grp%d Getting the DF Axis\n', GR)
    Mdl_F = fitcdiscr(Score_F(:, 1:NPC_opt_F),SexDeafMic_F, 'CrossVal','off', 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
    % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
    % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
    %  (default) — Matrix of size p-by-p, where p is the number of predictors
    % The axis of the DFA is the largest eigenvector of inv(Sigma) *
    % BetweenSigma which in the case of 2 classes...
    % is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
    % each category
    D_mean = mean(Score_F(contains(SexDeafMic_F, 'D'),1:NPC_opt_F));
    H_mean = mean(Score_F(contains(SexDeafMic_F, 'H'),1:NPC_opt_F));
    EigenVec_DFA_F = Mdl_F.Sigma \ (D_mean - H_mean)';
    PC_DF_females = PC_M(:,1:NPC_opt_F) * EigenVec_DFA_F(:,1);
    PC_DF_females = reshape(PC_DF_females, sum(IndWf), length(MPS_mic_wt));

fprintf(1, 'Classification performance between female K and female S calls in Acoustic Group %d: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\n#PC:%d \n', GR, PCC_females_mean_std(NPC_opt_ind,1), PCC_females_mean_std(NPC_opt_ind,2),PCCperm_females_mean_std(NPC_opt_ind,1), PCCperm_females_mean_std(NPC_opt_ind,2),PCCrand_females_mean_std(NPC_opt_ind,1), PCCrand_females_mean_std(NPC_opt_ind,2), NPC_opt_F)
save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFA_FemaleAcGroup%d.mat', GR)),  'PC_val_F', 'NPC_opt_F','L_females', 'PCC_females_mean_std', 'Lperm_females', 'PCCperm_females_mean_std', 'Lrand_females', 'PCCrand_females_mean_std', 'PC_DF_females','EigenVec_DFA_F', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append' )

% plot the positive direction of the DF1 axis in the MPS space
FIG;
subplot(1,3,3)
plot_mps(PC_DF_females, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
title('DF1 K vs S axis for females calls')
Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
    flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
colormap(Cmap)
Axis = caxis();
Lim = max(abs(Axis));
caxis([-Lim Lim])
suplabel(sprintf('Female Calls DFA performance Saline vs Kanamycin in Ac Group%d', GR), 't')
print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFA_FemaleAcGroup%d.png', GR)) , '-dpng')

clear Score_F PCC_F Mdl_F Mdlrand

    
end
%% Calculate the genelarized variance for males and females calls in the PCA
% space (determinant of the variance co-variance matrix). We only take the
% # PC that explains 90% of the variance for this calculation.
% 
% Calculations per bat
load(fullfile(LocalDataDir, 'DeafBats_PCA_AllCalls.mat'), 'Score', 'NPC90var')
BatIDMic = str2double(BatID(MicAudioGood01));
DValues = nan(length(BatName),1);
for bb=1:length(BatName)
    VCoVMat = cov(Score(BatIDMic==BatName(bb),1:NPC90var));
    DValues(bb) = (det(VCoVMat))^(1/(2*NPC90var));
end
DMax = (det(cov(Score(:,1:NPC90var))))^(1/(2*NPC90var));

CData = nan(length(BatName),3);
for bb=1:length(BatName)
    CData(bb,:) = UCSexDeaf(strcmp(USexDeaf,BatSexDeaf{bb}),:);
end

figure()
B=bar(1:10, DValues,'FaceColor', 'flat');
B.CData = CData;
xlabel('BatID')
B.Parent.XTickLabel = BatName;
ylabel('Vocalization dispersion')
for ct=1:length(USexDeaf)
    text(1,B.Parent.YLim(2)-diff(B.Parent.YLim)*(2+ct)/40, USexDeaf(ct), 'Color',UCSexDeaf(ct,:), 'FontWeight','bold' )
end
hold on
box off
hline(DMax, '--k')
hold off
fprintf(1,'Hearing female Vocalization dispersion: %.1f +- %.1f\n', mean(DValues(contains(BatSexDeaf, 'Hearing Female'))), std(DValues(contains(BatSexDeaf, 'Hearing Female'))));
fprintf(1,'Hearing female Vocalization dispersion: %.1f +- %.1f\n', mean(DValues(contains(BatSexDeaf, 'Hearing Male'))), std(DValues(contains(BatSexDeaf, 'Hearing Male'))));
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