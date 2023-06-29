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
UCT = UCT([1:4 9 5:8 10:14])
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
%% 
% First restrict the MPS spectral frequency to 3 cycles/kHz and normalized all 
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
%% 
% For AllCalls run a PCA  and then a permutation DFA first for both male and 
% females and then per sex


    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------All Calls----------------------------</strong>\n')
    [PC,Score,~, ~, VarExpl,~] = pca(MPS_mic_norm);
    
%% 
% Plot the % variance explained by the PC

    figure()
    CSVarExpl = cumsum(VarExpl);
    plot(CSVarExpl, 'Linewidth',2)
    xlabel('# PC')
    ylabel('%variance explained for All Calls')
    NPC = find(cumsum(VarExpl)>90,1);
    text(100,95, sprintf('90%% variance explained with %d PC',NPC ))
    text(100,90, sprintf('%.1f%% variance explained with 100 PC', CSVarExpl(100)))
    save(fullfile(LocalDataDir, 'Data4_DeafBats_PCA_AllCalls.mat'), 'PC', 'Score', 'VarExpl', 'NPC')
    save(fullfile(LocalDataDir, 'Data4_DeafBats_PermDFA_AllCalls.mat'),  'NPC')
%%
    NRandPerm=10;
%% 
% Permutation DFA for males only

    SexDeafMic = SexDeaf(IndMicAudioGood);
    MaleLogical = contains(SexDeafMic, 'M');
    SexDeafMic_M=SexDeafMic(MaleLogical);
    
    % Find the optimal regularization parameters Delta and Gamma with a PC
    % input that explains 90% variance
    fprintf(1, 'All male calls Finding the optimal regularization parameters\n')
    NPC = find(cumsum(VarExpl)>90,1);
    Mdl_DG_males = fitcdiscr(Score(MaleLogical, 1:NPC),SexDeafMic(MaleLogical), 'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', struct('Kfold',10), 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');

    
%%
    % Run the model with the best Delta and Gamma in cross validation and
    % get the Loss values accross validation sets
    fprintf(1, 'All male calls Calculating the classifier performance\n')
    Mdl = fitcdiscr(Score(MaleLogical, 1:NPC),SexDeafMic(MaleLogical), 'Delta', Mdl_DG_males.Delta, 'Gamma', Mdl_DG_males.Gamma, 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_males = kfoldLoss(Mdl, 'Mode', 'Individual');

    % Permutation of HD irrespective of ID
    fprintf(1, 'All male calls Permutation test irrespective of ID\n')
    Lrand_males = nan(NRandPerm,10);
    for bb=1:NRandPerm
        RandInd = randperm(length(SexDeafMic_M));
        Mdlrand = fitcdiscr(Score(MaleLogical,1:NPC),SexDeafMic_M(RandInd), 'Delta', Mdl_DG_males.Delta, 'Gamma', Mdl_DG_males.Gamma,'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lrand_males(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end

    % permutation of HD respecting ID
    fprintf(1, 'All male calls Permutation test respecting ID\n')
    Lperm_males = nan(3,10);
    Perm_males = [1:6 7 9 8 10; 1:6 9 8 7 10; 1:6 7 9 10 8];
    BatIDMicGR = BatID(IndMicAudioGood);
    BatIDMicGR_M = BatIDMicGR(MaleLogical);
    for bb=1:3
        SexDeafMic_M_temp = SexDeafMic_M;
        BatSexDeaf_local = BatSexDeaf(Perm_males(bb,:));
        for bat=1:length(BatName)
            if sum(contains(BatIDMicGR_M,num2str(BatName(bat))))
                Ind = find(contains(BatIDMicGR_M,num2str(BatName(bat))));
                for ii=1:length(Ind)
                    SexDeafMic_M_temp{Ind(ii)} = BatSexDeaf_local{bat};
                end
            end
        end
        Mdlrand = fitcdiscr(Score(MaleLogical,1:NPC),SexDeafMic_M_temp, 'Delta', Mdl_DG_males.Delta, 'Gamma', Mdl_DG_males.Gamma, 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lperm_males(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end

    % Rerun the DFA with all the data to obtain the filter
    fprintf(1, 'All male calls Getting the DF Axis\n')
    Mdl = fitcdiscr(Score(MaleLogical, 1:NPC),SexDeafMic(MaleLogical), 'Delta', Mdl_DG_males.Delta, 'Gamma', Mdl_DG_males.Gamma, 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
    % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
    % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
    %  (default) — Matrix of size p-by-p, where p is the number of predictors
    % The axis of the DFA is the largest eigenvector of inv(Sigma) *
    % BetweenSigma which in the case of 2 classes...
    % is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
    % each category
    MD_mean = mean(Score(contains(SexDeafMic, 'DM'),1:NPC));
    MH_mean = mean(Score(contains(SexDeafMic, 'HM'),1:NPC));
    PC_DFA = Mdl.Sigma \ (MD_mean - MH_mean)';
    PC_DF_males = PC(:,1:NPC) * PC_DFA(:,1);
    PC_DF_males = reshape(PC_DF_males, sum(IndWf), length(MPS_mic_wt));
    
    fprintf(1, 'Classification performance between K and S male calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\nGamma:%.2f and Delta:%.1e\n', mean((1-L_males)*100), std((1-L_males)*100),mean((1-Lperm_males)*100), std((1-Lperm_males)*100),mean((1-Lrand_males)*100), std((1-Lrand_males)*100), Mdl_DG_males.Gamma, Mdl_DG_males.Delta)
    save(fullfile(LocalDataDir, 'Data4_DeafBats_PermDFA_AllCalls.mat'), 'Mdl_DG_males', 'L_males', 'Lrand_males', 'Lperm_males', 'Perm_males', 'PC_DF_males', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append')
    
    
%%
% plot the positive direction of the DF1 axis in the MPS space
    figure()
    plot_mps(PC_DF_males, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
    title('DF1 K vs S axis for male calls')
Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
colormap(Cmap)
Axis = caxis();
Lim = max(abs(Axis));
caxis([-Lim Lim])
%% 
% Permutation DFA for females only

    SexDeafMic = SexDeaf(IndMicAudioGood);
    FemaleLogical = contains(SexDeafMic, 'F');
    SexDeafMic_F=SexDeafMic(FemaleLogical);
    
    % Find the optimal regularization parameters Delta and Gamma with a PC
    fprintf(1, 'All female calls Finding the optimal regularization parameters\n')
    % input that explains 90% variance
    Mdl_DG_females = fitcdiscr(Score(FemaleLogical, 1:NPC),SexDeafMic(FemaleLogical), 'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', struct('Kfold',10), 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    
    % Run the model with the best Delta and Gamma in cross validation and
    % get the Loss values accross validation sets
    fprintf(1, 'All female calls Calculating classifier performance\n')
    Mdl = fitcdiscr(Score(FemaleLogical, 1:NPC),SexDeafMic(FemaleLogical), 'Delta', Mdl_DG_females.Delta, 'Gamma', Mdl_DG_females.Gamma, 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_females = kfoldLoss(Mdl, 'Mode', 'Individual');

    % Permutation of HD irrespective of ID
    fprintf(1, 'All female calls permutation testNOT  respecting ID\n')
    Lrand_females = nan(NRandPerm,10);
    for bb=1:NRandPerm
        RandInd = randperm(length(SexDeafMic_F));
        Mdlrand = fitcdiscr(Score(FemaleLogical,1:NPC),SexDeafMic_F(RandInd), 'Delta', Mdl_DG_females.Delta, 'Gamma', Mdl_DG_females.Gamma,'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lrand_females(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end

    % permutation of HD respecting ID
    fprintf(1, 'All female calls permutation test respecting ID\n')
    C=nchoosek(1:6,3);
    Perm_females = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:)), repmat(7:10,size(C,1)/2,1)];
    Lperm_females = nan(size(Perm_females,1),10);
    BatIDMicGR = BatID(IndMicAudioGood);
    BatIDMicGR_F = BatIDMicGR(FemaleLogical);
    for bb=1:size(Perm_females,1)
        SexDeafMic_F_temp = SexDeafMic_F;
        BatSexDeaf_local = BatSexDeaf(Perm_females(bb,:));
        for bat=1:length(BatName)
            if sum(contains(BatIDMicGR_F,num2str(BatName(bat))))
                Ind = find(contains(BatIDMicGR_F,num2str(BatName(bat))));
                for ii=1:length(Ind)
                    SexDeafMic_F_temp{Ind(ii)} = BatSexDeaf_local{bat};
                end
            end
        end
        Mdlrand = fitcdiscr(Score(FemaleLogical,1:NPC),SexDeafMic_F_temp, 'Delta', Mdl_DG_females.Delta, 'Gamma', Mdl_DG_females.Gamma, 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lperm_females(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
    end

    % Rerun the DFA with all the data to obtain the filter
    fprintf(1, 'All female calls DF axis\n')
    Mdl = fitcdiscr(Score(FemaleLogical, 1:NPC),SexDeafMic(FemaleLogical), 'Delta', Mdl_DG_females.Delta, 'Gamma', Mdl_DG_females.Gamma, 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
    % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
    % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
    %  (default) — Matrix of size p-by-p, where p is the number of predictors
    % The axis of the DFA is the largest eigenvector of inv(Sigma) *
    % BetweenSigma which in the case of 2 classes...
    % is inv(Sigma) * (Mu1 - Mu2) where Mu1 and Mu2 are the sample means of
    % each category
    FD_mean = mean(Score(contains(SexDeafMic, 'DF'),1:NPC));
    FH_mean = mean(Score(contains(SexDeafMic, 'HF'),1:NPC));
    PC_DFA = Mdl.Sigma \ (FD_mean - FH_mean)';
    PC_DF_females = PC(:,1:NPC) * PC_DFA(:,1);
    PC_DF_females = reshape(PC_DF_females, sum(IndWf), length(MPS_mic_wt));
    
   
%%
 fprintf(1, 'Classification performance between K and S female calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\nGamma:%.2f and Delta:%.1e\n', mean((1-L_females)*100), std((1-L_females)*100),mean((1-Lperm_females)*100), std((1-Lperm_females)*100),mean((1-Lrand_females)*100), std((1-Lrand_females)*100), Mdl_DG_females.Gamma, Mdl_DG_females.Delta)
    save(fullfile(LocalDataDir, 'Data4_DeafBats_PermDFA_AllCalls.mat'), 'Mdl_DG_females', 'L_females', 'Lrand_females', 'Lperm_females', 'Perm_females', 'PC_DF_females', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append')
    
    % plot the positive direction of the DF1 axis in the MPS space
    figure()
    plot_mps(PC_DF_females, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
    title('DF1 K vs S axis for female calls')
    Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
    flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
    colormap(Cmap)
    Axis = caxis();
    Lim = max(abs(Axis));
    caxis([-Lim Lim])
%% 
% Now for each acoustic group run a permutation DFA per sex

UGroup = unique(TmicAll7);
for gg=1:length(UGroup)
    GR = UGroup(gg);
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------Acoustic Group %d----------------------------</strong>\n', GR)
    [PC,Score,~, ~, VarExpl,~] = pca(MPS_mic_norm(TmicAll7==GR,:));
    save(fullfile(LocalDataDir, sprintf('Data4_DeafBats_PermDFA_AcGroup%d.mat', GR)), 'PC', 'Score', 'VarExpl')
%% 
% (Plot the % variance explained by the PC) Not anymore we use the same PCA 
% space as the one used for all calls

%     figure()
%     CSVarExpl = cumsum(VarExpl);
%     plot(CSVarExpl, 'Linewidth',2)
%     xlabel('# PC')
%     ylabel(sprintf('%variance explained for Acoustic group %d', GR))
%     NPC = find(cumsum(VarExpl)>90,1);
%     text(100,95, sprintf('80%% variance explained with %d PC',NPC ))
%     text(100,90, sprintf('%.1f%% variance explained with 40 PC', CSVarExpl(40)))
%% 
% Permutation DFA for males only

    SexDeafMic = SexDeaf(MicAudioGood01);
    SexDeafMicGR = SexDeafMic(TmicAll7==GR);
    MaleLogical = contains(SexDeafMicGR, 'M');
    SexDeafMicGR_M=SexDeafMicGR(MaleLogical);
    Ind_GR_Males = logical((TmicAll7==GR).*contains(SexDeafMic, 'M'));
    
    % Find the optimal regularization parameters Delta and Gamma with a PC
    % input that explains 90% variance
    fprintf(1, 'Group %d male calls Finding optimal regularization parameters\n', GR)
    Mdl_DG_males = fitcdiscr(Score(Ind_GR_Males, 1:NPC),SexDeafMicGR(MaleLogical), 'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', struct('Kfold',10), 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    

     % Run the model with the best Delta and Gamma in cross validation and
    % get the Loss values accross validation sets
    fprintf(1, 'Group %d male calls calculating the performance of the classifier\n', GR)
    Mdl = fitcdiscr(Score(Ind_GR_Males, 1:NPC),SexDeafMicGR(MaleLogical), 'Delta', Mdl_DG_males.Delta, 'Gamma', Mdl_DG_males.Gamma, 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_males = kfoldLoss(Mdl, 'Mode', 'Individual');

    % Permutation of HD irrespective of ID
    fprintf(1, 'Group %d male calls permutation test irrespective of ID\n', GR)
    Lrand_males = nan(NRandPerm,1);
    for bb=1:NRandPerm
        RandInd = randperm(length(SexDeafMicGR_M));
        Mdlrand = fitcdiscr(Score(Ind_GR_Males,1:NPC),SexDeafMicGR_M(RandInd), 'Delta', Mdl_DG_males.Delta, 'Gamma', Mdl_DG_males.Gamma,'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lrand_males(bb) = kfoldLoss(Mdlrand);
    end

    % permutation of HD respecting ID
    fprintf(1, 'Group %d male calls permutation test respecting ID\n', GR)
    Lperm_males = nan(3,1);
    Perm_males = [1:6 7 9 8 10; 1:6 9 8 7 10; 1:6 7 9 10 8];
    BatIDMicGR = BatID(IndMicAudioGood(TmicAll7==GR));
    BatIDMicGR_M = BatIDMicGR(MaleLogical);
    for bb=1:3
        SexDeafMicGR_M_temp = SexDeafMicGR_M;
        BatSexDeaf_local = BatSexDeaf(Perm_males(bb,:));
        for bat=1:length(BatName)
            if sum(contains(BatIDMicGR_M,num2str(BatName(bat))))
                Ind = find(contains(BatIDMicGR_M,num2str(BatName(bat))));
                for ii=1:length(Ind)
                    SexDeafMicGR_M_temp{Ind(ii)} = BatSexDeaf_local{bat};
                end
            end
        end
        Mdlperm = fitcdiscr(Score(Ind_GR_Males,1:NPC),SexDeafMicGR_M_temp, 'Delta', Mdl_DG_males.Delta, 'Gamma', Mdl_DG_males.Gamma, 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lperm_males(bb) = kfoldLoss(Mdlperm);
    end

    % Rerun the DFA with all the data to obtain the filter
    fprintf(1, 'Group %d male calls getting DF axis\n', GR)
    Mdl = fitcdiscr(Score(Ind_GR_Males, 1:NPC),SexDeafMicGR_M, 'Delta', Mdl_DG_males.Delta, 'Gamma', Mdl_DG_males.Gamma, 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
    % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
    % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
    %  (default) — Matrix of size p-by-p, where p is the number of predictors
    % The firt PC of the pca of the between class covariance matrix divided
    % by within class covariance matrix is the first filter of the DFA in
    % DFA input space (first PCA above)
    %     [EigenVect,~] = eig(Mdl.BetweenSigma,Mdl.Sigma);
    %     PC_DF_males = PC(:,1:NPC) * EigenVect(:,1);
    [PC_DFA,~] = pca(Mdl.Sigma \ Mdl.BetweenSigma);
    PC_DF_males = PC(:,1:NPC) * PC_DFA(:,1);
    PC_DF_males = reshape(PC_DF_males, sum(IndWf), length(MPS_mic_wt));

    fprintf(1, 'Acoustic Group %d Classification performance between K and S male calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\nGamma:%.2f and Delta:%.1e\n', GR, mean((1-L_males)*100), std((1-L_males)*100),mean((1-Lperm_males)*100), std((1-Lperm_males)*100),mean((1-Lrand_males)*100), std((1-Lrand_males)*100), Mdl_DG_males.Gamma, Mdl_DG_males.Delta)
    save(fullfile(LocalDataDir, sprintf('Data4_DeafBats_PermDFA_AcGroup%d.mat', GR)), 'Mdl_DG_males', 'L_males', 'Lrand_males', 'Lperm_males', 'Perm_males', 'PC_DF_males', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt')


    % plot the positive direction of the DF1 axis in the MPS space
    figure()
    plot_mps(PC_DF_males, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
    title(sprintf('Ac Group %d DF1 K vs S axis for male calls', GR))
    Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
        flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
    colormap(Cmap)
    Axis = caxis();
    Lim = max(abs(Axis));
    caxis([-Lim Lim])
%% 
% Permutation DFA for females only

    SexDeafMicGR = SexDeaf(IndMicAudioGood(TmicAll7==GR));
    FemaleLogical = contains(SexDeafMicGR, 'F');
    SexDeafMicGR_F=SexDeafMicGR(FemaleLogical);
    Ind_GR_Females = logical((TmicAll7==GR).*contains(SexDeafMic, 'F'));

    % Find the optimal regularization parameters Delta and Gamma with a PC
    % input that explains 90% variance
    fprintf(1, 'Group %d female calls Finding the optimal regularization parameters\n', GR)
    Mdl_DG_females = fitcdiscr(Score(Ind_GR_Females, 1:NPC),SexDeafMicGR_F, 'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', struct('Kfold',10), 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    

     % Run the model with the best Delta and Gamma in cross validation and
    % get the Loss values accross validation sets
    fprintf(1, 'Group %d female calls Calculating performance of the classifier\n', GR)
    Mdl = fitcdiscr(Score(Ind_GR_Females, 1:NPC),SexDeafMicGR_F, 'Delta', Mdl_DG_females.Delta, 'Gamma', Mdl_DG_females.Gamma, 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
    L_females = kfoldLoss(Mdl, 'Mode', 'Individual');

    
    % Permutation of HD irrespective of ID
    fprintf(1, 'Group %d female calls Permutation test irrespective of ID\n', GR)
    Lrand_females = nan(NRandPerm,1);
    for bb=1:NRandPerm
        RandInd = randperm(length(SexDeafMicGR_F));
        Mdlrand = fitcdiscr(Score(Ind_GR_Females,1:NPC),SexDeafMicGR_F(RandInd), 'Delta', Mdl_DG_females.Delta, 'Gamma', Mdl_DG_females.Gamma,'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lrand_females(bb) = kfoldLoss(Mdlrand);
    end

    % permutation of HD respecting ID
    fprintf(1, 'Group %d female calls Permuattion test respecting ID\n', GR)
    C=nchoosek(1:6,3);
    Perm_females = [C(1:size(C,1)/2,:), flip(C((size(C,1)/2 +1) : end,:)), repmat(7:10,size(C,1)/2,1)];
    Lperm_females = nan(size(Perm_females,1),1);
    BatIDMicGR = BatID(IndMicAudioGood(TmicAll7==GR));
    BatIDMicGR_F = BatIDMicGR(FemaleLogical);
    for bb=1:size(Perm_females,1)
        SexDeafMicGR_F_temp = SexDeafMicGR_F;
        BatSexDeaf_local = BatSexDeaf(Perm_females(bb,:));
        for bat=1:length(BatName)
            if sum(contains(BatIDMicGR_F,num2str(BatName(bat))))
                Ind = find(contains(BatIDMicGR_F,num2str(BatName(bat))));
                for ii=1:length(Ind)
                    SexDeafMicGR_F_temp{Ind(ii)} = BatSexDeaf_local{bat};
                end
            end
        end
        Mdlperm = fitcdiscr(Score(Ind_GR_Females,1:NPC),SexDeafMicGR_F_temp, 'Delta', Mdl_DG_females.Delta, 'Gamma', Mdl_DG_females.Gamma, 'CrossVal','on', 'KFold', 10, 'Prior', [0.5 0.5], 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lperm_females(bb) = kfoldLoss(Mdlperm);
    end
    
    % Rerun the DFA with all the data to obtain the filter
    fprintf(1, 'Group %d female calls getting DF axis\n', GR)
    Mdl = fitcdiscr(Score(Ind_GR_Females, 1:NPC),SexDeafMicGR_F, 'Delta', Mdl_DG_females.Delta, 'Gamma', Mdl_DG_females.Gamma, 'Prior', [0.5 0.5], 'SaveMemory', 'off', 'FillCoeffs', 'on');
    % Mdl.BetweenSigma is the p-by-p matrix, the between-class covariance, where p is the number of predictors.
    % Mdl.Sigma % Within-class covariance matrix or matrices. The dimensions depend on DiscrimType: 'linear'...
    %  (default) — Matrix of size p-by-p, where p is the number of predictors
    % The firt PC of the pca of the between class covariance matrix divided
    % by within class covariance matrix is the first filter of the DFA in
    % DFA input space (first PCA above)
    %     [EigenVect,~] = eig(Mdl.BetweenSigma,Mdl.Sigma);
    %     PC_DF_males = PC(:,1:NPC) * EigenVect(:,1);
    [PC_DFA,~] = pca(Mdl.Sigma \ Mdl.BetweenSigma);
    PC_DF_females = PC(:,1:NPC) * PC_DFA(:,1);
    PC_DF_females = reshape(PC_DF_females, sum(IndWf), length(MPS_mic_wt));

    fprintf(1, 'Acoustic Group %d Classification performance between K and S female calls: %.1f%% +/-%.1f%%\nPermutation value (respecting ID): %.1f%% +/-%.1f%%\nPermutation value (irrespective of ID): %.1f%% +/-%.1f%%\nGamma:%.2f and Delta:%.1e\n', GR, mean((1-L_females)*100), std((1-L_females)*100),mean((1-Lperm_females)*100), std((1-Lperm_females)*100),mean((1-Lrand_females)*100), std((1-Lrand_females)*100), Mdl_DG_females.Gamma, Mdl_DG_females.Delta)
    save(fullfile(LocalDataDir, sprintf('Data4_DeafBats_PermDFA_AcGroup%d.mat', GR)), 'Mdl_DG_females', 'L_females', 'Lrand_females', 'Lperm_females', 'Perm_females', 'PC_DF_females', 'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append')


    % plot the positive direction of the DF1 axis in the MPS space
    figure()
    plot_mps(PC_DF_females, MPS_mic_wf(IndWf),MPS_mic_wt, 60,nan,0,[0 max(MPS_mic_wf(IndWf).*10^3)], [-150 150]);
    title(sprintf('Ac Group %d DF1 K vs S axis for female calls', GR))
    Cmap = flip([ones(128,1) (0:1/127:1)' (0:1/127:1)';
      flip(0:1/127:1)' flip(0:1/127:1)' ones(128,1)]);
    colormap(Cmap)
    Axis = caxis();
    Lim = max(abs(Axis));
    caxis([-Lim Lim])
end
%% 
% Calculate the genelarized variance for males and females calls in the PCA 
% space (determinant of the variance co-variance matrix). We only take the first 
% 100 PC for this calculation, similar to what is done for the DFA.
% 
% Caluclations per bat

NPC =100;
BatIDMic = str2double(BatID(MicAudioGood01));
DValues = nan(length(BatName),1);
for bb=1:length(BatName)
    VCoVMat = cov(Score(BatIDMic==BatName(bb),1:NPC));
    DValues(bb) = (det(VCoVMat))^(1/(2*NPC));
end
DMax = (det(cov(Score(:,1:NPC))))^(1/(2*NPC));

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