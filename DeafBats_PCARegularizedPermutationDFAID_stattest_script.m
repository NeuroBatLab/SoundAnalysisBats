%% Permutation DFA of bat ID for 7 acoustic groups

LocalDataDir = '/Users/elie/Documents/DeafBats/Data';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
SalineOnly = 1; % Set to 1 to only run calculations with saline Bats, 0 otherwise
%%
% Loading previous data

load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls2.mat'), 'CallType', 'BatID')
load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls.mat'),'MicAudioGood','TmicAll7');
MicAudioGood01 = MicAudioGood;
MicAudioGood01(isnan(MicAudioGood01)) = 0;
MicAudioGood01 = logical(MicAudioGood01);
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

% Number of bootsrtaps
NRandPerm=10^4;
%% For All Male Calls statistical test of permutation DFA
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------All Male Calls----------------------------</strong>\n')
SexDeafMic = SexDeaf(MicAudioGood01);
if SalineOnly
    MaleLogical = contains(SexDeafMic, 'HM');
else
    MaleLogical = contains(SexDeafMic, 'M');
end
BatIDMic = BatID(MicAudioGood01);
BatIDMic_M=BatIDMic(MaleLogical);
NumInd = length(unique(BatIDMic_M));
[~,Score_M,~] = pca(MPS_mic_norm(MaleLogical,:));
if SalineOnly
    load(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAIDS_MaleCalls.mat'), 'PC_val_M', 'NPC_opt_M','L_males')
else
    load(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAID_MaleCalls.mat'),  'PC_val_M', 'NPC_opt_M','L_males')
end
%% Permutation tests for male calls
% Permutation of ID irrespective of Sex
fprintf(1, 'Male calls Permutation test')
MeanLrand_males_bootstrap = cell(NRandPerm,1);
parfor bb=1:NRandPerm
    if ~rem(bb,NRandPerm/10)
        fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
    end
    RandInd = randperm(length(BatIDMic_M));
    Mdlrand = fitcdiscr(Score_M(:,1:NPC_opt_M),BatIDMic_M(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', 1/NumInd .* ones(1,NumInd), 'SaveMemory', 'on', 'FillCoeffs', 'off');
    MeanLrand_males_bootstrap{bb} = mean(kfoldLoss(Mdlrand, 'Mode', 'Individual'));
end
MeanPCCrand_males_bootstrap = 100-100.*cell2mat(MeanLrand_males_bootstrap);
if SalineOnly
    save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAIDS_MaleCalls.mat'),  'MeanPCCrand_males_bootstrap', '-append')
else
    save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAID_MaleCalls.mat'),  'MeanPCCrand_males_bootstrap', '-append')
end

FIG = figure;
histogram(MeanPCCrand_males_bootstrap)
hold on
V=vline(mean(100-100*L_males{PC_val_M==NPC_opt_M}), 'r--', 'observed');
V.LineWidth = 2;
xlabel('Correct classification (%)')
L=legend('random', 'Location', 'north');
L.Box='off';
title(sprintf('p = %.4f', sum(MeanPCCrand_males_bootstrap>=mean(100-100*L_males{PC_val_M==NPC_opt_M}))/NRandPerm))
suplabel('Male Calls DFA performance ID', 't');


if SalineOnly
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAIDS_MaleCalls_bootstrap.png') , '-dpng')
else
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAID_MaleCalls_bootstrap.png') , '-dpng')
end

clear Score_M Mdlrand
%% For All female Calls statistical test of permutation DFA
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
fprintf(1,'<strong>-------------------------------All Female Calls----------------------------</strong>\n')
SexDeafMic = SexDeaf(MicAudioGood01);
if SalineOnly
    FemaleLogical = contains(SexDeafMic, 'HF');
else
    FemaleLogical = contains(SexDeafMic, 'F');
end
BatIDMic = BatID(MicAudioGood01);
BatIDMic_F=BatIDMic(FemaleLogical);
NumInd = length(unique(BatIDMic_F));
[~,Score_F,~] = pca(MPS_mic_norm(FemaleLogical,:));
if SalineOnly
    load(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAIDS_FemaleCalls.mat'), 'PC_val_F', 'NPC_opt_F','L_females')
else
    load(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAID_FemaleCalls.mat'),  'PC_val_F', 'NPC_opt_F','L_females')
end
%% Permutation tests for female calls
% Permutation of ID
fprintf(1, 'Female calls Permutation test')
MeanLrand_females_bootstrap = cell(NRandPerm,1);
parfor bb=1:NRandPerm
    if ~rem(bb,NRandPerm/10)
        fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
    end
    RandInd = randperm(length(BatIDMic_F));
    Mdlrand = fitcdiscr(Score_F(:,1:NPC_opt_F),BatIDMic_F(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', 1/NumInd .* ones(1,NumInd), 'SaveMemory', 'on', 'FillCoeffs', 'off');
    MeanLrand_females_bootstrap{bb} = mean(kfoldLoss(Mdlrand, 'Mode', 'Individual'));
end
MeanPCCrand_females_bootstrap = 100-100.*cell2mat(MeanLrand_females_bootstrap);
if SalineOnly
    save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAIDS_FemaleCalls.mat'),  'MeanPCCrand_females_bootstrap', '-append')
else
    save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAID_FemaleCalls.mat'),  'MeanPCCrand_females_bootstrap', '-append')
end

FIG = figure;
histogram(MeanPCCrand_females_bootstrap)
hold on
V=vline(mean(100-100*L_females{PC_val_F==NPC_opt_F}), 'r--', 'observed');
V.LineWidth = 2;
xlabel('Correct classification (%)')
L=legend('random', 'Location', 'north');
L.Box='off';
title(sprintf('p = %.4f', sum(MeanPCCrand_females_bootstrap>=mean(100-100*L_females{PC_val_F==NPC_opt_F}))/NRandPerm))
suplabel('Female Calls DFA performance ID', 't');


if SalineOnly
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAIDS_FemaleCalls_bootstrap.png') , '-dpng')
else
    print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAID_FemaleCalls_bootstrap.png') , '-dpng')
end

clear Score_F Mdlrand

%% Now same analysis for each acoustic group
SexDeafMic = SexDeaf(MicAudioGood01);
UGroup = unique(TmicAll7);
for gg=1:length(UGroup)
    GR = UGroup(gg);
    %% For All Male Calls statistical test of permutation DFA
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-----------------------Acoustic Group %d Male Calls----------------------------</strong>\n', GR)
    
    if SalineOnly
        MaleLogical = logical(contains(SexDeafMic, 'HM').*(TmicAll7==GR));
    else
        MaleLogical = logical(contains(SexDeafMic, 'M').*(TmicAll7==GR));
    end
    BatIDMic = BatID(MicAudioGood01);
    BatIDMic_M=BatIDMic(MaleLogical);
    NumInd = length(unique(BatIDMic_M));
    [~,Score_M,~, ~, ~,~] = pca(MPS_mic_norm(MaleLogical,:));
    if SalineOnly
        load(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAIDS_MaleAcGroup%d.mat', GR)), 'PC_val_M', 'NPC_opt_M','L_males')
    else
        load(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAID_MaleAcGroup%d.mat', GR)),  'PC_val_M', 'NPC_opt_M','L_males')
    end
    
    %% Permutation tests for male calls
    % Permutation of ID irrespective of Sex
    fprintf(1, 'Male calls Permutation test')
    MeanLrand_males_bootstrap = cell(NRandPerm,1);
    parfor bb=1:NRandPerm
        if ~rem(bb,NRandPerm/10)
            fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
        end
        RandInd = randperm(length(BatIDMic_M));
        Mdlrand = fitcdiscr(Score_M(:,1:NPC_opt_M),BatIDMic_M(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', 1/NumInd .* ones(1,NumInd), 'SaveMemory', 'on', 'FillCoeffs', 'off');
        MeanLrand_males_bootstrap{bb} = mean(kfoldLoss(Mdlrand, 'Mode', 'Individual'));
    end
    MeanPCCrand_males_bootstrap = 100-100.*cell2mat(MeanLrand_males_bootstrap);
    if SalineOnly
        save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAIDS_MaleAcGroup%d.mat', GR)),  'MeanPCCrand_males_bootstrap', '-append')
    else
        save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAID_MaleAcGroup%d.mat', GR)),  'MeanPCCrand_males_bootstrap', '-append')
    end

    FIG = figure;
    histogram(MeanPCCrand_males_bootstrap)
    hold on
    V=vline(mean(100-100*L_males{PC_val_M==NPC_opt_M}), 'r--', 'observed');
    V.LineWidth = 2;
    xlabel('Correct classification (%)')
    L=legend('random', 'Location', 'north');
    L.Box='off';
    title(sprintf('p = %.4f', sum(MeanPCCrand_males_bootstrap>=mean(100-100*L_males{PC_val_M==NPC_opt_M}))/NRandPerm))
    suplabel(sprintf('Male Calls Ac Grp%d DFA performance ID', GR), 't');


    if SalineOnly
        print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFAIDS_MaleCallsAcGrp%d_bootstrap.png',GR)) , '-dpng')
    else
        print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAID_MaleCallsAcGrp%d_bootstrap.png',GR) , '-dpng')
    end

    clear Score_M Mdlrand

    %% For All female Calls run a PCA  and then a permutation DFA for females
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>-------------------------------------------------------------------------------</strong>\n')
    fprintf(1,'<strong>---------------------Female Calls Acoustic Group %d----------------------------</strong>\n', GR)
    if SalineOnly
        FemaleLogical = logical(contains(SexDeafMic, 'HF').*(TmicAll7==GR));
    else
        FemaleLogical = logical(contains(SexDeafMic, 'F').*(TmicAll7==GR));
    end
    BatIDMic = BatID(MicAudioGood01);
    BatIDMic_F=BatIDMic(FemaleLogical);
    NumInd = length(unique(BatIDMic_F));
    [~,Score_F,~, ~, ~,~] = pca(MPS_mic_norm(FemaleLogical,:));
    if SalineOnly
        load(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAIDS_FemaleAcGroup%d.mat', GR)), 'PC_val_F', 'NPC_opt_F','L_females')
    else
        load(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAID_FemaleAcGroup%d.mat', GR)),  'PC_val_F', 'NPC_opt_F','L_females')
    end
    %% Permutation tests for female calls
    % Permutation of ID
    fprintf(1, 'Male fecalls Permutation test')
    MeanLrand_females_bootstrap = cell(NRandPerm,1);
    parfor bb=1:NRandPerm
        if ~rem(bb,NRandPerm/10)
            fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
        end
        RandInd = randperm(length(BatIDMic_F));
        Mdlrand = fitcdiscr(Score_F(:,1:NPC_opt_F),BatIDMic_F(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', 1/NumInd .* ones(1,NumInd), 'SaveMemory', 'on', 'FillCoeffs', 'off');
        MeanLrand_females_bootstrap{bb} = mean(kfoldLoss(Mdlrand, 'Mode', 'Individual'));
    end
    MeanPCCrand_females_bootstrap = 100-100.*cell2mat(MeanLrand_females_bootstrap);
    if SalineOnly
        save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAIDS_FemaleAcGroup%d.mat', GR)),  'MeanPCCrand_females_bootstrap', '-append')
    else
        save(fullfile(LocalDataDir, sprintf('DeafBats_RegularizedPermDFAID_FemaleAcGroup%d.mat', GR)),  'MeanPCCrand_females_bootstrap', '-append')
    end

    FIG = figure;
    histogram(MeanPCCrand_females_bootstrap)
    hold on
    V=vline(mean(100-100*L_females{PC_val_F==NPC_opt_F}), 'r--', 'observed');
    V.LineWidth = 2;
    xlabel('Correct classification (%)')
    L=legend('random', 'Location', 'north');
    L.Box='off';
    title(sprintf('p = %.4f', sum(MeanPCCrand_females_bootstrap>=mean(100-100*L_females{PC_val_F==NPC_opt_F}))/NRandPerm))
    suplabel(sprintf('Female Calls Ac Grp%d DFA performance ID', GR), 't');


    if SalineOnly
        print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', sprintf('RegPermDFAIDS_FemaleCallsAcGrp%d_bootstrap.png',GR)) , '-dpng')
    else
        print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAID_FemaleCallsAcGrp%d_bootstrap.png',GR) , '-dpng')
    end

    clear Score_F Mdlrand

    
end
