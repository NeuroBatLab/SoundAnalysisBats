%% Regularized DFA of identity for 4 acoustic groups
% here we are using the clustering obtained with the 3 dimensions UMAP
% calculated in Notebbok 8ter. The optimization of the number
% of dimensions for the UMAP clustering performance was done by looking at
% median values of Silhouette (notebook 8ter)

LocalDataDir = '/Users/elie/Documents/DeafBats/Data';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
BoxPath = '/Users/elie/Box';
Path2Paper = fullfile(BoxPath, 'BatmanData', 'Deaf Paper');
%%
% Loading previous data

load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls2.mat'), 'CallType', 'BatID','MicAudioGood','MPS_mic', 'MPS_mic_wf', 'MPS_mic_wt')
%load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls.mat'),'MicAudioGood','TmicAll7');
load(fullfile(BaseDataDir, 'Data4_DeafBats_CatCalls_UMAPMic.mat'), 'TmicAll4');
MicAudioGood01 = MicAudioGood;
MicAudioGood01(isnan(MicAudioGood01)) = 0;
MicAudioGood01 = logical(MicAudioGood01);
%%
% Get the color vector ready for BatName, Sex and Deafness
ColorCode = [get(groot, 'DefaultAxesColorOrder'); 0 1 1; 0.5 0.5 0.5; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0 1; 0 0 0];
Path2RecordingTable = fullfile(BoxPath,'JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx');
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

%% For AllCalls run a PCA  and then a permutation DFA for both male and
% females
NRandPerm=100;
PLim = 0.01;
NumPCsup = 10;

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
% save(fullfile(LocalDataDir, 'DeafBats_PCA_AllCalls.mat'), 'PC',  'VarExpl', 'NPC90var')
% save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFA_AllCalls.mat'),  'NPC90var')

%% DFA for both male and female calls
DeafMic = Deaf(IndMicAudioGood);
BatIDMic = BatID(MicAudioGood01);
NumInd = length(unique(BatIDMic));

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
    Mdl = fitcdiscr(Score(:, 1:PC_val(npc)),BatIDMic, 'CrossVal','on','KFold',10, 'Prior', 1/NumInd.*ones(1,NumInd), 'SaveMemory', 'on', 'FillCoeffs', 'off');
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
% Contingency of correct and wrong classification for the optimal NPC for
% Fisher exact test with call permutation DFA
% CorrectObs = round((1-mean(L_all{PC_val==NPC_opt_all})).*size(Score,1));
% ErrorObs = round(mean(L_all{PC_val==NPC_opt_all}).*size(Score,1));

% Plot the figure
FIG=figure();
FIG.Position(3) = 2*FIG.Position(3);
FIG.PaperPosition(3) = 2*FIG.PaperPosition(3);
subplot(1,2,1)
shadedErrorBar(PC_val, PCC_all_mean_std(:,1),PCC_all_mean_std(:,2),{'-','color',ColorCode(5,:), 'LineWidth',2})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [0 100])
hold on
H=hline(100/NumInd, 'k--','Chance level');
H.LineWidth = 2;
hold on
V = vline(NPC_opt_all, 'b--', sprintf('Optimal #PC = %d', NPC_opt_all));
V.LineWidth = 2;
%% Permutation tests for both male and female calls
% Permutation of calls
fprintf(1, 'All calls Permutation test')
PC_val_rand = PC_val(1:min(length(PC_val),(NPC_opt_ind+NumPCsup)));
% PC_val_rand = PC_val;
Lrand_all = cell(length(PC_val_rand),1);
% PFisher = nan(1,NRandPerm);
for bb=1:NRandPerm
    fprintf(1, '\n Permutation %d/%d', bb, NRandPerm)
    RandInd = randperm(length(BatIDMic));
    for npc=1:length(PC_val_rand)
        fprintf(1, '  #PC = %d (%d/%d)', PC_val_rand(npc), npc, length(PC_val_rand))
        if bb==1
            Lrand_all{npc} = nan(NRandPerm,10);
        end
        Mdlrand = fitcdiscr(Score(:,1:PC_val_rand(npc)),BatIDMic(RandInd), 'CrossVal','on', 'KFold', 10, 'Prior', 1/NumInd.*ones(1,NumInd), 'SaveMemory', 'on', 'FillCoeffs', 'off');
        Lrand_all{npc}(bb,:) = kfoldLoss(Mdlrand, 'Mode', 'Individual');
        % if npc==NPC_opt_ind
        %     % perform a right tail Fisher's exact test to determine if the
        %     % odds of correct classification is higher in the observed case
        %     % than the permutation case
        %     CorrectRand = round((1-mean(Lrand_all{npc}(bb,:))).*size(Score,1));
        %     ErrorRand = round(mean(Lrand_all{npc}(bb,:)).*size(Score,1));
        %     Tbl = table([CorrectObs ; CorrectRand],[ErrorObs ; ErrorRand], 'VariableNames', {'Correct', 'Error'}, 'RowNames', {'Observed', 'Permutation'});
        %     [~,PFisher(bb),~] = fishertest(Tbl, 'Tail','right');
        % end
    end
end
PCCrand_all_mean_std = nan(length(PC_val_rand),2);
for npc=1:length(PC_val_rand)
    PCCrand_all_mean_std(npc,1) = mean(reshape(100*(1-Lrand_all{npc}),numel(Lrand_all{npc}),1));
    PCCrand_all_mean_std(npc,2) = std(reshape(100*(1-Lrand_all{npc}),numel(Lrand_all{npc}),1));
end


% Add the permutation values to the figure
FIG;
subplot(1,2,2)
shadedErrorBar(PC_val_rand, PCC_all_mean_std(1:length(PC_val_rand),1),PCC_all_mean_std(1:length(PC_val_rand),2),{'-','color',ColorCode(5,:), 'LineWidth',2})
xlabel('# PC')
ylabel('% Classification correct')
set(gca, 'YLim', [40 100])
hold on
V = vline(NPC_opt_all, 'b--', sprintf('Optimal #PC = %d', NPC_opt_all));
V.LineWidth = 2;
hold on
shadedErrorBar(PC_val_rand, PCCrand_all_mean_std(:,1),PCCrand_all_mean_std(:,2),{'-','color','k', 'LineWidth',2})
hold off
XLim = get(gca,'XLim');
text(XLim(2)/2,97, 'Observed data', 'Color', ColorCode(5,:))
text(XLim(2)/2,93, 'Full permutation', 'Color', 'k')
% title(sprintf('# Significant exact Fisher test (obs vs perm) at optimal #PC: %d/%d', sum(PFisher<PLim), length(PFisher)))


fprintf(1, 'Classification performance of calls along ID: %.1f%% +/-%.1f%%\nPermutation value: %.1f%% +/-%.1f%%\n#PC:%d \n', PCC_all_mean_std(NPC_opt_ind,1), PCC_all_mean_std(NPC_opt_ind,2),PCCrand_all_mean_std(NPC_opt_ind,1), PCCrand_all_mean_std(NPC_opt_ind,2), NPC_opt_all)
save(fullfile(LocalDataDir, 'DeafBats_RegularizedPermDFAID_AllCalls.mat'),  'PC_val', 'NPC_opt_all','L_all', 'PCC_all_mean_std', 'Lrand_all', 'PCCrand_all_mean_std', 'PFisher',  'MPS_mic_wf', 'IndWf', 'MPS_mic_wt', '-append' )

suplabel('All Calls DFA performance ID classification', 't');
for cc=1:length(FIG.Children)
    FIG.Children(cc).FontSize=12;
end
print(FIG, fullfile(GGPath.folder, GGPath.name,'My Drive/BatmanData/FigureLabMeeting/DeafBatsProject', 'RegPermDFAID_AllCalls.png') , '-dpng')
clear PC Mdl_all Mdlrand



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