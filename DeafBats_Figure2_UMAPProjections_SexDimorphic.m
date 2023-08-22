%% Set paths
LocalDataDir = '/Users/elie/Documents/DeafBats/Data';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
GGPath = dir('/Users/elie/Google Drive*');
Path2Paper = fullfile(GGPath.folder, GGPath.name, 'My Drive', 'BatmanData', 'Deaf Paper');

%% Load previous data
load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls2.mat'), 'BatID')
load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls.mat'),'MicAudioGood')
MicAudioGood01 = MicAudioGood;
MicAudioGood01(isnan(MicAudioGood01)) = 0;
MicAudioGood01 = logical(MicAudioGood01);

% load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_SexDimorph2.mat'), 'KSDistance','BatID_Pair_KSDistance')
load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_SexDimorph2.mat'), 'KSDistance_100PC','BatID_Pair_KSDistance_100PC')
load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_SexDimorph2.mat'), 'KSDistance_3PC','BatID_Pair_KSDistance_3PC')
% load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_SexDimorph2.mat'), 'KSDistance_UMAP','BatID_Pair_KSDistance_UMAP')
% load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_SexDimorph.mat'), 'SexDistance', 'BatID_Pair_SexDistance','NullDistance')
load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_SexDimorph.mat'), 'SexDistance_100PC','BatID_Pair_SexDistance_100PC','NullDistance_100PC')
load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_SexDimorph.mat'), 'SexDistance_3PC', 'BatID_Pair_SexDistance_3PC','NullDistance_3PC')
% load(fullfile(LocalDataDir, 'Data4_DeafBats_CatCalls_SexDimorph.mat'),'SexDistance_UMAP','BatID_Pair_SexDistance_UMAP')

%% Get the color vector ready for BatName, Sex and Deafness
ColorCode = [get(groot, 'DefaultAxesColorOrder'); 0 1 1; 0.5 0.5 0.5; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0 1; 0 0 0];
GGPath = dir('/Users/elie/Google Drive*');
try
    Path2RecordingTable = fullfile(GGPath.folder, GGPath.name,'My Drive/JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx');
    [~,~,RecTableData]=xlsread(Path2RecordingTable,2,'A1:k3','basic');
catch ME
    Path2Paper = fullfile(GGPath.folder, GGPath.name, 'Mon Drive', 'BatmanData', 'Deaf Paper');
    Path2RecordingTable = fullfile(GGPath.folder, GGPath.name,'Mon Drive/JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx');
    [~,~,RecTableData]=xlsread(Path2RecordingTable,2,'A1:k3','basic');
end

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
Sex_mic = Sex(MicAudioGood01);
SexDeaf_mic = SexDeaf(MicAudioGood01);
BatID_mic = str2double(BatID(MicAudioGood01));
Sex_mic_sign = 2.*contains(Sex_mic, 'M') -1;
NVoc = sum(MicAudioGood01);

%% Plot all data of distance
% UPair = cell(1,6);
% UPair{1} = unique(BatID_Pair_KSDistance_100PC)';
% UPair{2} = unique(BatID_Pair_SexDistance_100PC)';
% UPair{3} = unique(BatID_Pair_KSDistance_3PC)';
% UPair{4} = unique(BatID_Pair_SexDistance_3PC)';
% UPair{5} = unique(BatID_Pair_KSDistance_UMAP)';
% UPair{6} = unique(BatID_Pair_SexDistance_UMAP)';
% UPairall = unique([UPair{:}]);
% [Colors, GRAD, MinColVal, MaxColVal, Slide]=cubehelix_ColorVal(1:length(UPairall),1);
% 
% Fig19 = figure(19);
% clf
% subplot(3,2,1)
% ColPair = nan(NVoc,3);
% for pp=1:length(UPair{1})
%     LocalCol = Colors(strcmp(UPairall, UPair{1}{pp}),:);
%     ColPair(strcmp(BatID_Pair_KSDistance_100PC,UPair{1}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_KSDistance_100PC,UPair{1}{pp})),1);
% end
% 
% scatter(abs(KSDistance_100PC - NullDistance_100PC), abs(SexDistance_100PC - NullDistance_100PC), 20, ColPair, 'filled');
% title('Color: Opposite Treatment Distance')
% ylabel('Distance Sex (100PC)'); xlabel('Distance treatment (100PC)')
% XLim = xlim();
% YLim = ylim();
% for pp=1:length(UPair{1})
%     LocalCol = Colors(strcmp(UPairall, UPair{1}{pp}),:);
%     hold on
%     text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{1})*(pp-1), UPair{1}{pp},'Color', LocalCol, 'FontSize',12)
% end
% hold off
% 
% subplot(3,2,2)
% ColPair = nan(NVoc,3);
% for pp=1:length(UPair{2})
%     LocalCol = Colors(strcmp(UPairall, UPair{2}{pp}),:);
%     ColPair(strcmp(BatID_Pair_SexDistance_100PC,UPair{2}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_SexDistance_100PC,UPair{2}{pp})),1);
% end
% scatter(abs(KSDistance_100PC - NullDistance_100PC), abs(SexDistance_100PC - NullDistance_100PC), 20, ColPair, 'filled');
% title('Color: Opposite Sex Distance')
% ylabel('Distance Sex (100PC)'); xlabel('Distance treatment (100PC)')
% XLim = xlim();
% YLim = ylim();
% for pp=1:length(UPair{2})
%     LocalCol = Colors(strcmp(UPairall, UPair{2}{pp}),:);
%     hold on
%     text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{2})*(pp-1), UPair{2}{pp},'Color', LocalCol, 'FontSize',8)
% end
% hold off
% 
% subplot(3,2,3)
% ColPair = nan(NVoc,3);
% for pp=1:length(UPair{3})
%     LocalCol = Colors(strcmp(UPairall, UPair{3}{pp}),:);
%     ColPair(strcmp(BatID_Pair_KSDistance_3PC,UPair{3}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_KSDistance_3PC,UPair{3}{pp})),1);
% end
% scatter(abs(KSDistance_3PC - NullDistance_3PC), abs(SexDistance_3PC - NullDistance_3PC), 20, ColPair, 'filled');
% title('Color: Opposite Treatment Distance')
% ylabel('Distance Sex (3PC)'); xlabel('Distance treatment (3PC)')
% XLim = xlim();
% YLim = ylim();
% for pp=1:length(UPair{3})
%     LocalCol = Colors(strcmp(UPairall, UPair{3}{pp}),:);
%     hold on
%     text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{3})*(pp-1), UPair{3}{pp},'Color', LocalCol, 'FontSize',12)
% end
% hold off
% 
% subplot(3,2,4)
% ColPair = nan(NVoc,3);
% for pp=1:length(UPair{4})
%     LocalCol = Colors(strcmp(UPairall, UPair{4}{pp}),:);
%     ColPair(strcmp(BatID_Pair_SexDistance_3PC,UPair{4}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_SexDistance_3PC,UPair{4}{pp})),1);
% end
% scatter(abs(KSDistance_3PC - NullDistance_3PC), abs(SexDistance_3PC - NullDistance_3PC), 20, ColPair, 'filled');
% title('Color: Opposite Sex Distance')
% ylabel('Distance Sex (3PC)'); xlabel('Distance treatment (3PC)')
% XLim = xlim();
% YLim = ylim();
% for pp=1:length(UPair{4})
%     LocalCol = Colors(strcmp(UPairall, UPair{4}{pp}),:);
%     hold on
%     text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{4})*(pp-1), UPair{4}{pp},'Color', LocalCol, 'FontSize',8)
% end
% hold off
% 
% subplot(3,2,5)
% ColPair = nan(NVoc,3);
% for pp=1:length(UPair{5})
%     LocalCol = Colors(strcmp(UPairall, UPair{5}{pp}),:);
%     ColPair(strcmp(BatID_Pair_KSDistance_UMAP,UPair{5}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_KSDistance_UMAP,UPair{5}{pp})),1);
% end
% scatter(KSDistance_UMAP, SexDistance_UMAP, 20, ColPair, 'filled');
% title('Color: Opposite Treatment Distance')
% ylabel('Distance Sex (UMAP)'); xlabel('Distance treatment (UMAP)')
% XLim = xlim();
% YLim = ylim();
% for pp=1:length(UPair{5})
%     LocalCol = Colors(strcmp(UPairall, UPair{5}{pp}),:);
%     hold on
%     text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{5})*(pp-1), UPair{5}{pp},'Color', LocalCol, 'FontSize',12)
% end
% hold off
% 
% subplot(3,2,6)
% ColPair = nan(NVoc,3);
% for pp=1:length(UPair{6})
%     LocalCol = Colors(strcmp(UPairall, UPair{6}{pp}),:);
%     ColPair(strcmp(BatID_Pair_SexDistance_UMAP,UPair{6}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_SexDistance_UMAP,UPair{6}{pp})),1);
% end
% scatter(KSDistance_UMAP, SexDistance_UMAP, 20, ColPair, 'filled');
% title('Color: Opposite Sex Distance')
% ylabel('Distance Sex (UMAP)'); xlabel('Distance treatment (UMAP)')
% XLim = xlim();
% YLim = ylim();
% for pp=1:length(UPair{6})
%     LocalCol = Colors(strcmp(UPairall, UPair{6}{pp}),:);
%     hold on
%     text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{6})*(pp-1), UPair{6}{pp},'Color', LocalCol, 'FontSize',8)
% end
% hold off
% suplabel('Distance beyond expected by point density in PCA', 't');

%% Look if for each subject the calls that are the most sexually dimorphics are also the one that are the most affected by treatment
% correct the values of distance by the distance based on point density
% 3 first PCs
figure(20)
clf
% set(gcf,'Visible','on')
CSD = [1:5 11:15];
CTD = [6:10 16:20];
UPair = cell(1,2);
for bb=1:length(BatName)
    IndBat = BatID_mic==BatName(bb);
    UPair{2} = unique(BatID_Pair_KSDistance_3PC(IndBat))';
    UPair{1} = unique(BatID_Pair_SexDistance_3PC(IndBat))';
    UPairall = [UPair{:}];

    subplot(4,5,CSD(bb))
    ColPair = nan(sum(IndBat),3);
    for pp=1:length(UPair{1})
        LocalCol = ColorCode(strcmp(UPairall, UPair{1}{pp}),:);
        ColPair(strcmp(BatID_Pair_SexDistance_3PC(IndBat),UPair{1}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_SexDistance_3PC(IndBat),UPair{1}{pp})),1);
    end
    scatter(KSDistance_3PC(IndBat) - NullDistance_3PC(IndBat), SexDistance_3PC(IndBat)- NullDistance_3PC(IndBat), 20, ColPair, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5);
    title(sprintf('%s %d - Opposite Sex Distance', BatSexDeaf{bb}, BatName(bb)))
    ylabel('Distance Sex (3PC)'); xlabel('Distance treatment (3PC)')
    XLim = xlim();
    YLim = ylim();
    for pp=1:length(UPair{1})
        LocalCol = ColorCode(strcmp(UPairall, UPair{1}{pp}),:);
        hold on
        text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{1})*(pp-1), UPair{1}{pp},'Color', LocalCol, 'FontSize',12)
    end
    plot([0 min(XLim(2), YLim(2))], [0 min(XLim(2), YLim(2))], '-k')
    hold off

    subplot(4,5,CTD(bb))
    ColPair = nan(sum(IndBat),3);
    for pp=1:length(UPair{2})
        LocalCol = ColorCode(strcmp(UPairall, UPair{2}{pp}),:);
        ColPair(strcmp(BatID_Pair_KSDistance_3PC(IndBat),UPair{2}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_KSDistance_3PC(IndBat),UPair{2}{pp})),1);
    end
    scatter(KSDistance_3PC(IndBat)- NullDistance_3PC(IndBat), SexDistance_3PC(IndBat)- NullDistance_3PC(IndBat), 20, ColPair, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5);
    title(sprintf('%s %d - Opposite Treatment Distance', BatSexDeaf{bb},BatName(bb)))
    ylabel('Distance Sex (3PC)'); xlabel('Distance treatment (3PC)')
    XLim = xlim();
    YLim = ylim();
    for pp=1:length(UPair{2})
        LocalCol = ColorCode(strcmp(UPairall, UPair{2}{pp}),:);
        hold on
        text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{2})*(pp-1), UPair{2}{pp},'Color', LocalCol, 'FontSize',12)
    end
    plot([0 min(XLim(2), YLim(2))], [0 min(XLim(2), YLim(2))], '-k')
    hold off
end
suplabel('Distance beyond expected by point density in PCA', 't');

%% Look if for each subject the calls that are the most sexually dimorphics are also the one that are the most affected by treatment
% correct the values of distance by the distance based on point density
% 100 first PCs
figure(20)
clf
% set(gcf,'Visible','on')
CSD = [1:5 11:15];
CTD = [6:10 16:20];
UPair = cell(1,2);
for bb=1:length(BatName)
    IndBat = BatID_mic==BatName(bb);
    UPair{2} = unique(BatID_Pair_KSDistance_100PC(IndBat))';
    UPair{1} = unique(BatID_Pair_SexDistance_100PC(IndBat))';
    UPairall = [UPair{:}];

    subplot(4,5,CSD(bb))
    ColPair = nan(sum(IndBat),3);
    for pp=1:length(UPair{1})
        LocalCol = ColorCode(strcmp(UPairall, UPair{1}{pp}),:);
        ColPair(strcmp(BatID_Pair_SexDistance_100PC(IndBat),UPair{1}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_SexDistance_100PC(IndBat),UPair{1}{pp})),1);
    end
    scatter(KSDistance_100PC(IndBat) - NullDistance_100PC(IndBat), SexDistance_100PC(IndBat)- NullDistance_100PC(IndBat), 20, ColPair, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5);
    title(sprintf('%s %d - Opposite Sex Distance', BatSexDeaf{bb}, BatName(bb)))
    ylabel('Distance Sex (100PC)'); xlabel('Distance treatment (100PC)')
    XLim = xlim();
    YLim = ylim();
    for pp=1:length(UPair{1})
        LocalCol = ColorCode(strcmp(UPairall, UPair{1}{pp}),:);
        hold on
        text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{1})*(pp-1), UPair{1}{pp},'Color', LocalCol, 'FontSize',12)
    end
    plot([0 min(XLim(2), YLim(2))], [0 min(XLim(2), YLim(2))], '-k')
    hold off

    subplot(4,5,CTD(bb))
    ColPair = nan(sum(IndBat),3);
    for pp=1:length(UPair{2})
        LocalCol = ColorCode(strcmp(UPairall, UPair{2}{pp}),:);
        ColPair(strcmp(BatID_Pair_KSDistance_100PC(IndBat),UPair{2}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_KSDistance_100PC(IndBat),UPair{2}{pp})),1);
    end
    scatter(KSDistance_100PC(IndBat)- NullDistance_100PC(IndBat), SexDistance_100PC(IndBat)- NullDistance_100PC(IndBat), 20, ColPair, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5);
    title(sprintf('%s %d - Opposite Treatment Distance', BatSexDeaf{bb},BatName(bb)))
    ylabel('Distance Sex (100PC)'); xlabel('Distance treatment (100PC)')
    XLim = xlim();
    YLim = ylim();
    for pp=1:length(UPair{2})
        LocalCol = ColorCode(strcmp(UPairall, UPair{2}{pp}),:);
        hold on
        text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{2})*(pp-1), UPair{2}{pp},'Color', LocalCol, 'FontSize',12)
    end
    plot([0 min(XLim(2), YLim(2))], [0 min(XLim(2), YLim(2))], '-k')
    hold off
end
suplabel('Distance beyond expected by point density in PCA', 't');
%% Plot all bats
figure(21)
clf
subplot(1,2,1)
UPair = cell(1,2);
UPair{1} = unique(BatID_Pair_KSDistance_3PC)';
UPair{2} = unique(BatID_Pair_SexDistance_3PC)';
UPairall = unique([UPair{:}]);
[Colors, GRAD, MinColVal, MaxColVal, Slide]=cubehelix_ColorVal(1:length(UPairall),1);

ColPair = nan(NVoc,3);
for pp=1:length(UPair{1})
    LocalCol = Colors(strcmp(UPairall, UPair{1}{pp}),:);
    ColPair(strcmp(BatID_Pair_KSDistance_3PC,UPair{1}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_KSDistance_3PC,UPair{1}{pp})),1);
end
scatter(KSDistance_3PC - NullDistance_3PC, SexDistance_3PC - NullDistance_3PC, 20, ColPair, 'filled');
title('Color: Opposite Treatment Distance')
ylabel('Distance Sex (3PC)'); xlabel('Distance treatment (3PC)')
XLim = xlim();
YLim = ylim();
for pp=1:length(UPair{1})
    LocalCol = Colors(strcmp(UPairall, UPair{1}{pp}),:);
    hold on
    text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{1})*(pp-1), UPair{1}{pp},'Color', LocalCol, 'FontSize',12)
end
hold off

subplot(1,2,2)
ColPair = nan(NVoc,3);
for pp=1:length(UPair{2})
    LocalCol = Colors(strcmp(UPairall, UPair{2}{pp}),:);
    ColPair(strcmp(BatID_Pair_SexDistance_3PC,UPair{2}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_SexDistance_3PC,UPair{2}{pp})),1);
end
scatter(KSDistance_3PC - NullDistance_3PC, SexDistance_3PC - NullDistance_3PC, 20, ColPair, 'filled');
title('Color: Opposite Sex Distance')
ylabel('Distance Sex (3PC)'); xlabel('Distance treatment (3PC)')
XLim = xlim();
YLim = ylim();
for pp=1:length(UPair{2})
    LocalCol = Colors(strcmp(UPairall, UPair{2}{pp}),:);
    hold on
    text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{2})*(pp-1), UPair{2}{pp},'Color', LocalCol, 'FontSize',8)
end
hold off

%% Look if for each subject the calls that are the most sexually dimorphics are also the one that are the most affected by treatment
% correct the values of distance by the distance based on point density
figure(22)
clf
% set(gcf,'Visible','on')
CSD = [1:5 11:15];
CTD = [6:10 16:20];
UPair = cell(1,2);
for bb=1:length(BatName)
    IndBat = BatID_mic==BatName(bb);
    UPair{2} = unique(BatID_Pair_KSDistance_3PC(IndBat))';
    UPair{1} = unique(BatID_Pair_SexDistance_3PC(IndBat))';
    UPairall = [UPair{:}];

    subplot(4,5,CSD(bb))
    ColPair = nan(sum(IndBat),3);
    for pp=1:length(UPair{1})
        LocalCol = ColorCode(strcmp(UPairall, UPair{1}{pp}),:);
        ColPair(strcmp(BatID_Pair_SexDistance_3PC(IndBat),UPair{1}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_SexDistance_3PC(IndBat),UPair{1}{pp})),1);
    end
    scatter(KSDistance_3PC(IndBat)./NullDistance_3PC(IndBat), SexDistance_3PC(IndBat)./NullDistance_3PC(IndBat), 20, ColPair, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5);
    title(sprintf('%s %d - Opposite Sex Distance', BatSexDeaf{bb}, BatName(bb)))
    ylabel('Distance Sex (3PC)'); xlabel('Distance treatment (3PC)')
    XLim = xlim();
    YLim = ylim();
    for pp=1:length(UPair{1})
        LocalCol = ColorCode(strcmp(UPairall, UPair{1}{pp}),:);
        hold on
        text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{1})*(pp-1), UPair{1}{pp},'Color', LocalCol, 'FontSize',12)
    end
    plot([min(XLim(1), YLim(1)) min(XLim(2), YLim(2))], [min(XLim(1), YLim(1)) min(XLim(2), YLim(2))], '-k')
    hold off

    subplot(4,5,CTD(bb))
    ColPair = nan(sum(IndBat),3);
    for pp=1:length(UPair{2})
        LocalCol = ColorCode(strcmp(UPairall, UPair{2}{pp}),:);
        ColPair(strcmp(BatID_Pair_KSDistance_3PC(IndBat),UPair{2}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_KSDistance_3PC(IndBat),UPair{2}{pp})),1);
    end
    scatter(KSDistance_3PC(IndBat)./NullDistance_3PC(IndBat), SexDistance_3PC(IndBat)./NullDistance_3PC(IndBat), 20, ColPair, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5);
    title(sprintf('%s %d - Opposite Treatment Distance', BatSexDeaf{bb},BatName(bb)))
    ylabel('Distance Sex (3PC)'); xlabel('Distance treatment (3PC)')
    XLim = xlim();
    YLim = ylim();
    for pp=1:length(UPair{2})
        LocalCol = ColorCode(strcmp(UPairall, UPair{2}{pp}),:);
        hold on
        text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{2})*(pp-1), UPair{2}{pp},'Color', LocalCol, 'FontSize',12)
    end
    plot([min(XLim(1), YLim(1)) min(XLim(2), YLim(2))], [min(XLim(1), YLim(1)) min(XLim(2), YLim(2))], '-k')
    hold off
end
suplabel('Distance normalized by point density in PCA', 't');

%% Look if for each subject the calls that are the most sexually dimorphics are also the one that are the most affected by treatment
% correct the values of distance by the distance based on point density
% 100 first PCs
figure(22)
clf
% set(gcf,'Visible','on')
CSD = [1:5 11:15];
CTD = [6:10 16:20];
UPair = cell(1,2);
for bb=1:length(BatName)
    IndBat = BatID_mic==BatName(bb);
    UPair{2} = unique(BatID_Pair_KSDistance_100PC(IndBat))';
    UPair{1} = unique(BatID_Pair_SexDistance_100PC(IndBat))';
    UPairall = [UPair{:}];

    subplot(4,5,CSD(bb))
    ColPair = nan(sum(IndBat),3);
    for pp=1:length(UPair{1})
        LocalCol = ColorCode(strcmp(UPairall, UPair{1}{pp}),:);
        ColPair(strcmp(BatID_Pair_SexDistance_100PC(IndBat),UPair{1}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_SexDistance_100PC(IndBat),UPair{1}{pp})),1);
    end
    scatter(KSDistance_100PC(IndBat)./NullDistance_100PC(IndBat), SexDistance_100PC(IndBat)./NullDistance_100PC(IndBat), 20, ColPair, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5);
    title(sprintf('%s %d - Opposite Sex Distance', BatSexDeaf{bb}, BatName(bb)))
    ylabel('Distance Sex (100PC)'); xlabel('Distance treatment (100PC)')
    XLim = xlim();
    YLim = ylim();
    for pp=1:length(UPair{1})
        LocalCol = ColorCode(strcmp(UPairall, UPair{1}{pp}),:);
        hold on
        text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{1})*(pp-1), UPair{1}{pp},'Color', LocalCol, 'FontSize',12)
    end
    plot([min(XLim(1), YLim(1)) min(XLim(2), YLim(2))], [min(XLim(1), YLim(1)) min(XLim(2), YLim(2))], '-k')
    hold off

    subplot(4,5,CTD(bb))
    ColPair = nan(sum(IndBat),3);
    for pp=1:length(UPair{2})
        LocalCol = ColorCode(strcmp(UPairall, UPair{2}{pp}),:);
        ColPair(strcmp(BatID_Pair_KSDistance_100PC(IndBat),UPair{2}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_KSDistance_100PC(IndBat),UPair{2}{pp})),1);
    end
    scatter(KSDistance_100PC(IndBat)./ NullDistance_100PC(IndBat), SexDistance_100PC(IndBat)./NullDistance_100PC(IndBat), 20, ColPair, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5);
    title(sprintf('%s %d - Opposite Treatment Distance', BatSexDeaf{bb},BatName(bb)))
    ylabel('Distance Sex (100PC)'); xlabel('Distance treatment (100PC)')
    XLim = xlim();
    YLim = ylim();
    for pp=1:length(UPair{2})
        LocalCol = ColorCode(strcmp(UPairall, UPair{2}{pp}),:);
        hold on
        text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{2})*(pp-1), UPair{2}{pp},'Color', LocalCol, 'FontSize',12)
    end
    plot([min(XLim(1), YLim(1)) min(XLim(2), YLim(2))], [min(XLim(1), YLim(1)) min(XLim(2), YLim(2))], '-k')
    hold off
end
suplabel('Distance normalized by point density in PCA', 't');
%% Plot all bats
figure(23)
clf
subplot(1,2,1)
UPair = cell(1,2);
UPair{1} = unique(BatID_Pair_KSDistance_3PC)';
UPair{2} = unique(BatID_Pair_SexDistance_3PC)';
UPairall = unique([UPair{:}]);
[Colors, GRAD, MinColVal, MaxColVal, Slide]=cubehelix_ColorVal(1:length(UPairall),1);

ColPair = nan(NVoc,3);
for pp=1:length(UPair{1})
    LocalCol = Colors(strcmp(UPairall, UPair{1}{pp}),:);
    ColPair(strcmp(BatID_Pair_KSDistance_3PC,UPair{1}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_KSDistance_3PC,UPair{1}{pp})),1);
end
scatter(KSDistance_3PC./NullDistance_3PC, SexDistance_3PC./NullDistance_3PC, 20, ColPair, 'filled');
title('Color: Opposite Treatment Distance')
ylabel('Distance Sex (3PC)'); xlabel('Distance treatment (3PC)')
XLim = xlim();
YLim = ylim();
for pp=1:length(UPair{1})
    LocalCol = Colors(strcmp(UPairall, UPair{1}{pp}),:);
    hold on
    text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{1})*(pp-1), UPair{1}{pp},'Color', LocalCol, 'FontSize',12)
end
hold off

subplot(1,2,2)
ColPair = nan(NVoc,3);
for pp=1:length(UPair{2})
    LocalCol = Colors(strcmp(UPairall, UPair{2}{pp}),:);
    ColPair(strcmp(BatID_Pair_SexDistance_3PC,UPair{2}{pp}),:) = repmat(LocalCol,sum(strcmp(BatID_Pair_SexDistance_3PC,UPair{2}{pp})),1);
end
scatter(KSDistance_3PC./NullDistance_3PC, SexDistance_3PC./NullDistance_3PC, 20, ColPair, 'filled');
title('Color: Opposite Sex Distance')
ylabel('Distance Sex (3PC)'); xlabel('Distance treatment (3PC)')
XLim = xlim();
YLim = ylim();
for pp=1:length(UPair{2})
    LocalCol = Colors(strcmp(UPairall, UPair{2}{pp}),:);
    hold on
    text(XLim(2)*0.75, YLim(2)*0.9 - YLim(2)*0.9/length(UPair{2})*(pp-1), UPair{2}{pp},'Color', LocalCol, 'FontSize',8)
end
suplabel('Normalized values by density', 't')
hold off

%% Look at the histograms of Corrected distance values
figure(24)
subplot(2,2,1)
histogram(SexDistance_3PC-NullDistance_3PC)
xlabel('Sex Distance beyond expected by point density')
subplot(2,2,2)
histogram(SexDistance_3PC./NullDistance_3PC)
xlabel('Sex Distance normalized by point density')
subplot(2,2,3)
histogram(KSDistance_3PC-NullDistance_3PC)
xlabel('Treatment Distance beyond expected by point density')
subplot(2,2,4)
histogram(KSDistance_3PC./NullDistance_3PC)
xlabel('Treatment Distance normalized by point density')

%% Now run a GLM to quantify the effect: abs(KSDistance_3PC - NullDistance_3PC) ~ abs(SexDistance_3PC - NullDistance_3PC) * Sex + (1|BatID)

%% INTERNAL FUNCTIONS

function [Colors, GRAD, MinColVal, MaxColVal, Slide]=cubehelix_ColorVal(VARcol,x, MinColVal, MaxColVal)
    % Getting VARZ ready for plotting.
    % If VARcol is always a non-null
    % integer then you can use it directly in cubehelix by specifying x=1
    % otherwise, the function first multiplies VARcol by 10^x to make sure its
    % variability is maximized after ceiling it to an integer
    if nargin<2
        x=3;% here I'm doing *1000 as an example
    end

    if nargin<3
        MinColVal=0;
    end

    if nargin<4
        MaxColVal = max(VARcol);
    end

%     if nargin<5
%         LogToggle=0;
%     end
    % Make sure you don't have nan values
    VARcol = VARcol(~isnan(VARcol));
%     if LogToggle==2
%         VARcol = log2(VARcol);
%     elseif LogToggle==10
%         VARcol = log10(VARcol);
%     end
    

    % Scale your minimum value as 0 or min Val
    Slide = 0;
    if ~MinColVal && sum(VARcol<=0)>0
        Slide = -min(VARcol)+ 10^(-x);
    elseif MinColVal==1
        Slide = -min(VARcol)+ 10^(-x);
    elseif MinColVal
        Slide = -MinColVal+ 10^(-x);
    end

    VARcol_cube = ceil((VARcol+Slide)*10^x);
    %GRAD=cubehelix(max(VARZ_cube));%Note that you can choose the values of...
    ...your gradient by doing cubehelix_view and then specifying the values in...
        ...the function, for instance: GRAD=cubehelix(max(VARZ_cube),0.5, -0.8, 1.5, 1.7, [1,0]);
%     GRAD=cubehelix(max(VARcol_cube),0.5, -1.1, 1.5, 0.5, [1,0]);
%     GRAD=cubehelix(max(VARcol_cube),0.5, -0.8, 1.5, 1.7, [1,0]);
    GRAD=cubehelix(ceil((MaxColVal+Slide)*10^x),0.5, -1.5, 1, 1, [1,0]);

    NU = length(VARcol_cube); % this is the number of points you have in your plot
    Colors = nan(NU, size(GRAD,2));
    for jj=1:NU
        Colors(jj,:) = GRAD(VARcol_cube(jj),:); % for each dot you call the line in the gradient of color that correspond to the value in Z
    end
    
% figure()
%     
%     cc=colorbar();
%     colormap(GRAD)
%     set(gca, 'Color', [0.7 0.7 0.7])
%     
%     YTL = get(cc, 'YTickLabel');
%     YTL_new = ceil((MaxColVal+Slide)*10^x) .* str2double(YTL)/10^x - Slide;
%     if LogToggle==2
%         YTL_new = 2.^YTL_new;
%     elseif LogToggle==10
%         YTL_new = 10.^YTL_new;
%     end
%     YTL_new = num2str(round(YTL_new.*10^(x-1))./10^(x-1));
%     set(cc, 'YTickLabel',YTL_new)%here you correct the value of the z axis that you artificially multiplied by 10^x

end