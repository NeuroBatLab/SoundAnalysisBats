% This script look at the performance of ResNet at predicting vocalizations
% and noise
FS= 50000;
Date = '200130';
ExpStartTime = '0949';
Loggers_dir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/20200130/audiologgers';
ResNetDataFile = dir(fullfile(Loggers_dir, sprintf('*%s.csv', Date)));
ResNetData = readtable(fullfile(ResNetDataFile.folder, ResNetDataFile.name), 'Delimiter',',');
WinSize  = 100/1000*FS;
PredTruth = cell(size(36:44));
Discrepency=[];
for df=36:44
    fprintf(1,'\nFile %d/44', df)
    Filename = sprintf('%s_%s_VocExtractData%d_200.mat', Date, ExpStartTime, df);
    DataFile1 = fullfile(Loggers_dir, Filename);
    load(DataFile1, 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartPiezo', 'IndVocStopPiezo');
    DataFile2 = fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData%d.mat', Date, ExpStartTime, df));
    load(DataFile2, 'Piezo_wave', 'VocFilename')
    Fns_AL = fieldnames(Piezo_wave);
    PredTruth{df} = cell(1,length(IndVocStartPiezo_merged));
    for vv=1:length(IndVocStartPiezo_merged)
        if ~rem(vv,10)
            fprintf(1,'Voc %d/%d...', vv,length(IndVocStartPiezo_merged))
        end
        PredTruth{df}{vv} = cell(1,length(Fns_AL));
        for ll=1:length(Fns_AL)
            RowsData = logical(contains(ResNetData.Var1(:), sprintf('15120%s%s_%d_%s_', Date, Filename(1:11), df,Fns_AL{ll})) .*contains(ResNetData.Var1(:), sprintf('_%d.mat', vv-1))); % ResNetData.Var1 referes to files with a zero indexing
            pVoclogger = ResNetData.Var3(RowsData)'; % all probability of a vocalization across all 100ms windows of that logger on that sequence
            VocVector = zeros(1,ceil(length(Piezo_wave.(Fns_AL{ll}){vv})/FS*10)); % This is a vector of the smae size as the number of 100ms window that could fit in that sequence
            if ~isempty(IndVocStartPiezo_merged{vv})
                for ii=1:length(IndVocStartPiezo_merged{vv}{ll})
                    IndVocStart = ceil(IndVocStartPiezo_merged{vv}{ll}(ii)/FS*10);
                    IndVocStop = floor(IndVocStopPiezo_merged{vv}{ll}(ii)/FS*10);
                    VocVector(IndVocStart:IndVocStop) = ones(size(IndVocStart:IndVocStop));
                end
            end
            Discrepency = [Discrepency length(VocVector)-length(pVoclogger)];
            if Discrepency(end)
                keyboard
            end
            PredTruth{df}{vv}{ll} = [VocVector ; pVoclogger(1:length(VocVector))];
        end
        PredTruth{df}{vv} = [PredTruth{df}{vv}{:}];
    end
    PredTruth{df} = [PredTruth{df}{:}];
end
PredTruth = [PredTruth{:}];
figure()
swarmchart(PredTruth(1,:), PredTruth(2,:))
set(gca, 'XTick', [0 1],'XTickLabels', {'Noise' 'Voc'})
xlabel('Manual annotation')
ylabel('ResNet Probability of a vocalization')
fprintf(1, 'Correct classification of noise with Posterior proba of 0.5: %d/%d, %.2f%%\n', sum(PredTruth(2,~PredTruth(1,:))<0.5),sum(~PredTruth(1,:)), sum(PredTruth(2,~PredTruth(1,:))<0.5)/sum(~PredTruth(1,:))*100)
fprintf(1, 'Correct classification of Vocalizations with Posterior proba of 0.5: %d/%d, %.2f%%\n', sum(PredTruth(2,logical(PredTruth(1,:)))>0.5),sum(PredTruth(1,:)), sum(PredTruth(2,logical(PredTruth(1,:)))>0.5)/sum(PredTruth(1,:))*100)

PostProba = 0:0.01:0.5;
PCC=nan(2,length(PostProba));
for pp=1:length(PostProba)
    PCC(1,pp) = sum(PredTruth(2,~PredTruth(1,:))<PostProba(pp))/sum(~PredTruth(1,:));
    PCC(2,pp) = sum(PredTruth(2,logical(PredTruth(1,:)))>PostProba(pp))/sum(logical(PredTruth(1,:)));
end
figure()
plot(PostProba,PCC')
legend({'Noise' 'Vocalization'})
xlabel('ResNet Posterior probability')
ylabel('Probability of correct classification of 100ms windows')
