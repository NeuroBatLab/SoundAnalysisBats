%% Load input data
% Path2Data = '/Volumes/server_home/users/JulieE/LMC_CoEd/logger/20190623';
% Path2Audio = '/Volumes/server_home/users/JulieE/LMC_CoEd/audio/20190623';
% % Load data manually extracted
% load(fullfile(Path2Data, '190623_1401_VocExtractData_200.mat'))
% load(fullfile(Path2Data, '190623_1401_VocExtractData.mat'))
% load(fullfile(Path2Audio, '190623_1401_VocExtractTimes.mat'))

% Path2Data1 = '/Users/ryanmoughan/Research/Debugging/20190927';
% Path2Data2=Path2Data1;
Path2Data1 = '/Volumes/server_home/users/JulieE/JuvenileRecordings155/20190927/audiologgers';
% Load data manually extracted
load(fullfile(Path2Data1,'190927_1014_VocExtractData_200.mat'))
load(fullfile(Path2Data1,'190927_1014_VocExtractData.mat'))
Path2Data2 = '/Volumes/server_home/users/JulieE/JuvenileRecordings155/20190927/audio';
load(fullfile(Path2Data2,'190927_1014_VocExtractTimes.mat'))

Fs_env=1000; % in Hertz, should have been saved in who_calls.m to correctly convert in time the starting and ending indices of vocalizations in IndVocStart
FS_Piezo = 50000; % could also be retrieved from Piezo_FS

% load results of the automatic detection
load(fullfile(Path2Data, 'allCallTimes.mat'))



% Duration of the data manually analysed
% ManDur = 10*60*10^3; %in ms
ManDur = 10*10*60*10^3; %in ms

%% Get in transceiver time the onset/offset of each vocalization manually
% extracted for each logger
AL_Id = fieldnames(Piezo_wave);
NLoggers = length(AudioLogs);
ManCallTimes = cell(1,NLoggers);
ManCallSamp = cell(1,NLoggers);
ManCallFile = cell(1,1400);
NFiles = length(IndVocStart_all);
voc_counter = 0;
for ff=1:NFiles
    [~,File]=fileparts(VocFilename{ff});
    Ind_ = strfind(File, '_');
    FileRaw = File(1:(Ind_(end)-1));
    % onset and offset of manually detected vocalization events for all
    % loggers
    IndVocStart_local = IndVocStart_all{ff};
    IndVocStop_local = IndVocStop_all{ff};
    
    % loop through loggers and append the 2 columns matrix of onset/offset
    % of vocalization events
    for ll=1:NLoggers
        if ~isempty(IndVocStart_local{ll})
            NewCalls = [IndVocStart_local{ll}' IndVocStop_local{ll}'] ./ Fs_env .* 10^3  + Voc_transc_time_refined(ff,1); % onset and offset time of the vocalization in ms in transceiver time
            ManCallTimes{ll} = [ManCallTimes{ll} ; NewCalls];
            
            
            NewCallsSample = [IndVocStart_local{ll}' IndVocStop_local{ll}'] ./ Fs_env .* FS  + Voc_samp_idx(ff,1); % onset and offset time of the vocalization in microphone (raw) samples from the beginning of each individual 10min file
            ManCallSamp{ll} = [ManCallSamp{ll} ; NewCallsSample];
            for voc=1:length(IndVocStart_local{ll})
                voc_counter = voc_counter+1;
                ManCallFile{voc_counter} = FileRaw;
            end
        end
    end
end
% ManCallSamples = cell2mat(ManCallSamp');
% ManCallFiles = ManCallFile(1:voc_counter);
% save('GroundTruthData.mat','ManCallFiles','ManCallSamples')
%% Massage the input of the automatic detection

% Find the boundaries in ms in transceiver time of the section of sound that was
% manually analysed
Delay2FirstFile = Voc_samp_idx(1,1)./FS.*10^3; %in ms (Voc_samp_idx(1,1) is the first sample number of the first isolated section of sound in the first 10 min file analyzed here
OnsetBoundary = Voc_transc_time(1,1) - Delay2FirstFile; % (Voc_Transc_time(1,1), is the absolute time onset in transceiver time of the first isolated section of sound in the first 10 min file analyzed here
OffsetBoundary = OnsetBoundary + ManDur;

% Select the events in the automatically detected results that are within
% the manually analysed window
for ll=1:NLoggers
    allCallTimes{ll} = cell2mat(allCallTimes{ll}').*10^(-3); % converting the cell array to a 2 column matrix and values from us to ms
    Idx_WBound = find(sum((allCallTimes{ll}>OnsetBoundary) .* (allCallTimes{ll}<OffsetBoundary), 2)>0); % any element that start or end within the boundaries is selected
    allCallTimes{ll} = allCallTimes{ll}(Idx_WBound,:);
end

% %% Plot (interesting but too few vocalizations to see nicely things)
% % order of loggers in auto detected events
% Ord = [4 1:3];
% 
% for ll=1:NLoggers
%     figure(ll)
%     for vv = 1:size(allCallTimes{Ord(ll)},1)
%         plot(allCallTimes{Ord(ll)}(vv,:), ones(2,1)*2, 'b-', 'LineWidth',2)
%         hold on
%     end
%     for vv = 1:size(ManCallTimes{ll},1)
%         plot(ManCallTimes{ll}(vv,:), ones(2,1), 'g-', 'LineWidth',1)
%         hold on
%     end
%     hold off
% end

%% Let's loop in the dataset

for ll=1:NLoggers
    for vv=1:size(ManCallTimes{ll},1)
        OnOffVoc = ManCallTimes{ll}(vv,:);
        
        Idx_WBound = find(sum((allCallTimes{ll}>OnOffVoc(1)-500) .* (allCallTimes{ll}<OnOffVoc(2)+500), 2)>0);
        
        if isempty(Idx_WBound)
            fprintf('No automatic call detected\n')
            disp(OnOffVoc)
        else
            fprintf('%d events correspond to that vocalization\n',length(Idx_WBound))
            figure(1)
            clf
            for ii=1:length(Idx_WBound)
                plot(allCallTimes{ll}(Idx_WBound(ii),:), ones(2,1), 'b-', 'LineWidth',2)
                hold on
            end
            plot(ManCallTimes{ll}(vv,:), ones(2,1) * 1.1, 'g-', 'LineWidth',1)
            hold off
            xlabel('Time in transceiver Time')
            legend({'Automatic', 'Manual'})
            title(sprintf('%s  Voc %d/%d', AL_Id{ll}, vv, size(ManCallTimes{ll},1)))
            ylim([0,2])
            pause()
        end
    end
end
    



%Call I'm looking at: 4.053667319626574e10   4.053670019626574e10 of logger 10
%closest file index is 23, corresponds to 184549377 index and
%4.054173326222046e10 start
%difference of starts is -5.060065954715014e+03 ms. (off by 40 ms in loop??)
%then the index difference should be -2.530032977357507e+05
%and the call should be at 1.842963737022642e+08 which will be 
%3.685927474045284e+06 on the graph
