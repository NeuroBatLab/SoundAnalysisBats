function who_calls_soundproofBooth(PupID)

Path2Audio='/Volumes/server_home-1/users/old lab members/Sandra/Lesion/Lesionrecordings/Lesion_Groups/';
Path2Data='/Volumes/server_home/users/JulieE/DeafSalineGroup151/';

% Gather dates
load(fullfile(Path2Audio, 'DataBooths.mat'), 'ListFiles');
IDDates_all = nan(size(ListFiles));
for ffile=1:length(ListFiles)
    [~, Name,~]=fileparts(ListFiles{ffile});
    Ind_ = strfind(Name, '_');
    IDDates_all(ffile) = str2double([Name(Ind_(1)-1) Name(Ind_(1) + (3:8))]); % combine both the pup ID and the date of recording sprintf('%d%d',PupID,yymmdd) 
end
UIDDate = unique(IDDates_all);

% for each pup listen to 2h of data randomly sampled through the recordings
% and annotate if there are vocalizations or not
Nhours = 2;
Nsamples = Nhours*60*60*2;
if nargin<1
    PupID = input('What is the ID of the pup? 4,6,7 or8?\n'); %4 6 7 8
end
BundleSize = 10;
Fhigh_power =50; % Frequency upper bound for calculating the envelope (time running RMS)
Fs_env = 1000; % Sample frequency of the enveloppe
FS_out = 50000; % Hz
FS = 192000;
DB_noise = 60; % Noise threshold for the spectrogram colormap
FS_high = 50000; % Max Frequency for spectrogram plot
% design filters of raw ambient recording, bandpass 100Hz - 20kHz
[z,p,k] = butter(6,[100 23500]/(FS/2),'bandpass');
sos_raw_band = zp2sos(z,p,k);
NVoc = 100; % number of triggers in output file
VolFactorMic = 0.5;

% sound detection param
Consecutive_binsMic = 10;
AmpThresh = 0.01;

% total number of triggers for this pup
NTrig = sum((IDDates_all>(PupID*1000000 + 160000)).*(IDDates_all<(PupID*1000000 + 170000)));
Proportion_manual = Nsamples/NTrig;
UIDate_pup = find((UIDDate>(PupID*1000000 + 160000)).*(UIDDate<(PupID*1000000 + 170000)));
NDates = length(UIDate_pup);
Dates2Cure = UIDDate(UIDate_pup(sort([2 1:round(1/Proportion_manual):NDates])));
fprintf(1, "%d dates need to be cured for this pup #%d\n", length(Dates2Cure), PupID);

% Loop through dates and annotate them
StartDate=1;
for dd=StartDate:length(Dates2Cure)
    Date2cure_str = num2str(Dates2Cure(dd));
    fprintf(1, 'Recording day 20%s Pup%d (%d/%d)\n', Date2cure_str(2:end), PupID, dd, length(Dates2Cure))
    NumFilesOut = ceil(sum(IDDates_all==Dates2Cure(dd))/NVoc);
    FileIndices = find(IDDates_all==Dates2Cure(dd));
    StartFile = 1; %1
    for ff=StartFile:NumFilesOut
        fprintf(1, 'File of 100 voc, %d/%d\n', ff, NumFilesOut)
        if ff<NumFilesOut
            NVocLocal = NVoc;
            SampleBundle = [1:BundleSize:NVocLocal; BundleSize:BundleSize:NVocLocal];
        elseif ff==NumFilesOut
            NVocLocal = length(FileIndices) - NVoc*(ff-1);
            if ~rem(NVocLocal,BundleSize)
                SampleBundle = [1:BundleSize:NVocLocal; BundleSize:BundleSize:NVocLocal];
            else
                SampleBundle = [1:BundleSize:NVocLocal; BundleSize:BundleSize:NVocLocal NVocLocal];
            end
        end
        IndVocStart_all = cell(NVocLocal,1);
        IndVocStop_all = cell(NVocLocal,1);
        DistantCall = nan(NVocLocal,1);
        ListFiles_local = ListFiles(FileIndices((ff-1)*NVoc + (1:NVocLocal)));
        
        for bb=1:size(SampleBundle,2)
            fprintf(1, 'Triggers %d to %d   ', SampleBundle(1,bb), SampleBundle(2,bb))
            % first treat the whole bundle to go faster
            SBi = SampleBundle(1,bb):SampleBundle(2,bb);
            BundleSound = nan(1,96000*length(SBi));
            Time=cell(1,length(SBi));
            for bbi=1:length(SBi)
                Sound = load(ListFiles_local{SBi(bbi)});
                [~,FName,~] = fileparts(ListFiles_local{SBi(bbi)});
                Ind_T = strfind(FName, 'T');
                Time{bbi} = FName(Ind_T + (1:6));
                BundleSound((bbi-1)*96000 + (1:96000)) = Sound.recsGroup;
                FS = Sound.fs;
            end
            fprintf(1, 'Time %s to %s    ', Time{1}, Time{bbi});
            % filter to 23500kHz
            SoundFilt = filtfilt(sos_raw_band,1,BundleSound);
            % resample at 48kHz readble by computer
            SampleMic = resample((SoundFilt - mean(SoundFilt))/(std(SoundFilt)/VolFactorMic),FS/4,FS);
            figure(2)
            [~] = spec_only_bats(SampleMic, FS/4, DB_noise, 24000);
            title(sprintf('Trigger %d-%d',SampleBundle(1,bb), SampleBundle(2,bb)))
            % Play and take decision if a call or not
            Call1Noise0Bundle = playuntildecision(SampleMic, FS/4);
            
            % record results
            if Call1Noise0Bundle==0
                for bbi=1:length(SBi)
                    IndVocStart_all{SBi(bbi)} = cell(2,1);
                    IndVocStop_all{SBi(bbi)} = cell(2,1);
                    IndVocStart_all{SBi(bbi)}{1} = []; % weird organization, but that is to be compatible with ResNet data input
                    IndVocStop_all{SBi(bbi)}{1} = [];
                    DistantCall(SBi(bbi)) = 0;
                end
            else % there is a vocalization, analyse trigger by trigger
                Call1Noise0Trigger = nan(length(SBi),1);
                for bbi=1:length(SBi)
                    IndVocStart_all{SBi(bbi)} = cell(2,1);
                    IndVocStop_all{SBi(bbi)} = cell(2,1);
                    fprintf(1, 'Trigger %d   ', SBi(bbi))
                    SampleTrig = SampleMic((bbi-1)*FS/8 + (1:(FS/8)));
                    SampleTrigger = (SampleTrig - mean(SampleTrig))/(std(SampleTrig)/VolFactorMic);
                    Call1Noise0Trigger(SBi) = playuntildecision(SampleTrigger, FS/4);
                    if Call1Noise0Trigger(SBi) == 0
                        IndVocStart_all{SBi(bbi)}{1} = [];
                        IndVocStop_all{SBi(bbi)}{1} = [];
                        DistantCall(SBi(bbi)) = 0;
                    elseif Call1Noise0Trigger(SBi) == 2 % this is a distant call
                        DistantCall(SBi(bbi)) = 1;
                    elseif Call1Noise0Trigger(SBi) == 1
                        DistantCall(SBi(bbi)) = 0;
                        % calculate ampenvelope
                        Amp_env_Mic = running_rms(SoundFilt((bbi-1)*FS/2 + (1:(FS/2))), FS, Fhigh_power, Fs_env);
                        % Plot the spectrogram of the ambient microphone
                        F1=figure(1);
                        clf(F1)
                        ColorCode = [get(groot,'DefaultAxesColorOrder');0 0 0; 0 1 1; 1 1 0];
                        [~] = spec_only_bats(BundleSound((bbi-1)*FS/2 + (1:(FS/2))), FS, DB_noise, FS_high,300);
                        title(sprintf('Trigger %d/%d',bbi, length(SBi)))
                        hold on
                        yyaxis right
                        plot((1:length(Amp_env_Mic))/Fs_env*1000, Amp_env_Mic, 'r-', 'LineWidth',2)
                        ylabel(sprintf('Amp\nMic'))
                        YLIM = get(gca, 'YLim');
                        ylim([-YLIM(2)/20 YLIM(2)])

                        % detect vocalization onset/offset
                        VocpMic = Amp_env_Mic>(AmpThresh*max(Amp_env_Mic)); % Time points above amplitude threshold
                        IndVocStart = strfind(VocpMic, ones(1,Consecutive_binsMic)); %find the first indices of every sequences of length "Consecutive_bins" higher than RMS threshold
                        if isempty(IndVocStart)
                            fprintf(1,'No vocalization detected on microphone\n');
                            keyboard

                        else% Some vocalizations are automatically detected
                            IndVocStart_diffind = find(diff(IndVocStart)>1);
                            IndVocStart = [IndVocStart(1) IndVocStart(IndVocStart_diffind +1)]; % these two lines get rid of overlapping sequences that werer detected several times
                            NV = length(IndVocStart); % This is the number of detected potential vocalization
                            IndVocStop = nan(1,NV);
                            Call1Noise0_temp = nan(NV,1);
                            for ii=1:NV
                                IVStop = find(VocpMic(IndVocStart(ii):end)==0, 1, 'first');
                                if ~isempty(IVStop)
                                    IndVocStop(ii) = IndVocStart(ii) + IVStop -1;
                                else
                                    IndVocStop(ii) = length(VocpMic);
                                end
                                % Plot the localization of that sound
                                % extract on figure 1
                                figure(1)
                                hold on
                                yyaxis left; ylim([-1000 FS_high])
                                yyaxis right
                                plot([IndVocStart(ii)/Fs_env IndVocStop(ii)/Fs_env]*1000, [-0.5*YLIM(2)*1000/FS_high -0.5*YLIM(2)*1000/FS_high], 'k-', 'LineWidth',2)
                                hold off
                                Call1Noise0_temp(ii) = input('Indicate your choice: new call (1);    noise (0)\n');
                                
                                % update figure(1) with the decision
                                figure(1)
                                yyaxis right
                                hold on
                                Dur = IndVocStop(ii) - IndVocStart(ii);
                                if Call1Noise0_temp(ii)
                                    plot([IndVocStart(ii)/Fs_env IndVocStop(ii)/Fs_env]*1000, [-0.5*YLIM(2)*1000/FS_high -0.5*YLIM(2)*1000/FS_high], 'Color',ColorCode(5,:), 'LineWidth',2, 'LineStyle','-')
                                    text((IndVocStart(ii)+Dur/3)/Fs_env*1000, -YLIM(2)/20*0.75, 'New call','Color','k','FontWeight','bold')
                                else
                                    plot([IndVocStart(ii)/Fs_env IndVocStop(ii)/Fs_env]*1000, [-0.5*YLIM(2)*1000/FS_high -0.5*YLIM(2)*1000/FS_high], 'k-', 'LineWidth',2)
                                    text((IndVocStart(ii)+Dur/3)/Fs_env*1000, -YLIM(2)/20*0.75, 'Noise','Color','k','FontWeight','bold')
                                end
                                hold off
                            end

                            % Stash noise and only keep vocalizations
                            IndVocStart_all{SBi(bbi)}{1} = round(IndVocStart(logical(Call1Noise0_temp))/Fs_env *FS_out); % onset of vocalizations in samples with FS=50000
                            IndVocStop_all{SBi(bbi)}{1} = round(IndVocStop(logical(Call1Noise0_temp))/Fs_env *FS_out);% offset of vocalizations in samples with FS=50000

                        end
                    end
                end
            end
        end
        save(fullfile(Path2Data, ['20' Date2cure_str(2:end)], sprintf('%s_%s_VocExtractData%d.mat', Date2cure_str(2:end), repmat(Date2cure_str(1),1,4),ff)), "IndVocStart_all",'IndVocStop_all','DistantCall', '-append')
        clear IndVocStart_all IndVocStop_all DistantCall
    end
end
end


function Call1Noise0 = playuntildecision(Sound, FS)
PlayerMic= audioplayer(Sound, FS,24);
play(PlayerMic)
pause(length(Sound)/FS +1)

TempIn = input('Indicate your choice: call (1); distant call (2);  noise (0);  listen again to mic(any other number)\n');
if isempty(TempIn)
    fprintf(1, 'NO ENTRY GIVEN, playing again the sound and asking the same question\n')
    Call1Noise0 = 3;
else
    Call1Noise0 = TempIn;
end
while Call1Noise0~=0 && Call1Noise0~=1 && Call1Noise0~=2
%     pause(length(Sound)/FS/2)
    TempIn = input('Indicate your choice: call (1);    distant call (2);    noise (0);    listen again to mic(any other number)\n');
    play(PlayerMic)
    if isempty(TempIn)
        fprintf(1, 'NO ENTRY GIVEN, playing again the sound and asking the same question\n')
        Call1Noise0 = 3;
    else
        Call1Noise0 = TempIn;
    end
end
end
