Path2Audio='/Volumes/server_home-1/users/old lab members/Sandra/Lesion/Lesionrecordings/Lesion_Groups/';
Path2Data='/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
Pups = [4 6 7 8]; % Bats ID: HF4 (65696) DM1 (71351) HF6 (71353) DF2 (71354)
DOB = [20160819 20160914 20160924 20160927];
% load(fullfile(Path2Audio, 'DataBooths.mat'), 'ListFiles');
% IDDates_all = nan(size(ListFiles));
% for ffile=1:length(ListFiles)
%     [~, Name,~]=fileparts(ListFiles{ffile});
%     Ind_ = strfind(Name, '_');
%     IDDates_all(ffile) = str2double([Name(Ind_(1)-1) Name(Ind_(1) + (3:8))]); % combine both the pup ID and the date of recording sprintf('%d%d',PupID,yymmdd) 
% end
% UIDDate = unique(IDDates_all);
Debug = 1;
% Import biosound library
py.importlib.import_module('soundsig')
NVoc = 100; % max number of voc per file
VolFactorMic = 0.5;
FS= 192000; % Hz
DB_noise = 60;
F_high = 50000;
HighFc = 20000;
% design filters of raw ambient recording, bandpass 100Hz - 20kHz
[z,p,k] = butter(6,[100 50000]/(FS/2),'bandpass');
sos_raw_band = zp2sos(z,p,k);

for dd=1:length(UIDDate)
    fprintf(1, 'Recording day %d/%d\n', dd, length(UIDDate))
    NumFilesOut = ceil(sum(IDDates_all==UIDDate(dd))/NVoc);
    FileIndices = find(IDDates_all==UIDDate(dd));
    UIDDate_str = num2str(UIDDate(dd));
    PupID = num2str(UIDDate_str(1));
    for ff=1:NumFilesOut
        fprintf(1, 'File of 100 voc %d/%d\n', ff, NumFilesOut)
        if ff<NumFilesOut
            NVocLocal = NVoc;
        elseif ff==NumFilesOut
            NVocLocal = length(FileIndices) - NVoc*(ff-1);
        end
        BioSound_all = cell(NVocLocal,1);
        ListFiles_local = ListFiles(FileIndices((ff-1)*NVoc + (1:NVocLocal)));
        for vv=1:NVocLocal
            Sound = load(ListFiles_local{vv});
            [~,NameIn,~]=fileparts(ListFiles_local{vv});
            IndT = strfind(NameIn,'T');
            Time_voc = NameIn(IndT + (1:6));
            if Sound.fs~=FS
                warning('the extract does not have the expected sampling frequency\n')
                
            end
            % filter to 50kHz
            SoundFilt = filtfilt(sos_raw_band,1,Sound.recsGroup);
            % Calculate Biosound features
            BioSoundObj = runBiosound(SoundFilt, FS, F_high, HighFc);
            if Debug
                figure(1); clf
                plotBiosound(BioSoundObj, F_high, SoundFilt);
                suplabel(sprintf('Pup%s %s %s voc %d/%d', PupID, UIDDate_str(2:end),Time_voc,(ff-1)*NVoc + vv, length(FileIndices)),'t');
            end
            %save
            BioSound_all{vv} = BioSoundObj;
        end
        BioSound_all.(sprintf('Logger%s',PupID)) = BioSound_all;
        save(fullfile(Path2Data, ['20' UIDDate_str(2:end)], sprintf('%s_%s_VocExtractData%d.mat', UIDDate_str(2:end), repmat(UIDDate_str(1),1,4),ff)), "BioSound_all",'-append')
    end
end