function audioQuality_correct_logger_data_voc(Loggers_dir, Date, ExpStartTime)

% Load data
DataFiles = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractDat*_*.mat', Date, ExpStartTime)));

% Loop through the datafiles
for df=1:length(DataFiles) 
    DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData%d_200.mat', Date, ExpStartTime,df)));
    fprintf(1,'Set %d/%d\nwith file %s\n', df, length(DataFiles), DataFile.name)
    load(fullfile(DataFile.folder, DataFile.name), 'ManualAnnotationOK', 'BioSoundFilenames','BioSoundCalls','IndVocStartPiezo','IndVocStopPiezo','IndVocStartRaw_merged','IndVocStartPiezo_merged','IndVocStopPiezo_merged','BatID','LoggerName');
    if ~exist('BatID', 'var')
        Ind_ = strfind(DataFile.name, '_');
        load(fullfile(DataFile.folder, sprintf('%s_%s_VocExtractData%d%s', Date, ExpStartTime,1, DataFile.name(Ind_(end):end))), 'BatID','LoggerName');
    end
    load(fullfile(DataFile.folder, [DataFile.name(1:end-8) '.mat']), 'VocFilename', 'Piezo_wave')
    IndVocStartPiezoAnnotationOk = cell(size(IndVocStartPiezo));
    IndVocStopPiezoAnnotationOk = cell(size(IndVocStopPiezo));

    % Number (NV) and indices (VocInd) of call sequences with identified vocalizations
    VocInd = find(~cellfun('isempty',IndVocStartRaw_merged));
    VocInd_true = nan(size(VocInd));
    for vvi=1:length(VocInd)
        VocInd_true(vvi) = any(~cellfun('isempty',IndVocStartRaw_merged{VocInd(vvi)}));
    end
    VocInd = VocInd(logical(VocInd_true));
    VocInd2 = find(~cellfun('isempty',IndVocStartPiezo));
    VocInd2_true = nan(size(VocInd2));
    for vvi=1:length(VocInd2)
        VocInd2_true(vvi) = any(~cellfun('isempty',IndVocStartPiezo{VocInd2(vvi)}));
    end
    VocInd2 = VocInd2(logical(VocInd2_true));
    if length(VocInd) ~= length(VocInd2) || any(VocInd~=VocInd2)
        %             keyboard % there should be the same sequences labelled as having vocalizationss
        VocInd = intersect(VocInd, VocInd2);
    end
    NV = length(VocInd);

    % Names of loggers
    Fns_AL = fieldnames(Piezo_wave);
    clear Piezo_wave
    

    NVocFile = 0; %counter that keeps track of the indices value in ManualAnnotationOK and BioSoundFilenames
    VocCall = size(BioSoundFilenames,1); % this is the expected number of vocal sequence for that file
    for vv=1:NV
        if strfind(VocFilename{VocInd(vv)}, '/')
            [~,FileVoc]=fileparts(VocFilename{VocInd(vv)});
        else
            ParseVocFile = strsplit(VocFilename{VocInd(vv)}, '\');
            FileVoc = ParseVocFile{end};
        end
        IndVocStartPiezoAnnotationOk{VocInd(vv)} = cell(size(IndVocStartPiezo{VocInd(vv)}));
        IndVocStopPiezoAnnotationOk{VocInd(vv)} = cell(size(IndVocStopPiezo{VocInd(vv)}));
        for ll=1:length(Fns_AL)
            % Logger number
            AL_local = Fns_AL{ll};
            ALNum = AL_local(7:end);
            % ID of the bat
            ALIndex = contains(LoggerName, 'AL') .* contains(LoggerName, ALNum);
            BatID_local =BatID{find(ALIndex)}; %#ok<FNDSB>
            Ncall = length(IndVocStartRaw_merged{VocInd(vv)}{ll});
            if Ncall ~= length(IndVocStartPiezo_merged{VocInd(vv)}{ll})
                keyboard
            end
            if Ncall
                LocalStart = cell(1,Ncall);
                LocalStop = cell(1,Ncall);
                for nn=1:Ncall
                        NVocFile = NVocFile +1;
                        fprintf(1,'\n\n%d/%d Vocalization\n',NVocFile,VocCall)
                        LocalBioSoundFN = sprintf('%s_Bat%d_AL%s_Elmt%d_Piezo',FileVoc, BatID_local,ALNum,nn);
                        if isempty(BioSoundFilenames{NVocFile,2})
                            warning('Issue of allignment that resulted in truncated piezo data, ignore')
                            continue
                        end
                        [~,BioSoundFilename,~] =fileparts(BioSoundFilenames{NVocFile,2});
%                         open([BioSoundFilenames{NVocFile,2}(1:end-4) '.pdf'])
                        if ~strcmp(BioSoundFilename, LocalBioSoundFN)
                            warning('Problem of indexing!')
                            keyboard
                        end
                        % Find out if that sound was taken into elmts
                        if isfield(BioSoundCalls{NVocFile,2}, 'OnOffSets_elmts')
                            IndVocStartPiezo_local = BioSoundCalls{NVocFile,2}.OnOffSets_elmts(:,1) + IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn) -1;
                            IndVocStopPiezo_local = BioSoundCalls{NVocFile,2}.OnOffSets_elmts(:,2) + IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn) -1;
                        else
                            IndVocStartPiezo_local = IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn);
                            IndVocStopPiezo_local = IndVocStopPiezo_merged{VocInd(vv)}{ll}(nn);
                        end
                        if isempty(ManualAnnotationOK{NVocFile})
                            ManualAnnotationOK{NVocFile} = zeros(size(IndVocStopPiezo_local));
                        end
                        if length(ManualAnnotationOK{NVocFile}) ~= length(IndVocStartPiezo_local)
                            warning('Problem of indexing!')
                            keyboard
                        end
                        LocalStart{nn} = IndVocStartPiezo_local(ManualAnnotationOK{NVocFile}==1)';
                        LocalStop{nn} = IndVocStopPiezo_local(ManualAnnotationOK{NVocFile}==1)';
                end
                IndVocStartPiezoAnnotationOk{VocInd(vv)}{ll} = [LocalStart{:}];
                IndVocStopPiezoAnnotationOk{VocInd(vv)}{ll} = [LocalStop{:}];

            end  
        end
    end
    save(fullfile(DataFile.folder, DataFile.name), 'IndVocStartPiezoAnnotationOk','IndVocStopPiezoAnnotationOk', '-append')
end