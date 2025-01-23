Path2Data='/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
ListDir = dir(fullfile(Path2Data, '2016*'));
StartDate = 1;
for dd=StartDate:length(ListDir)
    if str2double(ListDir(dd).name)>20160800
        fprintf(1, '%s   date %d/%d\n', ListDir(dd).name, dd, length(ListDir))
        ListFiles = dir(fullfile(Path2Data, ListDir(dd).name, '16*'));
        StartFile = 1;
        for ff=StartFile:length(ListFiles)
            fprintf(1, 'File %d/%d\n', ff,length(ListFiles))
            Data =load(fullfile(Path2Data, ListDir(dd).name, ListFiles(ff).name));
            if isfield(Data, 'IndVocStart_all')
                PupID = str2double(ListFiles(ff).name(8));
                IndVocStart_all = cell(length(Data.IndVocStart_all),1);
                IndVocStop_all = cell(length(Data.IndVocStop_all),1);
                for vv=1:length(Data.IndVocStart_all)
                    IndVocStart_all{vv} = cell(2,1);
                    IndVocStop_all{vv} = cell(2,1);
                    if ~iscell(Data.IndVocStart_all{vv})
                        if ~isempty(Data.IndVocStart_all{vv})
                            IndVocStart_all{vv}{1} = Data.IndVocStart_all{vv};
                            IndVocStop_all{vv}{1} = Data.IndVocStop_all{vv};
                        end
                    elseif any(~cellfun('isempty', Data.IndVocStart_all{vv})) || any(~cellfun('isempty', Data.IndVocStop_all{vv})) % there were vocalizations
                        IndVocStart = find(~cellfun('isempty', Data.IndVocStart_all{vv}));
                        IndVocStop = find(~cellfun('isempty', Data.IndVocStop_all{vv}));
                        if any(IndVocStart~=IndVocStop)
                            error('issue!!')
                        else
                            IndVoc = IndVocStart;
                        end
                        if length(IndVoc)>1
                            error("that should be 1")
                        end

%                     if ~iscell(Data.IndVocStart_all{vv})
%                         if isempty(Data.IndVocStart_all{vv})
%                             Data.IndVocStop_all{vv} = cell(1,8);
%                             Data.IndVocStart_all{vv} = cell(1,8);
%                         else
%                             keyboard
%                         end
%                     end
                        IndVocStart_all{vv}{1} = Data.IndVocStart_all{vv}{IndVoc};
                        IndVocStop_all{vv}{1} = Data.IndVocStop_all{vv}{IndVoc};
                    end
                end
                save(fullfile(Path2Data, ListDir(dd).name, ListFiles(ff).name), 'IndVocStop_all', 'IndVocStart_all', '-append')
            end
            clear Data
        end
    end
end