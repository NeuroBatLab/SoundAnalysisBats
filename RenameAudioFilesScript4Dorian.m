DirIn = 'C:\Users\tobias\Documents\BioSoundTutorial-master\BatCallsOld';
DirOut = 'C:\Users\tobias\Documents\BioSoundTutorial-master\BatCalls';
AllFiles = dir(fullfile(DirIn, '*.wav'));
NFiles = length(AllFiles);
BatStatus.name = [11648 14461 14463 14464 65696 71043 71047 71351 71353 71354];
BatStatus.sex = {'F' 'M' 'F' 'M' 'F' 'M' 'F' 'M' 'F' 'F'};
BatStatus.deaf = [0 0 1 0 0 1 1 1 0 1];
for ff=1:NFiles
    fprintf(1, 'Copying File %d/%d\n', ff, NFiles)
    BatInd = strfind(AllFiles(ff).name, 'Bat');
    BatID = str2double(AllFiles(ff).name(BatInd+(3:7)));
    BatSex = BatStatus.sex{BatStatus.name == BatID};
    BatDeaf = BatStatus.deaf(BatStatus.name == BatID);
    if BatDeaf
        NewName = sprintf('%s%d_Kanamycin_%s', BatSex, BatID, AllFiles(ff).name);
    else
        NewName = sprintf('%s%d_Saline_%s', BatSex, BatID, AllFiles(ff).name);
    end
    copyfile(fullfile(DirIn, AllFiles(ff).name), fullfile(DirOut, NewName), 'f')
end