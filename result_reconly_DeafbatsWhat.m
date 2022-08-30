% This script is about running What_calls on Deaf bats data manually curated by
% CallCura

% Indicate the paths to the code and to the data
BaseDataDir = 'Z:\users\JulieE\DeafSalineGroup151\';

% Create the log file that will keep track of experiments for which What
% calls has been run
WhatLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWhat.txt');
if ~exist(WhatLog, 'file')
    FidWhat = fopen(WhatLog, 'a');
    fprintf(FidWhat, 'Subject\tDate\tTime\n');
    DoneListWhat = [];
else
    FidWhat = fopen(WhatLog, 'r');
    Header = textscan(FidWhat,'%s\t%s\t%s\n',1);
    DoneListWhat = textscan(FidWhat,'%s\t%s\t%s\t');
    fclose(FidWhat);
    FidWhat = fopen(WhatLog, 'a');
end

% Load the file that contains the list of experiments manually curated
WhoLog = fullfile(BaseDataDir, 'RecOnlyLogDeafSalWho.txt');
%how to do this w variable BaseDataDir
if ~exist(WhoLog, 'file')
    error('Cannot find the list of file to run in: %s \n',WhoLog);
else
    FidWho = fopen(WhoLog, 'r');
    Header = textscan(FidWho,'%s\t%s\t%s\t%s\t%s\n',1);
    DoneListDetect = textscan(FidWho,'%s\t%s\t%s\t%.1f\t%d');
    fclose(FidWho);
end

% This is the loop that goes through all manually manually curated exp
NExpe = length(DoneListDetect{1});

for exp=1:NExpe
    % Create the path to the data folder
    BatsID = DoneListDetect{1}{exp};
    Date = DoneListDetect{2}{exp};
    ExpStartTime = DoneListDetect{3}{exp};

    Loggerdir = fullfile(BaseDataDir,['20', Date],'audiologgers');
    % Apply what calls to that experiment
    fprintf(1,'Now running %s\n', Loggerdir)
    what_calls(Loggerdir, Date, ExpStartTime);
    
    % Append to Logfile to indicate this data is processed
    fprintf(FidWhat, '%s\t%s\t%s\n',BatsID,Date,ExpStartTime);
end
fclose(FidWhat);