%% 
% Once you have your .wav files, use this to analyze and create h5 files with 
% BioSound
% 
% Assumes you have BioSound1_Script.py and is in python search path
%% 
% Variables

%Enter path to wav files to be analyzed
wavfile_path = '/home/mtb/Redwood/Data/TestFolder/Test4';

%Input variables as strings to be concatenated.
plotMe = 'plotMe = True'; 
normalize = 'normalize = False'; 
spec_sample_rate = 'spec_sample_rate = 1000';
freq_spacing = 'freq_spacing = 50';
min_freq = 'min_freq = 300';
max_freq = 'max_freq = 50000';
cutoff_freq = 'cutoff_freq = 75';
amp_sample_rate = 'amp_sample_rate = 1000';
f_high = 'f_high = 50000';
maxFund = 'maxFund = 4000';
minFund = 'minFund = 300';
lowFc = 'lowFc = 500';
highFc = 'highFc = 10000';
minSaliency = 'minSaliency = 0.5';
debugFig = 'debugFig = 0';
minFormantFreq = 'minFormantFreq = 500'; 
maxFormantBW = 'maxFormantBW = 500';
method= strcat('method =',"""", 'Stack',"""");
N = 'N = 2014'; % number of files to analyze and save
%% 
% Command Concatenation

pt1 = "python -c 'import BioSound1_Script; BioSound1_Script.BioSound1(";
pt2 = strcat("dir_path=","""",wavfile_path,"""", ...
    ',',plotMe,',',normalize,',',spec_sample_rate,',',freq_spacing, ...
    ',',min_freq,',',max_freq,',',cutoff_freq,',',amp_sample_rate, ...
    ',',maxFund,',',minFund,',',lowFc,',',highFc,',',minSaliency, ...
    ',',debugFig,',',minFormantFreq,',',maxFormantBW,',',method, ',', N);
    
pt3 = ")'";

command = strcat(pt1,pt2,pt3)

%%
[status,cmd_out] = system(command)