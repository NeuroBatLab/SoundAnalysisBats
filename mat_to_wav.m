%% 
% Take mat files with raw audio in file.cut and convert to wav file

% Go to directory with analyzed auto 
% and extract numnber of files
cd /home/mtb/Redwood/Data/180907_Analyzed_auto/

mat_files = dir('*.mat');
num_mat_files = length(mat_files); 

%% 
% Extract all cuts to see distribution

% Get 'cut' data for statistics
cut_length = nan(num_mat_files,1);
for ii = 1:num_mat_files
    temp = load(mat_files(ii).name);
    cut_length(ii) = length(temp.cut);
end

%% 
% Plot Histogram, get percentiles and longest vocalization etc.

% fs = 250000;
% histogram(cut_lengths/fs)
% prctl90 = prctile(cut_lengths, 90);
% prctl95 = prctile(cut_lengths, 95);
% prctl99 = prctile(cut_lengths, 99);
% m = max(cut_lengths);
%% 
% Write the cuts to folder of choice as .wav files

% Enter the path to the folder for wav files to be saved
Path2Data='/home/mtb/Redwood/Data/180907_wavFiles';

win = 0.4; %standard window length in s
fs = 250000; %assuming sample frequency is uniform
strd_length = win*fs;

% Start to loop through matfiles
for ii = 1:num_mat_files
    name = mat_files(ii).name;
    temp = load(name);
    cut = temp.cut;
    silence = round((strd_length - length(cut))/2);
    if silence < 0 %if v is longer than standard length
        New_cut = cut(1:strd_length); %start at beginning and cutoff end
    else
        New_cut = zeros(strd_length,1);
        New_cut(silence:(silence+length(cut)-1)) = cut;
    end
    
    %save each file as a .wav file
    loc = strfind(name, '.');
    filename = [fullfile(Path2Data,name(1:loc-1)) '.wav'];
    audiowrite(filename,New_cut, fs);

end



%% 
%