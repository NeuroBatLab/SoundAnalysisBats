function extract_error_audio(Loggers_dir, Date, ExpStartTime)
% Converts vocal extraction data in .mat files into relevant .wav files,
% separating the audio by type. Writes into folder Z:\tobias\vocOperant\error_clips

OutputDataPath = 'Z:\tobias\vocOperant\error_clips';
DataFile = fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime));
load(DataFile, 'Piezo_wave', 'Piezo_FS',  'Raw_wave','FS', 'RatioRMS', 'DiffRMS','BandPassFilter', 'AudioLogs', 'RMSHigh', 'RMSLow','VocFilename');

Nvoc = length(Raw_wave);
IndVocStart_all = cell(1,Nvoc);
IndVocStop_all = cell(1,Nvoc);
IndNoiseStart_all = cell(1,Nvoc);
IndNoiseStop_all = cell(1,Nvoc);
IndVocStartRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc)
% a cell array of the size the number of loggers + 1 in case only one bat without a logger
% or +2 incase no identification possible but you want to keep onset/offset of each voc and
% for each logger the index onset of when the animal start vocalizing in the raw recording before merge
IndVocStartPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording before merge
IndVocStopRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording before merge
IndVocStopPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording before merge
IndNoiseStartRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the raw recording before merge
IndNoiseStartPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording before merge
IndNoiseStopRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording before merge
IndNoiseStopPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording before merge
IndVocStartRaw_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the raw recording
IndVocStopRaw_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording
IndVocStartPiezo_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording
IndVocStopPiezo_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording

% convert with / 1000 * FS
clips = cell(Nvoc);
clip_names = cell(Nvoc);
for vv=vv:Nvoc
  clips(vv) = Raw_wave(vv, IndVocStart_all / 1000 * FS, IndVocStop_all / 1000 * FS);
  clip_names(vv) = VocFilename(vv);

  audiowrite(append(OutputDataPath, clip_names(vv)), clips(vv), FS)
end
