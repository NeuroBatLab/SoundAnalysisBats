% Converts vocal extraction data in .mat files into relevant .wav files,
% separating the audio by type. Writes into folder Z:\tobias\vocOperant\error_clips

OutputDataPath = 'Z:\tobias\vocOperant\error_clips';
BoxOfInterest = [3 4 6 8];

for bb=1:length(BoxOfInterest)
    DatesDir = dir(fullfile('Z:\tobias\vocOperant\box3\piezo',sprintf('box%d',BoxOfInterest(bb)),'piezo', '1*'));
    for dd=1:length(DatesDir)
        indsrc = dir(fullfile(DatesDir(dd).folder, DatesDir(dd).name,'audiologgers', '*_VocExtractData_*'));
        wavsrc = dir(fullfile(DatesDir(dd).folder, DatesDir(dd).name,'audiologgers', '*_VocExtractData'));
        for ff=1:length(wavsrc)
            if (~isempty(wavsrc)) && (~isempty(indsrc))
                load(fullfile(wavsrc(ff).folder,wavsrc(ff).name), 'Raw_wave','FS');
                load(fullfile(indsrc(ff).folder,indsrc(ff).name),  'IndVocStartRaw', 'IndVocStopRaw', 'IndNoiseStartRaw', 'IndNoiseStopRaw');

                % Filter for the Mic signal
                [z,p,k] = butter(3,100/(FS/2),'high');
                sos_high_raw = zp2sos(z,p,k);

                % Do we even need to save the clip names? Use this if yes:
                % clip_names = cell(sum(cellfun(@length,IndVocStartRaw)));
                % clip_names(vv * (length(IndVocStartRaw(vv)) - 1) + logger) = ...

                Nvoc = length(Raw_wave);
                for vv=1:Nvoc
                  WL = Raw_wave(vv);
                    for logger=1:length(IndVocStartRaw(vv))
                        %[WL, FS] = audioread(fullfile(WaveFiles(vv).folder, WaveFiles(vv).name));
                        % convert with / 1000 * FS ??

                        indices = IndVocStartRaw(vv);
                        vocSnippet = WL(indices(logger):indices(logger));
                        % filter and center the data
                        vocFiltWL = filtfilt(sos_high_raw, 1, vocSnippet);
                        vocFiltWL = vocFiltWL-mean(vocFiltWL);
                        % get name and write to file
                        voc_name = sprintf('%s_%d.wav', wavsrc(ff).name(1:end-4), indices(logger));
                        audiowrite(fullfile(OutputDataPath, voc_name), vocSnippet, FS)

                        %repeat lines 48-54 for IndNoiseStart/Stop
                        indices = IndNoiseStartRaw(vv);
                        noiseSnippet = WL(indices(logger):indices(logger));
                        noiseFiltWL = filtfilt(sos_high_raw, 1, noiseSnippet);
                        noiseFiltWL = noiseFiltWL-mean(noiseFiltWL);
                        noise_name = sprintf('%s_%d.wav', wavsrc(ff).name(1:end-4), indices(logger));
                        audiowrite(fullfile(OutputDataPath, noise_name), noiseSnippet, FS)
                    end
                end
            end
        end
    end
end
