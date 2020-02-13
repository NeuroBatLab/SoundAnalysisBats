
function [DAS]=cal_duration_acoustic_state(BioSoundPath)

% 1. Calculate acoustic features in time i (one FFT frame): pitch, FM, Wiener entropy, Goodness of pitch
% 2. Calculate same features in itime i+1
% 3. scale the features according to the table below (subtract median and divide by MAD)
% 4. calculate Euclidean distance across features, e.i., sqrt (pitch i - pitch i+1)^2 + sqrt(FM....
% 5. if distance < 2 (MADs) repeat for i+2, i+3... each time using the same origin , i, as reference, until distance from origin is =>2
% 6. Repeat for i-1, ... i-2... 
% 7 DAS is the interval around i, where MAD<2

%% List the data
AllBioS = dir(fullfile(BioSoundPath, '*.mat'));
NFiles = length(AllBioS);

%% Calculate median and MAD for all time varying features
AllSal = cell(1,NFiles);
AllSpecMeanSlope = cell(1,NFiles);
AllAmpSlope = cell(1,NFiles);
%AllFo = cell(1,NFiles);
AllAmp =cell(1,NFiles);
AllSpecMean = cell(1,NFiles);

for ff=1:NFiles
    load(fullfile(AllBioS(ff).folder,AllBioS(ff).name),'BioSoundCall')
    AllSal{ff} = BioSoundCall.sal;
    AllSpecMeanSlope{ff} = diff(BioSoundCall.SpectralMean);
    AllAmpSlope{ff} = diff(BioSoundCall.amp);
%     AllFo{ff} = BioSoundCall.Fo;
    AllAmp{ff} = BioSoundCall.amp;
    AllSpecMean{ff} = BioSoundCall.SpectralMean;
end
AllSal = cell2mat(AllSal);
AllSpecMeanSlope = cell2mat(AllSpecMeanSlope);
AllAmpSlope = cell2mat(AllAmpSlope);
% AllFo = cell2mat(AllFo);
AllAmp = cell2mat(AllAmp);
AllSpecMean = cell2mat(AllSpecMean);
MedianSal = median(AllSal,'omitnan');
MADSal = median(abs(AllSal - MedianSal),'omitnan');
MedianSpecMeanSlope = median(AllSpecMeanSlope,'omitnan');
MADSpecMeanSlope = median(abs(AllSpecMeanSlope - MedianSpecMeanSlope),'omitnan');
MedianAmpSlope = median(AllAmpSlope,'omitnan');
MADAmpSlope = median(abs(AllAmpSlope - MedianAmpSlope),'omitnan');
% MedianFo = median(AllFo,'omitnan');
% MADFo = median(abs(AllFo - MedianFo),'omitnan');
MedianAmp = median(AllAmp,'omitnan');
MADAmp = median(abs(AllAmp - MedianAmp),'omitnan');
MedianSpecMean = median(AllSpecMean,'omitnan');
MADSpecMean = median(abs(AllSpecMean - MedianSpecMean),'omitnan');

%% Loop again through files correct input values and calculate DAS
for ff=1:NFiles
    % Initialize the output variable
    DAS = cell(NFiles,1);
    % load the data for that bat/file
    load(fullfile(AllBioS(ff).folder,AllBioS(ff).name),'BioSoundCall')
    % correct values by median and MAD
    Sal = (BioSoundCall.sal(1:end-1)-MedianSal)/MADSal;
    SpecMeanSlope = (diff(BioSoundCall.SpectralMean) - MedianSpecMeanSlope)/MADSpecMeanSlope;
    AmpSlope = (diff(BioSoundCall.amp) - MedianAmpSlope)/MADAmpSlope;
    Amp = (BioSoundCall.amp(1:end-1)-MeadianAmp)/MADAmp;
    SpecMean = (BioSoundCall.SpectralMean(1:end-1) - MedianSpecMean)/MADSpecMean;
    % Now iterate through time point and calculate DAS
    NP = length(Sal);
    DAS{ff} = nan(NP,1);
    for pp=1:NP
        % iterate through future time point and find when euclidian
        % distance >2
        for  ifuture=1:NP-pp
            ED = ((Sal(pp)-Sal(pp+ifuture))^2)^0.5 + ((SpecMean(pp)-SpecMean(pp+ifuture))^2)^0.5 + ((SpecMeanSlope(pp)-SpecMeanSlope(pp+ifuture))^2)^0.5 +((Amp(pp)-Amp(pp+ifuture))^2)^0.5 + ((AmpSlope(pp)-AmpSlope(pp+ifuture))^2)^0.5;
            if ED<2
                break
            end
        end
            
        % iterate through past time point and find when euclidian
        % distance >2
        for  ipast=1:pp-1
            ED = ((Sal(pp)-Sal(pp-ipast))^2)^0.5 + ((SpecMean(pp)-SpecMean(pp-ipast))^2)^0.5 + ((SpecMeanSlope(pp)-SpecMeanSlope(pp-ipast))^2)^0.5 +((Amp(pp)-Amp(pp-ipast))^2)^0.5 + ((AmpSlope(pp)-AmpSlope(pp-ipast))^2)^0.5;
            if ED<2
                break
            end
        end
        DAS{ff}(pp) = ipast + ifuture -2;
    end
end

