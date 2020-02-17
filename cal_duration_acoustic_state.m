
function cal_duration_acoustic_state(BioSoundPath)

% 1. Import from Biosound time varying acoustic features: pitch saliency,
% spectral mean and amplitude
% 2. calculate the median and the MAD of each features accross all times
% points and all files
% 3. scale the time varying features by subtracting median and dividing by MAD
% 4. at each time point i calculate Euclidean distance between i and i+1 across features, e.i., sqrt (Sal i - Sal i+1)^2 + sqrt(SM....
% if distance < 2 (MADs) repeat for i+2, i+3... each time using the same origin , i, as reference, until distance from origin is =>2
% Repeat for i-1, ... i-2... 
% 7. DAS is the interval around i, where euclidian distance in unit of MAD<2

%% List the data
AllBioS = dir(fullfile(BioSoundPath,'*.mat'));
NFiles = length(AllBioS);

%% Calculate median and MAD for all time varying features
AllSal = cell(1,NFiles);
AllSpecMeanSlope = cell(1,NFiles);
AllAmpSlope = cell(1,NFiles);
%AllFo = cell(1,NFiles);
AllAmp =cell(1,NFiles);
AllSpecMean = cell(1,NFiles);

for ff=1:NFiles
    fprintf('Calculating Median MAD: File %d/%d\n',ff, NFiles);
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
    fprintf('Calculating DAS: File %d/%d\n',ff, NFiles);
    % load the data for that bat/file
    load(fullfile(AllBioS(ff).folder,AllBioS(ff).name),'BioSoundCall')
    % correct values by median and MAD
    Sal = (BioSoundCall.sal(1:end-1)-MedianSal)/MADSal;
    Sal(isnan(Sal)) = 0;
    SpecMeanSlope = (diff(BioSoundCall.SpectralMean) - MedianSpecMeanSlope)/MADSpecMeanSlope;
    AmpSlope = (diff(BioSoundCall.amp) - MedianAmpSlope)/MADAmpSlope;
    Amp = (BioSoundCall.amp(1:end-1)-MedianAmp)/MADAmp;
    SpecMean = (BioSoundCall.SpectralMean(1:end-1) - MedianSpecMean)/MADSpecMean;
    % Now iterate through time point and calculate DAS
    NP = length(Sal);
    DAS = nan(NP,1);
    for pp=1:NP
        % iterate through future time point and find when euclidian
        % distance >2
        for  ifuture=1:NP-pp
            ED = ((Sal(pp)-Sal(pp+ifuture))^2)^0.5 + ((SpecMean(pp)-SpecMean(pp+ifuture))^2)^0.5 + ((SpecMeanSlope(pp)-SpecMeanSlope(pp+ifuture))^2)^0.5 +((Amp(pp)-Amp(pp+ifuture))^2)^0.5 + ((AmpSlope(pp)-AmpSlope(pp+ifuture))^2)^0.5;
            if ED>2
                break
            end
        end
            
        % iterate through past time point and find when euclidian
        % distance >2
        for  ipast=1:pp-1
            ED = ((Sal(pp)-Sal(pp-ipast))^2)^0.5 + ((SpecMean(pp)-SpecMean(pp-ipast))^2)^0.5 + ((SpecMeanSlope(pp)-SpecMeanSlope(pp-ipast))^2)^0.5 +((Amp(pp)-Amp(pp-ipast))^2)^0.5 + ((AmpSlope(pp)-AmpSlope(pp-ipast))^2)^0.5;
            if ED>2
                break
            end
        end
        if isempty(ifuture)
            ifuture =0;
        end
        if isempty(ipast)
            ipast = 0;
        end
        DAS(pp) = ipast + ifuture -2;
    end
    save(fullfile(AllBioS(ff).folder,AllBioS(ff).name),'DAS','-append')
end

