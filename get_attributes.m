%% 
% Get attributes from h5 files for analysis using h5utils.m

import h5utils.*
%%
% Go to directory with h5 files
cd '/home/mtb/Redwood/Data/180907_wavFiles/h5files'

h5files = dir('*.h5');
perm = 'r'; % read only

h5 = h5utils();
%%
% Get attributes from first file
fid = h5.open(h5files(2).name, perm);
fields = h5.get_subgroups(fid, '/');

%%
Feature_matrix = nan(length(h5files),length(fields));
%%
% For debugging and reference

temp_array = num2cell(zeros(1,length(fields)));
features = cell2struct(temp_array,fields,2);
for k = 1:length(fields)
    features.(fields{k}) = h5.get_ds(fid, fields{k});
end
%% 
% Loop through each file, and for each file get means of attributes and store 
% them in a nan array


for ii = 1:length(h5files)
    fid = h5.open(h5files(ii).name, perm);
    for jj = 1:length(fields)
        attr = h5.get_ds(fid, fields{jj});
        if isnumeric(attr) == true
            Feature_matrix(ii,jj) = nanmean(attr, 'all');
        end
    end
end

%% 
% PCA

% standardize matrix
mu = nanmean(Feature_matrix);
sigma = nanstd(Feature_matrix);
Stand_Feats = (Feature_matrix-repmat(mu,length(Feature_matrix),1))./repmat(sigma,length(Feature_matrix),1);

%%
% Complete method with nan fields
% Use vec to ignore irrelevant fields 
vec1 = 1:length(fields);
vec1([3,6,12,15,30,34,37,38,39,41,42]) = []; % Picked out by hand

[coeff1,score1,latent1,tsquared1,explained1] = pca(Stand_Feats(:,vec1));
%%
% Plot explained values
bar(1:length(explained1), explained1)
title("Explained Percentages for each PC")
%%
%Plot data on first three/two PC's

biplot(coeff1(:,1:2),'Scores',score1(:,1:2),'VarLabels',fields(vec1));
title("Biplot of First Two PCs")
%% 
% 
% Complete method with no nan fields
vec2 = 1:length(fields);

temp = [];
for k = 1:length(fields)
    if any(isnan(Stand_Feats(:,k)))
        temp = [temp,k];
    end
end
vec2([3,6,12,15,30,34,37,38,39,41,42,temp]) = [];
%%
[coeff2,score2,latent2,tsquared2,explained2] = pca(Stand_Feats(:,vec2));
%%
bar(1:length(explained2), explained2);
title("Explained Percentages for each PC on All Obs");
%%
biplot(coeff2(:,1:2),'Scores',score2(:,1:2),'VarLabels',fields(vec2));
title("Biplot of First Two PCs with All Observations");
%%
% 
%% 
% Get some stats on each field

% RMS
fn = 'rms';
fn_loc = find(contains(fields, fn));

RMS = Feature_matrix(:,fn_loc);
mean_RMS = mean(RMS);
std_RMS = std(RMS);
histogram(RMS, 10)
title('RMS')

%%
% entropy time
fn = 'entropytime';
fn_loc = find(contains(fields, fn));

enTime = Feature_matrix(:,fn_loc);
mean_et = mean(enTime);
std_et = std(enTime);
histogram(enTime)
title('Entropy Time')

%%
%saliency
fn = 'sal';
fn_loc = find(contains(fields, fn));

saliency = Feature_matrix(:,fn_loc);
mean_sal = mean(saliency);
histogram(saliency, 10)
title('Saliency')
%%
%standard deviation of time
fn = 'stdtime';
fn_loc = find(contains(fields, fn));

std_time = Feature_matrix(:,fn_loc);
histogram(std_time, 10)
title('Standard Deviation of Time')
%%
%Mean spectrum
fn = 'meanspect';
fn_loc = find(contains(fields, fn));

mean_spect = Feature_matrix(:,fn_loc);
histogram(mean_spect,10 )
title('Mean Spectrum')

%%
scatter(saliency,RMS)
xlabel('Saliency')
ylabel('RMS')
%%
scatter(saliency,std_time)
xlabel('Saliency')
ylabel('STD of Time')
%%
scatter(saliency,enTime)
xlabel('Saliency')
ylabel('Entropy Time')
%%
scatter(enTime,std_time)
xlabel('Entropy Time')
ylabel('STD of Time')
%%
scatter(saliency,mean_spect)
xlabel('Saliency')
ylabel('Mean Spectrum')
%%
scatter(mean_spect,std_time)
xlabel('Mean Spectrum')
ylabel('std_time')