%% 
% This notebook investigate various parameters for UMAP 
% 
% it depends on scriptTestLoggerDetection.m

addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'))
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'))

addpath /Users/elie/Documents/CODE/umapDistribution/umap
addpath /Users/elie/Documents/CODE/umapDistribution/util
javaaddpath('/Users/elie/Documents/CODE/umapDistribution/umap/umap.jar')
%%
Path2Results1 = '/Volumes/Julie4T/JuvenileRecordings151/20190927/audiologgers/GroundTruthResultsPipelineCheck';
Data190927 = load(fullfile(Path2Results1, 'SoundEvent2.mat'));
DataSet190927 = ~isnan(Data190927.BioSoundUniqParam(:,21));
Data190927.BioSoundUniqParam = Data190927.BioSoundUniqParam(DataSet190927,:);
Path2Results1 = '/Volumes/Julie4T/LMC_CoEd/logger/20190603/GroundTruthResultsRecOnly';
Data190603 = load(fullfile(Path2Results1, 'SoundEvent2.mat'));
DataSet190603 = ~isnan(Data190603.BioSoundUniqParam(:,21));
Data190603.BioSoundUniqParam = Data190603.BioSoundUniqParam(DataSet190603,:);

BioSoundUniqParam = [Data190603.BioSoundUniqParam; Data190927.BioSoundUniqParam];
BioSoundParamNames = Data190927.BioSoundParamNames;
%% Change the code of Noise/call in BioSoundParam
Call1Noise0 = BioSoundUniqParam(:,21);
BioSoundUniqParam(:,21) = BioSoundUniqParam(:,21)+1; %1=Noises; 2 =Calls
%%
UsefulParams = [1:20 21];
NIterations = 10;
rng shuffle
pdfRatioScoresThresh = [0:0.01:0.09 0.1:1.5 2:20];
NoiseContamination = nan(NIterations,length(pdfRatioScoresThresh));
LostCalls = nan(NIterations,length(pdfRatioScoresThresh));
ProbaThresh = [0:0.001:0.04 0.05:0.05:1];
PercMissVoc = nan(NIterations,length(ProbaThresh));
PercNoiseContam = nan(NIterations,length(ProbaThresh));

for ii=1:NIterations
    
    fprintf(1,'\n\n ***** Creating training/testing sets %d/%d *****\n',ii,NIterations)
    RandomSet = randperm(size(BioSoundUniqParam,1));
    TrainingSet = RandomSet(1:round(0.8*size(BioSoundUniqParam,1)));
    TestingSet = RandomSet((round(0.8*size(BioSoundUniqParam,1))+1) :end);
    
    %% Reduction in 3D with target_weight =0.6,
    N_component =3;
    Target_weight = 0.55;
    
    % Applying supervised UMAP on the training set
    [ReductionTrain,UMAPTrain]= run_umap(BioSoundUniqParam(TrainingSet,UsefulParams),'parameter_names',BioSoundParamNames(UsefulParams(1:(end-1))),'label_column',length(UsefulParams),'save_template_file',fullfile(Path2Results1,'UMAP_templateNoiseVoc.mat'), 'metric','euclidean','label_file', fullfile(Path2Results1,'NoiseCallLabels.properties.txt'), 'target_weight',Target_weight,'n_component',N_component);
    
    
    %% Applying UMAP template obtained above to testing set
    [ReductionTest,UMAPTest]= run_umap(BioSoundUniqParam(TestingSet,UsefulParams(1:end-1)),'parameter_names',BioSoundParamNames(UsefulParams(1:(end-1))),'template_file',fullfile(Path2Results1,'UMAP_templateNoiseVoc.mat'), 'metric','euclidean','match_supervisors',1,'qf_tree',true,'qf_dissimilarity', true,'target_weight',Target_weight, 'n_component',N_component);
    
    %% Plot the results of the projection with ground truth
%     figure(10)
%     subplot(1,2,1)
%     scatter3(ReductionTrain(:,1), ReductionTrain(:,2),ReductionTrain(:,3),5,[Call1Noise0(TrainingSet) zeros(size(ReductionTrain,1),2)],'filled')
%     xlabel('Dim1')
%     ylabel('Dim2')
%     zlabel('Dim3')
%     title('TRAINING set UMAP projection')
%     hold off
%     % Plot the projection and the ground truth
%     subplot(1,2,2)
%     scatter3(ReductionTest(:,1), ReductionTest(:,2),ReductionTest(:,3),5,[Call1Noise0(TestingSet) zeros(size(ReductionTest,1),2)],'filled')
%     xlabel('Dim1')
%     ylabel('Dim2')
%     zlabel('Dim3')
%     title('TESTING set UMAP projection')
    
    %% Modeling the projection of the training set with Gaussians and testing different pdf ratio score for the classification noise/calls in the testing set
    
    
    %Gaussian Model
    CallsIdx = logical(Call1Noise0(TrainingSet));
    NGaussCalls = floor(sum(CallsIdx)/90)-1; % the best estimate of the maximum of gaussian we can reasonnably fit given that in 3D we need 9 parameters to define a gaussian and estimate that 10 data points per parameter is a good number
    NGaussNoise = min(floor(sum(~CallsIdx)/90)-1, 9); % the best estimate of the maximum of gaussian we can reasonnably fit given that in 3D we need 9 parameters to define a gaussian and estimate that 10 data points per parameter is a good number
    MGMC = 0;
    while MGMC==0
        try
            MixtureGaussianModelCalls = fitgmdist(ReductionTrain(CallsIdx,:),NGaussCalls);
            MGMC=1;
        end
    end
    MGMN = 0;
    while MGMN==0
        try
            MixtureGaussianModelNoise = fitgmdist(ReductionTrain(~CallsIdx,:),NGaussNoise);
            MGMN=1;
        end
    end
    
    pdfCalls = pdf(MixtureGaussianModelCalls, ReductionTest);
    pdfNoise = pdf(MixtureGaussianModelNoise, ReductionTest);
    figure(); scatter(log(pdfNoise), log(pdfCalls),5,[Call1Noise0(TestingSet) zeros(size(ReductionTest,1),2)],'filled')
    hold on
    plot(log([-1.5e-10 2.28]), log([-1.5e-10 2.28]), 'b-', 'LineWidth',2)
    ylabel('Log(pdfCalls)')
    xlabel('Log(pdfNoise)')
    hold off
    
    pdfRatioScores = pdfCalls./pdfNoise;
    
    for pp=1:length(pdfRatioScoresThresh)
        LostCalls(ii,pp) = sum((pdfRatioScores<pdfRatioScoresThresh(pp)) .* Call1Noise0(TestingSet))./sum(Call1Noise0(TestingSet))*100;
        NoiseContamination(ii,pp) = sum((pdfRatioScores>pdfRatioScoresThresh(pp)) .* ~Call1Noise0(TestingSet))./sum(pdfRatioScores>pdfRatioScoresThresh(pp))*100;
    end
    
    
    
    
    
    %% Comparing with SVM classifier
    
    % Let's try to predict data using a Binary SVM
    X_train = BioSoundUniqParam(TrainingSet,UsefulParams(1:end-1));
    Y_train = Call1Noise0(TrainingSet);
    X_test = BioSoundUniqParam(TestingSet,UsefulParams(1:end-1));
    Y_test = Call1Noise0(TestingSet);
    %Train an SVM classifier. Standardize the data . Conserve memory by reducing the size of the trained SVM classifier.
    SVMModelSplit = fitcsvm(X_train,Y_train,'Standardize',true,'KernelFunction','RBF',...
        'KernelScale','auto','Prior','Uniform');
    CompactSVMModelSplit = compact(SVMModelSplit);
    whos('SVMModel','CompactSVMModel')
    
    % The CompactClassificationSVM classifier (CompactSVMModel) uses less space than the ClassificationSVM classifier (SVMModel) because SVMModel stores the data.
    % Estimate the optimal score-to-posterior-probability transformation function.
    CompactSVMModelSplit = fitPosterior(CompactSVMModelSplit,...
        X_train,Y_train);
    
    %The optimal score transformation function (CompactSVMModel.ScoreTransform)
    %is a sigmoid  function because the classes are inseparable.
    % Predict the out-of-sample labels and class posterior probabilities. Because true labels are available, compare them with the predicted labels.
    [labels,PostProbs] = predict(CompactSVMModelSplit,X_test);
    
%     figure();
%     scatter(PostProbs(:,1), PostProbs(:,2), 40, [Y_test zeros(size(Y_test)) zeros(size(Y_test))], 'filled')
%     hold on
%     scatter(PostProbs(:,1), PostProbs(:,2), 42, [labels zeros(size(Y_test)) zeros(size(Y_test))])
%     hold off
%     xlabel(sprintf('Posterior probability class %d', CompactSVMModelSplit.ClassNames(1)))
%     ylabel(sprintf('Posterior probability class %d', CompactSVMModelSplit.ClassNames(2)))
%     
%     fprintf(1,'Percentage of misses (vocalizations detected as noise): %.1f or %d/%d\n', sum(~labels.*Y_test)/sum(Y_test)*100, sum(~labels.*Y_test), sum(Y_test)) % 2.5 %
%     fprintf(1,'Percentage of noise contamination (noise detected as vocalizations): %.1f or %d/%d\n', sum(labels.*~Y_test)/sum(labels)*100, sum(labels.*~Y_test), sum(labels)) % 0.5%
%     
    % Now choose a less restrictive label attribution -> any sound with a
    % probability of being a vocalization (class 1) above 0.1 is labeled
    
    for pp=1:length(ProbaThresh)
        NewLabels = PostProbs(:,2)>=ProbaThresh(pp);
        PercMissVoc(ii,pp) = sum(~NewLabels.*Y_test)/sum(Y_test)*100;
        PercNoiseContam(ii,pp) = sum(NewLabels.*~Y_test)/sum(NewLabels)*100;
    end
%     pause()
    close all
    
end
%% Plot results
figure()
shadedErrorBar(pdfRatioScoresThresh, mean(LostCalls,1),std(LostCalls,[],1), {'r-', 'LineWidth',2})
hold on
shadedErrorBar(pdfRatioScoresThresh, mean(NoiseContamination,1), std(NoiseContamination,[],1), {'k-', 'LineWidth',2})
hold off

ylabel('Percentage of error (r:lost voc; k:noise contam)')
xlabel('pCalls/pNoise thresholds')
title('UMAP + Mixture of Gaussian CV')
hold off
ylim([0 100])

figure()
shadedErrorBar(ProbaThresh, mean(PercMissVoc,1),std(PercMissVoc,[],1), {'r','LineWidth',2})
hold on
shadedErrorBar(ProbaThresh, mean(PercNoiseContam,1), std(PercNoiseContam,[],1),{ 'k','LineWidth',2})
hold off
xlabel('Threshold on Posterior probability of class 1 (vocalization)')
ylabel('Percentage of error (r:lost voc; k:noise contam)')
%     legend('Missed Vocalizations', 'Noise contamination')
title('SVM classifier CV')


%% Run UMAP on all the dataset and calculate Mixtures of Gaussian models to predict classification
[Reduction,UMAP]= run_umap(BioSoundUniqParam(:,UsefulParams),'parameter_names',BioSoundParamNames(UsefulParams(1:(end-1))),'label_column',length(UsefulParams),'save_template_file',fullfile(Path2Results1,'UMAP_templateNoiseVoc.mat'), 'metric','euclidean','n_component',3,'target_weight',0.55);
figure()
scatter3(Reduction(:,1), Reduction(:,2),Reduction(:,3),5,[Call1Noise0 zeros(size(Reduction,1),2)],'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
% IndCallCluster = (Reduction(:,1)>-4) .* (Reduction(:,1)<4) .* (Reduction(:,2)<4) .* (Reduction(:,2)>-2)  .* (Reduction(:,3)<20) .* (Reduction(:,3)>10);
IndCallCluster =  (Reduction(:,3)<-10);
fprintf('Percent of noise detected as vocalizations: %d/%d or %.2f%%\n', sum(IndCallCluster .* ~Call1Noise0), sum(IndCallCluster),sum(IndCallCluster .* ~Call1Noise0)*100/sum(IndCallCluster))
fprintf('Percent of lost calls in the noise: %d/%d or %.2f%%\n', sum(~IndCallCluster .* Call1Noise0), sum(Call1Noise0), sum(~IndCallCluster .* Call1Noise0)*100/sum(Call1Noise0))
LostCalls = find(~IndCallCluster .* Call1Noise0);

%Gaussian Model
CallsIdx = logical(Call1Noise0);
MixtureGaussianModelCalls = fitgmdist(Reduction(CallsIdx,:),5);
MixtureGaussianModelNoise = fitgmdist(Reduction(~CallsIdx,:),10);
figure()
subplot(1,2,1)
scatter3(Reduction(CallsIdx,1), Reduction(CallsIdx,2),Reduction(CallsIdx,3),5,'r','filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')
subplot(1,2,2)
scatter3(Reduction(~CallsIdx,1), Reduction(~CallsIdx,2),Reduction(~CallsIdx,3),5,'k','filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')

pdfCalls = pdf(MixtureGaussianModelCalls, Reduction);
pdfNoise = pdf(MixtureGaussianModelNoise, Reduction);
figure(); scatter(log(pdfNoise), log(pdfCalls),5,[Call1Noise0 zeros(size(Reduction,1),2)],'filled')
hold on
plot(log([-1.5e-10 2.28]), log([-1.5e-10 2.28]), 'b-', 'LineWidth',2)

% Get a feel of what the Gaussian looks like
Npoints = 1000;
XRange = [-5 10];
YRange = [-5 10];
ZRange = [-20 10];
np = 0;
DensityGaussian = nan(Npoints,4);
while np<Npoints
    PointReduction = [rand(1)*diff(XRange) + XRange(1) rand(1)*diff(YRange) + YRange(1) rand(1)*diff(ZRange) + ZRange(1)];
    pdfPoint = pdf(MixtureGaussianModelCalls, PointReduction);
    if pdfPoint>rand(1)*max(pdfCalls)
        np=np+1;
        DensityGaussian(np, :) = [PointReduction pdfPoint];
    end
end


% Get a feel of what the Gaussian looks like
Npoints = 10000;
XRange = [-5 10];
YRange = [-10 15];
ZRange = [-20 15];
np = 0;
DensityGaussianNoise = nan(Npoints,4);
while np<Npoints
    PointReduction = [rand(1)*diff(XRange) + XRange(1) rand(1)*diff(YRange) + YRange(1) rand(1)*diff(ZRange) + ZRange(1)];
    pdfPoint = pdf(MixtureGaussianModelNoise, PointReduction);
    if pdfPoint>rand(1)*max(pdfNoise)
        np=np+1;
        DensityGaussianNoise(np, :) = [PointReduction pdfPoint];
    end
end
        
figure()
subplot(1,2,1)
scatter3(Reduction(CallsIdx,1), Reduction(CallsIdx,2),Reduction(CallsIdx,3),20,'r','filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')  
subplot(1,2,2)
colormap hot
scatter3(DensityGaussian(:,1), DensityGaussian(:,2),DensityGaussian(:,3),20,DensityGaussian(:,4),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3') 
    
figure()
subplot(1,2,1)
scatter3(Reduction(~CallsIdx,1), Reduction(~CallsIdx,2),Reduction(~CallsIdx,3),5,'k','filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')  
subplot(1,2,2)
colormap hot
scatter3(DensityGaussianNoise(:,1), DensityGaussianNoise(:,2),DensityGaussianNoise(:,3),5,DensityGaussianNoise(:,4),'filled')
xlabel('Dim1')
ylabel('Dim2')
zlabel('Dim3')


%% Calculate UMAP on dynamic time warping distance of spectrograms
WorkingDirSpectro = '/Users/elie/Documents/GroundTruthWorkDir';
AllEventFiles = dir(WorkingDirSpectro);


%%

[Reduction,UMAP,ClustID]= run_umap(BioSoundUniqParam(TrainingSet,UsefulParams),'parameter_names',BioSoundParamNames(UsefulParams(1:(end-1))),'label_column',length(UsefulParams),'save_template_file',fullfile(Path2Results1,'UMAP_templateNoiseVoc.mat'), 'metric','euclidean','target_weight',0.6,'n_component',3);
figure()
scatter3(Reduction(:,1), Reduction(:,2),Reduction(:,3),5,[BioSoundUniqParam(TrainingSet,21) zeros(size(Reduction,1),2)],'filled')
title('Euclidean distance')