%
% This m-file runs the G-ETMV method as described in:
%
% Robust point-process Granger causality analysis in presence of exogenous
% temporal modulations and trial-by-trial variability in spike trains.
%
% by Casile A., Faghih R. T. & Brown E. N.
%
%
% Code tested in Matlab R2019B
%
% author:   Antonino Casile
% toninocasile@gmail.com
%

% add path of G-ETM to use the plotResults_G_ETM function
path(path, '../G_ETM');

% load spike trains
D = load('../../Results/ExampleDataSet_TrialByTrialVar.mat');
Topology = D.Topology;
coeffsSpkRate = D.coeffsSpkRate;
time_Topology_ms = D.time_Topology_ms;
SpikeTrains = D.SpikeTrains;
sample_Hz = D.sample_Hz;
sampleTime_ms = 1000 / sample_Hz;
clear D

% get information about spike trains
[nNeurons, lenTrial_samples, nTrials] = size(SpikeTrains);
lenTrial_ms = lenTrial_samples * sampleTime_ms;
disp('Done with the spikes trains! Now I fit GLM models');

% ---------------- define global regressor ------------------------
% number of windows used to divide each trial
globalRegressor.nBins = [1, 5, 10, 15, 20, 25, 30, 40];

% temporal duration of each bin of the global regressor
% in BOTH MILLISECONDS and SAMPLES
globalRegressor.binDuration_samples = round(lenTrial_samples ./ globalRegressor.nBins);
globalRegressor.binDuration_ms = globalRegressor.binDuration_samples * sampleTime_ms;

% ---------------- define history regressor -----------------------
% here is the duration of each bin used for the history of the neuron
historyRegressor.binDuration_samples = 3;
historyRegressor.binDuration_ms = historyRegressor.binDuration_samples * sampleTime_ms;

% define maximum number of bins
historyRegressor.maxNBins = 20;
historyRegressor.winHistory_samples = ones(1, historyRegressor.binDuration_samples);
historyRegressor.winHistory_ms = ones(1, historyRegressor.binDuration_ms);

% number of steps for the history regressor that we test
historyRegressorNBins = [2:2:historyRegressor.maxNBins];

% now run G-ETMV on the loaded spike trains
OutStruct = runGranger_G_ETMV(SpikeTrains, globalRegressor, historyRegressor, historyRegressorNBins);
% comment the line above and uncomment the line below to run G-ETM instead
% OutStruct = runGranger_G_ETM(SpikeTrains, globalRegressor, historyRegressor, historyRegressorNBins);

% now save all the results
fName = ['./Out.mat'];
fprintf('Saving %s \n', fName);
save(fName);

% ... and now let's plot the results
plotResults_Granger(fName);



