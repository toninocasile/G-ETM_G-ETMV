%
%
% This M-file fits a GLM that takes into account the spiking history
% of a neuron as well as global component that takes care of
% non-stationarity in the firing of that neuron
%
% Code tested in Matlab R2019B
%
% author:   Antonino Casile
% toninocasile@gmail.com
%
%
% USAGE :
%
% function [beta_Global, beta_History, dev, beta_Global_pVals, beta_History_pVals] = ...
%     fitGLM_G_ETM(SpikesData_in, globalRegressor, historyRegressor, targetNeuron, nBinsHistory, triggerNeuron)
%
%      beta_History will have zeros line for the beta coefficients of the
%          trigger neuron
%
function [beta_Global, beta_History, dev, beta_Global_pVals, beta_History_pVals] = ...
    fitGLM_G_ETM(SpikesData_in, globalRegressor, historyRegressor, targetNeuron, nBinsHistory, triggerNeuron)

% triggerNeuron is defined
if nargin == 6
    % Dimension of SpikesData_in
    [nNeurons, lenTrial_samples, nTrials] = size(SpikesData_in);
    
    % Remove trigger neuron
    if triggerNeuron == 1
        SpikesData = SpikesData_in(triggerNeuron+1:end,:,:);
    elseif triggerNeuron == nNeurons
        SpikesData = SpikesData_in(1:triggerNeuron-1,:,:);
    else
        SpikesData = [SpikesData_in(1:triggerNeuron-1,:,:); SpikesData_in(triggerNeuron+1:end,:,:)];
    end
else
    SpikesData = SpikesData_in;
end

% get information about our spike data
[nNeurons, lenTrial_samples, nTrials] = size(SpikesData);

% compute total duration of the history IN SAMPLES
totalHistory_samples = historyRegressor.binDuration_samples*nBinsHistory;

% here we compute the final size of the Xneuron and Yneuron matrices
dim2 = globalRegressor.nBins + nNeurons*nBinsHistory;
dim1 = (lenTrial_samples-nBinsHistory*historyRegressor.binDuration_samples)*nTrials;

Xneuron = zeros(dim1, dim2);
Yneuron = zeros(dim1, 1);

% now start filling in the matrices X and Y
currRowInd = 1;
for currTrialInd = 1:nTrials
    for currTimeInd = totalHistory_samples+1:lenTrial_samples
        
        % find the index which is non-zero in the vector of parameters
        if globalRegressor.nBins
            parmsGlobal = zeros(1, globalRegressor.nBins);
            nonZeroInd = floor((currTimeInd-1) / globalRegressor.binDuration_samples) + 1;
            parmsGlobal(nonZeroInd) = 1;
        else
            parmsGlobal = 1;
        end
        
        % extract values for the history of the neuron
        parmsHistory = zeros(nBinsHistory, nNeurons);
        for currInd = 1:nNeurons
            TMP = squeeze(SpikesData(currInd, currTimeInd-totalHistory_samples:currTimeInd-1, currTrialInd));
            TMP = reshape(TMP, [historyRegressor.binDuration_samples, nBinsHistory]);
            parmsHistory(:, currInd) = (historyRegressor.winHistory_samples * TMP)';
        end;
        
        % ... and join all of them together in the final vector
        Xneuron(currRowInd, :) = [parmsGlobal, parmsHistory(:)'];
        
        currRowInd = currRowInd + 1;
    end;
end;

% build output vector
currRowInd = 1;
for currTrialInd = 1:nTrials
    for currTimeInd = totalHistory_samples+1:lenTrial_samples
        % save data into the Y matrix
        Yneuron(currRowInd) = SpikesData_in(targetNeuron, currTimeInd, currTrialInd);
        currRowInd = currRowInd + 1;
    end
end

% perform some sanity check ... just in case
if currRowInd-1 ~= dim1
	error('More elements than I was expecting');
end;

if currRowInd > dim1 + 1
    error('More rows than I can chew')
end;

% if a 1 in a column in Xneuron is always followed by a Y=0 then
% the likelihood diverges ... and we thus have to remove these cases
nBetas = size(Xneuron, 2);
indsToDel = numel(1, nBetas);
currIndToDel = 1;
for currRegrInd = 1:size(Xneuron, 2)
	inds1 = find(squeeze(Xneuron(:, currRegrInd)) == 1);
	if numel(find(Yneuron(inds1) == 1)) == 0
		indsToDel(currIndToDel) = currRegrInd;
		currIndToDel = currIndToDel + 1;
	end
end
if currIndToDel == 1
	indsToDel = [];
else
	indsToDel = indsToDel(1:currIndToDel-1);
end
% indices of the betas that we can keep
indsToKeep = setdiff(1:nBetas, indsToDel);
X = Xneuron(:, indsToKeep);

% run Matlab's glmfit routine
[par_est1 dev stats] = glmfit(X, Yneuron, 'binomial', 'link', 'logit', 'estdisp', 'off', 'constant', 'off');

% now we put the betas together and we set to -20 those that
% are not predictive of the Y
par_est = zeros(nBetas, 1);
par_est_pVals = zeros(nBetas, 1);
par_est(indsToKeep) = par_est1;
par_est_pVals(indsToKeep) = stats.p;
if numel(indsToDel) > 0
	% we set the non-predictive zeros to -20
	par_est(indsToDel) = -20;
	par_est_pVals(indsToDel) = 1;
end

% ... and we split the parameter vector into local and global component
beta_Global = par_est(1:globalRegressor.nBins);
beta_Global_pVals = par_est_pVals(1:globalRegressor.nBins);

% we have to flip the vector of history because time 0 is represented
% on the rightmost element of that vector
% TMP = fliplr(par_est(globalRegressor.nBins+1:end));
% and reshape such that we can see the regressor for each neuron
% for each row time 0 is at index 1!
TMP = reshape(par_est(globalRegressor.nBins+1:end), nBinsHistory, nNeurons)';
TMP1 = reshape(par_est_pVals(globalRegressor.nBins+1:end), nBinsHistory, nNeurons)';

% now reshape the TMP matrix in case we have a target neuron
if nargin == 6
    % add a zero line for the coefficients of the trigger neuron
    if triggerNeuron == 1
        TMP = [zeros(1, nBinsHistory); TMP];
        TMP1 = [ones(1, nBinsHistory); TMP1];
    elseif triggerNeuron == nNeurons
        TMP = [TMP; zeros(1, nBinsHistory)];
        TMP1 = [TMP1; ones(1, nBinsHistory)];
    else
        TMP = [TMP(1:triggerNeuron-1, :); zeros(1, nBinsHistory); ...
            TMP(triggerNeuron:end, :)];
        TMP1 = [TMP1(1:triggerNeuron-1, :); ones(1, nBinsHistory); ...
            TMP1(triggerNeuron:end, :)];
    end
end;
beta_History = fliplr(TMP);
beta_History_pVals = fliplr(TMP1);



