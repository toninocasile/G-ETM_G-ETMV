%
% This M-file computes the log likelihood given the regressors for the
% method G-ETM described in:
%
% Robust point-process Granger causality analysis in presence of exogenous
% temporal modulations and trial-by-trial variability in spike trains.
%
% by Casile A., Faghih R. T. & Brown E. N.
%
% Code tested in Matlab R2019B
%
% author:   Antonino Casile
% toninocasile@gmail.com
%

function llk = computeLogLikelihood_SpikeTrains(...
	Spikes, beta_Global, ...
	beta_History, globalRegressor, historyRegressor, ...
	targetNeuron)

% get information about our data
[nNeurons, lenTrial_samples, nTrials] = size(Spikes);

% compute parameters of the history regressor
nBinsHistory = size(beta_History, 2);
totalHistory_samples = historyRegressor.binDuration_samples*nBinsHistory;

% expand global regressor to the size of a trial
TMP = repmat(beta_Global, 1, globalRegressor.binDuration_samples)';
globalLambda = TMP(:)';

% initialize the value of the variable llk
llk = 0;

% here we go with the computations
for currTrialInd = 1:nTrials
	
	% compute the effect of local regressors
	historyLambda = zeros(1, lenTrial_samples);
	
	% compute the effect of other and same neuron on own spiking history
	for currNeuronInd = 1:nNeurons
		% get regressor for this neuron. We have to "replicate" the samples as
		% the bins for the local regressor are bigger than the time bins
		currentRegressor = repmat(beta_History(currNeuronInd, :), historyRegressor.binDuration_samples, 1);
		currentRegressor = currentRegressor(:);
		
		% the bins of the local regressor are bigger than the time bins
		TMP = conv(squeeze(Spikes(currNeuronInd, :, currTrialInd)), currentRegressor);
		% the first elements must be 0 as there can be no influence of past
		% history on the first bin (as there is actually NO past history at
		% that point!)
		historyLambda = historyLambda + [0, TMP(1:lenTrial_samples-1)];
	end;
	
	% compute probability
	totalLambda = globalLambda + historyLambda;
	% P in the case of a binomial distribution
	P = exp(totalLambda) ./ (1+exp(totalLambda));
	
	% we consider only from some point on to make results homogeneous
	% across conditions with different history durations
	startInd = totalHistory_samples + 1;
	TMP = Spikes(targetNeuron, startInd:end, currTrialInd);
	
	% now compute log-likelihood and store it in the llk variable
	llk = llk + sum(TMP .* log(P(startInd:end)) + (1-TMP).*log(1-P(startInd:end)));
end;





