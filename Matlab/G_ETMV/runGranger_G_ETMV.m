%
% This m-file runs the G-ETMV method as described in:
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
function GrangerRes = runGranger_G_ETMV(SpikeTrains, globalRegressor_In, historyRegressor, historyRegressorNBins)

[nNeurons, nSamples, nTrials] = size(SpikeTrains);

% get the number of bins
nhistoryRegressorNBins = numel(historyRegressorNBins);
nglobalRegressorNBins = numel(globalRegressor_In.nBins);

% here we save the fitted parameters
beta_Trials = cell(nNeurons, nhistoryRegressorNBins, nglobalRegressorNBins);
beta_Trials_pVals = cell(nNeurons, nhistoryRegressorNBins, nglobalRegressorNBins);
beta_Global = cell(nNeurons, nhistoryRegressorNBins, nglobalRegressorNBins);
beta_Global_pVals = cell(nNeurons, nhistoryRegressorNBins, nglobalRegressorNBins);
beta_History = cell(nNeurons, nhistoryRegressorNBins, nglobalRegressorNBins);
beta_History_pVals = cell(nNeurons, nhistoryRegressorNBins, nglobalRegressorNBins);

% here we go with the first round of the causal test
glm_dev = zeros(nNeurons, nhistoryRegressorNBins, nglobalRegressorNBins);
% log-likelihood of the data
spikes_LL = zeros(nNeurons, nhistoryRegressorNBins, nglobalRegressorNBins);
% Akaike's information criterion
aic = zeros(nNeurons, nhistoryRegressorNBins, nglobalRegressorNBins);

for currGlobalRegressorInd = 1:nglobalRegressorNBins
	currNGlobalBins = globalRegressor_In.nBins(currGlobalRegressorInd);
	
	globalRegressor.nBins = currNGlobalBins;
	% temporal duration of each bin of the global regressor
	% in SAMPLES
	globalRegressor.binDuration_samples = round(nSamples / globalRegressor.nBins);
	
    % here we use Matlab's parallel toolbox
	% if that is not available, substitute the parfor with a simple for
	parfor currHistoryRegressorInd = 1:nhistoryRegressorNBins
		currNRegressorSteps = historyRegressorNBins(currHistoryRegressorInd);
		
		for currNeuronInd = 1:nNeurons
			fprintf('Processing neuron #%d -- history regressor steps = %d -- global regressor steps = %d\n', ...
				currNeuronInd, currNRegressorSteps, currNGlobalBins);
			
			% here we call function for fitting a GLM model to our data
			[tmp_Trials, tmp_Global, tmp_History, dev, tmp_Trials_pVals, tmp_Global_pVals, tmp_History_pVals] = ...
				fitGLM_G_ETMV(SpikeTrains, globalRegressor, ...
				historyRegressor, currNeuronInd, currNRegressorSteps);
			
			% save results of the GLM fitting
			beta_Trials{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd} = tmp_Trials;
			beta_Trials_pVals{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd} = tmp_Trials_pVals;
			beta_Global{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd} = tmp_Global;
			beta_History{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd} = tmp_History;
			beta_Global_pVals{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd} = tmp_Global_pVals;
			beta_History_pVals{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd} = tmp_History_pVals;
			
			% save deviance, likelihood and Akaike criterion for this neuron
			glm_dev(currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd) = dev;
			
			% compute log-likelihood of the spike trains given the fitted
			% parameters. We need that to compute the order of the regressors
			spikes_LL(currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd) = ...
				computeLogLikelihood_SpikeTrains_trialByTrialVar(SpikeTrains, ...
				beta_Trials{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd}, ...
				beta_Global{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd}, ...
				beta_History{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd}, ...
				globalRegressor, historyRegressor, currNeuronInd);
			
			% compute Akaike's information criterion
			aic(currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd) = ...
				-2*spikes_LL(currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd) + ...
				2*(numel(beta_History{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd}) + ...
				numel(beta_Global{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd}) + ...
				numel(beta_Trials{currNeuronInd, currHistoryRegressorInd, currGlobalRegressorInd}) ...
				+ 1);
		end;
	end;
end;

% ----------------------------------------------------------------------
% here we re-fit the model by excluding the effect of the trigger neuron
% ----------------------------------------------------------------------
beta_Trials_causal = cell(nNeurons, nNeurons);
beta_Global_causal = cell(nNeurons, nNeurons);
beta_History_causal = cell(nNeurons, nNeurons);
glm_dev_ratio = zeros(nNeurons, nNeurons);
% compute number of bins we have to use for the history component
indsNBinsHistory_neuron = zeros(1, nNeurons);
indsNBinsGlobal_neuron = zeros(1, nNeurons);
for currNeuronInd = 1:nNeurons
	TMPAIC = squeeze(aic(currNeuronInd, :, :));
	% The following lines of code are for compensating for a peculiar
	% behavior of Matlab.
	% If size(A) = [n m j] with j ~= 1 then size(squeeze(A(1, :, :)) = [m j]
	% HOWEVER
	% If size(A) = [n m j] with j = 1 then size(squeeze(A(1, :, :)) = [1 m]
	% NOTE THAT this "hack" will work only if aic has size 3 so, if in the
	% future I will change the dimensionality of aic I have to "re-hack"
	% this code for it to work!
	if size(TMPAIC, 1) ~= size(aic, 2)
		TMPAIC = TMPAIC';
	end
	
	[TT, ind] = min(TMPAIC(:));
	[indsNBinsHistory_neuron(currNeuronInd), ...
		indsNBinsGlobal_neuron(currNeuronInd)] = ...
		ind2sub(size(TMPAIC), ind);
end

disp(['----- nBinsHistory_neuron = ', ...
	num2str(historyRegressorNBins(indsNBinsHistory_neuron))]);
disp(['----- nBinsGlobal_neuron = ', ...
	num2str(globalRegressor_In.nBins(indsNBinsGlobal_neuron))]);

disp('Refitting model excluding the effect of the triggering neuron');
for targetNeuronInd = 1:nNeurons
	
	nhistoryRegressorNBins = historyRegressorNBins(indsNBinsHistory_neuron(targetNeuronInd));
	% define the characteristics of the global regressor for this neurons
	globalRegressor.nBins = globalRegressor_In.nBins(indsNBinsGlobal_neuron(targetNeuronInd));
	% temporal duration of each bin of the global regressor
	% in SAMPLES
	globalRegressor.binDuration_samples = round(nSamples / globalRegressor.nBins);
	
	% here we use Matlab's parallel toolbox
	% if that is not available, substitute the parfor with a simple for
	parfor triggerNeuronInd = 1:nNeurons
		fprintf('Causal Step - Processing target neuron #%d - -- trigger neuron #%d\n', ...
			targetNeuronInd, triggerNeuronInd);
		fprintf('Target Neuron: nHistoryBins %d -- nGlobalBins\n', ...
			nhistoryRegressorNBins, globalRegressor.nBins);
		
		% here we call function for fitting a GLM model to our data
		[tmp_Trials, tmp_Global, tmp_History, dev, tmp_Trials_pVals, tmp_Global_pVals, tmp_History_pVals] = ...
			fitGLM_G_ETMV(SpikeTrains, ...
			globalRegressor, historyRegressor, targetNeuronInd, nhistoryRegressorNBins, triggerNeuronInd);
		
		% save results of the GLM fitting
		beta_Trials_causal{targetNeuronInd, triggerNeuronInd} = tmp_Trials;
		beta_Global_causal{targetNeuronInd, triggerNeuronInd} = tmp_Global;
		beta_History_causal{targetNeuronInd, triggerNeuronInd} = tmp_History;
		glm_dev_ratio(targetNeuronInd, triggerNeuronInd) = ...
			dev - glm_dev(targetNeuronInd, indsNBinsHistory_neuron(targetNeuronInd), indsNBinsGlobal_neuron(targetNeuronInd));
	end;
end;

% ==== Significance Testing ====
% Causal connectivity matrix, Psi, w/o FDR
D = glm_dev_ratio;
alpha = 0.05;
temp1 = zeros(size(D));
neuronHistoryRegressorNBins = zeros(1, nNeurons);
for targetNeuronInd = 1:nNeurons
	neuronHistoryRegressorNBins(targetNeuronInd) = historyRegressorNBins(indsNBinsHistory_neuron(targetNeuronInd));
	temp1(targetNeuronInd,:) = D(targetNeuronInd,:) > chi2inv(1-alpha, neuronHistoryRegressorNBins(targetNeuronInd));
end
Psi1 = temp1;

% Causal connectivity matrix, Psi, w/ FDR
fdrv = 0.05;
temp2 = FDR(D, fdrv, neuronHistoryRegressorNBins);
Psi2 = temp2;

% here we save the results in a structure that we return as output
GrangerRes.Psi1 = Psi1;
GrangerRes.Psi2 = Psi2;
GrangerRes.alpha = alpha;
GrangerRes.beta_Trials = beta_Trials;
GrangerRes.beta_Trials_pVals = beta_Trials_pVals;
GrangerRes.beta_Global = beta_Global;
GrangerRes.beta_Global_pVals = beta_Global_pVals;
GrangerRes.beta_History = beta_History;
GrangerRes.beta_History_pVals = beta_History_pVals;
GrangerRes.glm_dev = glm_dev;
GrangerRes.glm_dev_ratio = glm_dev_ratio;
GrangerRes.spikes_LL = spikes_LL;
GrangerRes.aic = aic;
GrangerRes.beta_Global_causal = beta_Global_causal;
GrangerRes.beta_History_causal = beta_History_causal;
GrangerRes.historyRegressorNBins = historyRegressorNBins;
GrangerRes.neuronHistoryRegressorNBins = neuronHistoryRegressorNBins;
GrangerRes.indsNBinsHistory_neuron = indsNBinsHistory_neuron;
GrangerRes.globalRegressorNBins = globalRegressor_In.nBins;
GrangerRes.indsNBinsGlobal_neuron = indsNBinsGlobal_neuron;



