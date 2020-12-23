%
% given the spike trains across trials (SpikeTrainIn) this function returns
% a matrix containing the raster plot and one curve containing the smoothed spike
% rate averaged across trials
%
% Code tested in Matlab R2019B
%
% author:   Antonino Casile
% toninocasile@gmail.com
%

function [SpikeTrains, SpikeTrainsSmoothed, SpikeRates] = ...
    computeRasterAndRate(SpikeTrainsIn, sampFreq_Hz, nPixels, FilterType)

if nargin == 3
    % by default we smooth with a triangular kernel with a sigma of 25ms (following Nawrot,
    % Aertsen and Rotter 1999 who say that a triangular kernel has the same smoothing
    % properties of a Gaussian kernel, but it has the advantage of a finite
    % support)
    % sigma of the smoothing kernel (IN S!)
    FilterType.type = 'triangle';
    FilterType.causal = false;
    FilterType.pars = [.1];
end

% select which smoothing we want to use
% sigma of the smoothing kernel (IN S!)
sampTime_s = 1/sampFreq_Hz;
switch FilterType.type
    case 'triangle'
        sig = FilterType.pars(1);
        T2 = [0:sampTime_s:sqrt(6)*sig];
        T1 = [0:-sampTime_s:-sqrt(6)*sig];
        T1 = fliplr(T1);
        T = [T1(1:end-1), T2];
        % ... and here we compute the triangular kernel
        KS = (1/(6*sig^2))*(sqrt(6)*sig - abs(T));
        if FilterType.causal == true
            [~, indMax] = max(KS);
            KS(1:indMax-1) = 0;
            % renormalize KS
            KS = (1/(sum(KS)*sampTime_s)) * KS;
        end
    case 'boxcar'
        sig = FilterType.pars(1);
        T2 = [0:sampTime_s:sig/2];
        T1 = [0:-sampTime_s:-sig/2];
        T1 = fliplr(T1);
        T = [T1(1:end-1), T2];
        KS = ones(size(T));
        if FilterType.causal == true
            [~, indZeroTime] = min(abs(T));
            KS(1:indZeroTime) = 0;
        end
        % normalize KS
        KS = (1/(sum(KS)*sampTime_s)) * KS;
    case 'gaussian'
        sig = FilterType.pars(1);
        T2 = [0:sampTime_s:3*sig];
        T1 = [0:-sampTime_s:-3*sig];
        T1 = fliplr(T1);
        T = [T1(1:end-1), T2];
        % ... and here we compute the triangular kernel
        KS = (1/(sig*sqrt(2*pi)))*exp(-(T.^2/(2*sig^2)));
        if FilterType.causal == true
            [~, indMax] = max(KS);
            KS(1:indMax-1) = 0;
            % renormalize KS
            KS = (1/(sum(KS)*sampTime_s)) * KS;
        end
    otherwise
        error('I have not idea of which filter you are talking about');
end

% get dimension of the spike trains
[nNeurons, nSamples, nTrials] = size(SpikeTrainsIn);

% contains the raster plots of the responses
SpikeTrains = zeros(nNeurons, nTrials, nPixels);
SpikeTrainsSmoothed = zeros(nNeurons, nTrials, nPixels);
SpikeRates = zeros(nNeurons, nTrials, nPixels);

% This is the kernel used for "replicating" the spikes in order to
% make it more visible in the plot
K = [1 1];

fromX = linspace(0, 1, nSamples);
toX = linspace(0, 1, nPixels);
for currNeuronInd = 1:nNeurons
    for currTrialInd = 1:nTrials
        % get the current trial
        fromSpikes = squeeze(SpikeTrainsIn(currNeuronInd, :, currTrialInd));
        
        % initialize vector of results
        toSpikes = zeros(1, nPixels);
        
        inds = find(fromSpikes);
        toSpikes(ceil(inds * (nPixels/nSamples))) = 1;
        SpikeTrains(currNeuronInd, currTrialInd, :) = toSpikes;
        
        % now add some more pixels in places where there is a spike to increase the visibility in the rasterplot
        TMP1 = filtfilt(K, 1, toSpikes);
        TMP1(TMP1>1) = 1;
        SpikeTrainsSmoothed(currNeuronInd, currTrialInd, :) = TMP1;
        
        % compute the spike rates ... and store them
        currRates = conv(fromSpikes, KS, 'same');
        TMP = interp1(fromX, currRates, toX);
        SpikeRates(currNeuronInd, currTrialInd, :) = TMP;
    end
end




