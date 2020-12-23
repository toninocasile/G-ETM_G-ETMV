%
% This M-file plots the results obtained from our G-ETM Method
%
% author:   Antonino Casile
% toninocasile@gmail.com
%

%
% The structure parsForPlotting is meant to "collect" different parameters
% used for plotting the results
%
function plotResults_G_ETM(dataFile, parsForPlotting, figuresDir)

if nargin < 2
    parsForPlotting = struct();
end

% get the offset for determining the zero of the temporal axis
if isfield(parsForPlotting, 'offsetTAxis_s')
    offsetTAxis_s = parsForPlotting.offsetTAxis_s;
else
    offsetTAxis_s = 0;
end

% load the data
D = load(dataFile);
inds = find(dataFile == '/');
if numel(inds) == 0
    startInd = 1;
else
    startInd = inds(end) + 1;
end
fName = dataFile(startInd:end-4);

% get information about our data to plot
[nNeurons, lenTrial_samples, nTrials] = size(D.SpikeTrains);
lenTrial_s = lenTrial_samples /D.sample_Hz;

% compute raster plots and firingRates
Filter.type = 'gaussian';
Filter.causal = 'false';
% variance of the Gaussian (IN MILLISECONDS)
Filter.pars = .05;
[SpikeTrains, SpikeTrainsToPlotSmoothed, SpikeRates] = ...
    computeRasterAndRate(D.SpikeTrains, D.sample_Hz, 1000, Filter);

% compute maximum for rescaling the plots
TMP = squeeze(mean(SpikeRates, 2));
maxRate = ceil(max(TMP(:))/5)*5;
maxMaxRate = ceil(max(SpikeRates(:))/5) * 5;

% p-level used as statistical threshold
pLevel = 0.05;

% compute the recoverd connectivity matrix with a given p-value
fdrv = pLevel;
Psi2 = FDR(D.OutStruct.glm_dev_ratio, fdrv, D.OutStruct.neuronHistoryRegressorNBins);


%% ------------------------------------------------------------------
% -------------- Plot spike trains of the neurons -------------------
% -------------------------------------------------------------------

% part of the spike rate curve that is not plotted given that it is not
% valid due to "border effect" of the filtering process
notPlotted_s = 1.5*Filter.pars;

for currNeuronInd = 1:size(SpikeTrainsToPlotSmoothed, 1)
    figure('Name', sprintf('spike train - neuron %d', currNeuronInd))
    % plot spike train
    subplot(2, 1, 1)
    TMP = squeeze(SpikeTrainsToPlotSmoothed(currNeuronInd, :, :));
    X = linspace(0, lenTrial_s, size(TMP, 2)) + offsetTAxis_s;
    imagesc(X, 1:nTrials, TMP);
    % set(gca, 'XTick', []);
    ylabel('trial #');
    colormap(1-gray);
    box off
    drawnow
    
    % set a couple of "cosmetic" things
    axis off
    title(['neuron ', num2str(currNeuronInd)], 'FontSize', 26);
    
    % plot spike rate
    subplot(2, 1, 2)
    TMP = squeeze(mean(SpikeRates(currNeuronInd, :, :), 2));
    notPlotted_inds = round((notPlotted_s / lenTrial_s) * numel(X));
    currentMaxRate = max(TMP(notPlotted_inds:end-notPlotted_inds));
    
    plot(X(notPlotted_inds:end-notPlotted_inds), ...
        TMP(notPlotted_inds:end-notPlotted_inds), 'LineWidth', 2, 'Color', 'k');
    box off;
    xlabel('time (s)');
    ylabel('rate (spikes/s)');
    % ylim([0 maxRate]);
    ylim([0 1.2*currentMaxRate]);
    xlim([min(X) max(X)]);
    
    % move axes a little bit to improve image quality
    TMP = get(gca, 'Position');
    set(gca, 'FontSize', 24, 'Position', TMP + [0 .1 0 0]);
    
    if exist('figuresDir', 'var') == true
        outName = [figuresDir, '/', fName, '_rasters_', num2str(currNeuronInd), '.eps'];
        fprintf('Printing %s\n', outName);
        print('-depsc', outName);
    end
end

if isfield(parsForPlotting, 'plotRateAllTrials')
    % -------------------------------------------------------------------
    % -------------- Plot spike trains of the neurons -------------------
    %				We plot the firing rates of ALL TRIALS
    % -------------------------------------------------------------------
    colIni = [1 0 0];
    colEnd = [1 0 0];
    lambda = linspace(0, 1, nTrials);
    for currNeuronInd = 1:size(SpikeTrainsToPlotSmoothed, 1)
        figure('Name', sprintf('spike train - neuron %d - rates ALL trials', currNeuronInd))
        
        % plot spike train
        subplot(2, 1, 1)
        TMP = squeeze(SpikeTrainsToPlotSmoothed(currNeuronInd, :, :));
        X = linspace(0, lenTrial_s, size(TMP, 2)) + offsetTAxis_s;
        imagesc(X, 1:nTrials, TMP);
        % set(gca, 'XTick', []);
        ylabel('trial #');
        colormap(1-gray);
        box off
        drawnow
        
        % set a couple of "cosmetic" things
        axis off
        title(['neuron ', num2str(currNeuronInd)], 'FontSize', 26);
        
        % plot spike rates for all trials
        subplot(2, 1, 2)
        hold
        TMP = squeeze(SpikeRates(currNeuronInd, :, :));
        for currTrialInd = 1:nTrials
            currCol = lambda(currTrialInd)*colIni + (1-lambda(currTrialInd))*colEnd;
            plot(X, TMP(currTrialInd, :), 'Color', currCol);
        end
        % plot average firing rate
        plot(X, mean(TMP), 'LineWidth', 2, 'Color', 'k');
        % do some "cosmetic" adjustments
        box off;
        xlabel('time (s)');
        ylabel('rate (spikes/s)');
        ylim([0 maxMaxRate]);
        drawnow
        
        % move axes a little bit to improve image quality
        TMP = get(gca, 'Position');
        set(gca, 'FontSize', 24, 'Position', TMP + [0 .1 0 0]);
        
        if exist('figuresDir', 'var') == true
            outName = [figuresDir, '/', fName, '_rastersAllTrials_', num2str(currNeuronInd), '.eps'];
            fprintf('Printing %s\n', outName);
            print('-depsc', outName);
        end
    end
end

% -----------------------------------------------------------
%		Plot ALL spike trains in a square grid
% -----------------------------------------------------------
figure('Name', 'spike trains - ALL neurons')

if  ~isfield(parsForPlotting, 'plotGridGlobal')
    % compute side of the square grid for plotting
    sideLenX = ceil(sqrt(nNeurons));
    sideLenY = sideLenX;
else
    sideLenX = parsForPlotting.plotGridGlobal(1);
    sideLenY = parsForPlotting.plotGridGlobal(2);
end

% enlarge figure so that panels look better
hF = gcf;
P = hF.Position;
P = P + P .* [-.3, -.3 .3 .3];
hF.Position = P;

currPlotInd = 1;
exitCycles = false;
for xInd = 1:sideLenX
    
    % get betas for that neurons
    for yInd = 1:sideLenY
        currNeuronInd = currPlotInd;
        
        % if we have way more neurons than panels then exit for cycles
        if currNeuronInd > nNeurons
            exitCycles = true;
            break;
        end
        
        % create the subplot
        subplot(sideLenX, sideLenY, currPlotInd);
        
        % plot also spike rates for comparison
        spikeRate = squeeze(mean(SpikeRates(currNeuronInd, :, :), 2));
        xRates_s = linspace(0, lenTrial_s, numel(spikeRate)) + offsetTAxis_s;
        plot(xRates_s, spikeRate, 'Color', 'k', 'LineWidth', 2);
        
        % set limits of the plot
        xlim([min(xRates_s) max(xRates_s)]);
        ylim([0 maxRate]);
        
        % set general font size
        set(gca, 'FontSize', 16);
        
        % set title of the plot
        title(['neuron ', num2str(currNeuronInd)], 'FontSize', 16);
        
        % set ylabel only on the leftmost panels
        if yInd == 1
            hy = ylabel('rate (spikes/s)', 'FontSize', 14);
        end
        
        if xInd == ceil(currNeuronInd / sideLenY)
            xlabel('time (s)');
        end
        
        currPlotInd = currPlotInd + 1;
    end
    if exitCycles == true
        break;
    end
    
    % export figure on file
    if exist('figuresDir', 'var') == true
        outName = [figuresDir, '/', fName, '_global_AllSquare'];
        fprintf('Printing %s\n', outName);
        print('-depsc', outName);
    end
end


%% ----------------------------------------------------------------------
%                  Plot results of Granger analysis
% ----------------------------------------------------------------------
% get information about our data to plot
[nNeurons, lenTrial_samples, nTrials] = size(D.SpikeTrains);
lenTrial_s = lenTrial_samples * (1000/D.sample_Hz);

%% ----------------------------------------------------------------------
%                Plot SIGNIFICANT interaction functions
% -----------------------------------------------------------------------
indsSel = find(Psi2);
[indsX, indsY] = ind2sub(size(Psi2) , indsSel);
if isfield(parsForPlotting, 'fontSizeAxisBetas')
    fontSizeAxis = parsForPlotting.fontSizeAxisBetas;
else
    fontSizeAxis = 18;
end
if isfield(parsForPlotting, 'fontSizeTitleBetas')
    fontSizeTitle = parsForPlotting.fontSizeTitleBetas;
else
    fontSizeTitle = 20;
end

if numel(indsSel) > 0
    
    if  ~isfield(parsForPlotting, 'plotGridBetas')
        % compute side of the square grid for plotting
        sideLenX = ceil(sqrt(numel(indsX)));
        sideLenY = sideLenX;
    else
        sideLenX = parsForPlotting.plotGridBetas(1);
        sideLenY = parsForPlotting.plotGridBetas(2);
    end
    
    figure('Name', 'SIGNIFICANT betas')
    % enlarge figure so that panels look better
    hF = gcf;
    P = hF.Position;
    P = P + P .* [-.3, -.3 .3 .3];
    hF.Position = P;
    
    if isfield(parsForPlotting, 'limBeta')
        limBeta = parsForPlotting.limBeta;
    else
        % first compute the maximum and minimum for the plot
        limBeta = 0;
        for xInd = indsX
            for yInd = indsY
                selInd = D.OutStruct.indsNBinsHistory_neuron(xInd);
                Betas_vals = D.OutStruct.beta_History{xInd, selInd};
                Betas_vals = Betas_vals(indsY, :);
                
                % -20 represents NaN in my Granger causality code
                inds = find(Betas_vals~= -20);
                TMP = abs(Betas_vals(inds));
                if max(TMP) > limBeta
                    limBeta = max(TMP);
                end
            end
        end
        limBeta = ceil(limBeta);
    end
    
    currPlotInd = 1;
    pLevel = 0.05;
    for subPlotInd = 1:numel(indsX)
        xInd = indsX(subPlotInd);
        yInd = indsY(subPlotInd);
        
        % compute current position in the subplot grid
        gridY = rem(currPlotInd-1, sideLenY+1) + 1;
        gridX = floor((currPlotInd-1) / sideLenY) + 1;
        
        % get the optimal history regressor bins
        optimalGlobalRegressorInd = D.OutStruct.indsNBinsGlobal_neuron(xInd);
        selInd = D.OutStruct.indsNBinsHistory_neuron(xInd);
        Betas_vals = D.OutStruct.beta_History{xInd, selInd, optimalGlobalRegressorInd};
        Betas_p = D.OutStruct.beta_History_pVals{xInd, selInd, optimalGlobalRegressorInd};
        X_ms = [0:D.historyRegressorNBins(selInd)-1] * D.historyRegressor.binDuration_ms;
        maxX_ms = (D.historyRegressorNBins(end)-1) * D.historyRegressor.binDuration_ms;
        
        % get betas for that neuron
        subplot(sideLenX, sideLenY, currPlotInd);
        if numel(X_ms) > 0
            plot(X_ms, Betas_vals(yInd, :), 'LineWidth', 3);
            inds = find(squeeze(Betas_p(yInd, :)) < pLevel);
            hold;
            plot(X_ms(inds), Betas_vals(yInd, inds), 'ro', 'MarkerFaceColor', 'r', ...
                'MarkerSize', 8);
            
        end
        ylim([-limBeta limBeta]);
        xlim([0 maxX_ms]);
        if isfield(D, 'Topology')
            plot(D.time_Topology_ms, squeeze(D.Topology(xInd, yInd, :)), 'k', 'LineWidth', 2);
        end
        
        % set some plot attributes
        nTicks = 3;
        deltaTick = round((maxX_ms/nTicks)/10)*10;
        set(gca, 'FontSize', fontSizeAxis, 'XTick', [0:deltaTick:maxX_ms]);
        title(sprintf('%d \\rightarrow %d', yInd, xInd), 'FontSize', fontSizeTitle);
        
        % add xlabel only on the bottom row panels
        if gridX == sideLenX
            xlabel('time (ms)', 'FontSize', fontSizeAxis);
        end
        
        currPlotInd = currPlotInd + 1;
    end
    
    if exist('figuresDir', 'var') == true
        outName = [figuresDir, '/', fName, '_selectedBetas'];
        fprintf('Printing %s\n', outName);
        print('-depsc', outName);
    end
end

% ------------------------------------------------------------------
%					plot global regressors
% ------------------------------------------------------------------

% compute the maximum to rescale curves
maxLambda = -D.sample_Hz;
for currNeuronInd = 1:nNeurons
    TMP = D.OutStruct.beta_Global{currNeuronInd, 1};
    if max(TMP(:)) > maxLambda
        maxLambda = max(TMP(:));
    end
end

% first compute the maximum of the global regressor to scale plots
maxBeta = -inf;
for currNeuronInd = 1:nNeurons
    % get the optimal global and history regressor bins
    optimalGlobalRegressorInd = D.OutStruct.indsNBinsGlobal_neuron(currNeuronInd);
    optimalHistoryRegressorInd = D.OutStruct.indsNBinsHistory_neuron(currNeuronInd);
    yTrial = D.sample_Hz * exp(D.OutStruct.beta_Global{currNeuronInd, optimalHistoryRegressorInd, optimalGlobalRegressorInd});
    maxBeta = max(maxBeta, max(yTrial));
end;
 

if numel(D.globalRegressor.nBins) > 1
    % ... and now plot
    for currNeuronInd = 1:nNeurons
        % squeeze a little bit the figure on the y axis
        hf = figure('Name', sprintf('global regressors - neuron %d', currNeuronInd));
        P = get(hf, 'Position');
        P(4) = P(4) * .8;
        set(hf, 'Position', P);
        
        % get the optimal global and history regressor bins
        optimalGlobalRegressorInd = D.OutStruct.indsNBinsGlobal_neuron(currNeuronInd);
        optimalHistoryRegressorInd = D.OutStruct.indsNBinsHistory_neuron(currNeuronInd);
        binDuration_ms = D.globalRegressor.binDuration_ms(optimalGlobalRegressorInd);
        xTrial = [1:binDuration_ms:lenTrial_s] + binDuration_ms/2 + (offsetTAxis_s * 1000);
        
        xTrial_s = xTrial / 1000;
        yTrial = D.sample_Hz * exp(D.OutStruct.beta_Global{currNeuronInd, optimalHistoryRegressorInd, optimalGlobalRegressorInd});
        plot(xTrial_s, yTrial, 'LineWidth', 3);
        hold;
        inds = find(D.OutStruct.beta_Global_pVals{currNeuronInd, optimalHistoryRegressorInd, optimalGlobalRegressorInd} < pLevel);
        plot(xTrial_s(inds), yTrial(inds), 'ro', 'MarkerFaceColor', 'r', ...
            'MarkerSize', 8);
        
        % plot also spike rates for comparison
        spikeRate = squeeze(mean(SpikeRates(currNeuronInd, :, :), 2));
        xRates_s = linspace(0, lenTrial_s/1000, numel(spikeRate)) + offsetTAxis_s;
        plot(xRates_s, spikeRate, 'Color', 'k', 'LineWidth', 2);
        
        % set limits of the plot
        xlim([min(xRates_s) max(xRates_s)]);
        ylim([0 max(maxRate, ceil(maxBeta/5)*5)]);
        xlabel('time (s)');
        ylabel('rate (spikes/s)');
        % set axis font
        set(gca, 'FontSize', 24);

        % set plot title
        title(sprintf('neuron %d', currNeuronInd), 'FontSize', 28);
   
        % export figure as .eps file
        if exist('figuresDir', 'var') == true
            outName = [figuresDir, '/', fName, '_global_', num2str(currNeuronInd)];
            fprintf('Printing %s\n', outName);
            print('-depsc', outName);
        end
    end
end

% -----------------------------------------------------------
% ... now replot everything as a square grid
% If there is no global regressor we plot nothing!
if numel(D.globalRegressor.nBins) > 1
    figure('Name', 'global regressor - ALL neurons')
    
    if  ~isfield(parsForPlotting, 'plotGridGlobal')
        % compute side of the square grid for plotting
        sideLenX = ceil(sqrt(nNeurons));
        sideLenY = sideLenX;
    else
        sideLenX = parsForPlotting.plotGridGlobal(1);
        sideLenY = parsForPlotting.plotGridGlobal(2);
    end
    
    % first compute the maximum of the global regressor to scale plots
    maxBeta = -inf;
    for currNeuronInd = 1:nNeurons
        % get the optimal global and history regressor bins
        optimalGlobalRegressorInd = D.OutStruct.indsNBinsGlobal_neuron(currNeuronInd);
        optimalHistoryRegressorInd = D.OutStruct.indsNBinsHistory_neuron(currNeuronInd);
        yTrial = D.sample_Hz * exp(D.OutStruct.beta_Global{currNeuronInd, optimalHistoryRegressorInd, optimalGlobalRegressorInd});
        maxBeta = max(maxBeta, max(yTrial));
    end;

    % enlarge figure so that things look better
    hF = gcf;
    P = hF.Position;
    P = P + P .* [-.3, -.3 .3 .3];
    hF.Position = P;
    
    currPlotInd = 1;
    pLevel = 0.05;
    exitCycles = false;
    for xInd = 1:sideLenX
        
        % get betas for that neurons
        for yInd = 1:sideLenY
            currNeuronInd = currPlotInd;
            
            % if we have way more neurons than panels then exit for cycles
            if currNeuronInd > nNeurons
                exitCycles = true;
                break;
            end
            
            % create the subplot
            subplot(sideLenX, sideLenY, currPlotInd);
            
            % get the optimal global and history regressor bins
            optimalGlobalRegressorInd = D.OutStruct.indsNBinsGlobal_neuron(currNeuronInd);
            optimalHistoryRegressorInd = D.OutStruct.indsNBinsHistory_neuron(currNeuronInd);
            binDuration_ms = D.globalRegressor.binDuration_ms(optimalGlobalRegressorInd);
            xTrial = [1:binDuration_ms:lenTrial_s] + binDuration_ms/2 + (offsetTAxis_s * 1000);
            
            xTrial_s = xTrial / 1000;
            yTrial = D.sample_Hz * exp(D.OutStruct.beta_Global{currNeuronInd, optimalHistoryRegressorInd, optimalGlobalRegressorInd});
            maxBeta = max(maxBeta, max(yTrial));
            plot(xTrial_s, yTrial, 'LineWidth', 2);
            hold;
            inds = find(D.OutStruct.beta_Global_pVals{currNeuronInd, optimalHistoryRegressorInd, optimalGlobalRegressorInd} < pLevel);
            plot(xTrial_s(inds), yTrial(inds), 'ro', 'MarkerFaceColor', 'r', ...
                'MarkerSize', 6);
            
            % plot also spike rates for comparison
            spikeRate = squeeze(mean(SpikeRates(currNeuronInd, :, :), 2));
            xRates_s = linspace(0, lenTrial_s/1000, numel(spikeRate)) + offsetTAxis_s;
            plot(xRates_s, spikeRate, 'Color', 'k', 'LineWidth', 2);
            
            % plot limits
            xlim([min(xRates_s) max(xRates_s)]);
            ylim([0 max(maxRate, ceil(maxBeta/5)*5)]);
            
            % place labels and title
            set(gca, 'FontSize', 18);
            title(['neuron ', num2str(currNeuronInd)], 'FontSize', 16);
            if yInd == 1
                ylabel('rate (spikes/s)', 'FontSize', 16);
            end
            
            if xInd == sideLenX
                xlabel('time (s)', 'FontSize', 18);
            end
            
            currPlotInd = currPlotInd + 1;
        end
        if exitCycles == true
            break;
        end
    end
    
    % export figure on file
    if exist('figuresDir', 'var') == true
        outName = [figuresDir, '/', fName, '_global_AllSquare'];
        fprintf('Printing %s\n', outName);
        print('-depsc', outName);
    end
end

%% ------------------------------------------------------------------
% if present we plot the estimates of the spike rates' coefficients
% -------------------------------------------------------------------
if isfield(D.OutStruct, 'beta_Trials')
    figure('Name', 'estimates of spike rates'' coefficients');
    % plot(D.coeffsSpkRate(1, :), 'k', 'LineWidth', 2);
    hold
    
    colIni = [1 0 0];
    colEnd = [0 0 1];
    lambda = linspace(0, 1, nNeurons);
    strLegend = cell(1, nNeurons+1);
    for currNeuronInd = 1:nNeurons
        TMP = D.OutStruct.beta_Trials{currNeuronInd, ...
            D.OutStruct.indsNBinsHistory_neuron(currNeuronInd), ...
            D.OutStruct.indsNBinsGlobal_neuron(currNeuronInd)};
        currCol =  lambda(currNeuronInd)*colIni + (1-lambda(currNeuronInd))*colEnd;
        strLegend{currNeuronInd} = sprintf('neuron %d', currNeuronInd);
        
        plot(exp(TMP), 'Color', currCol, 'LineWidth', 2);
    end;
    % plot ground-truth of the trial-by-trial variability
    plot(D.coeffsSpkRate(1, :), 'k', 'LineWidth', 2);
    strLegend{end} = 'ground-truth';
    
    set(gca, 'FontSize', 24);
    ylim([0 4]);
    xlabel('trial #');
    ylabel('resp. magnitude');
    legend(strLegend);
    
    % In case plot the figure
    if exist('figuresDir', 'var') == true
        outName = [figuresDir, '/', fName, '_estimatedCorrVar'];
        print('-depsc', outName);
    end
end;


%% ------------------------------------------------------------------
%  if present we plot the coefficient of the spike rates on each trial
% -------------------------------------------------------------------
if isfield(D, 'coeffsSpkRate')
    figure('Name', 'spike rates'' coefficients');
    plot(D.coeffsSpkRate', 'Color', 'k', 'LineWidth', 3);
    xlabel('trial #');
    ylabel('resp. magnitude');
    set(gca, 'FontSize', 24);
    ylim([0, 3])
    
    % title('Recovered connectivity');
    if exist('figuresDir', 'var') == true
        outName = [figuresDir, '/', fName, '_groundTruthCorrVar'];
        fprintf('Printing %s\n', outName);
        print('-depsc', outName);
    end 
end


%% ------------------------------------------------------------------
%			Plot the topology as recovered by Granger causality
% ------------------------------------------------------------------

figure('Name', 'Recovered Topology');
imagesc(Psi2, [0 1]);
colormap([0 0 0; 0 1 0]);
set(gca, 'FontSize', 26, 'XTick', [1:nNeurons], 'YTick', [1:nNeurons]);
xlabel('source', 'FontSize', 42);
ylabel('target', 'FontSize', 42);

% title('Recovered connectivity');
if exist('figuresDir', 'var') == true
    outName = [figuresDir, '/', fName, '_recoveredConnectivity'];
    fprintf('Printing %s\n', outName);
    print('-depsc', outName);
end


%% ------------------------------------------------------------------
%					Plot the GROUND-TRUTH topology
% ------------------------------------------------------------------
if isfield(D, 'Topology')
    TMP = squeeze(sum(abs(D.Topology), 3));
    figure('Name', 'Ground-truth Topology')
    imagesc(TMP ~= 0, [0 1]);
    colormap([0 0 0; 0 1 0]);
    set(gca, 'FontSize', 26, 'XTick', [1:size(D.Topology, 1)], ...
        'YTick', [1:size(D.Topology, 1)]);
    xlabel('source', 'FontSize', 42);
    ylabel('target', 'FontSize', 42);
    % title('Real connectivity');
    if exist('figuresDir', 'var') == true
        outName = [figuresDir, '/', fName, '_realConnectivity'];
        print('-depsc', outName);
    end
end

