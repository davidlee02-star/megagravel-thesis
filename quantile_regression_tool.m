%% Quantile Regression Tool for Megagravel Transport Analysis
%
% Fits quantile regression models to the 2013-2014 Irish boulder dataset
% at five quantile levels (tau = 0.10, 0.25, 0.50, 0.75, 0.90), with
% optional bootstrap confidence intervals and holdout validation.
%
% Produces coefficient-path plots (natural and standardised units),
% effect plots for each predictor, and observed-vs-predicted plots for
% the median model.
%
% Requires: qreg_lp.m (quantile regression solver, linear programming)
%
% Usage:
%   1. Set responseName to either the horizontal or vertical displacement
%      variable.
%   2. Adjust predictor list, quantile levels, and run switches as needed.
%   3. Run the script.

clear; clc; close all;

%% User inputs
filename = '2014_Boulder_Movements_dataset_of_record.xlsx';

% Set to either 'MovementVertical_m_' or
% 'TotalHorizontalMovement_m__absoluteValue_'
responseName = 'MovementVertical_m_';

predictorNames = { ...
    'DistanceInlandAHW_m_', ...
    'Vol_m3_', ...
    'HeightAHW_m_', ...
    'Mass_t_', ...
    'Steepness' ...
};

taus = [0.10 0.25 0.50 0.75 0.90];

% Run switches
runMainModel   = true;
runBootstrap   = true;
runValidation  = true;
runCoefPlot    = true;
runEffectPlot  = true;
runObsPredPlot = true;

saveResults = false;

% Validation settings
testFraction   = 0.20;
seedValidation = 0;

% Bootstrap settings
nBoot         = 500;
seedBootstrap = 1;

% Output folder
outputFolder = 'quantile_regression_results';

%% Axis labels and titles
axisLabel_response = 'Total vertical movement (m)';

axisLabel_predictors = { ...
    'Distance inland above high water (m)', ...
    'Volume (m^3)', ...
    'Elevation above high water (m)', ...
    'Mass (t)', ...
    'Coastal steepness (-)' ...
};

axisLabel_tau           = 'Quantile (\tau)';
axisLabel_coeffPath     = 'Coefficient';
axisLabel_coeffPath_std = 'Standardised coefficient';
axisLabel_obs           = 'Observed total vertical movement (m)';
axisLabel_pred          = 'Predicted total vertical movement (m)';

title_coeffPath     = 'Quantile regression coefficient paths';
title_coeffPath_std = 'Standardised quantile regression coefficient paths';
title_effect        = 'Quantile regression response curves';
title_obsPred       = 'Observed vs predicted (median model)';

%% Setup
runTag = datestr(now, 'yyyymmdd_HHMMSS');
safeResponse = matlab.lang.makeValidName(responseName);

if saveResults && ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Load and validate data
if ~isfile(filename)
    error('File not found: %s', filename);
end

data = readtable(filename);

if ~ismember(responseName, data.Properties.VariableNames)
    error('Response variable "%s" not found.', responseName);
end

for i = 1:numel(predictorNames)
    if ~ismember(predictorNames{i}, data.Properties.VariableNames)
        error('Predictor "%s" not found.', predictorNames{i});
    end
end

if any(taus <= 0 | taus >= 1)
    error('All quantiles must be strictly between 0 and 1.');
end

%% Build analysis table
varsNeeded = [{responseName}, predictorNames];
T = rmmissing(data(:, varsNeeded));

y = T.(responseName);
X = table2array(T(:, predictorNames));

n = size(X, 1);
p = size(X, 2);

fprintf('Cleaned dataset: %d rows, %d predictors\n', n, p);
fprintf('Response: %s\n', responseName);
fprintf('Quantiles: %s\n', num2str(taus));

%% Check for collinearity
rankX = rank(X);
if rankX < p
    fprintf('Warning: predictor matrix is rank deficient (rank = %d, p = %d).\n', rankX, p);
    C = corrcoef(X);
    for ii = 1:p
        for jj = ii+1:p
            if abs(C(ii, jj)) > 0.9999
                fprintf('  Near-perfect collinearity: %s and %s (r = %.4f)\n', ...
                    predictorNames{ii}, predictorNames{jj}, C(ii, jj));
            end
        end
    end
end

%% Preallocate
Bfull = [];
Bfull_std = [];
Btrain = [];
yPred = [];
yTest = [];
rowNames = ["Intercept"; string(predictorNames(:))];
tauNames = matlab.lang.makeValidName(compose("tau_%0.2f", taus));

%% Fit full-data model (natural units)
if runMainModel || runBootstrap || runCoefPlot || runEffectPlot
    fprintf('\nFitting full-data quantile regression (natural units)...\n');

    Bfull = zeros(p+1, numel(taus));
    for i = 1:numel(taus)
        Bfull(:, i) = qreg_lp(X, y, taus(i));
    end

    coefTable = array2table(Bfull, ...
        'RowNames', cellstr(rowNames), ...
        'VariableNames', tauNames);

    disp('Full-data coefficients (natural units):');
    disp(coefTable);

    if saveResults
        writetable(coefTable, ...
            fullfile(outputFolder, ['coef_full_' safeResponse '_' runTag '.csv']), ...
            'WriteRowNames', true);
    end
end

%% Fit standardised model (for coefficient-path plot)
if runCoefPlot
    Xmean = mean(X, 1);
    Xstd = std(X, 0, 1);
    Xstd(Xstd == 0) = 1;
    Xz = (X - Xmean) ./ Xstd;

    Bfull_std = zeros(p+1, numel(taus));
    for i = 1:numel(taus)
        Bfull_std(:, i) = qreg_lp(Xz, y, taus(i));
    end

    coefTable_std = array2table(Bfull_std, ...
        'RowNames', cellstr(rowNames), ...
        'VariableNames', tauNames);

    if saveResults
        writetable(coefTable_std, ...
            fullfile(outputFolder, ['coef_full_std_' safeResponse '_' runTag '.csv']), ...
            'WriteRowNames', true);
    end
end

%% Bootstrap confidence intervals
if runBootstrap
    if isempty(Bfull)
        error('Bootstrap requires the main model. Set runMainModel = true.');
    end

    fprintf('\nRunning bootstrap (%d resamples)...\n', nBoot);

    rng(seedBootstrap, 'twister');
    Bboot = nan(p+1, numel(taus), nBoot);
    failCount = zeros(numel(taus), 1);

    for b = 1:nBoot
        idx = randsample(n, n, true);
        Xb = X(idx, :);
        yb = y(idx);

        for i = 1:numel(taus)
            try
                Bboot(:, i, b) = qreg_lp(Xb, yb, taus(i));
            catch
                failCount(i) = failCount(i) + 1;
            end
        end

        if mod(b, 50) == 0 || b == nBoot
            fprintf('  %d / %d\n', b, nBoot);
        end
    end

    for i = 1:numel(taus)
        if failCount(i) > 0
            fprintf('  tau = %.2f: %d solver failures of %d\n', ...
                taus(i), failCount(i), nBoot);
        end
    end

    B_lo  = prctile(Bboot, 2.5, 3);
    B_med = prctile(Bboot, 50, 3);
    B_hi  = prctile(Bboot, 97.5, 3);

    % Build long-format CI table
    nRows = (p+1) * numel(taus);
    coefNameCol = strings(nRows, 1);
    tauCol = zeros(nRows, 1);
    estCol = zeros(nRows, 1);
    medCol = zeros(nRows, 1);
    loCol = zeros(nRows, 1);
    hiCol = zeros(nRows, 1);

    k = 1;
    for i = 1:numel(taus)
        for j = 1:(p+1)
            coefNameCol(k) = rowNames(j);
            tauCol(k) = taus(i);
            estCol(k) = Bfull(j, i);
            medCol(k) = B_med(j, i);
            loCol(k) = B_lo(j, i);
            hiCol(k) = B_hi(j, i);
            k = k + 1;
        end
    end

    ciTable = table(coefNameCol, tauCol, estCol, medCol, loCol, hiCol, ...
        'VariableNames', {'Coefficient', 'tau', 'Estimate', ...
                          'BootstrapMedian', 'CI95_Lower', 'CI95_Upper'});

    disp('Bootstrap 95% confidence intervals:');
    disp(ciTable);

    if saveResults
        writetable(ciTable, ...
            fullfile(outputFolder, ['bootstrapCI_' safeResponse '_' runTag '.csv']));

        failTable = table(taus(:), failCount, ...
            'VariableNames', {'tau', 'FailedResamples'});
        writetable(failTable, ...
            fullfile(outputFolder, ['bootstrapFailures_' safeResponse '_' runTag '.csv']));
    end
end

%% Holdout validation
if runValidation
    fprintf('\nRunning holdout validation (%.0f%% test)...\n', 100*testFraction);

    rng(seedValidation, 'twister');
    idx = randperm(n);
    nTest = round(testFraction * n);

    testIdx = idx(1:nTest);
    trainIdx = idx(nTest+1:end);

    XTrain = X(trainIdx, :);
    yTrain = y(trainIdx);
    XTest = X(testIdx, :);
    yTest = y(testIdx);

    fprintf('  Training: %d, test: %d\n', size(XTrain, 1), size(XTest, 1));

    Btrain = zeros(p+1, numel(taus));
    for i = 1:numel(taus)
        Btrain(:, i) = qreg_lp(XTrain, yTrain, taus(i));
    end

    XTestI = [ones(size(XTest, 1), 1), XTest];
    yPred = XTestI * Btrain;

    pinballLoss = zeros(numel(taus), 1);
    MAE = zeros(numel(taus), 1);

    for i = 1:numel(taus)
        tau = taus(i);
        e = yTest - yPred(:, i);
        pinballLoss(i) = mean(max(tau*e, (tau-1)*e));
        MAE(i) = mean(abs(e));
    end

    metricsTable = table(taus(:), pinballLoss, MAE, ...
        'VariableNames', {'tau', 'PinballLoss', 'MAE'});

    disp('Holdout validation metrics:');
    disp(metricsTable);

    if saveResults
        writetable(metricsTable, ...
            fullfile(outputFolder, ['validationMetrics_' safeResponse '_' runTag '.csv']));
    end
end

%% Coefficient path plots
if runCoefPlot
    % Natural units
    figure;
    plot(taus, Bfull(2:end, :).', '-o', 'LineWidth', 1.5);
    xlabel(axisLabel_tau);
    ylabel(axisLabel_coeffPath);
    legend(axisLabel_predictors, 'Location', 'best');
    title(title_coeffPath);
    grid on;
    set(gcf, 'Color', 'w');
    box on;
    set(gca, 'LineWidth', 1);

    % Standardised
    figure;
    plot(taus, Bfull_std(2:end, :).', '-o', 'LineWidth', 1.5);
    xlabel(axisLabel_tau);
    ylabel(axisLabel_coeffPath_std);
    legend(axisLabel_predictors, 'Location', 'best');
    title(title_coeffPath_std);
    grid on;
    set(gcf, 'Color', 'w');
    box on;
    set(gca, 'LineWidth', 1);
end

%% Effect plots (one per predictor)
if runEffectPlot
    gridN = 150;

    for predictorToPlot = 1:numel(predictorNames)
        gridX = repmat(median(X, 1), gridN, 1);
        gridX(:, predictorToPlot) = linspace( ...
            min(X(:, predictorToPlot)), ...
            max(X(:, predictorToPlot)), gridN)';

        gridXI = [ones(gridN, 1), gridX];
        gridY = gridXI * Bfull;

        figure;
        plot(X(:, predictorToPlot), y, '.', 'MarkerSize', 10);
        hold on;
        for i = 1:numel(taus)
            plot(gridX(:, predictorToPlot), gridY(:, i), 'LineWidth', 2);
        end
        hold off;

        xlabel(axisLabel_predictors{predictorToPlot});
        ylabel(axisLabel_response);
        legendLabels = ["Observed", compose('\\tau = %.2f', taus)];
        legend(legendLabels, 'Location', 'best');
        title([title_effect ' vs ' axisLabel_predictors{predictorToPlot}]);
        grid on;
        set(gcf, 'Color', 'w');
        box on;
        set(gca, 'LineWidth', 1);
    end
end

%% Observed vs predicted plot (median model)
if runObsPredPlot
    medianIdx = find(abs(taus - 0.50) < 1e-9, 1);

    figure;
    plot(yTest, yPred(:, medianIdx), '.', 'MarkerSize', 12);
    hold on;
    upper = prctile([yTest; yPred(:, medianIdx)], 95);
    plot([0 upper], [0 upper], '--', 'LineWidth', 1.5);
    xlim([0 upper]);
    ylim([0 upper]);
    hold off;

    xlabel(axisLabel_obs);
    ylabel(axisLabel_pred);
    title(title_obsPred);
    grid on;
    set(gcf, 'Color', 'w');
    box on;
    set(gca, 'LineWidth', 1);
end

%% Save run config for reproducibility
if saveResults
    configTable = table( ...
        string(runTag), ...
        string(responseName), ...
        string(strjoin(predictorNames, '; ')), ...
        string(num2str(taus)), ...
        n, p, rankX, ...
        nBoot, seedBootstrap, ...
        testFraction, seedValidation, ...
        'VariableNames', {'RunTag', 'Response', 'Predictors', 'Taus', ...
                          'N', 'P', 'RankX', 'nBoot', 'SeedBootstrap', ...
                          'TestFraction', 'SeedValidation'});
    writetable(configTable, ...
        fullfile(outputFolder, ['runConfig_' safeResponse '_' runTag '.csv']));
end

fprintf('\nDone.\n');
if saveResults
    fprintf('Outputs saved to: %s\n', outputFolder);
end
