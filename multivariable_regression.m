%% Multivariable Regression Analysis for Megagravel Transport
%
% Fits three regression models to total horizontal boulder movement:
%   1. Multivariable ordinary least squares (natural units)
%   2. Log-transformed OLS
%   3. Principal component regression (first three PCs)
%
% Produces residual diagnostics, coefficient tables, and PCA output used
% in Chapter 4 of the thesis (Tables 4-1 to 4-5, Figures 4-4 to 4-7).
%
% Requires: 2014_Boulder_Movements_dataset_of_record.xlsx in the working
% directory.

clear; clc; close all;

%% Load data
data = readtable('2014_Boulder_Movements_dataset_of_record.xlsx');

%% Response variable
Y = data.TotalHorizontalMovement_m__absoluteValue_;
responseLabel = 'Total horizontal movement (m)';

fprintf('Response: %s\n', responseLabel);
fprintf('  n = %d, missing = %d\n', numel(Y), sum(isnan(Y)));
fprintf('  range = [%.2f, %.2f], mean = %.2f, median = %.2f\n', ...
    min(Y), max(Y), mean(Y, 'omitnan'), median(Y, 'omitnan'));

% Plot response distribution
figure;
histogram(Y);
xlabel(responseLabel);
ylabel('Frequency');
title('Distribution of total horizontal movement');
set(gcf, 'Color', 'w');
box on;
set(gca, 'LineWidth', 1, 'FontSize', 11);

%% Predictors
predictorNames = {'Mass_t_', 'Vol_m3_', 'HeightAHW_m_', ...
                  'DistanceInlandAHW_m_', 'Steepness'};

X = data(:, predictorNames);

% Build model table and remove rows with missing values
T = [table(Y), X];
T = rmmissing(T);

fprintf('\nCleaned model table: %d rows, %d predictors\n', ...
    height(T), numel(predictorNames));

%% Multivariable OLS regression (natural units)
Yreg = T.Y;
Xreg = [ones(height(T), 1), table2array(T(:, 2:end))];

b = Xreg \ Yreg;
Yhat = Xreg * b;
resid = Yreg - Yhat;

R2 = 1 - sum(resid.^2) / sum((Yreg - mean(Yreg)).^2);

fprintf('\nOLS (natural units):\n');
fprintf('  R^2 = %.4f\n', R2);

% Residuals vs fitted
figure;
scatter(Yhat, resid, 'filled');
xlabel('Fitted values');
ylabel('Residuals');
title('Residuals vs fitted values');
yline(0, '--');
set(gcf, 'Color', 'w');
box on;
set(gca, 'LineWidth', 1, 'FontSize', 11);

% Residual histogram
figure;
histogram(resid);
xlabel('Residual');
ylabel('Frequency');
title('Residual distribution');
set(gcf, 'Color', 'w');
box on;
set(gca, 'LineWidth', 1, 'FontSize', 11);

%% Log-transformed OLS regression
% Exclude non-positive movements before log transform
Tpos = T(T.Y > 0, :);
Ylog = log(Tpos.Y);
Xlog = [ones(height(Tpos), 1), table2array(Tpos(:, 2:end))];

b_log = Xlog \ Ylog;
Yhat_log = Xlog * b_log;
resid_log = Ylog - Yhat_log;

R2_log = 1 - sum(resid_log.^2) / sum((Ylog - mean(Ylog)).^2);

fprintf('\nOLS (log-transformed):\n');
fprintf('  n = %d (after excluding non-positive movement)\n', height(Tpos));
fprintf('  R^2 = %.4f\n', R2_log);

% Residuals vs fitted (log)
figure;
scatter(Yhat_log, resid_log, 'filled');
xlabel('Fitted values (log scale)');
ylabel('Residuals');
title('Residuals vs fitted values (log-transformed model)');
yline(0, '--');
set(gcf, 'Color', 'w');
box on;
set(gca, 'LineWidth', 1, 'FontSize', 11);

% Residual histogram (log)
figure;
histogram(resid_log);
xlabel('Residual');
ylabel('Frequency');
title('Residual distribution (log-transformed model)');
set(gcf, 'Color', 'w');
box on;
set(gca, 'LineWidth', 1, 'FontSize', 11);

%% Principal component analysis
Xmat = table2array(T(:, 2:end));
Xz = (Xmat - mean(Xmat)) ./ std(Xmat);

% PCA via SVD (no toolbox required)
[U, S, V] = svd(Xz, 'econ');

latent = diag(S).^2 / (size(Xz, 1) - 1);
explained = 100 * latent / sum(latent);

score = U * S;   % component scores
coeff = V;       % component loadings

fprintf('\nPCA variance explained (%%):\n');
for i = 1:numel(explained)
    fprintf('  PC%d: %.2f  (cumulative: %.2f)\n', ...
        i, explained(i), sum(explained(1:i)));
end

%% Principal component regression (first 3 PCs)
PC = score(:, 1:3);
PCreg = [ones(size(PC, 1), 1), PC];

b_pcr = PCreg \ T.Y;
Yhat_pcr = PCreg * b_pcr;
resid_pcr = T.Y - Yhat_pcr;

R2_pcr = 1 - sum(resid_pcr.^2) / sum((T.Y - mean(T.Y)).^2);

fprintf('\nPrincipal component regression (3 PCs):\n');
fprintf('  R^2 = %.4f\n', R2_pcr);

%% Save results
saveResults = false;   % set to true to write outputs to disk
outputFolder = 'multivariable_regression_results';

if saveResults
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    runTag = datestr(now, 'yyyymmdd_HHMMSS');

    coefNames = ["Intercept", string(predictorNames)]';

    % OLS coefficients
    coefTable = table(coefNames, b, 'VariableNames', {'Variable', 'Coefficient'});
    writetable(coefTable, fullfile(outputFolder, ...
        ['multivar_coef_' runTag '.csv']));

    % Log-OLS coefficients
    coefTable_log = table(coefNames, b_log, 'VariableNames', {'Variable', 'Coefficient'});
    writetable(coefTable_log, fullfile(outputFolder, ...
        ['log_coef_' runTag '.csv']));

    % R^2 summary
    R2Table = table(R2, R2_log, R2_pcr, ...
        'VariableNames', {'R2_OLS', 'R2_log', 'R2_PCR'});
    writetable(R2Table, fullfile(outputFolder, ...
        ['R2_summary_' runTag '.csv']));

    % PCA variance
    pcaTable = table((1:numel(explained))', explained, cumsum(explained), ...
        'VariableNames', {'PC', 'VarianceExplained', 'CumulativeVariance'});
    writetable(pcaTable, fullfile(outputFolder, ...
        ['pca_variance_' runTag '.csv']));

    % PCA loadings
    loadingsTable = array2table(coeff, ...
        'VariableNames', compose("PC%d", 1:size(coeff, 2)), ...
        'RowNames', cellstr(string(predictorNames)));
    writetable(loadingsTable, fullfile(outputFolder, ...
        ['pca_loadings_' runTag '.csv']), 'WriteRowNames', true);

    fprintf('\nResults saved to: %s\n', outputFolder);
end

fprintf('\nDone.\n');
