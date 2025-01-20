%% Signal Analysis Toolkit
% Author: Yousif Shwetar
% GitHub: @yshwetar
% Publication: Signal Processing in Electroretinography
% Version: 1.0.0
% Description: A comprehensive signal analysis toolkit that provides multiple
% analysis methods including Power Spectrum, STFT, CWT, and DWT analysis.
% For documentation and usage examples, please see the README.md file.

%% Initialize workspace
close all;
clear all;
clc;

%% Load and Process Input Signal
% Load the MAT file
load('0001.mat');
signal = s0001.RE_1;

% Pad signal by repeating the final value once
signal = [signal; signal(end)];

% Sampling Parameters
fs = 1700;                  % Sampling frequency (Hz)
N = length(signal);         % Number of samples including padding
t = (0:N-1)/fs;            % Time vector calculated from padded signal length

%% Analysis Parameters

% 1. Plot Y-Axis Frequency Limits
freq_min = 0;              % Minimum frequency to show in plots (Hz)
freq_max = fs/2;           % Maximum frequency to show in plots (Hz)

% 2. STFT Parameters
window_type = 'hann';      % Options: 'hann', 'hamming', 'blackman', 'rectwin'
window_length_ms = 128;    % Window length in milliseconds
overlap_percent = 50;      % Overlap percentage (0-100)

% 3. CWT Parameters
cwt_wavelet_type = 'amor'; % Options: 'amor', 'bump', 'morl', 'morse', 'mexh'

% 4. DWT Parameters
dwt_wavelet_type = 'db2';  % Options: 'haar', 'db4', 'sym4', etc.

% 5. Visualization Parameters
fontname = 'Arial';        % Font for all plots
fontsize_labels = 10;      % Font size for axis labels
fontsize_title = 12;       % Font size for titles

%% Power Spectrum Analysis
% Calculate frequency axis
freq = linspace(0, fs/2, floor(N/2) + 1);

% Compute power spectrum
power = abs(fft(signal)).^2 / N;
power = power(1:length(freq));
power = power / max(power);

% Create visualization with custom frequency limits
createPowerSpectrumPlot(freq, power, freq_min, freq_max, fontname, fontsize_labels, fontsize_title);

%% STFT Analysis
% Parameters matching the reference code
window_length = 32;  % For signal 2 (use 128 for signal 1)
overlap = 16;       % For signal 2 (use 64 for signal 1)
nfft = 128;        % For signal 2 (use 256 for signal 1)

% Create window
window = createWindow(window_type, window_length);

% Compute STFT
[s, f, t] = spectrogram(signal, window, overlap, nfft, fs);
power_spectrogram = abs(s).^2;

% Filter frequency range and normalize
freq_range = (f >= 0 & f <= 100);
frequencies_to_plot = f(freq_range);
power_spectrogram = power_spectrogram(freq_range, :);
normalized_spectrogram = power_spectrogram / max(power_spectrogram(:));

% Create visualization with custom frequency limits
createSTFTPlot(t, frequencies_to_plot, normalized_spectrogram, window_type, freq_min, freq_max, ...
    fontname, fontsize_labels, fontsize_title);

%% CWT Analysis
% Perform CWT
[wt, freq] = cwt(signal, cwt_wavelet_type, fs);

% Normalize CWT coefficients
wt = abs(wt);
wt = wt / max(wt(:));

% Create visualization with custom frequency limits
createCWTPlot(t, freq, wt, freq_min, freq_max, fontname, fontsize_labels, fontsize_title);

%% DWT Analysis
% Calculate maximum decomposition level
max_possible = wmaxlev(length(signal), dwt_wavelet_type);
max_level = determineDecompositionLevel(dwt_wavelet_type, max_possible);

% Perform wavelet decomposition
[c, l] = wavedec(signal, max_level, dwt_wavelet_type);

% Process coefficients and create heatmap matrix
[heatmap_matrix_combined, time_vector] = processDWTCoefficients(c, l, signal, ...
    fs, max_level, dwt_wavelet_type);

% Create visualization
createDWTPlot(time_vector, max_level, heatmap_matrix_combined, ...
    fontname, fontsize_labels, fontsize_title);

%% Helper Functions

function window = createWindow(window_type, window_length)
    switch window_type
        case 'hann'
            window = hann(window_length);
        case 'hamming'
            window = hamming(window_length);
        case 'blackman'
            window = blackman(window_length);
        case 'rectwin'
            window = rectwin(window_length);
        otherwise
            error('Unsupported window type');
    end
end

function max_level = determineDecompositionLevel(wavelet_type, max_possible)
    if strcmpi(wavelet_type, 'haar')
        max_level = max_possible - 1;  % Use n-1 for Haar
    else
        max_level = max_possible;      % Use n for all others
    end
end

function [heatmap_matrix_combined, time_vector] = processDWTCoefficients(c, l, signal, fs, max_level, wavelet_type)
    duration = (length(signal) - 1)/fs;
    time_vector = linspace(0, duration, length(signal));
    
    % Process detail coefficients
    heatmap_matrix_detail = zeros(max_level, length(signal));
    for level = 1:max_level
        detail_coeffs = detcoef(c, l, level);
        detail_energy = detail_coeffs .^ 2;
        
        % Calculate step size for each coefficient
        step_size = floor(length(signal) / length(detail_energy));
        detail_extended = [];
        
        % Repeat each value step_size times
        for i = 1:length(detail_energy)
            detail_extended = [detail_extended, repmat(detail_energy(i), 1, step_size)];
        end
        
        % Pad or trim to match signal length
        if length(detail_extended) < length(signal)
            detail_extended = [detail_extended, repmat(detail_energy(end), 1, length(signal) - length(detail_extended))];
        elseif length(detail_extended) > length(signal)
            detail_extended = detail_extended(1:length(signal));
        end
        
        heatmap_matrix_detail(level, :) = detail_extended;
    end
    
    % Process approximation coefficients
    approx_coeffs = appcoef(c, l, wavelet_type, max_level);
    approx_energy = approx_coeffs .^ 2;
    
    % Apply the same step-wise extension for approximation coefficients
    step_size = floor(length(signal) / length(approx_energy));
    approx_extended = [];
    
    % Repeat each value step_size times
    for i = 1:length(approx_energy)
        approx_extended = [approx_extended, repmat(approx_energy(i), 1, step_size)];
    end
    
    % Pad or trim to match signal length
    if length(approx_extended) < length(signal)
        approx_extended = [approx_extended, repmat(approx_energy(end), 1, length(signal) - length(approx_extended))];
    elseif length(approx_extended) > length(signal)
        approx_extended = approx_extended(1:length(signal));
    end
    
    approx_interp = approx_extended;
    
    % Normalize and combine
    heatmap_matrix_detail = normalizeMatrix(heatmap_matrix_detail);
    approx_row_norm = normalizeVector(approx_interp);
    
    % Ensure all matrices are the correct size
    if size(heatmap_matrix_detail, 2) ~= length(approx_row_norm)
        approx_row_norm = approx_row_norm(1:size(heatmap_matrix_detail, 2));
    end
    
    heatmap_matrix_combined = [heatmap_matrix_detail; approx_row_norm];
end

function normalized = normalizeMatrix(matrix)
    min_val = min(matrix(:));
    max_val = max(matrix(:));
    if max_val == min_val
        normalized = zeros(size(matrix));
    else
        normalized = (matrix - min_val) / (max_val - min_val);
    end
end

function normalized = normalizeVector(vector)
    min_val = min(vector);
    max_val = max(vector);
    if max_val == min_val
        normalized = zeros(size(vector));
    else
        normalized = (vector - min_val) / (max_val - min_val);
    end
end

function createPowerSpectrumPlot(freq, power, freq_min, freq_max, fontname, fontsize_labels, fontsize_title)
    fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 10, 6]);
    set(fig, 'Color', 'w');
    
    % Filter data within frequency range
    freq_mask = freq >= freq_min & freq <= freq_max;
    freq_plot = freq(freq_mask);
    power_plot = power(freq_mask);
    
    plot(freq_plot, power_plot, 'b', 'LineWidth', 1.5);
    title('Power Spectrum', 'FontWeight', 'bold', 'FontSize', fontsize_title, 'FontName', fontname);
    ylabel('Normalized Power', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
    xlabel('Frequency (Hz)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
    grid on;
    xlim([freq_min, freq_max]);
    ylim([0, 1.05]);
    
    set(gca, 'LineWidth', 1.5, ...
        'FontName', fontname, ...
        'FontSize', fontsize_labels, ...
        'Position', [0.15, 0.15, 0.75, 0.75]);
end

function createSTFTPlot(t, f, normalized_spectrogram, window_type, freq_min, freq_max, ...
    fontname, fontsize_labels, fontsize_title)
    fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 10, 6]);
    set(fig, 'Color', 'w');
    
    imagesc(t, f, normalized_spectrogram);
    colormap('jet');
    h = colorbar;
    set(h, 'Ticks', 0:0.2:1);
    set(h, 'TickLabels', arrayfun(@num2str, 0:0.2:1, 'UniformOutput', false));
    set(h, 'LineWidth', 1.5, 'FontName', fontname, 'FontSize', fontsize_labels);
    
    title(['STFT Analysis (', window_type, ' window)'], 'FontWeight', 'bold', ...
        'FontSize', fontsize_title, 'FontName', fontname);
    xlabel('Time (s)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
    ylabel('Frequency (Hz)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
    
    axis xy;
    set(gca, 'YDir', 'reverse', 'LineWidth', 1.5, ...
        'FontName', fontname, ...
        'FontSize', fontsize_labels, ...
        'Position', [0.15, 0.15, 0.7, 0.75]);
end

function createCWTPlot(t, freq, wt, freq_min, freq_max, fontname, fontsize_labels, fontsize_title)
    figure('Units', 'inches', 'Position', [0.1, 0.1, 10, 6], 'Color', 'w');
    
    % Filter data within frequency range
    freq_mask = freq >= freq_min & freq <= freq_max;
    freq_plot = freq(freq_mask);
    wt_plot = wt(freq_mask, :);
    
    imagesc(t, freq_plot, wt_plot);
    colormap('jet');
    cbar = colorbar('FontSize', fontsize_labels, 'LineWidth', 1.5);
    caxis([0, 1]);
    
    xlabel('Time (s)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
    ylabel('Frequency (Hz)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
    title('Continuous Wavelet Transform', 'FontSize', fontsize_title, 'FontWeight', 'bold', 'FontName', fontname);
    
    set(gca, 'FontSize', fontsize_labels, 'LineWidth', 1.5, 'FontWeight', 'bold', ...
        'FontName', fontname, 'Position', [0.1, 0.15, 0.75, 0.75]);
    axis ij;
end

function createDWTPlot(time_vector, max_level, heatmap_matrix_combined, ...
    fontname, fontsize_labels, fontsize_title)
    figure('Units', 'inches', 'Position', [0.1, 0.1, 10, 6], 'Color', 'w');
    
    imagesc(time_vector, 1:(max_level+1), heatmap_matrix_combined);
    colormap('jet');
    cb = colorbar;
    caxis([0 1]);
    set(cb, 'LineWidth', 1.5, 'FontWeight', 'bold', 'FontName', fontname);
    
    titleStr = sprintf('Wavelet Transform (Level %d decomposition)', max_level);
    title(titleStr, 'FontWeight', 'bold', 'FontSize', fontsize_title, 'FontName', fontname);
    xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', fontsize_labels, 'FontName', fontname);
    ylabel('Decomposition Level', 'FontWeight', 'bold', 'FontSize', fontsize_labels, 'FontName', fontname);
    
    % Y-axis labels
    yticks(1:(max_level+1));
    ylabels = cell(max_level+1, 1);
    for i = 1:max_level
        ylabels{i} = ['D' num2str(i)];
    end
    ylabels{end} = ['A' num2str(max_level)];
    yticklabels(ylabels);
    
    set(gca, 'FontWeight', 'bold', 'LineWidth', 1.5, 'YDir', 'normal', ...
        'FontName', fontname, 'FontSize', fontsize_labels);
    box on;
end