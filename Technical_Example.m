% Step 1: Clear Environment
close all; clear all; clc;

% Sampling rate
sampling_rate = 1700;

% Design a 4th-order Butterworth low-pass filter with a 100 Hz cutoff
[b, a] = butter(4, 100 / (sampling_rate / 2), 'low');

% Define participant IDs and corresponding ERG titles
participants = {'0001', '0019'};
titles_erg = {'Healthy', 'Macular Dystrophy'};

% Initialize figure with portrait dimensions
figure('Units', 'inches', 'Position', [0.1, 0.1, 7.5, 10], 'PaperPositionMode', 'auto');

% Define plot layout parameters
row_count = 5;
col_count = length(participants);

% Loop over participants
for col = 1:length(participants)
    % Load participant data from .mat file
    participant_id = participants{col};
    load([participant_id '.mat']); % This will load the table with variable name s0001 or s0019
    
    % Get the variable name dynamically
    var_name = ['s' participant_id];
    erg_data = eval(var_name);

    % Identify the first `RE` recording
    re_recordings = regexpi(erg_data.Properties.VariableNames, '^RE_\d+$', 'match');
    re_recordings = [re_recordings{:}];

    if isempty(re_recordings)
        error('No `RE` recordings found for participant %s.', participant_id);
    end

    % Load the first `RE` recording
    signal = erg_data.(re_recordings{1});

    % Step 1: Pad the signal by repeating the last value
    if length(signal) == 255
        signal = [signal; signal(end)]; % Add the last value as the 256th sample
    end

    % Create a separate variable for the power spectrum
    signal_for_ps = detrend(signal, 'linear');
    signal_for_ps = filtfilt(b, a, signal_for_ps);

    % Time axis for original signal
    time_axis = (0:length(signal)-1) / sampling_rate;

    % Row 1: ERG Signal (Original Signal)
    ax1 = subplot(row_count, col_count, col); % Row 1
    text(0.02, 0.98, char('a' + (col-1)), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    plot(time_axis, signal, 'LineWidth', 1.2);
    ylim([-7, 7]); % Standardize scale for ERG plots
    xlim([0, 0.1499]);
    title(titles_erg{col}, 'FontSize', 10, 'Interpreter', 'none');
    if col == 1
        yl = ylabel('Amplitude (\muV)', 'FontSize', 8, 'FontWeight', 'bold');
        set(yl, 'Units', 'normalized');
        set(yl, 'Position', [-0.15, 0.5, 0]);
    end
    xlabel('Time (sec)', 'FontSize', 8, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 8);

    % Row 2: Power Spectrum (Filtered Signal)
    ax2 = subplot(row_count, col_count, col + col_count); % Row 2
    text(0.02, 0.98, char('d' + (col-1)), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    N = length(signal_for_ps);
    fft_result = fft(signal_for_ps);
    freqs = (0:floor(N/2)) * (sampling_rate / N);
    power_spectrum = abs(fft_result(1:floor(N/2)+1)).^2;
    normalized_power_spectrum = (power_spectrum / sum(power_spectrum)) * 100 / 100; % Normalize to 0-1 scale
    valid_indices = freqs <= 100;
    plot(freqs(valid_indices), normalized_power_spectrum(valid_indices), ...
        'LineWidth', 1.2, 'Color', 'r');
    xlim([0, 100]);
    ylim([0, 1]);  % Changed to 0-1 range
    if col == 1
        yl = ylabel('Normalized Power', 'FontSize', 8, 'FontWeight', 'bold');
        set(yl, 'Units', 'normalized');
        set(yl, 'Position', [-0.15, 0.5, 0]);
    end
    xlabel('Frequency (Hz)', 'FontSize', 8, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 8);

    % Save the position of the first two rows to match later rows
    pos1 = get(ax1, 'Position');
    pos2 = get(ax2, 'Position');

    % Row 3: STFT Spectrogram (Filtered Signal)
    ax3 = subplot(row_count, col_count, col + 2 * col_count); % Row 3
    text(0.02, 0.98, char('g' + (col-1)), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    window_size = 32; % Bartlett window size
    overlap = 4; % Overlap between windows
    nfft = max(128, 2^nextpow2(window_size)); % FFT points
    [s, f, t] = spectrogram(signal_for_ps, bartlett(window_size), overlap, nfft, sampling_rate);
    power_spectrogram = abs(s).^2; % Compute power
    freq_range = (f >= 0 & f <= 100); % Logical index for 0–100 Hz
    f_filtered = f(freq_range); % Filtered frequencies
    power_spectrogram_filtered = power_spectrogram(freq_range, :);
    normalized_power_spectrogram = power_spectrogram_filtered / max(power_spectrogram_filtered(:)); % Normalize by max value
    imagesc(t, flip(f_filtered), flipud(normalized_power_spectrogram));
    colormap('jet');
    % if col == 2  % Changed from 3 to 2
    %     colorbar('FontSize', 8);
    % end
    xlim([0, 0.1499]);
    if col == 1
        yl = ylabel('Frequency (Hz)', 'FontSize', 8, 'FontWeight', 'bold');
        set(yl, 'Units', 'normalized');
        set(yl, 'Position', [-0.15, 0.5, 0]);
    end
    xlabel('Time (sec)', 'FontSize', 8, 'FontWeight', 'bold');
    set(gca, 'XTick', 0:0.05:0.1);
    set(gca, 'FontSize', 8);

    % Adjust width and alignment
    pos3 = get(ax3, 'Position');
    pos3(1) = pos1(1);
    pos3(3) = pos1(3);
    set(ax3, 'Position', pos3);

    % Row 4: CWT (Original Signal)
    ax4 = subplot(row_count, col_count, col + 3 * col_count); % Row 4
    text(0.02, 0.98, char('j' + (col-1)), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    [wt, freq] = cwt(signal, 'amor', sampling_rate);
    wt = abs(wt) / max(abs(wt(:))); % Normalize wavelet transform
    valid_indices = freq <= 100;
    imagesc(time_axis, freq(valid_indices), wt(valid_indices, :));
    colormap('jet');
    % if col == 2  % Changed from 3 to 2
    %     colorbar('FontSize', 8);
    % end
    xlim([0, 0.1499]);
    if col == 1
        yl = ylabel('Frequency (Hz)', 'FontSize', 8, 'FontWeight', 'bold');
        set(yl, 'Units', 'normalized');
        set(yl, 'Position', [-0.15, 0.5, 0]);
    end
    xlabel('Time (sec)', 'FontSize', 8, 'FontWeight', 'bold');
    set(gca, 'FontSize', 8);

    % Adjust width and alignment
    pos4 = get(ax4, 'Position');
    pos4(1) = pos1(1);
    pos4(3) = pos1(3);
    set(ax4, 'Position', pos4);

    % Row 5: DWT Heatmap (Original Signal)
    ax5 = subplot(row_count, col_count, col + 4 * col_count); % Row 5
    text(0.02, 0.98, char('m' + (col-1)), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    waveletType = 'haar';
    max_level = 7;
    [c, l] = wavedec(signal, max_level, waveletType);

    % Get approximation coefficients at level 7
    a7 = appcoef(c, l, waveletType, max_level);
    a7_energy = a7 .^ 2;

    % Create heatmap matrix for detail coefficients
    selected_levels = 1:7; % Only detail coefficients
    heatmap_matrix = zeros(length(selected_levels), length(signal));

    % Fill in detail coefficients (levels 1-7)
    for level_idx = 1:length(selected_levels)
        level = selected_levels(level_idx);
        coeff = detcoef(c, l, level);
        coeff_energy = coeff .^ 2;
        heatmap_matrix(level_idx, :) = ...
            repelem(coeff_energy, ceil(length(signal) / length(coeff_energy)));
    end

    % Normalize detail coefficients separately
    heatmap_matrix_normalized = heatmap_matrix / max(heatmap_matrix(:));

    % Normalize approximation coefficients separately and reshape to match matrix width
    a7_extended = repelem(a7_energy, ceil(length(signal) / length(a7_energy)));
    % Trim or pad to match signal length
    a7_extended = a7_extended(1:size(heatmap_matrix_normalized, 2));
    a7_normalized = a7_extended / max(a7_extended);
    % Reshape to row matrix
    a7_normalized = reshape(a7_normalized, 1, []);

    % Combine normalized matrices
    final_heatmap = [heatmap_matrix_normalized; a7_normalized];

    % Display the heatmap
    imagesc(linspace(0, length(signal) / sampling_rate, size(final_heatmap, 2)), ...
        1:8, final_heatmap);
    colormap('jet');
    if col == 2  % Changed from 3 to 2
        colorbar('FontSize', 8);
    end
    xlabel('Time (sec)', 'FontSize', 8, 'FontWeight', 'bold');
    if col == 1
        yl = ylabel('Decomposition Level', 'FontSize', 8, 'FontWeight', 'bold');
        set(yl, 'Units', 'normalized');
        set(yl, 'Position', [-0.15, 0.5, 0]);
    end
    yticks(1:8);
    yticklabels({'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'A7'});
    xlim([0, 0.1499]);
    set(gca, 'FontSize', 8, 'YDir', 'normal');

    % Adjust width and alignment
    pos5 = get(ax5, 'Position');
    pos5(1) = pos1(1);
    pos5(3) = pos1(3);
    set(ax5, 'Position', pos5);
end

% Beautify the Figure
set(gcf, 'Color', 'w');