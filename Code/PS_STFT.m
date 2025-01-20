clear all; close all; clc;

% Load PERG data for Signal 2
load('0001.mat');  % Load the MAT file from current directory
signal_2 = s0001.RE_1;  % Get the RE_1 column from s0001 table

% Signal padding
signal_2 = [signal_2; signal_2(end)];

% Detrend Signal 2 (PERG signal)
signal_2_original = signal_2;
signal_2 = detrend(signal_2, 'linear');

% Apply a 4th order 100 Hz low-pass Butterworth filter
fs_signal_2 = 1700;
[b, a] = butter(4, 100 / (fs_signal_2 / 2), 'low');
signal_2 = filtfilt(b, a, signal_2);

% Sampling parameters
t_signal_2 = (0:length(signal_2)-1) / fs_signal_2;

% Create Signal 1: 35 Hz sine wave
fs_signal_1 = 2000;
signal_length_1 = 1000;
t_signal_1 = (0:signal_length_1-1) / fs_signal_1;
signal_1 = sin(2 * pi * 35 * t_signal_1);

% Function to compute power spectrum
compute_power_spectrum = @(signal, fs) ...
    deal(linspace(0, fs/2, floor(length(signal)/2) + 1), ...
         abs(fft(signal)).^2 / length(signal));

% Compute and normalize power spectra
[freq_1, power_1] = compute_power_spectrum(signal_1, fs_signal_1);
[freq_2, power_2] = compute_power_spectrum(signal_2, fs_signal_2);
power_1 = power_1 / max(power_1);
power_2 = power_2 / max(power_2);

% Truncate to positive frequencies and limit to 100 Hz
freq_limit = 100;
power_1 = power_1(freq_1 <= freq_limit);
freq_1 = freq_1(freq_1 <= freq_limit);
power_2 = power_2(freq_2 <= freq_limit);
freq_2 = freq_2(freq_2 <= freq_limit);

% Create figure with publication-ready settings
fig = figure('Units', 'inches', 'Position', [0, 0, 8.5, 7.5]);
set(fig, 'Color', 'w', 'PaperPositionMode', 'auto');

% Set common font settings
fontname = 'Arial';
fontsize_labels = 10;
fontsize_title = 12;

% Define figure size and margins
figure_width = 8.5;
figure_height = 11;
left_margin = 0.15;
right_margin = 0.85;
top_margin = 0.95;
bottom_margin = 0.1;
vertical_spacing = 0.1;
horizontal_spacing = 0.1;

% Add figure labels explicitly to ensure visibility
% Column 1: Signal 1
ax1 = subplot(3, 2, 1);
plot(t_signal_1, signal_1, 'b', 'LineWidth', 1.5);
title('35 Hz Sine Wave', 'FontWeight', 'bold', 'FontSize', fontsize_title, 'FontName', fontname);
ylabel('Amplitude (\muV)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
xlabel('Time (s)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
grid on;
ylim([-1.5, 1.5]);
xlim([0, 0.5]);
set(gca, 'LineWidth', 1.5, 'FontName', fontname, 'FontSize', fontsize_labels);
annotation('textbox', [0.0825, 0.9725, 0.03, 0.03], 'String', '(a)', 'FontSize', fontsize_title, ...
           'FontWeight', 'bold', 'EdgeColor', 'none', 'FontName', fontname);

ax3 = subplot(3, 2, 3);
plot(freq_1, power_1, 'b', 'LineWidth', 1.5);
title('Power Spectrum', 'FontWeight', 'bold', 'FontSize', fontsize_title, 'FontName', fontname);
ylabel('Normalized Power', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
xlabel('Frequency (Hz)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
grid on;
xlim([0, 100]);
set(gca, 'LineWidth', 1.5, 'FontName', fontname, 'FontSize', fontsize_labels);
annotation('textbox', [0.0875, 0.6525, 0.03, 0.03], 'String', '(c)', 'FontSize', fontsize_title, ...
           'FontWeight', 'bold', 'EdgeColor', 'none', 'FontName', fontname);

% For Signal 1 STFT (subplot ax5)
ax5 = subplot(3, 2, 5);
[s, f, t_stft] = spectrogram(signal_1, 128, 64, 256, fs_signal_1);  % Use full signal
power_spectrogram = abs(s).^2;
freq_range = (f >= 0 & f <= 100);
normalized_spectrogram = power_spectrogram(freq_range, :);
normalized_spectrogram = normalized_spectrogram / max(normalized_spectrogram(:));
imagesc([0 0.5], f(freq_range), normalized_spectrogram);  % Set time range to 0.5s
colormap('jet');
title('STFT', 'FontWeight', 'bold', 'FontSize', fontsize_title, 'FontName', fontname);
xlabel('Time (s)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
ylabel('Frequency (Hz)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
axis xy;
set(gca, 'YDir', 'reverse', 'LineWidth', 1.5, 'FontName', fontname, 'FontSize', fontsize_labels);
annotation('textbox', [0.09, 0.33, 0.03, 0.03], 'String', '(e)', 'FontSize', fontsize_title, ...
           'FontWeight', 'bold', 'EdgeColor', 'none', 'FontName', fontname);

% Column 2: Signal 2
ax2 = subplot(3, 2, 2);
plot(t_signal_2, signal_2_original, 'r', 'LineWidth', 1.5);
title('PERG', 'FontWeight', 'bold', 'FontSize', fontsize_title, 'FontName', fontname);
ylabel('Amplitude (\muV)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
xlabel('Time (s)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
grid on;
ylim([-7, 7]);
xlim([0, max(t_signal_2)]);
set(gca, 'LineWidth', 1.5, 'FontName', fontname, 'FontSize', fontsize_labels);
annotation('textbox', [0.525, 0.9725, 0.03, 0.03], 'String', '(b)', 'FontSize', fontsize_title, ...
           'FontWeight', 'bold', 'EdgeColor', 'none', 'FontName', fontname);

ax4 = subplot(3, 2, 4);
plot(freq_2, power_2, 'r', 'LineWidth', 1.5);
title('Power Spectrum', 'FontWeight', 'bold', 'FontSize', fontsize_title, 'FontName', fontname);
ylabel('Normalized Power', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
xlabel('Frequency (Hz)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
grid on;
xlim([0, 100]);
set(gca, 'LineWidth', 1.5, 'FontName', fontname, 'FontSize', fontsize_labels);
annotation('textbox', [0.515, 0.655, 0.03, 0.03], 'String', '(d)', 'FontSize', fontsize_title, ...
           'FontWeight', 'bold', 'EdgeColor', 'none', 'FontName', fontname);

% For Signal 2 STFT (subplot ax6)
ax6 = subplot(3, 2, 6);
[s, f, t_stft] = spectrogram(signal_2, 32, 16, 128, fs_signal_2);
power_spectrogram = abs(s).^2;
freq_range = (f >= 0 & f <= 100);
frequencies_to_plot = f(freq_range);
power_spectrogram = power_spectrogram(freq_range, :);
normalized_spectrogram = power_spectrogram / max(power_spectrogram(:));
imagesc([0 0.1506], frequencies_to_plot, normalized_spectrogram);  % Keep PERG at 0.1506s
colormap('jet');
h = colorbar;
set(h, 'Ticks', 0:0.2:1);
set(h, 'TickLabels', arrayfun(@num2str, 0:0.2:1, 'UniformOutput', false));
set(h, 'LineWidth', 1.5, 'FontName', fontname, 'FontSize', fontsize_labels);
title('STFT', 'FontWeight', 'bold', 'FontSize', fontsize_title, 'FontName', fontname);
xlabel('Time (s)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
ylabel('Frequency (Hz)', 'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname);
axis xy;
set(gca, 'YDir', 'reverse', 'LineWidth', 1.5, 'FontName', fontname, 'FontSize', fontsize_labels);
annotation('textbox', [0.5225, 0.33, 0.03, 0.03], 'String', '(f)', 'FontSize', fontsize_title, ...
           'FontWeight', 'bold', 'EdgeColor', 'none', 'FontName', fontname);

% Adjust subplot positions with more spacing
subplot_width = 0.35;
subplot_height = 0.22;

% Left column positions (adjusted for better spacing)
set(ax1, 'Position', [0.15, 0.75, subplot_width, subplot_height]);
set(ax3, 'Position', [0.15, 0.43, subplot_width, subplot_height]);
set(ax5, 'Position', [0.15, 0.11, subplot_width, subplot_height]);

% Right column positions (adjusted for better spacing)
set(ax2, 'Position', [0.58, 0.75, subplot_width, subplot_height]);
% After signal_2_original plot but before the title in ax2 section
% Find P50 (maximum) and N95 (minimum) points in the signal
[p50_amp, p50_idx] = max(signal_2_original);
[n95_amp, n95_idx] = min(signal_2_original);
p50_time = t_signal_2(p50_idx);
n95_time = t_signal_2(n95_idx);

% Convert data points to normalized figure coordinates for annotation
ax_pos = get(ax2, 'Position');
p50_x = ax_pos(1) + (p50_time/max(t_signal_2))*ax_pos(3);
p50_y = ax_pos(2) + ((p50_amp-(-7))/(7-(-7)))*ax_pos(4);
n95_x = ax_pos(1) + (n95_time/max(t_signal_2))*ax_pos(3);
n95_y = ax_pos(2) + ((n95_amp-(-7))/(7-(-7)))*ax_pos(4);

% Add P50 and N95 annotations with smaller arrow heads
annotation('textarrow', [p50_x-0.03, p50_x], [p50_y-0.01, p50_y], 'String', 'P50', ...
    'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname, ...
    'HeadStyle', 'vback2', 'HeadWidth', 6, 'HeadLength', 4);

annotation('textarrow', [n95_x-0.03, n95_x], [n95_y+0.005, n95_y], 'String', 'N95', ...
    'FontSize', fontsize_labels, 'FontWeight', 'bold', 'FontName', fontname, ...
    'HeadStyle', 'vback2', 'HeadWidth', 6, 'HeadLength', 4);
set(ax4, 'Position', [0.58, 0.43, subplot_width, subplot_height]);
set(ax6, 'Position', [0.58, 0.11, subplot_width, subplot_height]);