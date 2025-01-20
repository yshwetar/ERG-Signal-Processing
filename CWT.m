clear all; close all; clc;

% Load PERG data from .mat file
load('0001.mat'); % This loads the table with variable name s0001
perg_data = s0001.RE_1; % Get the RE_1 column from the loaded table

% Signal padding
perg_data = [perg_data; perg_data(end)]; % Add the extra padded data point

% Perform Continuous Wavelet Transform (CWT)
sampling_rate = 1700; % Sampling rate in Hz
time_axis = (0:length(perg_data)-1) / sampling_rate; % Time axis in seconds
[wt, freq] = cwt(perg_data, 'amor', sampling_rate); % Perform CWT using 'amor' (Morlet wavelet)

% Normalize CWT coefficients
wt = abs(wt); % Take the absolute value of wavelet coefficients
wt = wt / max(wt(:)); % Normalize the coefficients
valid_indices = freq <= 100; % Filter for frequencies <= 100 Hz

% Generate the first 64-sample Morlet wavelet for overlay
fc1 = 0.5; % Center frequency of the first wavelet
waveletLength1 = 64; % Length of the wavelet (64 samples)
t_wavelet1 = linspace(-2, 2, waveletLength1); % Symmetric time range for wavelet

% Morlet wavelet formula with proper decay
morlet_wavelet1 = exp(-t_wavelet1.^2 / 2) .* cos(2 * pi * fc1 * t_wavelet1);

% Rescale the wavelet to match the PERG signal's amplitude range
pergMin = min(perg_data);
pergMax = max(perg_data);
scaling_factor = 0.3; % Set the scaling factor (less than 1 to shrink amplitude)
morlet_wavelet1 = scaling_factor * ((morlet_wavelet1 - min(morlet_wavelet1)) / (max(morlet_wavelet1) - min(morlet_wavelet1)) * (pergMax - pergMin) + pergMin);

% Generate the second 64-sample Morlet wavelet for overlay
fc2 = 2; % Center frequency of the second wavelet (increase this for higher oscillations)
waveletLength2 = 64; % Length of the wavelet (64 samples)
t_wavelet2 = linspace(-2, 2, waveletLength2); % Symmetric time range for wavelet

% Generate the second Morlet wavelet with higher oscillations
morlet_wavelet2 = exp(-t_wavelet2.^2 / 2) .* cos(2 * pi * fc2 * t_wavelet2);

% Scale the amplitude of the wavelet
scaling_factor2 = 3; % Set the scaling factor for larger amplitude
morlet_wavelet2 = scaling_factor2 * morlet_wavelet2; % Scale the amplitude

% Align the vertical middle of the wavelet with the waveform in the region
waveform_segment = perg_data(129:192); % Extract the region of the waveform
center_offset = mean(waveform_segment) - mean(morlet_wavelet2); % Calculate the vertical offset
morlet_wavelet2 = morlet_wavelet2 + center_offset; % Shift the wavelet vertically

% Create a figure with adjusted dimensions
figure('Units', 'inches', 'Position', [1, 1, 7.5, 6]); % Adjust figure size for heatmap scale visibility
set(gcf, 'Color', 'w'); % Set background color to white

% Top subplot: PERG signal with overlaid Morlet wavelets
ax1 = subplot(2, 1, 1);
hold on;
plot(time_axis, perg_data, 'k', 'LineWidth', 2, 'DisplayName', 'PERG');
wavelet_indices1 = time_axis(1:64);
plot(wavelet_indices1, morlet_wavelet1, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Wavelet 1');
wavelet_indices2 = time_axis(129:192);
plot(wavelet_indices2, morlet_wavelet2, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Wavelet 2');
ylabel('Amplitude (\muV)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('show', 'FontSize', 10);
title('PERG', 'FontSize', 11, 'FontWeight', 'bold'); % Reduced font size for the top plot title

% Add figure label (a)
text(-0.08, 1.1, '(a)', 'FontSize', 10, 'FontWeight', 'bold', 'Units', 'normalized', 'VerticalAlignment', 'top');

% Set border thickness for the PERG plot
set(gca, 'LineWidth', 1.5, 'Box', 'on', 'FontWeight', 'bold');
hold off;

% Bottom subplot: CWT plot
ax2 = subplot(2, 1, 2);
imagesc(time_axis, freq(valid_indices), wt(valid_indices, :));
colormap('jet'); % Revert to the original colormap
cbar = colorbar('FontSize', 10, 'LineWidth', 1.5); % Add colorbar and capture its handle
caxis([0, 1]); % Set color range for better contrast
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('Continuous Wavelet Transform', 'FontSize', 16, 'FontWeight', 'bold'); % Keep font size for the bottom plot
set(gca, 'FontSize', 10, 'LineWidth', 1.5, 'FontWeight', 'bold');
axis ij; % Flip y-axis for CWT

% Add a large red star at the center of wavelet 1 and 15 Hz
time_wavelet1_center = mean(wavelet_indices1); % Calculate the center time of wavelet 1
freq_star1 = 22; % Frequency at 15 Hz
hold on;
plot(time_wavelet1_center, freq_star1, 'w*', 'MarkerSize', 20, 'LineWidth', 2); % Larger red star for wavelet 1

% Highlight around the star with a circle for wavelet 1
% radius1 = 5; % Radius of the circle
% theta = linspace(0, 2*pi, 100); % Angle for the circle
% x_circle1 = time_wavelet1_center + radius1 * cos(theta) / sampling_rate; % X-coordinates of the circle
% y_circle1 = freq_star1 + radius1 * sin(theta); % Y-coordinates of the circle
% plot(x_circle1, y_circle1, 'g-', 'LineWidth', 1.5); % Green circular highlight

% Add a large red star at the center of wavelet 2 and 85 Hz
time_wavelet2_center = mean(wavelet_indices2); % Calculate the center time of wavelet 2
freq_star2 = 85; % Frequency at 85 Hz
plot(time_wavelet2_center, freq_star2, 'w*', 'MarkerSize', 20, 'LineWidth', 2); % Larger red star for wavelet 2

% Highlight around the star with a circle for wavelet 2
% radius2 = 5; % Radius of the circle
% x_circle2 = time_wavelet2_center + radius2 * cos(theta) / sampling_rate; % X-coordinates of the circle
% y_circle2 = freq_star2 + radius2 * sin(theta); % Y-coordinates of the circle
% plot(x_circle2, y_circle2, 'g-', 'LineWidth', 1.5); % Green circular highlight

hold off;


% Add figure label (b)
text(-0.08, 1.1, '(b)', 'FontSize', 10, 'FontWeight', 'bold', 'Units', 'normalized', 'VerticalAlignment', 'top');

% Adjust subplot spacing
set(ax1, 'Position', [0.1, 0.58, 0.75, 0.37]); % Adjust position of the top subplot
set(ax2, 'Position', [0.1, 0.1, 0.75, 0.37]); % Adjust position of the bottom subplot

% Save the figure
% set(gcf, 'PaperPositionMode', 'auto'); % Fit figure to page
% print('-dpng', '-r00', 'CWT.png'); % Save as PNG
% saveas(gcf, 'CWT.pdf'); % Save as PDF
