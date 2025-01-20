clear all; close all; clc;

% Create time vector
fs = 100; % Sampling frequency (Hz)
t = 0:1/fs:0.5-1/fs;  % This creates 50 points (excludes 0.5)
N = length(t);

% Create frequency components in order of increasing frequency
f1 = 2;  % Low frequency (Hz)
f2 = 10; % Medium frequency (Hz)
f3 = 25; % High frequency (Hz)

% Generate individual signals
low_freq = 1.5 * sin(2*pi*f1*t);
med_freq = 0.75 * sin(2*pi*f2*t);
high_freq = 0.5 * sin(2*pi*f3*t);
combined = low_freq + med_freq + high_freq;

% Calculate frequency axis for FFT plots
nPoints = floor(N/2) + 1;
freq = fs * (0:(nPoints-1))/N;

% Calculate FFTs
fft_low = abs(fft(low_freq)/N);
fft_med = abs(fft(med_freq)/N);
fft_high = abs(fft(high_freq)/N);
fft_combined = abs(fft(combined)/N);

% Only keep first half of FFT
fft_low = fft_low(1:nPoints);
fft_med = fft_med(1:nPoints);
fft_high = fft_high(1:nPoints);
fft_combined = fft_combined(1:nPoints);

% Normalize each spectrum with respect to its maximum value
fft_low = fft_low / max(fft_low);
fft_med = fft_med / max(fft_med);
fft_high = fft_high / max(fft_high);
fft_combined = fft_combined / max(fft_combined);

% Update y-axis limit for frequency plots
fft_max = 1.1; % Since all spectra are normalized to 1

% Create figure
inches_to_pixels = 300;
width_inches = 8.5;
height_inches = 11;
figure('Position', [100 100 width_inches*inches_to_pixels height_inches*inches_to_pixels]);

% Calculate axis limits
y_max = max(abs([low_freq med_freq high_freq combined])) * 1.1;
fft_max = max(abs([fft_low fft_med fft_high fft_combined])) * 1.1;

% Store signals in order from top to bottom to match subplot ordering
time_signals = {combined, high_freq, med_freq, low_freq};
freq_signals = {fft_combined, fft_high, fft_med, fft_low};
colors = {[0.5 0.2 0.7], [0.3 0.8 0.3], [0.8 0.4 0.2], [0.2 0.6 0.8]};

% Create subplot indices map for correct ordering
row_indices = [1 2 3 4];  % Top to bottom
subplot_indices = [row_indices*2-1; row_indices*2]';  % Creates pairs: [1 2; 3 4; 5 6; 7 8]
letters = 'h':-1:'a';

titles_time = {
    ['Combined Signal: (a) + (c) + (e)'],
    ['High Frequency Signal (' num2str(f3) ' Hz)'],
    ['Medium Frequency Signal (' num2str(f2) ' Hz)'],
    ['Low Frequency Signal (' num2str(f1) ' Hz)']
};

titles_freq = {
    'Combined Signal Spectrum',
    ['High Frequency Spectrum (' num2str(f3) ' Hz)'],
    ['Medium Frequency Spectrum (' num2str(f2) ' Hz)'],
    ['Low Frequency Spectrum (' num2str(f1) ' Hz)']
};

% Plot all subplots
for i = 1:4
    % Get subplot indices for this row
    idx_freq = subplot_indices(i,1);  % Left column index
    idx_time = subplot_indices(i,2);  % Right column index
    
    % Frequency domain plots (left column)
    subplot(4,2,idx_freq)
    plot(freq, freq_signals{i}, 'LineWidth', 1.5, 'Color', colors{i}, ...
         'Marker', '.', 'MarkerSize', 15)
    title(titles_freq{i})
    ylabel('Normalized Power', 'FontWeight', 'bold')
    if i == 1
        xlabel('Frequency (Hz)', 'FontWeight', 'bold')
    end
    grid on
    xlim([0 30])
    ylim([0 fft_max])
    text(-6.5, fft_max*1.1, ['(' letters(idx_freq) ')'], ...
         'FontSize', 10, 'FontWeight', 'bold')
    
    % Time domain plots (right column)
    subplot(4,2,idx_time)
    plot(t, time_signals{i}, 'LineWidth', 1.5, 'Color', colors{i}, ...
         'Marker', '.', 'MarkerSize', 15)
    title(titles_time{i})
    ylabel('Amplitude', 'FontWeight', 'bold')
    if i ==1
        xlabel('Time (s)', 'FontWeight', 'bold')
    end
    grid on
    xlim([0 0.5])
    ylim([-y_max y_max])
    text(-0.085, y_max*1.2, ['(' letters(idx_time) ')'], ...
         'FontSize', 10, 'FontWeight', 'bold')
end

% Adjust subplot spacing
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
pos(3:4) = [width_inches height_inches];
set(gcf, 'Position', pos);

% Add space between subplots
h = findobj(gcf, 'Type', 'axes');
set(h, 'Units', 'normalized');
for i = 1:length(h)
    pos = get(h(i), 'Position');
    if mod(i,2) == 1  % Left column (frequency domain)
        pos(1) = 0.1;  % Left edge
        pos(3) = 0.35; % Width
    else  % Right column (time domain)
        pos(1) = 0.6; % Left edge
        pos(3) = 0.35; % Width
    end
    % Adjust vertical spacing
    row = ceil(i/2);
    pos(2) = 0.08 + (4-row)*0.24;  % Increased vertical spacing
    pos(4) = 0.16;  % Height
    set(h(i), 'Position', pos);
end

% Style adjustments
ax = findall(gcf, 'type', 'axes');
for i = 1:length(ax)
    set(ax(i), 'LineWidth', 1.5);
    set(ax(i), 'Color', 'white');
    set(ax(i), 'XColor', 'black', 'YColor', 'black');
end

% White background
set(gcf, 'Color', 'white')

% Save figure
set(gcf, 'PaperPositionMode', 'auto');
% print('-dpng', '-r1000', 'FT.png')
% saveas(gcf, 'FT.pdf')