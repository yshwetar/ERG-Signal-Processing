clear all; close all; clc;

% Load and prepare PERG data
load('0001.mat');  % Load the MAT file from current directory
signal = s0001.RE_1;  % Get the RE_1 column from s0001 table
signal = [signal; signal(end)]; % Pad to 256 samples
fs = 1700; % Sampling frequency

% Preprocess signal
signal = detrend(signal, 'linear');
[b, a] = butter(4, 100 / (fs / 2), 'low');
signal = filtfilt(b, a, signal);

% Create figure with adjusted dimensions for paper
width_inches = 8.5;  % Standard paper width
height_inches = 11;  % Standard paper height
figure('Units', 'inches', 'Position', [1 1 width_inches height_inches]);

% Define analysis parameters
window_size = 64; % Fixed window size for left column
overlap_ratio = 0.5; % 50% overlap
window_types = {@rectwin, @hamming, @blackmanharris};
window_names = {'Rectangular', 'Hamming', 'Blackman-Harris'};

% Define subplot indices map
row_indices = [1 2 3];  % Top to bottom
subplot_indices = [row_indices*2; row_indices*2-1]';  % Creates pairs: [1 2; 3 4; 5 6]
letters = 'a':'f';

% Store parameters for right column
window_sizes = [32, 64, 128];
overlap_ratios = [0.25, 0.50, 0.75];

% Plot all subplots
for i = 3:-1:1
    % Get subplot indices for this row
    idx_left = subplot_indices(i, 1);  % Left column index
    idx_right = subplot_indices(i, 2); % Right column index
    
    % Left column: Hanning window with different sizes
    subplot(3, 2, idx_left)
    window = hanning(window_sizes(i));
    overlap = round(window_sizes(i) * overlap_ratios(i));
    [s, f, t] = spectrogram(signal, window, overlap, [], fs);
    
    % Plot and style spectrogram
    plot_spectrogram(s, f, t, fs, ...
        sprintf('Hanning Window\nSize: %d, Overlap: %.0f%%', ...
        window_sizes(i), overlap_ratios(i) * 100), ...
        true, letters(idx_left));
    
    % Right column: Different window types
    subplot(3, 2, idx_right)
    window = window_types{i}(window_size);
    overlap = round(window_size * overlap_ratio);
    [s, f, t] = spectrogram(signal, window, overlap, [], fs);
    
    % Plot and style spectrogram
    plot_spectrogram(s, f, t, fs, ...
        sprintf('%s Window\nSize: %d, Overlap: %.0f%%', ...
        window_names{i}, window_size, overlap_ratio * 100), ...
        false, letters(idx_right));
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
    if mod(i, 2) == 1  % Left column
        pos(1) = 0.10;  % Closer to center
        pos(3) = 0.35;  % Width
    else  % Right column
        pos(1) = 0.55;  % Reduced horizontal distance
        pos(3) = 0.35;  % Width
    end
    % Adjust vertical spacing
    row = ceil(i / 2);
    pos(2) = 0.08 + (3 - row) * 0.30;  % Vertical adjustment
    pos(4) = 0.20;  % Height
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

% Set paper size to match figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width_inches height_inches]);
set(gcf, 'PaperPosition', [0 0 width_inches height_inches]);

% % Save figure
% print('-dpng', '-r1000', 'Windowing.png')
% saveas(gcf, 'Windowing.pdf')

% Helper function for spectrogram plotting
function plot_spectrogram(s, f, t, fs, title_text, show_colorbar, letter)
    % Calculate power and normalize individually
    power_spectrogram = abs(s).^2;
    freq_range = (f >= 0 & f <= 100);
    % Normalize to the maximum of this specific spectrogram
    normalized_spectrogram = power_spectrogram(freq_range, :) / max(power_spectrogram(:));
    
    % Plot spectrogram
    imagesc(t, f(freq_range), normalized_spectrogram);
    axis xy
    colormap('jet')
    
    % Labels and title
    xlabel('\bf{Time (s)}')
    ylabel('\bf{Frequency (Hz)}')
    title(title_text)
    
    % Axes settings
    ylim([0 100])
    set(gca, 'YDir', 'reverse')
    
    % Colorbar with linewidth
    if show_colorbar
        c = colorbar;
        % c.Label.String = '\bf{Normalized Power}';
        c.LineWidth = 1.5;  % Added linewidth to colorbar
        % set(c, 'FontWeight', 'bold');  % Make colorbar text bold
    end
    caxis([0 1])
    
    % Add subplot letter in the top-left corner with custom positions for (f) and (b)
    if letter == 'f'
        text(min(t) - 0.025, -18, ...
             sprintf('(%s)', letter), ...
             'FontSize', 10, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    elseif letter == 'b'
        text(min(t) - 0.0305, -18, ...
             sprintf('(%s)', letter), ...
             'FontSize', 10, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    else
        text(min(t) - 0.0315, -18, ...
             sprintf('(%s)', letter), ...
             'FontSize', 10, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end
end