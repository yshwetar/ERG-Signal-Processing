% Clear workspace, close all figures, and clear command window
clear all; close all; clc;

% Load data from .mat file
load('0001.mat');
signal = s0001.RE_1;

% Signal padding
signal = [signal; signal(end)];

% Create time vector
fs = 1700; % Sampling frequency
t = (0:length(signal)-1) / fs;

% Plot original signal first
fig_orig = figure('Color', 'white', 'Position', [100 100 1000 800]);
plot(t, signal, 'k-', 'LineWidth', 1);
set(gca, 'Visible', 'off');

% Perform DWT using Haar wavelet
level = 7;
[c, l] = wavedec(signal, level, 'haar');

% Extract and process coefficients for detailed visualization
approx = {};
detail = {};
for i = 1:level
   % Get coefficients
   d = detcoef(c, l, i);
   a = appcoef(c, l, 'haar', i);
   
   % Create zero vectors
   detail{i} = zeros(size(signal));
   approx{i} = zeros(size(signal));
   
   % Place coefficients at correct positions
   step = 2^i;
   for j = 1:length(d)
       detail{i}((j-1)*step + 1:j*step) = d(j);
   end
   for j = 1:length(a)
       approx{i}((j-1)*step + 1:j*step) = a(j);
   end
end

% Scale coefficients
scaled_detail = {};
scaled_detail_sq = {};
scaled_approx = {};
scaled_approx_sq = {};
for i = 1:level
   % Scale detail coefficients
   d_max = max(abs(detail{i}));
   if d_max ~= 0
       scaled_detail{i} = detail{i} / d_max;
   else
       scaled_detail{i} = detail{i};
   end
   d_sq = detail{i}.^2;
   d_sq_max = max(d_sq);
   if d_sq_max ~= 0
       scaled_detail_sq{i} = d_sq / d_sq_max;
   else
       scaled_detail_sq{i} = d_sq;
   end
   
   % Scale approximation coefficients
   a_max = max(abs(approx{i}));
   if a_max ~= 0
       scaled_approx{i} = approx{i} / a_max;
   else
       scaled_approx{i} = approx{i};
   end
   a_sq = approx{i}.^2;
   a_sq_max = max(a_sq);
   if a_sq_max ~= 0
       scaled_approx_sq{i} = a_sq / a_sq_max;
   else
       scaled_approx_sq{i} = a_sq;
   end
end

% Plot scaled detail coefficients
fig_d = figure('Color', 'white', 'Position', [100 100 1000 800]);
subplot(8, 2, 1);
text(0.5, 0.5, 'Detail Coefficients', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
axis off;
subplot(8, 2, 2);
text(0.5, 0.5, 'Detail Energy', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
axis off;

for i = 1:level
    % Detail coefficients
    subplot(8, 2, 2*i+1)
    hold on
    plot(t, zeros(size(t)), 'k-', 'LineWidth', 0.2, 'Color', [0.7 0.7 0.7]);
    plot(t, scaled_detail{i}, 'k-', 'LineWidth', 1);
    ylim([-1 1])
    text(-0.02*max(t), 0, sprintf('D%d', i), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'Visible', 'off');
    hold off
    
    % Detail energy
    subplot(8, 2, 2*i+2)
    hold on
    plot(t, zeros(size(t)), 'k-', 'LineWidth', 0.2, 'Color', [0.7 0.7 0.7]);
    plot(t, scaled_detail_sq{i}, 'k-', 'LineWidth', 1);
    ylim([0 1])
    set(gca, 'Visible', 'off');
    hold off
end

% Plot approximation coefficients
fig_a = figure('Color', 'white', 'Position', [100 100 1000 800]);
subplot(8, 2, 1);
text(0.5, 0.5, 'Approximation Coefficients', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
axis off;
subplot(8, 2, 2);
text(0.5, 0.5, 'Approximation Energy', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
axis off;

for i = 1:level
    % Approximation coefficients
    subplot(8, 2, 2*i+1)
    hold on
    plot(t, zeros(size(t)), 'k-', 'LineWidth', 0.2, 'Color', [0.7 0.7 0.7]);
    plot(t, scaled_approx{i}, 'k-', 'LineWidth', 1);
    ylim([-1 1])
    text(-0.02*max(t), 0, sprintf('A%d', i), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'Visible', 'off');
    hold off
    
    % Approximation energy
    subplot(8, 2, 2*i+2)
    hold on
    plot(t, zeros(size(t)), 'k-', 'LineWidth', 0.2, 'Color', [0.7 0.7 0.7]);
    plot(t, scaled_approx_sq{i}, 'k-', 'LineWidth', 1);
    ylim([0 1])
    set(gca, 'Visible', 'off');
    hold off
end

% Create heatmap visualization
% Initialize matrices for detail coefficients
heatmap_matrix_detail = zeros(level, length(signal));
approx_coeffs = appcoef(c, l, 'haar', level);

% Fill matrix with detail coefficients
for i = 1:level
    detail_coeffs = detcoef(c, l, i);
    detail_energy = detail_coeffs .^ 2;
    heatmap_matrix_detail(i, :) = repelem(detail_energy, ceil(length(signal) / length(detail_energy)));
end

% Normalize detail coefficients matrix
heatmap_matrix_detail = (heatmap_matrix_detail - min(heatmap_matrix_detail(:))) / ...
    (max(heatmap_matrix_detail(:)) - min(heatmap_matrix_detail(:)));

% Process approximation coefficients
approx_energy = approx_coeffs .^ 2;
approx_row = zeros(1, length(signal));
half_length = length(signal)/2;
approx_row(1:half_length) = approx_energy(1);
approx_row(half_length+1:end) = approx_energy(2);
approx_row_norm = (approx_row - min(approx_row)) / (max(approx_row) - min(approx_row));

% Combine matrices
heatmap_matrix_combined = [heatmap_matrix_detail; approx_row_norm];

% Create heatmap figure
figure('Units', 'inches', 'Position', [0.1, 0.1, 6.5, 4], 'Color', 'w');
ax = axes('Position', [0.15 0.15 0.7 0.7]);

% Plot heatmap
imagesc(linspace(0, length(signal)/fs, size(heatmap_matrix_combined, 2)), ...
    1:(level+1), heatmap_matrix_combined);

% Style the heatmap
cb = colorbar;
caxis([0 1]);
set(cb, 'LineWidth', 1.5);
set(cb, 'FontWeight', 'bold');
set(cb, 'Ticks', 0:0.2:1);
set(gca, 'FontWeight', 'bold', 'LineWidth', 1.5);
xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
yticks(1:(level+1));
yticklabels({'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'A7'});
ylabel('Decomposition Level', 'FontWeight', 'bold', 'FontSize', 10);
set(gca, 'YDir', 'normal');
box on;
colormap('jet');

% Save all figures if needed
% print('original_signal', '-dpng', '-r1000');
% print('detail_coefficients', '-dpng', '-r1000');
% print('approximation_coefficients', '-dpng', '-r1000');
% print('dwt_heatmap', '-dpng', '-r1000');