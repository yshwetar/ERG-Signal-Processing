% ERG_WT_Analysis.m
% -----------------------------------------------------------
% Visualizes one or two ERG traces and builds numeric matrices
% Written by Yousif Shwetar (2025) | MIT License
% Modified to support multiple file formats (.xlsx, .mat, .csv, .asc)
% -----------------------------------------------------------

%% 1) RESET
close all; clearvars; clc;

%% 1.5) FILE LOADING
% Initialize variables
loadedSignal = [];
Fs = []; % Initialize Fs to empty

% First file selection
[fileName, filePath] = uigetfile({'*.xlsx;*.mat;*.csv;*.asc', 'All Supported Files (*.xlsx, *.mat, *.csv, *.asc)';...
                                  '*.xlsx', 'Excel Files (*.xlsx)';...
                                  '*.mat', 'MATLAB Files (*.mat)';...
                                  '*.csv', 'CSV Files (*.csv)';...
                                  '*.asc', 'ASCII Files (*.asc)'},...
                                  'Select first data file', 'MultiSelect', 'off');

% Process selection
if isequal(fileName, 0) || isequal(filePath, 0)
    % No file selected, will use demo data
    disp('No file selected. Using demo data.');
else
    loadedSignal = cell(1,2);
    sigLabel = cell(1,2);
    
    % Process first file
    fullPath = fullfile(filePath, fileName);
    [~, name, ext] = fileparts(fullPath);
    sigLabel{1} = name;
    
    switch lower(ext)
        case '.xlsx'
            % Load Excel file
            data = readmatrix(fullPath);
            % Assume first column is time if it exists and starts with 0 or small value
            if size(data, 2) > 1 && (data(1,1) == 0 || data(1,1) < 0.001)
                time = data(:,1);
                loadedSignal{1} = data(:,2);
                % Estimate Fs from time column
                Fs = 1/mean(diff(time));
            else
                loadedSignal{1} = data(:,1);
            end
            
        case '.mat'
            % Load MATLAB file
            matData = load(fullPath);
            fields = fieldnames(matData);
            
            % Check for common variable names for signals and sampling rate
            if any(strcmp(fields, 'signal'))
                loadedSignal{1} = matData.signal;
            elseif any(strcmp(fields, 'data'))
                loadedSignal{1} = matData.data;
            elseif any(strcmp(fields, 'y'))
                loadedSignal{1} = matData.y;
            elseif length(fields) >= 1
                % Just take the first variable if no common names found
                loadedSignal{1} = matData.(fields{1});
            end
            
            % Look for sampling rate in the MAT file
            if any(strcmp(fields, 'Fs'))
                Fs = matData.Fs;
            elseif any(strcmp(fields, 'fs'))
                Fs = matData.fs;
            elseif any(strcmp(fields, 'sampleRate'))
                Fs = matData.sampleRate;
            end
            
        case '.csv'
            % Load CSV file
            data = readmatrix(fullPath);
            % Assume first column is time if it exists and starts with 0 or small value
            if size(data, 2) > 1 && (data(1,1) == 0 || data(1,1) < 0.001)
                time = data(:,1);
                loadedSignal{1} = data(:,2);
                % Estimate Fs from time column
                Fs = 1/mean(diff(time));
            else
                loadedSignal{1} = data(:,1);
            end
            
        case '.asc'
            % Load ASCII file
            data = readmatrix(fullPath);
            % Assume first column is time if it exists and starts with 0 or small value
            if size(data, 2) > 1 && (data(1,1) == 0 || data(1,1) < 0.001)
                time = data(:,1);
                loadedSignal{1} = data(:,2);
                % Estimate Fs from time column
                Fs = 1/mean(diff(time));
            else
                loadedSignal{1} = data(:,1);
            end
    end
    
    % Ask if user wants to load a second file for comparison
    choice = questdlg('Would you like to load a second file to compare?', ...
                      'Load Second Signal', ...
                      'Yes', 'No', 'Yes');
    
    if strcmp(choice, 'Yes')
        % Second file selection
        [fileName2, filePath2] = uigetfile({'*.xlsx;*.mat;*.csv;*.asc', 'All Supported Files (*.xlsx, *.mat, *.csv, *.asc)';...
                                      '*.xlsx', 'Excel Files (*.xlsx)';...
                                      '*.mat', 'MATLAB Files (*.mat)';...
                                      '*.csv', 'CSV Files (*.csv)';...
                                      '*.asc', 'ASCII Files (*.asc)'},...
                                      'Select second data file', 'MultiSelect', 'off');
        
        if ~isequal(fileName2, 0) && ~isequal(filePath2, 0)
            % Process second file
            fullPath2 = fullfile(filePath2, fileName2);
            [~, name2, ext2] = fileparts(fullPath2);
            sigLabel{2} = name2;
            
            switch lower(ext2)
                case '.xlsx'
                    % Load Excel file
                    data = readmatrix(fullPath2);
                    % Assume first column is time if it exists and starts with 0 or small value
                    if size(data, 2) > 1 && (data(1,1) == 0 || data(1,1) < 0.001)
                        loadedSignal{2} = data(:,2);
                    else
                        loadedSignal{2} = data(:,1);
                    end
                    
                case '.mat'
                    % Load MATLAB file
                    matData = load(fullPath2);
                    fields = fieldnames(matData);
                    
                    % Check for common variable names for signals
                    if any(strcmp(fields, 'signal'))
                        loadedSignal{2} = matData.signal;
                    elseif any(strcmp(fields, 'data'))
                        loadedSignal{2} = matData.data;
                    elseif any(strcmp(fields, 'y'))
                        loadedSignal{2} = matData.y;
                    elseif length(fields) >= 1
                        % Just take the first variable if no common names found
                        loadedSignal{2} = matData.(fields{1});
                    end
                    
                case '.csv'
                    % Load CSV file
                    data = readmatrix(fullPath2);
                    % Assume first column is time if it exists and starts with 0 or small value
                    if size(data, 2) > 1 && (data(1,1) == 0 || data(1,1) < 0.001)
                        loadedSignal{2} = data(:,2);
                    else
                        loadedSignal{2} = data(:,1);
                    end
                    
                case '.asc'
                    % Load ASCII file
                    data = readmatrix(fullPath2);
                    % Assume first column is time if it exists and starts with 0 or small value
                    if size(data, 2) > 1 && (data(1,1) == 0 || data(1,1) < 0.001)
                        loadedSignal{2} = data(:,2);
                    else
                        loadedSignal{2} = data(:,1);
                    end
            end
        else
            % User canceled the second file selection
            loadedSignal = loadedSignal{1}; % Only use the first signal
        end
    else
        % User chose not to load a second file
        loadedSignal = loadedSignal{1}; % Only use the first signal
    end
end

%% 2) USER INPUTS
signal = loadedSignal;   % Use loaded signal or empty for demo
if isempty(signal)
    % If no file was loaded or the signal was not properly extracted, use demo data
    Fs = 1e3;  t = 0:1/Fs:0.2;
    demoA = 90*sin(2*pi*30*t)+35*sin(2*pi*80*t+.9)+8*randn(size(t));
    demoB = 40*sin(2*pi*20*t)+55*sin(2*pi*70*t   )+8*randn(size(t));
    signal   = {demoA , demoB};
    sigLabel = {"Demo A","Demo B"};
else
    % Make sure we have a valid sampling rate
    if isempty(Fs)
        % Ask user for sampling rate if not set during file loading
        answer = inputdlg('Sampling rate (Hz) not found. Please enter manually:', 'Sampling Rate', 1, {'1000'});
        if ~isempty(answer)
            Fs = str2double(answer{1});
        else
            Fs = 1000; % Default to 1kHz if user cancels
            disp('Using default sampling rate of 1000 Hz.');
        end
    end
    
    % If signals have different lengths, ask user how to handle it
    if iscell(signal) && length(signal) == 2 && ~isempty(signal{1}) && ~isempty(signal{2})
        len1 = length(signal{1});
        len2 = length(signal{2});
        
        if len1 ~= len2
            choice = questdlg(sprintf('Signal lengths differ (%d vs %d samples). How would you like to handle this?', len1, len2), ...
                             'Signal Length Mismatch', ...
                             'Truncate to shorter', 'Pad shorter with zeros', 'Truncate to shorter');
            
            if strcmp(choice, 'Truncate to shorter')
                minLen = min(len1, len2);
                signal{1} = signal{1}(1:minLen);
                signal{2} = signal{2}(1:minLen);
                disp(['Signals truncated to ', num2str(minLen), ' samples.']);
            elseif strcmp(choice, 'Pad shorter with zeros')
                maxLen = max(len1, len2);
                if len1 < maxLen
                    signal{1} = [signal{1}, zeros(1, maxLen - len1)];
                else
                    signal{2} = [signal{2}, zeros(1, maxLen - len2)];
                end
                disp(['Shorter signal padded with zeros to ', num2str(maxLen), ' samples.']);
            end
        end
    end
    
    % Ensure signal is in cell format for consistent processing
    if ~iscell(signal)
        signal = {signal};
    end
    
    % Make sure sigLabel is set properly
    if ~exist('sigLabel', 'var') || isempty(sigLabel) || length(sigLabel) ~= length(signal)
        sigLabel = cell(1, length(signal));
        for i = 1:length(signal)
            sigLabel{i} = ['Signal ', char('A' + i - 1)];
        end
    end
end

cwtWavelet = "amor";
dwtWavelet = "haar";

%% 3) PROCESS + COLLECT MATRICES
if ~iscell(signal); signal = {signal}; end
nSig = numel(signal);
N    = numel(signal{1});
for k=1:nSig; signal{k} = detrend(signal{k}(:).'); end
tVec = (0:N-1)/Fs;

AmpMat = cell(1,nSig);   % { [f , amplitude] }
CWTMat = cell(1,nSig);   % full matrix with f (col-1) & t (last-row)
DWTMat = cell(1,nSig);   % table with labeled columns

[~,fPSD] = periodogram(signal{1},[],2^nextpow2(N),Fs);
maxAmp = 0;  cwtMax = 0;

for k=1:nSig
    %% 3-A) Amplitude spectrum
    Px          = periodogram(signal{k},[],2^nextpow2(N),Fs);
    Ax          = sqrt(Px);
    maxAmp      = max(maxAmp,max(Ax));
    AmpMat{k} = array2table([fPSD(:) Ax(:)],'VariableNames',{'Frequency','Amplitude'});
    
    %% 3-B) CWT  (store f in first col, t in last row)
    [wt,fCWT]   = cwt(signal{k},cwtWavelet,Fs);     % size: nF × nT
    wtAbs       = abs(wt);
    cwtMax      = max(cwtMax,max(wtAbs(:)));
    
    % Pad matrix: (nF+1) × (nT+1)
    cwtFull          = zeros(numel(fCWT)+1 , numel(tVec)+1);
    cwtFull(1:end-1,1)  = fCWT(:);           % frequency column
    cwtFull(1:end-1,2:end) = wtAbs;          % energies
    cwtFull(end,2:end)    = tVec;            % time row
    CWTMat{k} = cwtFull;                     % raw; normalize later
    
    %% 3-C) DWT  (store as table with headers)
    maxLvl      = wmaxlev(N,dwtWavelet);
    [C,L]       = wavedec(signal{k},maxLvl,dwtWavelet);
    
    rows = [];
    idx  = 1;
    for lvl = 1:maxLvl
        dCoef  = detcoef(C,L,lvl);
        segLen = N / numel(dCoef);
        fStart = Fs/(2^(lvl+1));
        fEnd   = Fs/(2^lvl);
        dCoef  = detcoef(C,L,lvl);
        for n = 1:numel(dCoef)
            tStart = (n-1)*segLen/Fs;
            tEnd   = n*segLen/Fs;
            rows(idx,:) = [tStart tEnd fStart fEnd dCoef(n).^2];
            idx = idx+1;
        end
    end
    aCoef  = appcoef(C,L,dwtWavelet,maxLvl);
    segLen = N / numel(aCoef);
    for n = 1:numel(aCoef)
        tStart = (n-1)*segLen/Fs;
        tEnd   = n*segLen/Fs;
        rows(idx,:) = [tStart tEnd 0 Fs/(2^(maxLvl+1)) aCoef(n).^2];
        idx = idx+1;
    end
    
    % Convert to table with labeled headers
    DWTMat{k} = array2table(rows,...
        'VariableNames',{'TimeStart','TimeEnd','FreqStart','FreqEnd','Energy'});
end

% normalize CWT now that global max is known
for k=1:nSig
    C      = CWTMat{k};
    C(1:end-1,2:end) = C(1:end-1,2:end)/cwtMax;
    CWTMat{k} = C;
end

%% 4) PLOTTING
nSig = numel(signal);
assert(nSig<=2,"Provide at most two signals.");

N = numel(signal{1});
for k = 1:nSig
    assert(numel(signal{k})==N,"All signals must be equal length.");
    signal{k} = detrend(signal{k}(:).');
end
tVec = (0:N-1)/Fs;

% Time-domain limits
yMin = min(cellfun(@min,signal));
yMax = max(cellfun(@max,signal));
pad  = 0.05*(yMax-yMin);  yMin = yMin-pad;  yMax = yMax+pad;

% Amplitude spectra
Ax  = cell(1,nSig);  maxAmp = 0;
[~,fPSD] = periodogram(signal{1},[],2^nextpow2(N),Fs);
for k = 1:nSig
    [Px,~] = periodogram(signal{k},[],2^nextpow2(N),Fs);
    Ax{k}  = sqrt(Px);
    maxAmp = max(maxAmp,max(Ax{k}));
end

% Continuous Wavelet Transform
wtNorm = cell(1,nSig);  cwtMax = 0;  wtAbs = cell(1,nSig);
for k = 1:nSig
    [wt,fCWT]  = cwt(signal{k},cwtWavelet,Fs);
    wtAbs{k}   = abs(wt);
    cwtMax     = max(cwtMax,max(wtAbs{k}(:)));
end
for k = 1:nSig;  wtNorm{k} = wtAbs{k}/cwtMax;  end

% Discrete Wavelet Transform (maximum level)
maxLvl  = wmaxlev(N,dwtWavelet);
ctrFreq = round((3*Fs)./2.^((1:maxLvl)'+2));
ctrFreq = [ctrFreq ; round(Fs/2^(maxLvl+1))];

dwtNorm = cell(1,nSig);  dwtMax = 0;  dwtCell = cell(1,nSig);
for k = 1:nSig
    [C,L] = wavedec(signal{k},maxLvl,dwtWavelet);
    dwtHM = zeros(maxLvl+1,N);
    for lvl = 1:maxLvl
        e   = detcoef(C,L,lvl).^2;
        rep = repelem(e,ceil(N/numel(e)));
        dwtHM(lvl,:) = rep(1:N);
    end
    aCoeff       = appcoef(C,L,dwtWavelet,maxLvl).^2;
    rep          = repelem(aCoeff,ceil(N/numel(aCoeff)));
    dwtHM(end,:) = rep(1:N);
    dwtCell{k}   = dwtHM;
    dwtMax       = max(dwtMax,max(dwtHM(:)));
end
for k = 1:nSig;  dwtNorm{k} = dwtCell{k}/dwtMax;  end

% Figure
baseFS = 10;
figH   = figure('Color','w',...
                'Position',[50 50 1600 250+200*(nSig-1)]);

for k = 1:nSig
    col = (k-1)*4;
    
    % Time-domain
    subplot(nSig,4,col+1);
    plot(tVec,signal{k},'k','LineWidth',1.5);
    title("Time — " + sigLabel{k},'FontSize',11,'FontWeight','bold');
    xlabel('Time (s)'); ylabel('Amplitude (\muV)');
    xlim([tVec(1) tVec(end)]); ylim([yMin yMax]);
    set(gca,'FontSize',baseFS); box on;
    
    % Amplitude Spectrum
    subplot(nSig,4,col+2);
    stem(fPSD,Ax{k},'Marker','none','LineWidth',1.5,...
         'Color',[0 0.4470 0.7410]);
    title("Amplitude Spectrum — " + sigLabel{k},...
          'FontSize',11,'FontWeight','bold');
    xlabel('Frequency (Hz)'); ylabel('Amplitude (\muV/√Hz)');
    xlim([0 Fs/2]); ylim([0 1.05*maxAmp]);
    set(gca,'FontSize',baseFS); box on;
    
    % CWT
    subplot(nSig,4,col+3);
    imagesc(tVec,fCWT,wtNorm{k},[0 1]); axis xy;
    ax = gca;  ax.YDir='reverse';  ylim([min(fCWT) max(fCWT)]);
    colormap('jet');
    title("CWT (" + cwtWavelet + ") — " + sigLabel{k},...
          'FontSize',11,'FontWeight','bold');
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    set(ax,'FontSize',baseFS); box on;
    
    % DWT
    hAx = subplot(nSig,4,col+4);
    imagesc(tVec,1:(maxLvl+1),dwtNorm{k},[0 1]); axis xy;
    colormap('jet');
    title("DWT (" + dwtWavelet + ") — " + sigLabel{k},...
          'FontSize',11,'FontWeight','bold');
    xlabel('Time (s)'); ylabel('Center Freq (Hz)');
    set(gca,'YTick',1:(maxLvl+1),'YTickLabel',ctrFreq,...
        'FontSize',baseFS,'YAxisLocation','left'); box on;
    
    if k==nSig
        axPos = get(hAx,'Position');
        cb    = colorbar(hAx,'eastoutside');
        cb.Ticks      = [0];
        cb.TickLabels = {'0'};
        cb.FontSize   = baseFS;
        cbWidth = 0.007;  gap = 0.01;
        set(cb,'Units','normalized');
        cbPos = get(cb,'Position');
        cbPos(1) = axPos(1) + axPos(3) + gap;
        cbPos(3) = cbWidth;  cbPos(4) = axPos(4);
        set(cb,'Position',cbPos);
        set(hAx,'Position',axPos);
    end
end

%% 5) DISPLAY LOADED FILE INFORMATION
fprintf('\n==== Analysis Summary ====\n');
fprintf('Number of signals: %d\n', nSig);

% Display information for each signal separately
for k = 1:nSig
    fprintf('\n--- %s ---\n', sigLabel{k});
    N_signal = numel(signal{k});
    tVec_signal = (0:N_signal-1)/Fs;
    
    fprintf('Signal length: %d samples\n', N_signal);
    fprintf('Duration: %.4f seconds\n', tVec_signal(end));
    fprintf('Min amplitude: %.4f\n', min(signal{k}));
    fprintf('Max amplitude: %.4f\n', max(signal{k}));
    
    % Calculate basic frequency statistics from amplitude spectrum
    [maxAmp, maxIdx] = max(Ax{k});
    dominantFreq = fPSD(maxIdx);
    fprintf('Dominant frequency: %.2f Hz (amplitude: %.4f)\n', dominantFreq, maxAmp);
end

fprintf('\nSampling rate: %.2f Hz\n', Fs);
fprintf('========================\n\n');

% Optional: save figure
% print(figH,'ERG_scalogram_plot','-dpng','-r1000');