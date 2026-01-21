%% ============================================================
%  fusi_video_soner9.m
%  fUSI VIDEO ANALYSIS GUI (MATLAB 2017b)
% by Soner Caner Cagun - MPI for Biological Cybernetics - 2026
% Adapted from Literature and AUCT 
% ============================================================

clearvars;
close all;
clc;

%% ===================== USER PARAMETERS =====================
% NOTE:
% TR and TotalTimeSec are now AUTO-DETECTED from the raw .mat file if present.
% The values below are ONLY FALLBACKS and will be used ONLY if:
%   - TR / time vector / Fs are missing in the file.
% You may keep them, but you do NOT need to edit them manually anymore.

TR = 0.3;                 % [sec] fallback only (ignored if found in MAT)
TotalTimeSec = [];      % [] = auto-computed from data length × TR


% ------------------------------------------------------------
% VIDEO PLAYBACK SPEED
% ------------------------------------------------------------
initialFPS = 10;   % default playback speed (frames per second)
maxFPS     = 120;  % maximum allowed FPS via slider


% ------------------------------------------------------------
% TEMPORAL INTERPOLATION
% ------------------------------------------------------------
% Inserts intermediate frames between volumes for smoother playback.
% Does NOT change the data or temporal resolution.
%
% Higher value:
%   - smoother video
%   - slower rendering
%   - more memory
%
% Lower value:
%   - faster playback
%   - more discrete jumps
%
% Typical values: 4–8
par.interpol = 8;

% ------------------------------------------------------------
% TEMPORAL FILTERING
% ------------------------------------------------------------

% Low-pass filter cutoff (Hz)
% Removes fast temporal fluctuations (noise, motion, physiology).
% Keeps slow hemodynamic responses.
%
% PACAP / pharmacology-friendly values:
%   0.2–0.4 Hz
%
% Higher values = noisier
% Lower values = smoother but delayed responses
par.LPF = 0.95; %must be between 0-1 : Higher value more LPF; lower value less/none 

% High-pass filter cutoff (Hz)
% Removes slow drift (baseline instability).
%
% Usually OFF for fUSI (recommended).
% Turn ON only if you see strong baseline drift.
par.HPF = 0;

% ------------------------------------------------------------
% SPATIAL FILTERING
% ------------------------------------------------------------

% Gaussian smoothing kernel size (pixels)
% Defines spatial extent of smoothing.
% Must be an ODD number for symmetry.
par.gaussSize = 7;

% Gaussian sigma (pixels)
% Controls smoothing strength:
%   small sigma → sharp vessels
%   large sigma → smoother but blurrier
par.gaussSig  = 0.2;

% ------------------------------------------------------------
% CONNECTIVITY FILTER (VERY IMPORTANT)
% ------------------------------------------------------------
% Suppresses isolated noisy pixels and keeps spatially coherent clusters.
%
% How it works:
%   1) Threshold PSC at conectLev (%)
%   2) Check local neighborhood (disk of radius conectSize)
%   3) Penalize isolated pixels
%
% Increasing conectSize:
%   → stricter spatial coherence
%
% Increasing conectLev:
%   → only strong vessels survive

par.conectSize = 5;     % neighborhood radius (pixels)
par.conectLev  = 20;    % PSC threshold (%) for connectivity

% ------------------------------------------------------------
% DISPLAY (VISUALIZATION ONLY)
% ------------------------------------------------------------

% Color scale range for PSC visualization (%)
% Does NOT affect data or masking.
par.previewCaxis = [0 100];

% ------------------------------------------------------------
% BASELINE WINDOW
% ------------------------------------------------------------
% Used for PSC computation:
% PSC = (I - baseline) / baseline * 100
%
% Defined in SECONDS (not volumes)
baseline.mode  = 'sec';
baseline.start = 0;
baseline.end   = 300;


%% ===================== FILE SELECTION =====================
startPath = 'Z:\fUS\Project_PACAP_AVATAR_SC';

[file,path] = uigetfile( ...
    {'*.nii;*.mat','fUS data (*.nii, *.mat)'}, ...
    'Select fUS data file', startPath);

if isequal(file,0)
    error('No file selected.');
end

dataFile = fullfile(path,file);
[~,~,ext] = fileparts(file);
fprintf('Loading data: %s\n', dataFile);

%% ===================== LOAD DATA =====================
loadedMask = [];
loadedMaskIsInclude = true;
loadedPar = [];
loadedBaseline = [];
loadedTR = [];
loadedTotalTimeSec = [];

switch lower(ext)

  case '.mat'
    S = load(dataFile);

    if ~isfield(S,'I')
        error('MAT file must contain variable I.');
    end
    I = single(S.I);

    % ======================================================
    % AUTO-DETECT TR AND TOTAL ACQUISITION TIME (FAIL-SAFE)
    % ======================================================
    TR_found = false;
    TotalTimeSec_found = false;
    timeVec = [];

    % ---------- 1) Direct fields ----------
    if isfield(S,'TR') && isnumeric(S.TR) && isscalar(S.TR) && isfinite(S.TR) && S.TR>0
        TR = double(S.TR);
        TR_found = true;
    end

    if isfield(S,'TotalTimeSec') && isnumeric(S.TotalTimeSec) && ...
            isscalar(S.TotalTimeSec) && isfinite(S.TotalTimeSec) && S.TotalTimeSec>0
        TotalTimeSec = double(S.TotalTimeSec);
        TotalTimeSec_found = true;
    end

    % ---------- 2) Time vector inference ----------
    if isfield(S,'t') && isnumeric(S.t) && isvector(S.t)
        timeVec = double(S.t(:));
    elseif isfield(S,'time') && isnumeric(S.time) && isvector(S.time)
        timeVec = double(S.time(:));
    elseif isfield(S,'timestamps') && isnumeric(S.timestamps) && isvector(S.timestamps)
        timeVec = double(S.timestamps(:));
    end

    if ~isempty(timeVec)
        dt = diff(timeVec);
        dt = dt(isfinite(dt) & dt>0);
        if ~isempty(dt)
            if ~TR_found
                TR = median(dt);
                TR_found = true;
            end
            if ~TotalTimeSec_found
                TotalTimeSec = (timeVec(end) - timeVec(1)) + TR;
                TotalTimeSec_found = true;
            end
        end
    end

    % ---------- 3) Sampling rate ----------
    if ~TR_found && isfield(S,'Fs') && isnumeric(S.Fs) && isscalar(S.Fs) && isfinite(S.Fs) && S.Fs>0
        TR = 1 / double(S.Fs);
        TR_found = true;
    end

    % ---------- 4) Final fallbacks ----------
    if ~TR_found
        warning('TR not found in MAT file — using default TR = 0.3 s');
        TR = 0.3;
    end

    if ~TotalTimeSec_found
        TotalTimeSec = size(I,3) * TR;
    end

    % ---- restore parameters if present ----
    if isfield(S,'par'),      loadedPar = S.par; end
    if isfield(S,'baseline'), loadedBaseline = S.baseline; end
    if isfield(S,'TR'),       loadedTR = S.TR; end
    if isfield(S,'TotalTimeSec'), loadedTotalTimeSec = S.TotalTimeSec; end

    % ---- restore mask if present ----
    if isfield(S,'mask')
        loadedMask = logical(S.mask);
    end

    if isfield(S,'maskIsInclude')
        loadedMaskIsInclude = logical(S.maskIsInclude);
    end

  case '.nii'
    V = niftiread(dataFile);
    I = convertNiftiToI(V);

  otherwise
    error('Unsupported file type for this viewer.');
end

% Apply restored params ONLY if present
if ~isempty(loadedPar), par = loadedPar; end
if ~isempty(loadedBaseline), baseline = loadedBaseline; end
if ~isempty(loadedTR), TR = loadedTR; end
if ~isempty(loadedTotalTimeSec), TotalTimeSec = loadedTotalTimeSec; end


%% ===================== VOLUME COUNT (FIX) =====================
% Compute volumes from TotalTimeSec / TR, but be fail-safe:
% If loaded data length doesn't match requested length, prefer the data.
nVols_data = size(I,3);
nVols_req  = round(TotalTimeSec / TR);

if nVols_req <= 0 || isnan(nVols_req) || isinf(nVols_req)
    warning('Invalid TotalTimeSec/TR; falling back to data-derived volumes.');
    nVols = nVols_data;
    TotalTimeSec = nVols * TR;
else
    % If mismatch is large, prefer data to avoid index errors
    if abs(nVols_req - nVols_data) > 1
        warning([ ...
            'Requested volumes (round(TotalTimeSec/TR)=%d) do not match data volumes (%d). ' ...
            'Using DATA volumes to remain fail-safe. TotalTimeSec updated accordingly.'], ...
            nVols_req, nVols_data);
        nVols = nVols_data;
        TotalTimeSec = nVols * TR;
    else
        % Close enough: use requested
        nVols = nVols_req;
        % Fail-safe: trim/pad to requested length if needed
        if nVols_data > nVols
            I = I(:,:,1:nVols);
        elseif nVols_data < nVols
            % pad by repeating last volume
            I(:,:,end+1:nVols) = repmat(I(:,:,end), [1 1 (nVols-nVols_data)]);
        end
    end
end

Tmax = (nVols - 1) * TR;


%% ===================== FRAME-RATE / PNS QC (DIAGNOSTIC ONLY) =====================
QC = runFrameRateQC(I, TR);   % fast, no figures

set(0,'DefaultFigureColormap',parula);
%% ===================== PROCESSING =====================
I1 = interpFilm(I, par.interpol);
[nz,nx,nFrames] = size(I1);
nVols = size(I,3);
Tmax = (nVols - 1) * TR; 

b0 = round((baseline.start/TR)*par.interpol) + 1;
b1 = round((baseline.end/TR)*par.interpol);

% Clamp to valid frame range
b0 = max(1, b0);
b1 = min(nFrames, b1);

% ---------- FIX 2: BASELINE SAFETY ----------
if b0 >= b1 || b1 < 1 || b0 > nFrames
    warning('Baseline window invalid — using first 10%% of data instead.');
    b0 = 1;
    b1 = max(1, round(0.1 * nFrames));
end
% -------------------------------------------

ab = mean(I1(:,:,b0:b1), 3);

% Extra safety (rare but cheap)
ab(~isfinite(ab) | ab==0) = eps;

PSC = (I1 - ab) ./ ab * 100;

if par.LPF > 0
    [B,A] = butter(4,par.LPF,'low');
    PSC = filter(B,A,PSC,[],3);
end

% (kept from original)
if par.conectSize > 0
    h = fspecial('disk',par.conectSize);
    for k = 1:nFrames
        PSC(:,:,k) = PSC(:,:,k) .* filter2(h,PSC(:,:,k)>par.conectLev).^2;
    end
end

if par.gaussSize > 0
    hG = fspecial('gaussian',par.gaussSize,par.gaussSig);
    for k = 1:nFrames
        PSC(:,:,k) = filter2(hG,PSC(:,:,k));
    end
end

bg = mean(I,3);

m = max(bg(:));
if m <= 0 || ~isfinite(m)
    bg = zeros(size(bg));   % safe fallback
else
    bg = 20*log10(bg ./ m);
end



% Wait for QC figures to be closed (robust)
if isfield(QC,'figIntensity') && ishandle(QC.figIntensity)
    waitfor(QC.figIntensity);
end

if isfield(QC,'figRejected') && ishandle(QC.figRejected)
    waitfor(QC.figRejected);
end


%% ===================== OPTIONAL FRAME-RATE REJECTION =====================
applyRejection = false;

if QC.rejPct > 0
    choice = questdlg( ...
        sprintf(['Frame-rate QC detected %.1f%% unstable volumes.\n\n' ...
                 'What do you want to do?'], QC.rejPct), ...
        'Frame-rate QC', ...
        'Show QC plots','Interpolate (video only)','Ignore','Ignore');

    switch choice
        case 'Show QC plots'
            QC = runFrameRateQC(I, TR);
            applyRejection = false;

        case 'Interpolate (video only)'
            applyRejection = true;

        otherwise
            applyRejection = false;
    end
end

% ==========================================================
% >>> ADD THIS BLOCK RIGHT HERE <<<
% ==========================================================
if applyRejection
    fprintf('[INFO] %d / %d volumes interpolated for VIDEO DISPLAY ONLY.\n', ...
        nnz(QC.outliers), numel(QC.outliers));
else
    fprintf('[INFO] No frame-rate interpolation applied.\n');
end
% ==========================================================
% ==========================================================
% APPLY INTERPOLATION + INTERPOLATED QC (VIDEO ONLY)
% ==========================================================
if applyRejection

    % ---- Apply interpolation to DATA USED FOR VIDEO ----
    I_interp = interpolateRejectedVolumes(I, QC.outliers);

    % ---- Run QC AGAIN on interpolated data (save PNGs) ----
    QC_interp = runFrameRateQC(I_interp, TR, 'INTERPOLATED', true);

    fprintf('[INFO] Interpolated QC completed. PNGs saved.\n');

    % ---- Replace I ONLY for visualization / video ----
    I = I_interp;

end
% ==========================================================


%% ===================== GUI =====================
play_fusi_video_final( ...
    I, PSC, bg, par, initialFPS, maxFPS, TR, Tmax, baseline, ...
    loadedMask, loadedMaskIsInclude, nVols);


%% ===================== FUNCTIONS =====================
function I = convertNiftiToI(V)
    if ndims(V)==4
        V = squeeze(mean(V,3));
    end
    I = single(permute(V,[2 1 3]));
end

function I1 = interpFilm(I,N)
    [nz,nx,nt] = size(I);
    I1 = zeros(nz,nx,nt*N,'single');
    for ix=1:nx
        for iz=1:nz
            I1(iz,ix,:) = single(interp(double(squeeze(I(iz,ix,:))),N));
        end
    end
end

%% ============================================================
%  MAIN GUI
%% ============================================================
function play_fusi_video_final( ...
    I, PSC, bg, par, fps, maxFPS, TR, Tmax, baseline, ...
    loadedMask, loadedMaskIsInclude, nVols)


[nz,nx,nFrames] = size(PSC);

fig = figure('Color','k','Position',[60 60 1500 900], ...
    'Name','fUSI Video Analysis — Soner (Auto Mask v3)', ...
    'NumberTitle','off');

% --- Layout ---
ax  = axes('Parent',fig,'Units','normalized','Position',[0.14 0.12 0.56 0.70]);
axis(ax,'off','image');
img = image(ax,zeros(nz,nx,3,'single'));
set(ax,'HitTest','on');
set(img,'HitTest','off');

% --- Colormaps ---
Nc = 128;
mapA = hot(Nc); mapA(1,:)=0;
bgN = (bg+68)/(68-5);
bgN = max(0,min(1,bgN));
bgRGB = ind2rgb(round(bgN*127),gray(128));

% Info
info = uicontrol('Style','text','Units','normalized',...
    'Position',[0.01 0.88 0.74 0.09],...
    'ForegroundColor','w','BackgroundColor','k',...
    'FontName','Courier','FontSize',13,...
    'HorizontalAlignment','left');

% Right panel placement
figPos = get(fig,'Position');
rightX = round(figPos(3) * 0.72);

uiFontName = 'Helvetica';
uiFontSize = 11;

% Control Panel title
uicontrol('Style','text','Units','pixels',...
    'Position',[rightX 835 360 18],...
    'String','CONTROL PANEL',...
    'ForegroundColor',[0.9 0.9 0.9],...
    'BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',12,...
    'FontWeight','bold',...
    'HorizontalAlignment','left');
% ================= PLAYBACK CONTROLS =================

% ================= PLAYBACK CONTROLS =================

rowH = 28;
gap  = 10;

% Place controls directly under "CONTROL PANEL" title
y0 = 800;   % ← FPS row anchor (tweak ±10 if desired)

% -------- FPS (PRIMARY PLAYBACK CONTROL) --------
uicontrol('Style','text','Units','pixels', ...
    'Position',[rightX y0 100 rowH], ...
    'String','FPS', ...
    'ForegroundColor','w','BackgroundColor','k', ...
    'FontName',uiFontName,'FontSize',uiFontSize, ...
    'HorizontalAlignment','right');

fpsValue = uicontrol('Style','text','Units','pixels', ...
    'Position',[rightX+105 y0 80 rowH], ...
    'ForegroundColor','w','BackgroundColor','k', ...
    'FontName',uiFontName,'FontSize',uiFontSize, ...
    'HorizontalAlignment','left');

fpsSlider = uicontrol('Style','slider', ...
    'Min',1,'Max',240,'Value',fps, ...
    'Units','pixels', ...
    'Position',[rightX y0-rowH 360 rowH], ...
    'Callback',@(s,~) setFPS(s.Value));


% -------- VOLUME --------
y1 = y0 - (rowH*2 + gap);

uicontrol('Style','text','Units','pixels', ...
    'Position',[rightX y1 100 rowH], ...
    'String','Volume', ...
    'ForegroundColor','w','BackgroundColor','k', ...
    'FontName',uiFontName,'FontSize',uiFontSize, ...
    'HorizontalAlignment','right');

volValue = uicontrol('Style','text','Units','pixels', ...
    'Position',[rightX+105 y1 120 rowH], ...
    'ForegroundColor','w','BackgroundColor','k', ...
    'FontName',uiFontName,'FontSize',uiFontSize, ...
    'HorizontalAlignment','left');

volSlider = uicontrol('Style','slider', ...
    'Min',1,'Max',nVols,'Value',1, ...
    'Units','pixels', ...
    'Position',[rightX y1-rowH 360 rowH], ...
    'Callback',@(s,~) scrubVol(round(s.Value)));



% Colorbar left (thin)
caxis(par.previewCaxis);
cbar = colorbar('Position',[0.06 0.18 0.014 0.58]);
colormap(cbar,'hot');
cbar.Color = 'w';
cbar.FontSize = 13;
cbar.Label.String = 'Percent Signal Change (%)';
cbar.Label.FontSize = 14;
cbar.Label.Color = 'w';

% Footer
uicontrol('Style','text','Units','pixels',...
    'Position',[10 10 700 24],...
    'String','fUSI Video Analysis — Soner C., MPI for Biological Cybernetics, 2026',...
    'ForegroundColor',[0.7 0.7 0.7],...
    'BackgroundColor','k',...
    'HorizontalAlignment','left',...
    'FontName',uiFontName,'FontSize',11);

applyAllMaskBtn = uicontrol('Style','pushbutton', ...
    'Units','pixels', ...
    'Position',[rightX-20 118 295 36], ...   % wider so text fits
    'String','Apply current mask to ALL frames', ...
    'FontName',uiFontName, ...
    'FontSize',12, ...
    'FontWeight','bold', ...
    'ForegroundColor','w', ...
    'BackgroundColor',[0.25 0.55 0.25], ...
    'Callback',@applyMaskToAllFrames);

% Help + Close
helpBtn = uicontrol('Style','pushbutton','String','HELP',...
    'Units','pixels','Position',[rightX-20 78 105 32],...
    'FontName',uiFontName,'FontSize',12,...
    'FontWeight','bold',...
    'ForegroundColor','w',...
    'BackgroundColor',[0.10 0.35 0.95],...
    'Callback',@showHelpDialog);

closeBtn = uicontrol('Style','pushbutton','String','CLOSE',...
    'Units','pixels','Position',[rightX+90 78 105 32],...
    'FontName',uiFontName,'FontSize',12,...
    'FontWeight','bold',...
    'ForegroundColor','w',...
    'BackgroundColor',[0.75 0.15 0.15],...
    'Callback',@(~,~) close(fig));

% Play/Replay/Save
playBtn = uicontrol('Style','togglebutton','String','Play',...
    'Units','pixels','Position',[rightX-20 20 95 50],...
    'FontName',uiFontName,'FontSize',13,'Callback',@playPause);

uicontrol('Style','pushbutton','String','Replay',...
    'Units','pixels','Position',[rightX+85 20 95 50],...
    'FontName',uiFontName,'FontSize',13,'Callback',@replayVid);

uicontrol('Style','pushbutton','String','Save MP4',...
    'Units','pixels','Position',[rightX+190 20 105 50],...
    'FontName',uiFontName,'FontSize',13,'Callback',@saveVideo);

%% ===================== STATE =====================
volPos  = 1.0;   % continuous volume position (DOUBLE)
volume  = 1;     % displayed / indexed volume (INTEGER)
frame   = 1;
playing = false;

%playSpeed = 1.0; % playback speed multiplier



applyToAllFrames = true;   % keep as requested (global operations)

% mask volume
mask = false(nz,nx,nVols);
% ---- restore mask from session file if available ----
if exist('loadedMask','var') && ~isempty(loadedMask)
    if isequal(size(loadedMask), size(mask))
        mask = loadedMask;
        maskIsInclude = loadedMaskIsInclude;
        statusLine = 'Mask restored from session file.';
    else
        statusLine = sprintf( ...
            'Saved mask size mismatch (%s) vs current (%s). Mask ignored.', ...
            mat2str(size(loadedMask)), mat2str(size(mask)));
    end
end


% editor state
editorMode = false;
viewMaskedOnly = false;
maskIsInclude = true;

% brush
brushRadius = 4;
maskAlpha = 0.35;
maskColor = [0 1 0];

% strict mode selection used for AUTO MASK
% 1=Off, 2=A, 3=B
strictMode = 2;
percentileKeep = 90;       % default for auto mask (can change)

% percentile range (lower allowed)
percentileMin = 60;
percentileMax = 99;

% fill parameters (still used)
fillWindowR = 18;
fillSigmaFactor = 1.8;
fillMaxPixels = 300000;

% mouse / keyboard
mouseIsDown = false;
paintMode = '';
lastMouseXY = [NaN NaN];

% status line
statusLine = '';

%% ===================== MASK EDITOR UI =====================
uicontrol('Style','text','Units','pixels',...
    'Position',[rightX 680 360 20],...
    'String','Mask / Auto Tools',...
    'ForegroundColor',[0.9 0.9 0.9],...
    'BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',12,...
    'FontWeight','bold',...
    'HorizontalAlignment','left');

editBtn = uicontrol('Style','togglebutton','Units','pixels',...
    'Position',[rightX 650 360 28],...
    'String','Editor OFF',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'Callback',@toggleEditor);

viewBtn = uicontrol('Style','togglebutton','Units','pixels',...
    'Position',[rightX 615 175 28],...
    'String','VIEW: FULL',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'Callback',@toggleViewMasked);

posView = get(viewBtn,'Position');
includeDrop = uicontrol('Style','popupmenu','Units','pixels',...
    'Position',[posView(1)+posView(3)+10 posView(2) posView(3) posView(4)],...
    'String',{'Include','Exclude'},...
    'Value',1,...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'Callback',@setIncludeExclude);

% AUTO tools: apply to current vs all frames
applyAllBtn = uicontrol('Style','togglebutton','Units','pixels',...
    'Position',[rightX 580 175 28],...
    'String','AUTO: ALL',...
    'Value',1,...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'Callback',@toggleApplyAll);

autoBtn = uicontrol('Style','pushbutton','Units','pixels',...
    'Position',[rightX+185 580 175 28],...
    'String','AUTO MASK (M)',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'FontWeight','bold',...
    'Callback',@autoMaskButton);

% Brush size
uicontrol('Style','text','Units','pixels',...
    'Position',[rightX 555 120 18],...
    'String','Brush size',...
    'ForegroundColor',[0.8 0.8 0.8],...
    'BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'HorizontalAlignment','left');

brushSlider = uicontrol('Style','slider','Units','pixels',...
    'Position',[rightX 535 360 18],...
    'Min',1,'Max',25,'Value',brushRadius,...
    'Callback',@(s,~) setBrush(round(s.Value)));

brushVal = uicontrol('Style','text','Units','pixels',...
    'Position',[rightX 515 360 16],...
    'String',sprintf('Radius: %d px',brushRadius),...
    'ForegroundColor',[0.7 0.7 0.7],...
    'BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'HorizontalAlignment','left');

% Fill / Color / Clear
colorBtn = uicontrol('Style','pushbutton','Units','pixels',...
    'Position',[rightX 480 110 28],...
    'String','Color...',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'Callback',@pickColor);

fillBtn = uicontrol('Style','pushbutton','Units','pixels',...
    'Position',[rightX+120 480 110 28],...
    'String','Fill (F)',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'Callback',@fillRegion);

clearBtn = uicontrol('Style','pushbutton','Units','pixels',...
    'Position',[rightX+240 480 120 28],...
    'String','Clear mask',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'Callback',@clearMaskAll);

% Strict selector for AUTO MASK
uicontrol('Style','text','Units','pixels',...
    'Position',[rightX 450 360 18],...
    'String','Auto method (A/B)',...
    'ForegroundColor',[0.8 0.8 0.8],...
    'BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'HorizontalAlignment','left');

strictDrop = uicontrol('Style','popupmenu','Units','pixels',...
    'Position',[rightX 425 360 26],...
    'String',{'Off','A: robust','B: percentile'},...
    'Value',strictMode,...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'Callback',@setStrictMode);

percSlider = uicontrol('Style','slider','Units','pixels',...
    'Position',[rightX 395 360 16],...
    'Min',percentileMin,'Max',percentileMax,'Value',percentileKeep,...
    'Callback',@setPercentileKeep);

percVal = uicontrol('Style','text','Units','pixels',...
    'Position',[rightX 375 360 18],...
    'String',sprintf('Percentile: %.0f (lower = more voxels)',percentileKeep),...
    'ForegroundColor',[0.7 0.7 0.7],...
    'BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'HorizontalAlignment','left');

saveMaskBtn = uicontrol('Style','pushbutton','Units','pixels',...
    'Position',[rightX 330 360 32],...
    'String','Save Mask (.mat)',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'Callback',@saveMaskMat);

saveMaskedBtn = uicontrol('Style','pushbutton','Units','pixels',...
    'Position',[rightX 290 360 32],...
    'String','Save Masked PSC (.mat)',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'Callback',@saveMaskedPSCMat);

hintTxt = uicontrol('Style','text','Units','pixels',...
    'Position',[rightX 190 360 90],...
    'String',sprintf([ ...
        'AUTO MASK: press M or button (uses A/B)\n' ...
        'FILL: press F (uses mouse position)\n' ...
        'Paint: L add / R remove\n' ...
        'VIEW:MASKED hides overlay but applies mask\n' ...
        'Include/Exclude affects MASKED output']),...
    'ForegroundColor',[0.65 0.65 0.65],...
    'BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'HorizontalAlignment','left');

% Callbacks
set(fig,'WindowButtonDownFcn',@mouseDown);
set(fig,'WindowButtonUpFcn',@mouseUp);
set(fig,'WindowButtonMotionFcn',@mouseMove);
set(fig,'KeyPressFcn',@keyPressHandler);

% Initial render
render();

while ishandle(fig)
    if playing
        volume = min(volume + 1, nVols);
        volSlider.Value = volume;

        frame = (volume - 1) * par.interpol + 1;
        frame = min(frame, nFrames);

        render();

        drawnow limitrate;
        pause(1 / max(fps, 0.1));   % ← FPS controls timing

    else
        pause(0.02);
    end
end



%% ====== CONTINUES IN PART 2/3 ======
%% ============================================================
%  PART 2/3 — Core callbacks + AUTO MASK
%% ============================================================

function render()

    % =========================================================
    % FIX 3a — FRAME SAFETY (prevents invalid indexing)
    % =========================================================
    if frame < 1 || frame > size(PSC,3)
        img.CData = zeros(nz,nx,3,'single');
        return;
    end

    % =========================================================
    % PSC → RGB (safe normalization)
    % =========================================================
    A = (PSC(:,:,frame) - par.previewCaxis(1)) / diff(par.previewCaxis);
    A = max(0, min(1, A));

    pscRGB = ind2rgb(round(A*(Nc-1)), mapA) + bgRGB;

    % =========================================================
    % FIX 3b — NaN / Inf / white-screen guard
    % =========================================================
    if ~any(isfinite(pscRGB(:)))
        img.CData = bgRGB;   % show anatomy only
        return;
    end

    % =========================================================
    % Current VOLUME mask
    % =========================================================
    M = mask(:,:,volume);

    % =========================================================
    % VIEW logic
    % =========================================================
    if viewMaskedOnly

        if maskIsInclude
            show = M;
        else
            show = ~M;
        end

        baseRGB = bgRGB;
        idx = repmat(show, [1 1 3]);
        baseRGB(idx) = pscRGB(idx);
        img.CData = baseRGB;

    else
        baseRGB = pscRGB;

        % Tint overlay ONLY when editor is ON (clean visualization)
        if editorMode && any(M(:))
            overlay = zeros(nz,nx,3,'single');
            overlay(:,:,1) = maskColor(1);
            overlay(:,:,2) = maskColor(2);
            overlay(:,:,3) = maskColor(3);

            baseRGB = baseRGB .* (~M) + ...
                      (1-maskAlpha) .* baseRGB .* M + ...
                      maskAlpha .* overlay .* M;
        end

        img.CData = baseRGB;
    end

    % =========================================================
    % INFO TEXT (volume-correct time)
    % =========================================================
    t = (volume - 1) * TR;

    if editorMode, em='ON'; else, em='OFF'; end
    if viewMaskedOnly, vm='MASKED'; else, vm='FULL'; end
    if maskIsInclude, ms='Include'; else, ms='Exclude'; end

    methodLabel = 'Off';
    if strictMode==2, methodLabel='A'; end
    if strictMode==3, methodLabel='B'; end

    extra = '';
    if ~isempty(statusLine)
        extra = sprintf('\n%s', statusLine);
    end

info.String = sprintf( ...
    ['t = %.1f / %.1f s   |   Volume %d / %d   |   View: %s (%s)\n' ...
     'Baseline: %g–%g %s   |   Editor: %s   |   AUTO method: %s%s'], ...
    t, Tmax, volume, nVols, vm, ms, ...
    baseline.start, baseline.end, baseline.mode, ...
    em, methodLabel, extra);



    fpsValue.String = sprintf('%.1f', fps);
    volValue.String = sprintf('%d / %d', volume, nVols);

    % Enable percentile slider only for Method B
    if strictMode == 3
        set(percSlider,'Enable','on');
    else
        set(percSlider,'Enable','off');
    end

end


%% ---------------- PLAYBACK ----------------
    
   function scrubVol(v)

    playing = false;
    playBtn.Value = 0;
    playBtn.String = 'Play';

    volume = min(max(1,v), nVols);
    volPos = volume;   % ← CRITICAL

    frame = (volume - 1)*par.interpol + 1;
    frame = min(frame, nFrames);

    render();
end


 function setFPS(v)
    fps = max(1, min(240, round(v)));
    fpsValue.String = sprintf('%.0f', fps);
end


    function scrub(v)
        playing=false;
        playBtn.Value=0;
        playBtn.String='Play';
        frame=v;
        render();
    end

function setSpeed(v)
    playSpeed = v;
    speedVal.String = sprintf('%.2f×', v);
end



    function playPause(src,~)
        playing = logical(src.Value);
        if playing
            src.String='Pause';
        else
            src.String='Play';
        end
    end

 function replayVid(~,~)

    playing = false;

    volPos = 1;
    volume = 1;
    volSlider.Value = 1;

    frame = 1;

    playing = true;
    playBtn.Value = 1;
    playBtn.String = 'Pause';

    render();
end


 function saveVideo(~,~)

    [f,p] = uiputfile('*.mp4','Save fUSI video');
    if isequal(f,0), return; end

    % Compute export frame rate from TR and playSpeed
exportFPS = fps;


    vid = VideoWriter(fullfile(p,f),'MPEG-4');
    vid.FrameRate = exportFPS;
    vid.Quality   = 95;
    open(vid);

    % Temporary text overlay on AXES (not whole GUI)
    txt = text(ax, 0.01, 0.99, '', ...
        'Units','normalized', ...
        'Color','w', ...
        'FontName','Courier', ...
        'FontSize',16, ...
        'FontWeight','bold', ...
        'VerticalAlignment','top', ...
        'BackgroundColor','k', ...
        'Margin',4);

    oldVolume = volume;
    oldPlaying = playing;
    playing = false;

    for v = 1:nVols

        volume = v;
        volSlider.Value = v;

        frame = (v - 1)*par.interpol + 1;
        frame = min(frame, nFrames);

        render();

        t = (v - 1) * TR;
        txt.String = sprintf('t = %.1f / %.1f s   |   Volume %d / %d', ...
                             t, Tmax, v, nVols);

        drawnow;
        writeVideo(vid, getframe(ax));
    end

    delete(txt);
    close(vid);

    volume = oldVolume;
    playing = oldPlaying;
    render();

    statusLine = 'Video saved (image-only).';
end


%% ---------------- HELP DIALOG (PRIORITIZED + COLORED) ----------------
function showHelpDialog(~,~)

    % Resizable help window
    hf = figure( ...
        'Name','Help — fUSI Auto Mask v3', ...
        'Color',[0.06 0.06 0.06], ...
        'MenuBar','none', ...
        'ToolBar','none', ...
        'NumberTitle','off', ...
        'Position',[250 120 920 740], ...
        'Resize','on', ...
        'WindowStyle','modal');

    %% ---------- COLORS ----------
    colTitle    = [0.98 0.98 0.98];
    colWorkflow = [0.55 0.85 0.55];
    colAuto     = [0.55 0.75 0.95];
    colFill     = [0.95 0.75 0.50];
    colNormal   = [0.90 0.90 0.90];
    colMuted    = [0.70 0.70 0.70];

    %% ---------- TITLE ----------
    titleTxt = uicontrol('Style','text','Parent',hf,...
        'Units','pixels',...
        'Position',[20 690 880 36],...
        'String','fUSI Video Analysis — HELP', ...
        'ForegroundColor',colTitle,...
        'BackgroundColor',[0.06 0.06 0.06],...
        'FontName','Arial',...
        'FontSize',20,...
        'FontWeight','bold',...
        'HorizontalAlignment','left');

    %% ---------- MAIN TEXT (EDIT BOX) ----------
    msg = [ ...
        'RECOMMENDED WORKFLOW\n' ...
        '============================================================\n' ...
        '1) Select AUTO MASK method (A or B)\n' ...
        '2) Press M or click AUTO MASK\n' ...
        '3) Use FILL (F) to expand coherent regions\n' ...
        '4) Refine mask manually with paint tools\n' ...
        '5) Switch VIEW: MASKED for final visualization\n\n' ...
        ...
        'AUTO MASK (METHOD A / B)\n' ...
        '============================================================\n' ...
        '- Press M or click AUTO MASK\n' ...
        '- Method A: Robust automatic detection (broad, permissive)\n' ...
        '- Method B: Percentile-based detection (lower = more voxels)\n' ...
        '- AUTO: ALL applies automatic mask to all frames\n' ...
        '- VIEW: FULL + AUTO: ALL  -> draw and display mask\n' ...
        '- VIEW: MASKED + AUTO: ALL -> apply mask without showing overlay\n\n' ...
        ...
        'FILL TOOL\n' ...
        '============================================================\n' ...
        '- Press F with mouse over image\n' ...
        '- Or click Fill button\n' ...
        '- Expands region based on local signal similarity\n\n' ...
        ...
        'MANUAL EDITING\n' ...
        '============================================================\n' ...
        '- Left mouse  -> add to mask\n' ...
        '- Right mouse -> remove from mask\n' ...
        '- Brush size controls spatial extent\n\n' ...
        ...
        'VIEW MODES\n' ...
        '============================================================\n' ...
        '- VIEW: FULL   -> show PSC everywhere\n' ...
        '- VIEW: MASKED -> hide overlay but apply mask to data\n' ...
        '- Include / Exclude affects MASKED output only\n\n' ...
        ...
        'KEYBOARD SHORTCUTS\n' ...
        '============================================================\n' ...
        '- M  -> Automatic mask (A / B)\n' ...
        '- F  -> Fill at cursor position\n' ...
        ];

    txtBox = uicontrol('Style','edit','Parent',hf,...
        'Units','pixels',...
        'Position',[20 90 880 580],...
        'String',sprintf(msg),...
        'ForegroundColor',colNormal,...
        'BackgroundColor',[0.12 0.12 0.12],...
        'FontName','Arial',...
        'FontSize',14,...
        'HorizontalAlignment','left',...
        'Max',2,'Min',0,...
        'Enable','inactive');



    %% ---------- CLOSE BUTTON ----------
    closeBtn = uicontrol('Style','pushbutton','Parent',hf,...
        'Units','pixels',...
        'Position',[770 25 130 42],...
        'String','Close',...
        'FontName','Arial',...
        'FontSize',13,...
        'FontWeight','bold',...
        'ForegroundColor','w',...
        'BackgroundColor',[0.75 0.15 0.15],...
        'Callback',@(~,~) close(hf));

    %% ---------- RESIZE HANDLER ----------
    hf.SizeChangedFcn = @onResize;
    onResize();

    function onResize(~,~)
        p = hf.Position;
        W = p(3); H = p(4);

        set(titleTxt,'Position',[20 H-50 W-40 36]);
        set(txtBox,  'Position',[20 90 W-40 H-150]);
        set(closeBtn,'Position',[W-150 25 130 42]);
    end
end

%% ---------------- UI CALLBACKS ----------------
    function toggleEditor(src,~)
        editorMode = logical(src.Value);
        if editorMode, src.String='Editor ON'; else, src.String='Editor OFF'; end
        statusLine='';
        render();
    end

    function toggleViewMasked(src,~)
        viewMaskedOnly = logical(src.Value);
        if viewMaskedOnly, src.String='VIEW: MASKED'; else, src.String='VIEW: FULL'; end
        statusLine='';
        render();
    end

    function setIncludeExclude(src,~)
        maskIsInclude = (src.Value==1);
        statusLine='';
        render();
    end

    function toggleApplyAll(src,~)
        applyToAllFrames = logical(src.Value);
        if applyToAllFrames
            src.String = 'AUTO: ALL';
        else
            src.String = 'AUTO: FRAME';
        end
        statusLine='';
        render();
    end

function applyMaskToAllFrames(~,~)

    refMask = mask(:,:,frame);

    if ~any(refMask(:))
        statusLine = 'Current frame mask is empty — nothing to apply.';
        render();
        return;
    end

for v = 1:nVols
    mask(:,:,v) = refMask;
end

    statusLine = sprintf('Mask from frame %d applied to ALL volumes.', volume);
    render();
end


    function autoMaskButton(~,~)
        autoMask();
    end

    function setBrush(v)
        brushRadius = max(1,round(v));
        brushVal.String = sprintf('Radius: %d px',brushRadius);
    end

    function pickColor(~,~)
        c = uisetcolor(maskColor,'Pick mask overlay color');
        if numel(c)==3
            maskColor=c;
        end
        render();
    end

    function clearMaskAll(~,~)
        mask(:)=false;
        statusLine='Mask cleared.';
        render();
    end

    function setStrictMode(src,~)
        strictMode = src.Value;
        statusLine='';
        render();
    end

    function setPercentileKeep(src,~)
        percentileKeep = round(src.Value);
        percVal.String = sprintf('Percentile: %.0f (lower = more voxels)', percentileKeep);
        statusLine='';
        render();
    end

%% ---------------- MOUSE PAINT ----------------
   function mouseDown(~,~)
    if playing, return; end
    if ~editorMode, return; end

    mouseIsDown = true;

    sel = get(fig,'SelectionType');
    if strcmp(sel,'normal')
        paintMode = 'add';        % LEFT click
    elseif strcmp(sel,'alt')
        paintMode = 'remove';     % RIGHT click
    else
        paintMode = '';
        mouseIsDown = false;
        return;
    end

    applyPaintAtCursor();
end

function mouseUp(~,~)
    mouseIsDown = false;
    paintMode   = '';
end

    function mouseMove(~,~)

    % Cache mouse position for Fill (F)
    cp = get(ax,'CurrentPoint');
    x = cp(1,1);
    y = cp(1,2);
    if x>=1 && x<=nx && y>=1 && y<=nz
        lastMouseXY = [x y];
    end

    if ~mouseIsDown, return; end
    if playing, return; end
    if ~editorMode, return; end
    if isempty(paintMode), return; end

    applyPaintAtCursor();
end

function applyPaintAtCursor()

    cp = get(ax,'CurrentPoint');
    x = round(cp(1,1));
    y = round(cp(1,2));
    if x<1||x>nx||y<1||y>nz
        return;
    end

    brush = makeBrushMask(x,y,brushRadius,nz,nx);

    if strcmp(paintMode,'add')
        if applyToAllFrames
            for v = 1:nVols
                mask(:,:,v) = mask(:,:,v) | brush;
            end
        else
            mask(:,:,volume) = mask(:,:,volume) | brush;
        end
    else
        if applyToAllFrames
            for v = 1:nVols
                mask(:,:,v) = mask(:,:,v) & ~brush;
            end
        else
            mask(:,:,volume) = mask(:,:,volume) & ~brush;
        end
    end

    statusLine = '';
    render();
end


%% ---------------- KEYBOARD SHORTCUTS ----------------
    function keyPressHandler(~,evt)
        if ~isfield(evt,'Key'), return; end
        key = evt.Key;

        % Fill with F (mouse over image)
        if strcmpi(key,'f')
            if playing
                statusLine='Fill disabled while playing.';
                render(); return;
            end
            if any(isnan(lastMouseXY))
                statusLine='Move mouse over image, then press F.';
                render(); return;
            end
            fillAtXY(lastMouseXY(1), lastMouseXY(2));
            return;
        end

        % Auto mask with M
        if strcmpi(key,'m')
            autoMask();
            return;
        end
    end

    function fillAtXY(xf,yf)
        % actual fill core in PART 3
        x0 = round(xf);
        y0 = round(yf);
        if x0<1||x0>nx||y0<1||y0>nz
            statusLine='Fill aborted: mouse outside image.';
            render();
            return;
        end
        fillRegionAtSeed(x0,y0);
    end

%% ---------------- AUTO MASK CORE ----------------
    function autoMask()

        if strictMode == 1
            statusLine = 'Auto mask: select method A or B (not Off).';
            render();
            return;
        end
if applyToAllFrames
    for v = 1:nVols
        frameIdx = (v-1)*par.interpol + 1;
        frameIdx = min(frameIdx, nFrames);

        sig = PSC(:,:,frameIdx);
        auto = autoThreshold(sig, strictMode, percentileKeep);
        auto = cleanupMask(auto);

        mask(:,:,v) = auto;
    end
else
    frameIdx = (volume-1)*par.interpol + 1;
    frameIdx = min(frameIdx, nFrames);

    sig = PSC(:,:,frameIdx);
    auto = autoThreshold(sig, strictMode, percentileKeep);
    auto = cleanupMask(auto);

    mask(:,:,volume) = auto;

            statusLine = sprintf('Auto mask applied to frame %d (method %s).', frame, autoMethodLabel());
        end

        render();
    end

    function s = autoMethodLabel()
        if strictMode==2, s='A'; elseif strictMode==3, s='B'; else, s='Off'; end
    end

%% ============================================================
% ===== Fill core + thresholds + helpers in PART 3/3 ==========
%% ============================================================

%% ============================================================
%  PART 3/3 — Auto-threshold A/B + Fill core + helpers
%% ============================================================

%% ---------------- AUTO THRESHOLD (THIS IS THE KEY FIX) ----------------
% This function is responsible for AUTOMATIC marking.
% It MUST be permissive and ALWAYS return something reasonable.

    function BW = autoThreshold(sig, mode, pKeep)

        % Flatten & clean
        v = sig(:);
        v = v(~isnan(v) & ~isinf(v));

        if isempty(v)
            BW = false(size(sig));
            return;
        end

        switch mode

            case 2
                % ================= METHOD A =================
                % Very permissive, robust threshold
                % Designed to ALWAYS give a vascular-like mask

                med  = median(v);
                madv = mad(v,1);

                % EXTREMELY permissive threshold
                thr = med + 0.3 * madv;

            case 3
                % ================= METHOD B =================
                % Percentile-based, softened heavily

                % Slider value is softened by 40 percent points
                pSoft = max(40, pKeep - 40);
                thr   = prctile(v, pSoft);

            otherwise
                BW = false(size(sig));
                return;
        end

        % Initial binary mask
        BW = sig > thr;

        % Fallback: if mask is too small, relax automatically
        if nnz(BW) < 20
            thr2 = prctile(v, 30);   % emergency fallback
            BW = sig > thr2;
        end
    end

%% ---------------- CLEANUP MASK ----------------
% Keeps only spatially meaningful regions.
% Does NOT kill everything if Image Processing Toolbox is absent.

    function BW = cleanupMask(BW)

        % Remove tiny speckles
        minPix = 15;

        try
            % If Image Processing Toolbox exists
            BW = bwareaopen(BW, minPix);
        catch
            % Toolbox-free fallback
            CC = bwconncomp(BW,4);
            BW(:) = false;
            for i = 1:CC.NumObjects
                if numel(CC.PixelIdxList{i}) >= minPix
                    BW(CC.PixelIdxList{i}) = true;
                end
            end
        end
    end

%% ---------------- FILL CORE ----------------
% Region grow based on similarity.
% NEVER blocked by strict A/B.

    function fillRegion(~,~)
        cp = get(ax,'CurrentPoint');
        x0 = round(cp(1,1));
        y0 = round(cp(1,2));
        if x0<1||x0>nx||y0<1||y0>nz
            statusLine = 'Fill: move mouse over image, then press F.';
            render();
            return;
        end
        fillRegionAtSeed(x0,y0);
    end

    function fillRegionAtSeed(x0, y0)

    % ---------------- SAFETY CHECKS ----------------
    if playing
        statusLine = 'Fill disabled while playing.';
        render(); 
        return;
    end

    if ~editorMode
        statusLine = 'Enable Editor before Fill.';
        render(); 
        return;
    end

    % ---------------- GET CURRENT FRAME SIGNAL ----------------
    frameIdx = (volume - 1) * par.interpol + 1;
    frameIdx = min(frameIdx, nFrames);

    sig = PSC(:,:,frameIdx);

    % ---------------- REGION GROW ----------------
    region = regionGrowSimilarity( ...
        sig, y0, x0, ...
        fillWindowR, fillSigmaFactor, fillMaxPixels);

    % ---------------- EMPTY RESULT CHECK ----------------
    if ~any(region(:))
        statusLine = 'Fill found no similar region.';
        render();
        return;
    end

    % ---------------- APPLY MASK ----------------
    if applyToAllFrames
        for v = 1:nVols
            mask(:,:,v) = mask(:,:,v) | region;
        end
    else
        mask(:,:,volume) = mask(:,:,volume) | region;
    end

    % ---------------- FINALIZE ----------------
    statusLine = sprintf('Fill added %d pixels.', nnz(region));
    render();
end

%% ---------------- REGION GROW (VERY PERMISSIVE) ----------------

    function region = regionGrowSimilarity(sig, sy, sx, winR, sigmaFactor, maxPix)

        [H,W] = size(sig);
        region = false(H,W);

        % Local stats around seed
        yMin = max(1, sy-winR);
        yMax = min(H, sy+winR);
        xMin = max(1, sx-winR);
        xMax = min(W, sx+winR);

        patch = sig(yMin:yMax, xMin:xMax);
        mu = mean(patch(:));
        sd = std(patch(:));
        if sd < 1e-6, sd = 1e-6; end

        % VERY permissive tolerance
        tol = max( sigmaFactor*sd , 0.25*abs(mu) );

        % BFS queue
        qy = zeros(maxPix,1);
        qx = zeros(maxPix,1);
        head=1; tail=1;
        qy(1)=sy; qx(1)=sx;
        region(sy,sx)=true;

        while head<=tail && tail<maxPix

            y=qy(head); x=qx(head); head=head+1;

            % 4-connected neighbors
            if y>1 && ~region(y-1,x)
                if abs(sig(y-1,x)-mu)<=tol
                    tail=tail+1; qy(tail)=y-1; qx(tail)=x;
                    region(y-1,x)=true;
                end
            end
            if y<H && ~region(y+1,x)
                if abs(sig(y+1,x)-mu)<=tol
                    tail=tail+1; qy(tail)=y+1; qx(tail)=x;
                    region(y+1,x)=true;
                end
            end
            if x>1 && ~region(y,x-1)
                if abs(sig(y,x-1)-mu)<=tol
                    tail=tail+1; qy(tail)=y; qx(tail)=x-1;
                    region(y,x-1)=true;
                end
            end
            if x<W && ~region(y,x+1)
                if abs(sig(y,x+1)-mu)<=tol
                    tail=tail+1; qy(tail)=y; qx(tail)=x+1;
                    region(y,x+1)=true;
                end
            end
        end
    end

%% ---------------- BRUSH MASK ----------------
    function B = makeBrushMask(x,y,r,nz,nx)
        B = false(nz,nx);
        xMin=max(1,x-r); xMax=min(nx,x+r);
        yMin=max(1,y-r); yMax=min(nz,y+r);
        for yy=yMin:yMax
            for xx=xMin:xMax
                if (xx-x)^2+(yy-y)^2 <= r^2
                    B(yy,xx)=true;
                end
            end
        end
    end

%% ---------------- SAVE FUNCTIONS ----------------
  function saveMaskMat(~,~)

    [f,p] = uiputfile('*.mat','Save mask + data');
    if isequal(f,0), return; end

    out.I          = I;          % ORIGINAL DATA
    out.mask       = mask;       % FULL MASK VOLUME
    out.maskIsInclude = maskIsInclude;
    out.par        = par;
    out.TR         = TR;
    out.baseline   = baseline;
    out.frameInfo  = struct( ...
        'nFrames', nFrames, ...
        'interpol', par.interpol );

    save(fullfile(p,f),'-struct','out','-v7.3');

    statusLine = 'Session (data + mask) saved.';
    render();
end


    function saveMaskedPSCMat(~,~)

    [f,p] = uiputfile('*.mat','Save masked PSC + data');
    if isequal(f,0), return; end

   
    PSCm = PSC;

for v = 1:nVols
    frame0 = (v-1)*par.interpol + 1;
    frame1 = min(nFrames, v*par.interpol);

    M = mask(:,:,v);

    for k = frame0:frame1
        if maskIsInclude
            PSCm(:,:,k) = PSC(:,:,k) .* single(M);
        else
            PSCm(:,:,k) = PSC(:,:,k) .* single(~M);
        end
    end
end


    out.I            = I;          % ORIGINAL DATA
    out.PSC_masked   = PSCm;       % RESULT
    out.mask         = mask;
    out.maskIsInclude = maskIsInclude;
    out.par          = par;
    out.TR           = TR;
    out.baseline     = baseline;

    save(fullfile(p,f),'-struct','out','-v7.3');

    statusLine = 'Masked PSC session saved.';
    render();
end


%% ---------------- CLOSE GUI FUNCTION ----------------
end   % play_fusi_video_final

%% ================= END OF FILE =================

function QC = runFrameRateQC(I, TR, tag, savePNG)
if nargin < 3 || isempty(tag)
    tag = 'ORIGINAL';
end
if nargin < 4
    savePNG = false;
end


% Diagnostic frame-rate / PNS QC 
% No data modification

[nz,nx,nVols] = size(I);

% Global mean per volume
g = squeeze(mean(mean(I,1),2));
g = g(:);

% Normalize
gNorm = g ./ median(g);

% Robust noise estimation (lower tail only)
gLow = gNorm(gNorm < 1);
sigma = sqrt(mean((gLow - 1).^2));

threshold = 1 + 3*sigma;
outliers  = gNorm > threshold;
rejPct    = 100 * mean(outliers);

%% ---- FIGURE 1: Intensity distribution ----
fig1 = figure( ...
    'Name','Frame-rate / PNS QC — Intensity distribution', ...
    'Color','w','Position',[200 200 900 450]);

nbins = 100;
[counts, edges] = histcounts(gNorm, nbins);
centers = edges(1:end-1) + diff(edges)/2;

bar(centers, counts, 'FaceColor',[0.7 0.7 0.7],'EdgeColor','none'); hold on

x = linspace(min(gNorm), max(gNorm), 500);
gauss = exp(-0.5*((x-1)/sigma).^2);
gauss = gauss * sum(counts) * (centers(2)-centers(1)) / ...
        (sigma*sqrt(2*pi));
plot(x, gauss,'k','LineWidth',2)

yl = ylim;
plot([threshold threshold], yl,'r','LineWidth',2)

xlabel('Normalized global intensity')
ylabel('Number of volumes')
title('Global signal stability (Urban / Montaldo)')

annotation(fig1,'textbox',[0.60 0.60 0.35 0.30], ...
    'String',sprintf([ ...
        'Threshold (3σ): %.3f\n' ...
        'Rejected volumes: %.1f %%\n\n' ...
        'Interpretation:\n' ...
        '• < 10 %%  → stable acquisition\n' ...
        '• > 10 %%  → mechanical / motion instability\n' ...
        '• Clustered rejections → task-related movement\n\n' ...
        'QC only — no data modified.' ], threshold, rejPct), ...
    'FitBoxToText','on','BackgroundColor','w','FontSize',11);

grid on
hold off

%% ---- FIGURE 2: Rejected volumes over time ----
fig2 = figure( ...
    'Name','Frame-rate / PNS QC — Rejected volumes', ...
    'Color','w','Position',[200 700 900 250]);

t = (0:nVols-1) * TR;
stem(t, outliers,'filled','MarkerSize',4);
ylim([-0.1 1.1])
xlabel('Time (s)')
ylabel('0 = accepted   1 = rejected')
title('Rejected volumes over time')
grid on

%% ---- OUTPUT ----
QC.globalNormSignal = gNorm;
QC.outliers  = outliers;
QC.threshold = threshold;
QC.sigma     = sigma;
QC.rejPct    = rejPct;
QC.figIntensity = fig1;
QC.figRejected  = fig2;
end
function Iout = interpolateRejectedVolumes(I, outliers)
% Interpolates rejected volumes along time (Urban/Montaldo style)
% INPUT:
%   I        : [nz × nx × nVols]
%   outliers : logical [nVols × 1]
% OUTPUT:
%   Iout     : same size as I

[nz,nx,nVols] = size(I);
Iout = I;

good = find(~outliers);
bad  = find(outliers);

% Safety
if numel(good) < 2
    warning('Too few valid volumes for interpolation. Data left unchanged.');
    return;
end

tAll  = 1:nVols;

for z = 1:nz
    for x = 1:nx
        sig = squeeze(I(z,x,good));
        Iout(z,x,:) = interp1(good, sig, tAll, 'linear', 'extrap');
    end
end
end
