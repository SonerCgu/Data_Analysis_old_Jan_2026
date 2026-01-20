%% ============================================================
%  fusi_video_soner4.m
%  fUSI VIDEO ANALYSIS GUI (MATLAB 2017b)
%
%  KEY FEATURE (requested):
%   - Automatic marking based on Method A/B:
%       * AUTO MASK button (and key 'M')
%       * Uses selected Strict mode A or B to create mask automatically
%       * Works on current frame or all frames (toggle)
%
%  Manual editing:
%   - Left click add, right click remove
%   - Fill button or key 'F' (uses current mouse position)
%
%  Visualization:
%   - VIEW:FULL  -> PSC everywhere, optional mask overlay tint (when Editor ON)
%   - VIEW:MASKED-> overlay hidden, mask applied (Include or Exclude semantics)
%
% ============================================================

clearvars; clc;

%% ===================== USER PARAMETERS =====================
% Temporal interpolation factor
% Higher = smoother video, more frames, slower rendering
% Lower  = faster playback, more discrete time
par.interpol = 6;
% par.interpol = 4;        % standard (faster)
% par.interpol = 8;        % very smooth preview (slow)

% Low-pass filter cutoff (Hz)
% Higher = preserve fast dynamics (more noise)
% Lower  = smoother PSC (PACAP-friendly)
par.LPF = 0.4;
% par.LPF = 0.3;           % smoother
% par.LPF = 0.2;           % very smooth long effects

% High-pass filter cutoff (Hz) - usually OFF for fUSI
par.HPF = 0;
% par.HPF = 0.01;          % only if strong drift exists

% Gaussian smoothing kernel size (pixels)
par.gaussSize = 5;
% par.gaussSize = 3;       % sharper vessels
% par.gaussSize = 7;       % heavy smoothing

% Gaussian sigma (pixels)
par.gaussSig  = 0.8;
% par.gaussSig = 0.6;      % sharper
% par.gaussSig = 1.0;      % smoother

% Connectivity filter radius (pixels)
par.conectSize = 5;
% par.conectSize = 3;      % permissive
% par.conectSize = 7;      % strict

% Connectivity threshold (% PSC)
par.conectLev  = 25;       % PACAP-safe default
% par.conectLev = 20;      % permissive
% par.conectLev = 40;      % strict

% Display PSC range (visualization only)
par.previewCaxis = [0 60];
% par.previewCaxis = [0 40]; % higher contrast
% par.previewCaxis = [0 80]; % conservative

% Baseline window (sec)
baseline.mode  = 'sec';
baseline.start = 0;
baseline.end   = 60;
% baseline.end = 30;       % shorter baseline

TR = 1;
initialFPS = 30;
maxFPS     = 120;

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
switch lower(ext)
    case '.mat'
        S = load(dataFile);
        if ~isfield(S,'I')
            error('MAT file must contain variable I.');
        end
        I = single(S.I);

    case '.nii'
        V = niftiread(dataFile);
        I = convertNiftiToI(V);

    otherwise
        error('Unsupported file type for this viewer.');
end

%% ===================== PROCESSING =====================
I1 = interpFilm(I, par.interpol);
[nz,nx,nFrames] = size(I1);
nVols = size(I,3);
Tmax  = nVols * TR;

b0 = round((baseline.start/TR)*par.interpol) + 1;
b1 = round((baseline.end/TR)*par.interpol);
b0 = max(1,b0);
b1 = min(nFrames,b1);

ab = mean(I1(:,:,b0:b1),3);
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
bg = 20*log10(bg ./ max(bg(:)));

%% ===================== GUI =====================
play_fusi_video_final(PSC,bg,par,initialFPS,maxFPS,TR,Tmax,baseline);

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
function play_fusi_video_final(PSC,bg,par,fps,maxFPS,TR,Tmax,baseline)

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

% FPS
fpsLabel = uicontrol('Style','text','Units','pixels',...
    'Position',[rightX 805 100 20],...
    'String','FPS:',...
    'ForegroundColor','w','BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'HorizontalAlignment','right');

fpsValue = uicontrol('Style','text','Units','pixels',...
    'Position',[rightX+105 805 80 20],...
    'ForegroundColor','w','BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',uiFontSize);

fpsSlider = uicontrol('Style','slider',...
    'Min',1,'Max',maxFPS,'Value',fps,...
    'Units','pixels',...
    'Position',[rightX 780 360 22],...
    'Callback',@(s,~) setFPS(s.Value));

% Frame slider
frameLabel = uicontrol('Style','text','Units','pixels',...
    'Position',[rightX 745 100 20],...
    'String','Frame:',...
    'ForegroundColor','w','BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',uiFontSize,...
    'HorizontalAlignment','right');

frameValue = uicontrol('Style','text','Units','pixels',...
    'Position',[rightX+105 745 160 20],...
    'ForegroundColor','w','BackgroundColor','k',...
    'FontName',uiFontName,'FontSize',uiFontSize);

frameSlider = uicontrol('Style','slider',...
    'Min',1,'Max',nFrames,'Value',1,...
    'Units','pixels',...
    'Position',[rightX 720 360 22],...
    'Callback',@(s,~) scrub(round(s.Value)));

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
frame = 1;
playing = false;

applyToAllFrames = true;   % keep as requested (global operations)

% mask volume
mask = false(nz,nx,nFrames);

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

% Main loop
while ishandle(fig)
    if playing
        frame = min(frame+1,nFrames);
        frameSlider.Value = frame;
        render();
        pause(1/max(fps,0.01));
    else
        pause(0.02);
    end
end

%% ====== CONTINUES IN PART 2/3 ======
%% ============================================================
%  PART 2/3 — Core callbacks + AUTO MASK
%% ============================================================

%% ---------------- RENDER ----------------
    function render()

        % PSC -> RGB
        A = (PSC(:,:,frame)-par.previewCaxis(1))/diff(par.previewCaxis);
        A = max(0,min(1,A));
        pscRGB = ind2rgb(round(A*(Nc-1)),mapA) + bgRGB;

        % current mask
        M = mask(:,:,frame);

        % VIEW logic
        if viewMaskedOnly
            if maskIsInclude
                show = M;
            else
                show = ~M;
            end
            baseRGB = bgRGB;
            idx = repmat(show,[1 1 3]);
            baseRGB(idx) = pscRGB(idx);
            img.CData = baseRGB;
        else
            baseRGB = pscRGB;

            % tint overlay only if editor ON (cleaner)
            if editorMode && any(M(:))
                overlay = zeros(nz,nx,3);
                overlay(:,:,1) = maskColor(1);
                overlay(:,:,2) = maskColor(2);
                overlay(:,:,3) = maskColor(3);
                baseRGB = baseRGB.*(~M) + (1-maskAlpha)*baseRGB.*M + maskAlpha*overlay.*M;
            end

            img.CData = baseRGB;
        end

        % Info text
        t = (frame/par.interpol)*TR;
        if editorMode, em='ON'; else, em='OFF'; end
        if viewMaskedOnly, vm='MASKED'; else, vm='FULL'; end
        if maskIsInclude, ms='Include'; else, ms='Exclude'; end

        methodLabel = 'Off';
        if strictMode==2, methodLabel='A'; end
        if strictMode==3, methodLabel='B'; end

        extra = '';
        if ~isempty(statusLine)
            extra = sprintf('\n%s',statusLine);
        end

        info.String = sprintf( ...
            't = %.1f / %.1f s   |   Frame %d / %d   |   View: %s (%s)\nBaseline: %g–%g %s   |   Editor: %s   |   AUTO method: %s%s', ...
            t,Tmax,frame,nFrames,vm,ms,...
            baseline.start,baseline.end,baseline.mode,...
            em,methodLabel,extra);

        fpsValue.String   = sprintf('%.1f',fps);
        frameValue.String = sprintf('%d / %d',frame,nFrames);

        % enable percentile slider only if method B
        if strictMode == 3
            set(percSlider,'Enable','on');
        else
            set(percSlider,'Enable','off');
        end
    end

%% ---------------- PLAYBACK ----------------
    function setFPS(v)
        fps = v;
    end

    function scrub(v)
        playing=false;
        playBtn.Value=0;
        playBtn.String='Play';
        frame=v;
        render();
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
        frame=1;
        frameSlider.Value=1;
        playing=true;
        playBtn.Value=1;
        playBtn.String='Pause';
    end

    function saveVideo(~,~)
        [f,p] = uiputfile('*.mp4','Save fUSI video');
        if isequal(f,0), return; end
        vid = VideoWriter(fullfile(p,f),'MPEG-4');
        vid.FrameRate = fps;
        vid.Quality   = 95;
        open(vid);
        for k = 1:nFrames
            frame = k;
            render(); drawnow;
            writeVideo(vid,getframe(ax));
        end
        close(vid);
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

    for k = 1:nFrames
        mask(:,:,k) = refMask;
    end

    statusLine = sprintf('Mask from frame %d applied to ALL frames.', frame);
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

        if applyToAllFrames
            if strcmp(paintMode,'add')
                for k=1:nFrames
                    mask(:,:,k) = mask(:,:,k) | brush;
                end
            else
                for k=1:nFrames
                    mask(:,:,k) = mask(:,:,k) & ~brush;
                end
            end
        else
            if strcmp(paintMode,'add')
                mask(:,:,frame) = mask(:,:,frame) | brush;
            else
                mask(:,:,frame) = mask(:,:,frame) & ~brush;
            end
        end

        statusLine='';
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
            % Apply automatic mask to ALL frames
            for k = 1:nFrames
                sig = PSC(:,:,k);
                auto = autoThreshold(sig, strictMode, percentileKeep);

                % simple cleanup: keep only sufficiently large regions
                auto = cleanupMask(auto);

                mask(:,:,k) = auto;
            end
            statusLine = sprintf('Auto mask applied to ALL frames (method %s).', autoMethodLabel());
        else
            % Apply only current frame
            sig = PSC(:,:,frame);
            auto = autoThreshold(sig, strictMode, percentileKeep);
            auto = cleanupMask(auto);
            mask(:,:,frame) = auto;
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

    function fillRegionAtSeed(x0,y0)

        if playing
            statusLine='Fill disabled while playing.';
            render(); return;
        end
        if ~editorMode
            statusLine='Enable Editor before Fill.';
            render(); return;
        end

        sig = PSC(:,:,frame);

        region = regionGrowSimilarity(sig, y0, x0, fillWindowR, fillSigmaFactor, fillMaxPixels);

        if ~any(region(:))
            statusLine='Fill found no similar region.';
            render();
            return;
        end

        if applyToAllFrames
            for k=1:nFrames
                mask(:,:,k) = mask(:,:,k) | region;
            end
        else
            mask(:,:,frame) = mask(:,:,frame) | region;
        end

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
        [f,p] = uiputfile('*.mat','Save mask');
        if isequal(f,0), return; end
        out.mask = mask;
        out.maskIsInclude = maskIsInclude;
        out.par = par;
        out.TR = TR;
        out.baseline = baseline;
        save(fullfile(p,f),'-struct','out','-v7.3');
        statusLine='Mask saved.';
        render();
    end

    function saveMaskedPSCMat(~,~)
        [f,p] = uiputfile('*.mat','Save masked PSC');
        if isequal(f,0), return; end
        PSCm = PSC;
        for k=1:nFrames
            if maskIsInclude
                PSCm(:,:,k)=PSC(:,:,k).*single(mask(:,:,k));
            else
                PSCm(:,:,k)=PSC(:,:,k).*single(~mask(:,:,k));
            end
        end
        out.PSC_masked = PSCm;
        out.mask = mask;
        out.par = par;
        save(fullfile(p,f),'-struct','out','-v7.3');
        statusLine='Masked PSC saved.';
        render();
    end

%% ---------------- CLOSE GUI FUNCTION ----------------
end   % play_fusi_video_final

%% ================= END OF FILE =================
