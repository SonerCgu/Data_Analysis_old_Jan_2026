%% ============================================================
%  fusi_video_soner.m
%  FINAL fUSI VIDEO ANALYSIS GUI (MATLAB 2017b)
%
%  FIX:
%   - Save MP4 now correctly uses uiputfile (save dialog)
%
%  NOTHING ELSE CHANGED
% ============================================================

clearvars; clc;

%% ===================== USER PARAMETERS =====================
par.interpol = 6;
par.LPF = 0.4;
par.HPF = 0;

par.gaussSize = 5;
par.gaussSig  = 0.8;

par.conectSize = 5;
par.conectLev  = 30;

par.previewCaxis = [0 60];

baseline.mode  = 'sec';
baseline.start = 0;
baseline.end   = 60;

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

function play_fusi_video_final(PSC,bg,par,fps,maxFPS,TR,Tmax,baseline)

[nz,nx,nFrames] = size(PSC);

fig = figure('Color','k','Position',[100 100 1100 850]);

ax  = axes('Parent',fig,'Position',[0.05 0.12 0.75 0.70]);
axis(ax,'off','image');
img = image(ax,zeros(nz,nx,3,'single'));

Nc = 128;
mapA = hot(Nc); mapA(1,:)=0;
bgN = (bg+68)/(68-5);
bgRGB = ind2rgb(round(bgN*127),gray(128));

%% ---------- TOP LEFT INFO ----------
info = uicontrol('Style','text','Units','normalized',...
    'Position',[0.01 0.88 0.55 0.08],...
    'ForegroundColor','w','BackgroundColor','k',...
    'FontName','Courier','FontSize',13,...
    'HorizontalAlignment','left');

%% ---------- TOP RIGHT CONTROLS ----------
fpsLabel = uicontrol('Style','text','Units','pixels',...
    'Position',[780 800 100 20],...
    'String','FPS:',...
    'ForegroundColor','w','BackgroundColor','k',...
    'FontSize',12,...
    'HorizontalAlignment','right');

fpsValue = uicontrol('Style','text','Units','pixels',...
    'Position',[885 800 80 20],...
    'ForegroundColor','w','BackgroundColor','k',...
    'FontSize',12);

fpsSlider = uicontrol('Style','slider',...
    'Min',1,'Max',maxFPS,'Value',fps,...
    'Units','pixels',...
    'Position',[780 775 300 22],...
    'Callback',@(s,~) setFPS(s.Value));

frameLabel = uicontrol('Style','text','Units','pixels',...
    'Position',[780 740 100 20],...
    'String','Frame:',...
    'ForegroundColor','w','BackgroundColor','k',...
    'FontSize',12,...
    'HorizontalAlignment','right');

frameValue = uicontrol('Style','text','Units','pixels',...
    'Position',[885 740 120 20],...
    'ForegroundColor','w','BackgroundColor','k',...
    'FontSize',12);

frameSlider = uicontrol('Style','slider',...
    'Min',1,'Max',nFrames,'Value',1,...
    'Units','pixels',...
    'Position',[780 715 300 22],...
    'Callback',@(s,~) scrub(round(s.Value)));

%% ---------- FOOTER ----------
uicontrol('Style','text','Units','pixels',...
    'Position',[10 10 520 24],...
    'String','fUSI Video Analysis — Soner C., MPI for Biological Cybernetics, 2026',...
    'ForegroundColor',[0.7 0.7 0.7],...
    'BackgroundColor','k',...
    'HorizontalAlignment','left',...
    'FontSize',11);

%% ---------- BUTTONS ----------
playBtn = uicontrol('Style','togglebutton','String','Play',...
    'Units','pixels','Position',[760 20 95 40],...
    'FontSize',13,'Callback',@playPause);

uicontrol('Style','pushbutton','String','Replay',...
    'Units','pixels','Position',[865 20 95 40],...
    'FontSize',13,'Callback',@replayVid);

uicontrol('Style','pushbutton','String','Save MP4',...
    'Units','pixels','Position',[970 20 105 40],...
    'FontSize',13,'Callback',@saveVideo);

%% ---------- COLORBAR ----------
caxis(par.previewCaxis);
cbar = colorbar('Position',[0.90 0.18 0.025 0.58]);
colormap(cbar,'hot');
cbar.Color = 'w';
cbar.FontSize = 13;
cbar.Label.String = 'Percent Signal Change (%)';
cbar.Label.FontSize = 14;
cbar.Label.Color = 'w';

frame = 1; 
playing = false;

render();

while ishandle(fig)
    if playing
        frame=min(frame+1,nFrames);
        frameSlider.Value = frame;
        render();
        pause(1/max(fps,0.01));
    else
        pause(0.02);
    end
end

    function render()
        A=(PSC(:,:,frame)-par.previewCaxis(1))/diff(par.previewCaxis);
        A=max(0,min(1,A));
        img.CData = ind2rgb(round(A*(Nc-1)),mapA)+bgRGB;

        t=(frame/par.interpol)*TR;
        info.String = sprintf( ...
            't = %.1f / %.1f s   |   Frame %d / %d\nBaseline: %g–%g %s', ...
            t,Tmax,frame,nFrames,...
            baseline.start,baseline.end,baseline.mode);

        fpsValue.String   = sprintf('%.1f',fps);
        frameValue.String = sprintf('%d / %d',frame,nFrames);
    end

    function setFPS(v), fps = v; end

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
end
