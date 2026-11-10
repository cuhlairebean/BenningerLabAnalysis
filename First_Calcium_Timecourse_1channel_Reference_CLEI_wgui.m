%% Extracting Single Cell Timecourse
% 04.04.2023 Claire H Levitt
% 10.17.26 Evelyn Ibarra

% This code is designed to extract timecourses of fluorescent signals
% For multiple channels: skip down to 'Masking Islet > For Colocalization'
% and change the parameters as needed, it also normalizes and generates
% matrices per islet


% NEED: filename.czi and filename_pre.czi for this code to work
%_pre is a high quality image acquire

%% FIRST TIME USING THIS CODE?
% Run this script and choose "New user" in the first dialog. It downloads both GitHub repos
% (BenningerLabUniversalFunctions + BenningerLabAnalysis) into ~/Documents/GitHub, the latest
% Bio-Formats jar, and the required MATLAB toolboxes (Image Processing, Statistics and Machine
% Learning). Requires git (https://git-scm.com/downloads) and a MATLAB license with those toolboxes.
% Then restart MATLAB, run the script again and choose "Other".
% To get bug fixes later: choose "New user" again (it pulls the latest code), or run
% `git pull origin main` in both repo folders.
%% SETTINGS (edit here)
% RERUN MODE: set rerun = true to skip all image processing/clicking and rebuild the
% normalized/frequency outputs from files saved by a previous run (no popups).
rerun             = true;
rerunFolder       = 'out';      % folder containing the previous outputs, e.g. '/Volumes/.../output'
rerunFilename     = '10-29-26_5mMPC_d1';      % the name typed at "Input Filename" in the original run
rerunMultiplyNeg1 = true;   % same choice as the original run
rerunUser         = 'Evelyn'; % 'Evelyn' or 'Other'
rerunDrug         = 'KCl';      % Evelyn extras (only used if rerunUser = 'Evelyn')
rerunSex          = 'F';
rerunAge          = NaN;
% --- WHAT TO CHANGE IN THE RERUN (all optional; leave empty = redo everything with the saved choices) ---
rerunIslets       = [];      % islets to reprocess, e.g. [2]. [] = all islets
rerunReviewIslets = [];      % islets whose traces open in the click-to-remove window, e.g. [2]
rerunRemoveCells  = {};      % remove specific cells without a window: {islet, [cell numbers]; ...}
                             %   e.g. {2,[5 7]; 3,12}  = cells 5 and 7 of islet 2, cell 12 of islet 3
rerunRestoreIslets = [];     % islets whose previously removed cells are put back first, e.g. [2]
% Removals are ADDED to the ones saved by the previous run and saved again (_RemovedCells.mat).
% Only the islets in rerunIslets are re-exported; every output file of those islets is overwritten.

%% HOUSEKEEPING
close all
clearvars -except rerun rerunFolder rerunFilename rerunMultiplyNeg1 rerunUser rerunDrug rerunSex rerunAge rerunIslets rerunReviewIslets rerunRemoveCells rerunRestoreIslets
clc
% Who is running the analysis? 'Evelyn' also exports the Analysis_dataset_detrended.xlsx workbook
if rerun
    user = rerunUser;
else
    user = questdlg('Who is running this analysis?','User','Evelyn','Other','New user','Evelyn');
    if isempty(user), user = 'Other'; end
end

%% FIRST-TIME SETUP (only when 'New user' was chosen)
% Downloads/updates both GitHub repos into ~/Documents/GitHub, downloads the Bio-Formats jar to
% ~/Documents/MATLAB if missing, and installs the required toolboxes if they are missing
% (needs git installed and a MATLAB license that includes the toolboxes).
if strcmp(user,'New user')
    homeDir = char(java.lang.System.getProperty('user.home'));
    % 1) repos (clone if missing, otherwise pull the latest bug fixes)
    gh = fullfile(homeDir,'Documents','GitHub');
    if ~exist(gh,'dir'), mkdir(gh); end
    for r = {'BenningerLabUniversalFunctions','BenningerLabAnalysis'}
        dest = fullfile(gh,r{1});
        if exist(dest,'dir')
            system(['git -C "' dest '" pull origin main']);
        else
            system(['git clone https://github.com/cuhlairebean/' r{1} '.git "' dest '"']);
        end
    end
    % 2) newest Bio-Formats jar (the one in the repo is older and fails on some newer .czi files)
    mlDir = fullfile(homeDir,'Documents','MATLAB');
    if ~exist(mlDir,'dir'), mkdir(mlDir); end
    jarPath = fullfile(mlDir,'bioformats_package.jar');
    if ~exist(jarPath,'file')
        try
            websave(jarPath,'https://downloads.openmicroscopy.org/bio-formats/latest/artifacts/bioformats_package.jar');
        catch ME
            warning('Could not download Bio-Formats (%s). The older jar from the repo will be used.',ME.message);
        end
    end
    % 3) toolboxes
    if ~license('test','Image_Toolbox') || ~license('test','Statistics_Toolbox')
        try
            arch = lower(computer('arch')); mpmExe = fullfile(tempdir,'mpm');
            if startsWith(arch,'win'), mpmExe = [mpmExe '.exe']; mpmUrl = ['https://www.mathworks.com/mpm/' arch '/mpm.exe'];
            else, mpmUrl = ['https://www.mathworks.com/mpm/' arch '/mpm']; end
            websave(mpmExe,mpmUrl);
            if ~startsWith(arch,'win'), system(['chmod +x "' mpmExe '"']); end
            system(sprintf(['"%s" install --release=R%s --destination="%s" --products ' ...
                'Image_Processing_Toolbox Statistics_and_Machine_Learning_Toolbox'], ...
                mpmExe,version('-release'),matlabroot));   % may ask for admin rights
        catch ME
            warning(['Toolbox install failed (%s). Install "Image Processing Toolbox" and ' ...
                '"Statistics and Machine Learning Toolbox" from Home > Add-Ons > Get Add-Ons.'],ME.message);
        end
    end
    msgbox(['Setup finished. Please RESTART MATLAB (so the new Bio-Formats jar and toolboxes load), ' ...
        'then run this script again and choose "Other".'],'Setup complete');
    return   % stop here: restart MATLAB, then run again
end
analysisDatasetPath = '/Users/ibarevel/Desktop/CaMPARI/Calcium Imaging/Analysis/November13_AllCaMPARI_Sept-Nov/Analysis_dataset_detrended.xlsx'; % Evelyn only
%% Paths for code functions and packages (found automatically - no editing needed)
% Looks for the two cloned repos next to each other (e.g. ~/Documents/GitHub/).
analysisDir = fileparts(mfilename('fullpath'));
if isempty(analysisDir), analysisDir = pwd; end
univDir = fullfile(fileparts(analysisDir),'BenningerLabUniversalFunctions');
if ~exist(univDir,'dir')
    univDir = fullfile(getenv('HOME'),'Documents','GitHub','BenningerLabUniversalFunctions');
end
if ~exist(univDir,'dir')
    error(['BenningerLabUniversalFunctions not found. Clone it next to BenningerLabAnalysis ' ...
        '(see ONE-TIME SETUP at the top of this file).']);
end
addpath(univDir);
addpath(analysisDir);
% Bio-Formats jar: a copy in ~/Documents/MATLAB is used if present, otherwise the (current) copy shipped
% in BenningerLabUniversalFunctions. Newest jar: https://www.openmicroscopy.org/bio-formats/downloads/
homeDir = char(java.lang.System.getProperty('user.home'));
bfJar = fullfile(homeDir,'Documents','MATLAB','bioformats_package.jar');
if ~exist(bfJar,'file'), bfJar = fullfile(univDir,'bioformats_package.jar'); end
% A different Bio-Formats/LOCI jar already loaded in this MATLAB session (e.g. the older repo copy from an
% earlier run) wins over the one added here and can't be unloaded -> MATLAB must be restarted.
cpAll = javaclasspath('-all');
oldJars = cpAll(contains(cpAll,'bioformats_package') | contains(cpAll,'loci_tools'));
oldJars = oldJars(~strcmp(oldJars,bfJar));
if ~isempty(oldJars)
    error(['An older Bio-Formats jar is already loaded in this MATLAB session:\n  %s\n' ...
        'It cannot be unloaded and fails on some .czi files. Please RESTART MATLAB and run this script again.'], ...
        strjoin(oldJars,sprintf('\n  ')));
end

if ~any(strcmp(cpAll,bfJar))
    javaaddpath(bfJar);
end
fprintf('Bio-Formats jar: %s\n',bfJar);
% Required toolboxes
if ~license('test','Image_Toolbox') || ~license('test','Statistics_Toolbox')
    warning('Image Processing Toolbox and/or Statistics and Machine Learning Toolbox not licensed - some steps will fail.');
end


%% RE-RUN: reload everything from a previous run (rerun = true at the top of the script)
% Skips loading the .czi, drawing islets/cells, time course extraction and all questions. Reads the saved
% outputs (CaWaveForm, CellIDs, PeriodIndices, IsletTC, RemovedCells, TimeVec) from 'rerunFolder' for the
% name 'rerunFilename', then goes straight to the detrend / normalize / plot / frequency / export steps.
if rerun
    output_dir = rerunFolder;
    filename = rerunFilename;
    user = rerunUser;
    if ~exist(output_dir,'dir'), error('rerunFolder does not exist: %s',output_dir); end
    % filename left blank -> infer it from the files in the folder
    if isempty(filename)
        wfAll = dir(fullfile(output_dir,'*_islet*_CaWaveForm.mat'));
        pre = unique(regexprep({wfAll.name},'_islet\d+_CaWaveForm\.mat$',''));
        if numel(pre) == 1
            filename = pre{1};
        else
            error('Re-run: set rerunFilename. Found these names in the folder: %s',strjoin(pre,', '));
        end
    end
    wfs = dir(fullfile(output_dir,[filename '_islet*_CaWaveForm.mat']));
    groupIDs = [];
    for q = 1:numel(wfs)
        tokQ = regexp(wfs(q).name,['^' regexptranslate('escape',filename) '_islet(\d+)_CaWaveForm\.mat$'],'tokens','once');
        if ~isempty(tokQ), groupIDs(end+1) = str2double(tokQ{1}); end %#ok<SAGROW>
    end
    groupIDs = sort(unique(groupIDs));
    if isempty(groupIDs)
        error(['Re-run: no "%s_islet#_CaWaveForm.mat" files found in:\n  %s\n' ...
            'Check rerunFolder and rerunFilename (the name you typed at "Input Filename" in the original run).'],filename,output_dir);
    end
    nIslets = max(groupIDs);
    groupSuffix = arrayfun(@(k) sprintf('_islet%d',k),groupIDs,'UniformOutput',false);
    IsletPeriods_all = cell(1,nIslets);
    IsletTC_all = cell(1,nIslets);
    for k = groupIDs
        base = fullfile(output_dir,sprintf('%s_islet%d',filename,k));
        % Periods: (1) _PeriodIndices.mat (any variable that is a struct with name/start/end),
        % (2) _PeriodRanges.csv, (3) the 'Conditions' sheet of _normalized.xlsx. Older runs may only have some of these.
        P = [];
        if exist([base '_PeriodIndices.mat'],'file')
            Sp = load([base '_PeriodIndices.mat']);
            fn = fieldnames(Sp);
            for q = 1:numel(fn)
                v = Sp.(fn{q});
                if isstruct(v) && all(isfield(v,{'name','start','end'})), P = v; break; end
            end
        end
        if isempty(P) && exist([base '_PeriodRanges.csv'],'file')
            try
                T = readtable([base '_PeriodRanges.csv']);
                P = struct('name',cellstr(string(T{:,1})),'start',num2cell(T{:,2}),'end',num2cell(T{:,3}));
                P = P(:)';
            catch
            end
        end
        if isempty(P) && exist([base '_normalized.xlsx'],'file')
            try
                T = readtable([base '_normalized.xlsx'],'Sheet','Conditions');
                nm = cellstr(string(T{:,1}));
                nm(startsWith(nm,'BASELINE')) = {'Baseline'};
                P = struct('name',nm,'start',num2cell(T{:,2}),'end',num2cell(T{:,3}));
                P = P(:)';
            catch
            end
        end
        if isempty(P)
            error(['Re-run: could not recover the selected periods for islet %d (looked for _PeriodIndices.mat, ' ...
                '_PeriodRanges.csv and the Conditions sheet of _normalized.xlsx in %s).'],k,output_dir);
        end
        IsletPeriods_all{k} = P;
        tcf = [base '_IsletTC.mat'];
        if exist(tcf,'file')
            Stc = load(tcf); fn = fieldnames(Stc);
            if isfield(Stc,'IsletTC'), IsletTC_all{k} = Stc.IsletTC; else, IsletTC_all{k} = Stc.(fn{1}); end
        else
            IsletTC_all{k} = [];     % only used for plotting in the first run
        end
    end
    OverlayStack = [];
    fprintf('Re-run of "%s" from %s (islets: %s)\n',filename,output_dir,num2str(groupIDs));
else

%% Paths for data to analyze                                                                     
% Data folder = the folder currently open in MATLAB (the "Current Folder" bar), added to the path.
% So just open your data folder in MATLAB and set 'filename' below - no path to type.
% To use a different folder instead, uncomment the next line and set it:
% filepath = '/Volumes/ESD-USB/Calcium_imaging/10_29_26/';
if ~exist('filepath','var') || isempty(filepath), filepath = pwd; end
addpath(filepath);


%% LOAD CALCIUM IMAGE FILE

% comment out whichever file upload method you aren't using

% 1. SELECT MAUALLY
% 1. SELECT MAUALLY  (commented out: bfopen cannot open the text 'Select Calcium Imaging File')
%  disp('Select Calcium Imaging File');
%  R = bfopen('Select Calcium Imaging File');

% 2. INPUT FILE NAME                                                                            %%% CHANGE HERE  %%%%%%%%%

filename = ['10-29-26_5mMPC__d1.czi'];
% raw time-series file; 'filename' is later replaced by the output name you type.
% The file is looked for in 'filepath' (above) first, then in the current folder / as a full path.
if exist(fullfile(filepath,filename),'file') == 2
    czi_file = fullfile(filepath,filename);
elseif exist(filename,'file') == 2
    czi_file = filename;
else
    error(['Could not find the .czi file "%s".\nLooked in: %s\nand the current folder: %s\n' ...
        'Check that the drive is connected and that "filepath" and "filename" above are correct.'], ...
        filename,filepath,pwd);
end
czi_info = dir(czi_file);
fprintf('Opening %s (%.1f MB)\n',czi_file,czi_info(1).bytes/1e6);
if czi_info(1).bytes < 1e6
    warning('This .czi is very small (%.1f MB) - it may be incomplete/corrupted (interrupted acquisition or copy).',czi_info(1).bytes/1e6);
end

%https://www.openmicroscopy.org/bio-formats/downloads/

try
    R = bfopen(czi_file); % Uses bfopen program to open .czi/.lsm image files
catch ME
    error(['Bio-Formats could not open:\n  %s\n\nMost common causes:\n' ...
        '  1) An older Bio-Formats jar is loaded in this MATLAB session -> RESTART MATLAB and run again.\n' ...
        '  2) The .czi is incomplete/corrupted (interrupted acquisition or copy) -> open it in ZEN/Fiji to check.\n' ...
        '  3) Not enough Java memory for a very large file -> MATLAB Preferences > General > Java Heap Memory.\n\n' ...
        'Original error: %s'],czi_file,ME.message);
end

% Number of channels in file - important if its multichannel or ref.channel
                                                                                                  %  %%% CHANGE HERE  %%%%%%%%%
howmanychannel =2;

pics=R{1};
pics=pics(:,1);

for i=1:length(pics) %making the image into an array 
    IMG(:,:,i)=pics{i};
end

images=double(IMG); % converts images to double precision


%%%%%%%%%%%%
% Uncomment out the channel you want to quantify, comment out the rest.                         %%% CHANGE HERE  %%%%%%%%%

%channel 1 Red/fura
%images = images(:,:,1:howmanychannel:end);

%Channel 2 Lux/Green
images = images(:,:,2:howmanychannel:end);

%Channel 3 Green
%images = images(:,:,3:howmanychannel:end);

% Channel 4 Fura
%images = images(:,:,3:howmanychannel:end);

%%%%%%%%%%%



RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable

clear pics R IMG;
clear omeMeta;

% SELECT LOCATION FOR SAVING FILES
disp('Select Folder for Output Directory');
output_dir = uigetdir;
filename = input('Input Filename  ','s');

%% CODE TO MODIFY TIME COURSES
% This block converts the 3rd dimension of array into timecourse for
% later reference

st=1; % STARTING FRAME: Change this per sample if needed                                
%st = 100;
[xx, yy, zz] = size(images(:,:,:));
% ENDING FRAME: Change this per sample if needed
ed=zz;
%ed = 20;

T=st:1:ed;

images=images(:,:,st:ed);
% Define the start and end frames for your reference image
start_frame =336;                                                                                         %%% CHANGE HERE  %%%%%%%%%
end_frame = 1000;

% Check to ensure the selected range is valid
if start_frame < 1 || end_frame > size(images, 3) || start_frame > end_frame
    error('Invalid frame range selected. Please choose frames within the timeseries bounds.');
end

% Select the specified range of frames
ref = images(:,:,start_frame:end_frame);

% Create the max intensity projection from the selected frames
MaxIntensityMap = max(ref, [], 3);

% Convert the map to a displayable grayscale image
ref_pic = (mat2gray(MaxIntensityMap));
%% DECLARING IMAGE PROPERTIES

sx=size(images,1);
sy=size(images,2);
sz=size(images,3);

for i=1:sz
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]); %applies filter to clean up images
end

ImAv = mean(images,3); %compresses all frames into single array of intensities
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
HSV(:,:,1) = 0.3333;%converts image to green image
RGB = hsv2rgb(HSV); %converts to rgb image

OGFig = figure(1);
imshow(RGB);

%%Getting rid of Background using JD code
GrayFig = figure('Name','Remove Background');
%ImGray = rgb2gray(RGB); %converts image to gray
ImGray = mat2gray(ImAv); % Use the normalized average intensity image directly
imagesc(ImGray) %displays in default colormap; easier to see cell borders
colormap(GrayFig,'jet')
disp('Select Background Signal')
bkgrndC=createMask(imfreehand); %allows user to select background area to be removed
close(GrayFig);
BackC=ImGray.*bkgrndC; %multiplies mask by the image
BackC(~logical(BackC))=nan; %takes inverse of masked image and sets it to nan

BackC_nonan = BackC(~isnan(BackC)); % Select only the non-NaN values from BackC
thresh = max(BackC_nonan);           % Find the maximum value among those
%thresh = nanmax(nanmax(BackC)); %sets threshold based on the max of values after nan's are removed
close(OGFig);

[ImgNon, ~] = RemovingAreaswNoSignal(ImGray,thresh); %Uses function from JD to remove background based on threshold
ImGray(ImgNon)=0; %sets the removed areas to 0
NoSigFig = figure('Name','Background Removed; Draw ROI around Entire Islet');
imagesc(ImGray);
colormap(NoSigFig, 'jet');
clear OGFig GrayFig;

%% MASKING ISLET

% High-resolution overlay: the matching '<name>_pre.czi' next to the time-series file.
% It is automatically stretched onto the image grid (e.g. 512x512).
OverlayStack = []; OverlayChannelNames = {};
[czi_dir, czi_name, czi_ext] = fileparts(czi_file);
% Accepts <name>_pre.czi, <name>pre.czi, <name>_pre1.czi, <name>pre1.czi ... in the same folder.
preFile = '';
for pc = {'_pre','pre','_pre1','pre1','_Pre','Pre'}
    cand = fullfile(czi_dir, [czi_name pc{1} czi_ext]);
    if exist(cand,'file') == 2, preFile = cand; break; end
end
if isempty(preFile)
    disp(['No "' czi_name '_pre.czi" (or "' czi_name 'pre.czi") found in ' czi_dir]);
    [ovFile, ovPath] = uigetfile('*.czi','Select high-res _pre .czi (Cancel = no overlay)',[czi_dir filesep]);
    if ~isequal(ovFile,0), preFile = fullfile(ovPath,ovFile); else, preFile = ''; end
end
if ~isempty(preFile)
    try
        OverlayStack = loadCziChannels(preFile);
        disp(['Loaded overlay: ' preFile ' (' num2str(size(OverlayStack,3)) ' channels)']);
    catch ME
        warning('Could not load overlay (%s). Continuing without it.',ME.message);
    end
end

% User draws an ROI around EACH islet. Islets are numbered in the order drawn
% (islet1, islet2, ...). Every islet gets its own output files with '_islet#'.
% To rerun with the updated code you can instead IMPORT a previous islet mask.
isletAnswer = questdlg('Import a previous islet mask or draw new islets?', ...
    'Islet masks','Import islet mask','Draw new islets','Draw new islets');
IsletMask = [];
if strcmp(isletAnswer,'Import islet mask')
    [imf, imp] = uigetfile('*_IsletMask.mat','Select previous islet mask (*_IsletMask.mat)',[output_dir '/']);
    if ~isequal(imf,0)
        S_imp = load(fullfile(imp,imf));
        if isfield(S_imp,'IsletMask') && isequal(size(S_imp.IsletMask),[xx yy])
            IsletMask = S_imp.IsletMask;
            disp(['Imported islet mask: ' imf]);
        else
            warning('That file has no IsletMask matching this image size - drawing new islets instead.');
        end
    end
end
close(NoSigFig); % old static background-removed figure; the GUI shows the same image
if isempty(IsletMask)
    disp('Draw ROI around each islet (press DONE when finished)');
    [IsletMask, isletInfo] = drawROIsGUI(ImGray, ...
        'Name','Islet', ...
        'Title','Background Removed; Draw ROI around EACH islet (numbered in order drawn)', ...
        'Colormap','jet', ...
        'Overlay',OverlayStack, ...
        'ChannelNames',OverlayChannelNames);
end
nIslets = max(IsletMask(:));
if nIslets == 0
    error('No islets were drawn.');
end
save([output_dir '/' filename '_IsletMask.mat'],'IsletMask');

sz_FR=ed-st+1;
IsletTC_all = cell(1,nIslets);
periodNames = {'Low glucose','High glucose','Drug','Baseline','First phase','Second phase'};
IsletPeriods_all = cell(1,nIslets);
sharedCond = [];   % [start end] of Low glucose / High glucose / Drug, shared by all islets of this file
for k = 1:nIslets
    isletTag = sprintf('_islet%d',k);
    StartMaskStack = images.*(IsletMask==k);
    IsletTC = shiftdim(mean(mean(StartMaskStack)),2);
    IsletTC_all{k} = IsletTC;
    save([output_dir '/' filename isletTag '_IsletTC.mat'],'IsletTC'); % saved first, so it exists even if you stop later
    IsletTCfig = figure('Name',sprintf('Whole Islet Time Course - islet%d',k));
    plot(IsletTC); %plots intensity timecourse of entire islet so user can identify first responder and wave origin ranges
    xlabel('Frame index'); ylabel('Mean intensity');
    title(sprintf('%s islet%d',filename,k),'Interpreter','none');

    % Highlight low glucose, high glucose, drug, baseline, first phase, second phase.
    % Click START and END for each; ranges may overlap. You can re-select any of them.
    % Or import a previous selection (.mat / .csv); it is shown on the plot and you can still re-select.
    % Low glucose, High glucose and Drug are the same for every islet of this file, so they are
    % selected (clicked or typed) ONCE, on islet 1, and reused for the other islets.
    % Baseline, First phase and Second phase are selected separately for every islet.
    initialPeriods = [];
    prevFile = [output_dir '/' filename isletTag '_PeriodIndices.mat'];
    if exist(prevFile,'file'), impHint = ' (a previous file for this islet exists)'; else, impHint = ''; end
    impChoice = questdlg(sprintf('Islet %d: import previous period selections?%s',k,impHint), ...
        'Import period indices','Import previous','Select new','Select new');
    if strcmp(impChoice,'Import previous')
        if exist(prevFile,'file'), startFile = prevFile; else, startFile = [output_dir '/']; end
        [pf, pp] = uigetfile({'*_PeriodIndices.mat;*_PeriodRanges.csv','Period files'}, ...
            'Select previous period file',startFile);
        if ~isequal(pf,0)
            try
                initialPeriods = loadPreviousPeriods(fullfile(pp,pf),periodNames);
            catch ME
                warning('Could not import periods (%s). Select them manually.',ME.message);
            end
        end
    end
    if k == 1 || isempty(sharedCond)
        pickMask = [];                                   % islet 1: all six periods
    else
        if isempty(initialPeriods), initialPeriods = nan(numel(periodNames),2); pickMask = false(1,numel(periodNames)); pickMask(4:6) = true;
        else, pickMask = false(1,numel(periodNames)); end
        initialPeriods(1:3,:) = sharedCond;              % same low/high/drug as islet 1
    end
    [IsletPeriods, PeriodIndices] = selectTimeRanges(IsletTC, periodNames, ...
        sprintf('%s islet%d',filename,k), IsletTCfig, initialPeriods, pickMask);
    IsletPeriods_all{k} = IsletPeriods;
    if k == 1, sharedCond = [[IsletPeriods(1:3).start]', [IsletPeriods(1:3).end]']; end

    % Export right away (one CSV per islet; one column per period holding its frame indices)
    colNames = matlab.lang.makeValidName(strrep(periodNames,' ','_'));
    writetable(array2table(PeriodIndices,'VariableNames',colNames), ...
        [output_dir '/' filename isletTag '_PeriodIndices.csv']);
    rangeTbl = table(periodNames(:),[IsletPeriods.start]',[IsletPeriods.end]', ...
        'VariableNames',{'Period','StartIndex','EndIndex'});
    writetable(rangeTbl,[output_dir '/' filename isletTag '_PeriodRanges.csv']);
    save([output_dir '/' filename isletTag '_PeriodIndices.mat'],'IsletPeriods','PeriodIndices','periodNames');
    saveas(IsletTCfig,[output_dir '/' filename isletTag '_IsletTC_periods.tif']);
    disp(['Saved period indices for ' filename isletTag]);
end

clear StartMaskStack


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For colocalization% - comment this out if you don't have multiple channels
 % disp('Select Colocalization File')
 % [file2, path2] = uigetfile('*.czi')
 % id2 = [path2 file2];
 % disp(id2);
 % HI_image = bfopen(id2); %Human Islet Image (HI)
 % 
 % red_cells = HI_image{1,1}{1,1}; 
 % red_cells_resize = imresize(red_cells,[xx yy]);
 % 
 % green_cells = HI_image{1,1}{2,1}; %calcium
 % green_cells_resize = imresize(green_cells,[xx yy]);
 % 
 % % ca_cells = HI_image{1,1}{2,1}; %green cells
 % % ca_cells_resize = imresize(ca_cells,[xx yy]);
 % 
 % % r = imfuse(ca_cells,red_cells); %red cells + calcium
 % % g = imfuse(ca_cells,green_cells); %green cells + calcium
 % 
 % b = imfuse(r,g);
 % 
 % figure, imshow(r)
 %  title ('RED overlay')
 % 
 %  figure, imshow(g)
 % title ('GREEN overlay')

%% Cell Masking
 %#User Draws ROIs around each cell within the islet
 %#ROIs are saved in "CellMask" array and called back to throughout analysis
 %#Follow the prompts in the Command Window

InitialRefFig = figure; % Give the figure a handle so we can apply colormap
imshow(ref_pic)
colormap(InitialRefFig, 'jet'); % Apply 'jet' colormap here
pause(3)
close(InitialRefFig) % Close the figure using its handle
clear InitialRefFig
%figure, imshow(ref_pic)
%pause(3)

%close all

CellMask = double(zeros(xx,yy));
numcells = 1;

answer = questdlg('Import Mask Set or Draw New Masks?', ...
    'Import Dialogue', ...
    'Import Mask Set','Draw New Masks','Draw New Masks');

switch answer
    case 'Import Mask Set'
        disp('Importing Mask Set')
        disp('Please Select Mask File');
        [Mask,mpath] = uigetfile('*.mat','Please select Mask file'); % User selects .mat mask file - only valid for MatLab generated files!!!
        addpath(mpath);
        dummyMask = struct2cell(load(Mask)); % Loads in Mask file
        CellMask = dummyMask{:,:}; % Assigns loaded mask file to CellMask Variable
        disp('Please Select Cell Number File');
        [num,npath] = uigetfile('*.mat','Please select Cell Number file',mpath); % User selects .mat cell number file - only valid for MatLab generated files!!!
        addpath(npath);
        dummyNum = struct2cell(load(num)); % Loads in Cell Number file
        numcells = dummyNum{:,:};


    case 'Draw New Masks'
        disp(answer)
        % Reusable GUI: zoom, brightness/contrast, translucent hi-res .czi overlay.
        % Each cell is checked against IsletMask and gets that islet's number.
        [CellMask, cellInfo] = drawROIsGUI(ref_pic, ...
            'Name','Cell', ...
            'Title','Draw ROIs Around Cells of Interest (islet outlines shown dashed)', ...
            'Colormap','jet', ...
            'Overlay',OverlayStack, ...
            'ChannelNames',OverlayChannelNames, ...
            'ParentMask',IsletMask, ...
            'ParentName','Islet');
        numcells = max(CellMask(:));

        save([output_dir '/' filename '_' 'Masks' '.mat'],'CellMask'); %Saves CellMask array
        save([output_dir '/' filename '_' 'CellNumber' '.mat'],'numcells'); %Saves numcells array
end

% Which islet does each cell fall in? (same number as the islet mask: 0 = none)
[CellIsletIdx, CellIsletOverlap] = roiParentAssignment(CellMask, IsletMask);
for c = 1:numcells
    if CellIsletIdx(c) > 0
        fprintf('Cell %d -> islet%d (%.0f%% inside)\n',c,CellIsletIdx(c),100*CellIsletOverlap(c));
    else
        fprintf('Cell %d is not inside any islet\n',c);
    end
end
CellMask_All = CellMask;
save([output_dir '/' filename '_CellIsletAssignment.mat'],'CellIsletIdx','CellIsletOverlap');

%% RELOAD ALL CHANNELS FOR QUANTIFICATION
% Reload the original data from the file to include all channels.

R_quant = bfopen(czi_file); % Re-open the file
pics_quant = R_quant{1};
pics_quant = pics_quant(:,1);
IMG_quant = zeros(size(R_quant{1}{1,1}, 1), size(R_quant{1}{1,1}, 2), length(pics_quant));
for i=1:length(pics_quant)
    IMG_quant(:,:,i)=pics_quant{i};
end
images_all_channels = single(IMG_quant); % This now holds all channel data
% Apply the start/end frame limits (st, ed) again
images_all_channels = images_all_channels(:,:,st:ed); 
frameTimes_all = getFrameTimes(R_quant, howmanychannel); % per-frame time stamps (s) from the .czi metadata
clear R_quant pics_quant IMG_quant;

% %% TIMECOURSE EXTRACTION FOR ALL CHANNELS
% % The masks were created on Channel 3 (or whichever was selected).
% % Now we quantify all channels using these masks.
% howmanychannel = 4; % Reconfirm the total number of channels
% sz_2 = ed-st+1; % Total number of frames after cropping
% PlotHandles = zeros(1,numcells);
% PlotLabels = cell(1,numcells);
% 
% % New structure to store timecourses for all channels
% CellTC_AllChannels = struct(); 
% 
% for channel_num = 1:howmanychannel
%     disp(['Processing Channel ' num2str(channel_num)]);
% 
%     % Select only the frames for the current channel
%     current_channel_images = images_all_channels(:,:,channel_num:howmanychannel:end);
% 
%     % Initialize storage for this channel's timecourses
%     % Select only the frames for the current channel
%     current_channel_images = images_all_channels(:,:,channel_num:howmanychannel:end);
% 
%     % FIX: Get the *actual* frame count for the current channel stack
%     [~, ~, current_sz] = size(current_channel_images); 
% 
%     % Initialize storage using the *actual* size
%     CellTC = zeros(current_sz, numcells);
%     TC = zeros(current_sz, 1); 
%     tic
%     for i = 1:numcells % Loop through each cell mask
%         TCMask = CellMask; 
%         TCMask(find(TCMask ~= i)) = 0; 
%         MaskedIMGstack = current_channel_images.*logical(TCMask); 
% 
%         % FIX: Use the actual size for the inner loop
%         for ii = 1:current_sz 
%             TCnoZero = MaskedIMGstack(:,:,ii);
%             TCnoZero = TCnoZero(TCnoZero>0);
%             TC(ii) = mean(TCnoZero);
%         end
%         CellTC(:,i) = TC; 
%     end
%     toc
% 
%     % Store the results for this channel
%     channel_name = ['Channel' num2str(channel_num) '_TC'];
%     CellTC_AllChannels.(channel_name) = CellTC;
% 
%     % Plotting traces for the current channel time course
%     TCFig = figure('Name',['Average Intensity Over Time - Channel ' num2str(channel_num)]);
%     plot(CellTC);
%     title(['Time Course for Channel ' num2str(channel_num)]);
%     saveas(TCFig,[output_dir '/' filename '_Channel' num2str(channel_num) '_TCPlot' '.tif']);
% 
%     % Save data for this channel
% %     save([output_dir '/' filename '_Channel' num2str(channel_num) '_CaWaveForm.mat'],'CellTC'); 
% end
% %
% The rest of the original code (Normalization, plotting masks) can remain,
% but remember to adjust the 'ANALYSIS: COMPUTE NORMALIZED RATIO' section
% to process the specific channel you want to analyze (e.g., Channel 2 for Lux/Green)
%% TIMECOURSE

sz_2=ed-st; % includes all frames for exporting of entire trace

PlotHandles = zeros(1,numcells);
PlotLabels = cell(1,numcells);
CellTC = zeros(sz_2,numcells);
TC = zeros(sz_2,1);

tic
for i = 1:numcells
    TCMask = CellMask; %Pulls in CellMask array
    TCMask(find(TCMask ~= i)) = 0; %Gets rid of all masks besides current one
    MaskedIMGstack = images.*logical(TCMask); %Applies current mask to all frames of timecourse
    %get rid of for loop
    for ii = 1:sz_2 
        TCnoZero = MaskedIMGstack(:,:,ii); %Pulls in current frame
        TCnoZero = TCnoZero(TCnoZero>0); %Accounts for any zeros from preallocation
        TC(ii) = mean(TCnoZero); %Calculates mean intensity of current frame
    end

    CellTC(:,i) = TC;   
    PlotLabels{i} = [' Cell ' num2str(i)]; 
end
toc
CellTC_All = CellTC;   % full matrix (all cells); the per-islet export below takes its columns from it



%% EXPORT PER ISLET
% Every islet gets its own output files with the suffix '_islet#'
% (cells that were not inside any islet get '_noIslet').
groupIDs = unique(CellIsletIdx(:))';   % only groups that actually contain cells
groupSuffix = cell(1,numel(groupIDs));

for g = 1:numel(groupIDs)
    k = groupIDs(g);
    if k > 0, sfx = sprintf('_islet%d',k); else, sfx = '_noIslet'; end
    groupSuffix{g} = sfx;

    cellIDs = find(CellIsletIdx == k)';      % global cell numbers in this islet (draw order)
    CellTC = CellTC_All(:,cellIDs);
    % (the optional x -1 is no longer applied here: it is done LAST, after detrend + normalize,
    %  in the PROCESSING section)

    % Masks for this islet only, renumbered 1..n (single variable -> importable)
    CellMask = zeros(size(CellMask_All));
    for j = 1:numel(cellIDs)
        CellMask(CellMask_All == cellIDs(j)) = j;
    end
    save([output_dir '/' filename sfx '_Masks.mat'],'CellMask');
    numcells = numel(cellIDs);
    save([output_dir '/' filename sfx '_CellNumber.mat'],'numcells');
    save([output_dir '/' filename sfx '_CellIDs.mat'],'cellIDs'); % local cell j = global cell cellIDs(j)

    % Plotting traces for this islet
    TCFig = figure('Name',['Average Intensity Over Time' sfx]);
    plot(CellTC);
    title([filename sfx],'Interpreter','none');
    saveas(TCFig,[output_dir '/' filename sfx '_TCPlot.tif']);
    save([output_dir '/' filename sfx '_CaWaveForm.mat'],'CellTC'); % CaWaveForm - Cell TC
end
CellMask = CellMask_All;   % restore full mask
numcells = size(CellTC_All,2);

% View Cells w labels (colored by islet)
Cent = nan(numcells,2);
for j = 1:numcells
    STATS = regionprops(CellMask == j,'Centroid'); % centroid of this mask
    if ~isempty(STATS), Cent(j,:) = STATS(1).Centroid; end
end
LabelFig = figure; imshow(ref_pic)
title(filename,'Interpreter','none')
hold on
isletCols = lines(max(nIslets,1));
for k = 1:nIslets
    Bk = bwboundaries(IsletMask == k);
    for b = 1:numel(Bk)
        plot(Bk{b}(:,2),Bk{b}(:,1),'--','Color',isletCols(k,:),'LineWidth',1.2);
    end
end
for c = 1:numcells
    if isnan(Cent(c,1)), continue; end
    if CellIsletIdx(c) > 0, tcol = isletCols(CellIsletIdx(c),:); else, tcol = [1 1 1]; end
    text(Cent(c,1),Cent(c,2),num2str(c),'Color',tcol); %Labels cells with their region number
end
saveas(LabelFig,[output_dir '/' filename '_CellLabels_byIslet.tif']);

disp('END!')

end   % if rerun ... else (everything from "Paths for data to analyze" to here is skipped when rerun = true)


%% PROCESSING: PER-CELL DETREND -> NORMALIZE TO DETRENDED BASELINE -> (x -1) AND EXPORT (one workbook per islet)
% On the WHOLE matrix (every frame of every cell):
%   1) detrend    per cell: trendline = slope*time + intercept (own regression), subtracted from that cell
%   2) normalize  the DETRENDED data: (D - mean_baseline)/mean_baseline  (baseline = the period you selected)
%   3) x -1       LAST, only if selected
% Output per islet:  <filename>_islet#_normalized.xlsx
% Low glucose, High glucose, Drug -> labels in the 'condition' column + per-period means.
% (First phase / Second phase are saved in the _PeriodIndices files for later analysis.)
if rerun
    runAns = 'Yes';
else
    runAns = questdlg('Detrend, normalize to baseline and export now?','Processing','Yes','No','Yes');
end
if ~strcmp(runAns,'Yes')
    disp('Processing skipped - no normalized export.');
else
    if rerun
        applyNeg1 = rerunMultiplyNeg1;     % set at the top of the script
    else
        negAns = questdlg('Multiply the normalized data by -1 (applied last)?','x -1','No','Yes','No');
        applyNeg1 = strcmp(negAns,'Yes');
    end
    % time stamps of the frames used (czi metadata, else frame interval, else frame index)
    timeVec = [];
    if rerun
        % re-run: time stamps saved by the previous run (<name>_TimeVec.mat), else the 'time' column of
        % that run's normalized workbook
        tvFile = fullfile(output_dir,[filename '_TimeVec.mat']);
        if exist(tvFile,'file')
            Stv = load(tvFile,'timeVec'); timeVec = Stv.timeVec(:);
        else
            xlTv = fullfile(output_dir,sprintf('%s%s_normalized.xlsx',filename,groupSuffix{find(groupIDs>0,1)}));
            if exist(xlTv,'file')
                try
                    Ttv = readtable(xlTv,'Sheet','data');
                    if ismember('time',Ttv.Properties.VariableNames), timeVec = Ttv.time(:); end
                catch
                end
            end
        end
        if isempty(timeVec)
            warning('No saved time stamps found for the re-run - using the frame index as time.');
        else
            fprintf('Re-run: loaded %d time stamps.\n',numel(timeVec));
        end
    elseif exist('frameTimes_all','var') && ~isempty(frameTimes_all) && numel(frameTimes_all) >= ed
        timeVec = frameTimes_all(st:ed);
    else
        % Paste the time stamps (e.g. a column copied from Excel) OR type one frame interval.
        % Separators can be new lines, tabs, spaces, commas or semicolons.
        while true
            ansT = inputdlg({sprintf(['Time stamps not found in the .czi.\n' ...
                'FROM FIJI: circle an ROI, then Image > Stacks > Plot Z-axis Profile >\n' ...
                'Data > Copy all Data, and PASTE it here (only the [sec] column is used).\n' ...
                'Or paste a column of time stamps in seconds from Excel,\n' ...
                'or type ONE number = frame interval in seconds. Blank = use frame index.\n' ...
                'This file needs at least %d rows (frames 1 to %d of the recording).'],ed,ed)}, ...
                'Time stamps',[14 70],{''});
            if isempty(ansT) || isempty(strtrim(strjoin(cellstr(ansT{1}),' '))), break; end
            rowsT = cellstr(ansT{1});
            numsT = [];
            numRows = 0;
            for rr = 1:numel(rowsT)
                % first number on the line only: for a Fiji "[sec]  Mean" table the Mean column
                % (and the header line) is ignored; any extra columns are ignored as well
                tok = regexp(strrep(rowsT{rr},',','.'),'[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?','match','once');
                if ~isempty(tok), numsT(end+1,1) = str2double(tok); numRows = numRows + 1; end %#ok<AGROW>
            end
            if numRows == 1                     % single line: a list separated by spaces/commas/semicolons, or one interval
                numsT = str2double(regexp(strrep(strjoin(rowsT,' '),',',' '),'[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?','match'));
                numsT = numsT(:);
            end
            nCzi = zz;                          % number of frames of this .czi (quantified channel)
            nPaste = numel(numsT);
            if isscalar(numsT)                  % one number = frame interval
                timeVec = ((st:ed)' - st) * numsT; break;
            elseif nPaste == nCzi || nPaste == ed-st+1
                if nPaste == nCzi, timeVec = numsT(st:ed); else, timeVec = numsT; end
                break;                          % counts agree
            end
            % ---- MISMATCH: the pasted time series does not have this .czi's number of frames ----
            mmVal = nPaste - nCzi;
            mmMsg = sprintf(['TIME SERIES MISMATCH for %s\n\n' ...
                '  .czi frames:             %d\n' ...
                '  pasted time stamps:   %d\n' ...
                '  mismatch (pasted - .czi): %+d\n\n' ...
                'You may have copied the time stamps from the wrong file.'],filename,nCzi,nPaste,mmVal);
            fprintf('%s\n',strrep(mmMsg,sprintf('\n\n'),sprintf('\n')));
            usable = nPaste >= ed;              % enough rows to cover every frame used
            if usable
                opts = {'Paste another time series','Use these anyway (extra rows ignored)','Skip (use frame index)'};
            else
                opts = {'Paste another time series','Skip (use frame index)'};
            end
            choiceT = questdlg(sprintf('%s\n\nWould you like to paste another time series?',mmMsg), ...
                'Time series mismatch',opts{:},opts{1});
            if strcmp(choiceT,'Paste another time series')
                continue;                       % ask for a new paste
            elseif strcmp(choiceT,'Use these anyway (extra rows ignored)')
                timeVec = numsT(st:ed); break;
            else
                timeVec = []; break;            % skip / closed dialog: frame index is used
            end
    end
    end
    % keep the time stamps so a re-run (rerun = true) can reload them without the .czi / pasting
    if ~rerun, save(fullfile(output_dir,[filename '_TimeVec.mat']),'timeVec'); end

    % Extra info for the Analysis_dataset workbook (Evelyn only), asked once per file
    % (re-run: taken from the rerunDrug / rerunSex / rerunAge settings at the top of the script)
    if strcmp(user,'Evelyn')
        if rerun
            drgChoice = rerunDrug; sexChoice = rerunSex; ageVal = rerunAge;
        else
            drgChoice = questdlg('Which drug was used in this file?','DRG_ID','KCl','DZX','None','KCl');
            if isempty(drgChoice) || strcmp(drgChoice,'None'), drgChoice = ''; end
            sexChoice = questdlg('Sex of the animal:','Sex','M','F','Unknown','M');
            if isempty(sexChoice) || strcmp(sexChoice,'Unknown'), sexChoice = ''; end
            ageAns = inputdlg('Age (leave blank if unknown):','Age',1,{''});
            ageVal = NaN; if ~isempty(ageAns), ageVal = str2double(ageAns{1}); end
        end
        gTok = regexp(filename,'(\d+(?:\.\d+)?)\s*mM','tokens','once');   % e.g. ...-2mMF4 -> 2, ..._11mMPC_d1 -> 11
        if isempty(gTok), glucoseVal = NaN; else, glucoseVal = str2double(gTok{1}); end
    end

    % Optional: frequency / amplitude / duty cycle of the normalized data (2mM, first phase, second phase, drug)
    if rerun
        runFreqAns = true;
    else
        faAns = questdlg('Also calculate frequency, amplitude and duty cycle (2mM / first phase / second phase / drug) and export to Excel?', ...
            'Frequency analysis','Yes','No','Yes');
        runFreqAns = strcmp(faAns,'Yes');
    end

    for g = 1:numel(groupIDs)
        k = groupIDs(g);
        if k == 0, continue; end            % cells outside any islet have no period selection
        sfx = groupSuffix{g};
        wfS = load(fullfile(output_dir,[filename sfx '_CaWaveForm.mat'])); % raw data
        if isfield(wfS,'CellTC'), wf.CellTC = wfS.CellTC; else, fnW = fieldnames(wfS); wf.CellTC = wfS.(fnW{1}); end
        idF = fullfile(output_dir,[filename sfx '_CellIDs.mat']);
        if exist(idF,'file')
            idS = load(idF);
            if isfield(idS,'cellIDs'), ids.cellIDs = idS.cellIDs; else, fnI = fieldnames(idS); ids.cellIDs = idS.(fnI{1}); end
        else
            ids.cellIDs = 1:size(wf.CellTC,2);      % older run without saved IDs
            warning('No %s found - labelling cells 1..N.',idF);
        end
        cellIDs_use = ids.cellIDs(:)';
        cellLabels = arrayfun(@(j) sprintf('cell%d',j),cellIDs_use,'UniformOutput',false);
        CellTC_use = wf.CellTC;
        if rerun && ~isempty(timeVec) && numel(timeVec) ~= size(CellTC_use,1)
            warning('Time stamps (%d) do not match the number of frames (%d) - %s.',numel(timeVec),size(CellTC_use,1), ...
                'using the first frames / frame index');
            if numel(timeVec) > size(CellTC_use,1), timeVec = timeVec(1:size(CellTC_use,1)); else, timeVec = []; end
        end
        if rerun
            % re-run: optionally restrict to some islets
            if ~isempty(rerunIslets) && ~ismember(k,rerunIslets), continue; end
            % start from the removal saved by the previous run (unless restored)
            rf = fullfile(output_dir,[filename sfx '_RemovedCells.mat']);
            removedCells = {};
            if exist(rf,'file') && ~ismember(k,rerunRestoreIslets)
                Sr = load(rf); fnR = fieldnames(Sr);
                if isfield(Sr,'removedCells'), removedCells = cellstr(Sr.removedCells); else, removedCells = cellstr(Sr.(fnR{1})); end
            elseif ismember(k,rerunRestoreIslets)
                fprintf('  Islet %d: previous removals cleared (rerunRestoreIslets).\n',k);
            end
            % new removals typed in rerunRemoveCells
            for rr = 1:size(rerunRemoveCells,1)
                if rerunRemoveCells{rr,1} == k
                    removedCells = [removedCells, arrayfun(@(j) sprintf('cell%d',j),rerunRemoveCells{rr,2}(:)','UniformOutput',false)]; %#ok<AGROW>
                end
            end
            removedCells = unique(removedCells,'stable');
            % optional click-to-remove window (traces shown with the saved removals already applied)
            if ismember(k,rerunReviewIslets)
                curMask = ~ismember(cellLabels,removedCells);
                lblNow = cellLabels(curMask);
                keepNow = removeCellsGUI(CellTC_use(:,curMask), lblNow, [filename sfx]);
                removedCells = unique([removedCells, lblNow(~keepNow)],'stable');
            end
            rmMask = ismember(cellLabels,removedCells);
            if any(rmMask)
                fprintf('  Islet %d: %d cell(s) removed: %s\n',k,nnz(rmMask),strjoin(cellLabels(rmMask),', '));
                CellTC_use = CellTC_use(:,~rmMask);
                cellLabels = cellLabels(~rmMask);
                cellIDs_use = cellIDs_use(~rmMask);
            end
            save(rf,'removedCells');
            fprintf('  Islet %d: saved removal list to %s\n',k,rf);

        else
        % Optional: remove outlier cells BEFORE detrending/normalizing (they would skew the trendline)
        % The raw traces are shown FIRST so you can judge whether any cell is an outlier, then you are asked.
        previewFig = figure('Name',['Raw traces ' filename sfx ' - look for outlier cells'],'NumberTitle','off');
        plot(CellTC_use,'Color',[0.3 0.3 0.3 0.6]);
        xlabel('Frame index'); ylabel('Intensity');
        title([filename sfx ' | ' num2str(size(CellTC_use,2)) ' cells'],'Interpreter','none');
        drawnow
        rmAns = questdlg(sprintf('Looking at the traces in the figure: remove any outlier cells from %s?\n(Yes opens the click-to-remove window.)',[filename sfx]), ...
            'Remove cells','Yes','No','No');
        if isvalid(previewFig), close(previewFig); end
        if strcmp(rmAns,'Yes')
            keepCells = removeCellsGUI(CellTC_use, cellLabels, [filename sfx]);
            if any(~keepCells)
                fprintf('  Removed %d cell(s): %s\n',nnz(~keepCells),strjoin(cellLabels(~keepCells),', '));
                removedCells = cellLabels(~keepCells);
                save(fullfile(output_dir,[filename sfx '_RemovedCells.mat']),'removedCells');
                CellTC_use = CellTC_use(:,keepCells);
                cellLabels = cellLabels(keepCells);
                cellIDs_use = cellIDs_use(keepCells);
            end
        end
        end
        if isempty(CellTC_use) || size(CellTC_use,2) == 0
            warning('All cells removed for %s - skipping.',[filename sfx]); continue;
        end
        % sanity check: exactly one valid text label per remaining cell column
        if ~iscellstr(cellLabels) || numel(cellLabels) ~= size(CellTC_use,2) || any(cellfun(@isempty,cellLabels))
            warning('Cell labels did not match the %d remaining cell column(s) - regenerating them.',size(CellTC_use,2));
            if numel(cellIDs_use) == size(CellTC_use,2)
                cellLabels = arrayfun(@(j) sprintf('cell%d',j),cellIDs_use,'UniformOutput',false);
            else
                cellLabels = arrayfun(@(j) sprintf('cell%d',j),1:size(CellTC_use,2),'UniformOutput',false);
            end
        end
        cellLabels = cellLabels(:)';
        fprintf('Processing %s%s (%d cells, %d frames)\n',filename,sfx,size(CellTC_use,2),size(CellTC_use,1));
        processAndExportPeriods(CellTC_use, IsletPeriods_all{k}, timeVec, ...
            fullfile(output_dir,[filename sfx]), cellLabels, applyNeg1);

        % ---- plot: normalized data, every cell dim grey + mean of cells bold black ----
        try
            [~,matsPl] = processAndExportPeriods(CellTC_use, IsletPeriods_all{k}, timeVec, '', ...
                cellLabels, applyNeg1, false);        % compute only, no files
            if isfield(matsPl,'normalized')
                tPl = timeVec(:);
                if isempty(tPl) || numel(tPl) < size(matsPl.normalized,1)
                    tPl = (1:size(matsPl.normalized,1))'; xlab = 'Frame index';
                else
                    tPl = tPl(1:size(matsPl.normalized,1)); xlab = 'Time (s)';
                end
                NormFig = figure('Name',['Normalized traces ' filename sfx]);
                hold on
                plot(tPl, matsPl.normalized, 'Color',[0.75 0.75 0.75], 'LineWidth',0.5);
                plot(tPl, mean(matsPl.normalized,2,'omitnan'), 'k', 'LineWidth',2.5);
                hold off
                xlabel(xlab); ylabel('\DeltaF/F');
                title([filename sfx],'Interpreter','none');
                box off
                saveas(NormFig, fullfile(output_dir,[filename sfx '_NormalizedTraces.tif']));
                savefig(NormFig, fullfile(output_dir,[filename sfx '_NormalizedTraces.fig']));
            end
        catch ME
            warning('Could not plot normalized traces for %s (%s).',[filename sfx],ME.message);
        end

        % ---- optional: frequency / amplitude / duty cycle of the normalized data ----
        % phases: 2mM (= 'Low glucose'), First phase, Second phase, Drug.
        % Fs comes from the time stamps (timeVec). Output: <filename>_islet#_FreqAmpDutyCycle.xlsx
        if runFreqAns
            try
                [~,matsFA] = processAndExportPeriods(CellTC_use, IsletPeriods_all{k}, timeVec, '', ...
                    cellLabels, applyNeg1, false);        % compute only, no files
                if isfield(matsFA,'normalized')
                    freqAmpDutyCycleExport(matsFA.normalized, IsletPeriods_all{k}, timeVec, cellLabels, ...
                        fullfile(output_dir,[filename sfx '_FreqAmpDutyCycle.xlsx']));
                end
            catch ME
                warning('Frequency analysis failed for %s (%s).',[filename sfx],ME.message);
            end
        end

        % ---- Analysis_dataset_detrended.xlsx (Evelyn): one row + one matrices sheet per islet ----
        if strcmp(user,'Evelyn')
            try
                [~,matsDS] = processAndExportPeriods(CellTC_use, IsletPeriods_all{k}, timeVec, '', ...
                    cellLabels, applyNeg1, false);        % compute only, no files
                nT  = size(CellTC_use,1);
                tvDS = timeVec(:);
                if isempty(tvDS), tvDS = (1:nT)'; end
                if numel(tvDS) < nT, tvDS(end+1:nT,1) = NaN; end
                dataDS = table((1:nT)',tvDS(1:nT),'VariableNames',{'index','time'});
                if isfield(matsDS,'detrended')
                    dataDS = [dataDS, array2table(matsDS.detrended,'VariableNames',strcat('detrended_',cellLabels))];
                end
                if isfield(matsDS,'normalized')
                    dataDS = [dataDS, array2table(matsDS.normalized,'VariableNames',strcat('normalized_',cellLabels))];
                end
                P = IsletPeriods_all{k};
                [blI,blT]  = formatPeriodRange(P,'Baseline',tvDS,nFr);
                [~,~,hgI,hgT] = formatPeriodRange(P,'High glucose',tvDS,nFr);
                [fpI,fpT]  = formatPeriodRange(P,'First phase',tvDS,nFr);
                [spI,spT]  = formatPeriodRange(P,'Second phase',tvDS,nFr);
                [drI,drT]  = formatPeriodRange(P,'Drug',tvDS,nFr);
                S = struct();
                S.(sprintf('Islet_%d',k)) = formatIndexRanges(cellIDs_use);
                S.BL_idx = blI; S.HG_idx = hgI; S.FP_idx = fpI; S.SP_idx = spI; S.DRG_idx = drI;
                S.DRG_ID = drgChoice; S.Glucose_Conditions = glucoseVal; S.Sex = sexChoice; S.Age = ageVal;
                S.BL = blT; S.HG = round(hgT,2); S.FP = fpT; S.SP = spT; S.DRG = drT;
                appendAnalysisDataset(analysisDatasetPath,[filename sfx],S,dataDS);
            catch ME
                warning('Could not update %s (%s). Is the file open in Excel?',analysisDatasetPath,ME.message);
            end
        end
    end
end

%% CaMPARI RED / GREEN CELL THRESHOLD (optional)
% Uses the red/green channels of the high-res <name>_pre.czi (OverlayStack), lets you draw a
% background ROI plus green and red cell ROIs, click the red/green ratio threshold, and saves it
% into an EXISTING Excel file (sheets 'Thresholds' and 'background_values'; row replaced if present).
if rerun
    thrAns = 'No';   % the threshold needs the _pre.czi + drawing: not part of a re-run
else
    thrAns = questdlg('Select the CaMPARI red/green threshold for this file?','CaMPARI threshold','Yes','No','No');
end
if strcmp(thrAns,'Yes')
    if isempty(OverlayStack)
        warning('No high-res _pre.czi was loaded - cannot select the threshold.');
    else
        nOv = size(OverlayStack,3);
        if nOv < 2
            warning('The _pre.czi has only %d channel - need a red and a green channel.',nOv);
        else
            % show every channel and let the user label each one Green / Red / Other
            [gCh, rCh, chLabels] = classifyChannelsGUI(OverlayStack, filename);
            if ~isnan(gCh) && ~isnan(rCh)
                fprintf('Channel %d = GREEN, channel %d = RED\n',gCh,rCh);
                greenIm = OverlayStack(:,:,gCh); redIm = OverlayStack(:,:,rCh);

                % --- split the _pre.czi by channel: one tiff per channel, named by its label ---
                % Green -> _green, Red -> _red, anything else -> _ch<number> (raw values, no scaling)
                to16 = @(im) uint16(min(max(double(im),0),65535));
                chSfx = cell(1,nOv);
                for cI = 1:nOv
                    if strcmp(chLabels{cI},'Green'),     chSfx{cI} = 'green';
                    elseif strcmp(chLabels{cI},'Red'),   chSfx{cI} = 'red';
                    else,                                chSfx{cI} = sprintf('ch%d',cI);
                    end
                    imwrite(to16(OverlayStack(:,:,cI)),fullfile(output_dir,[filename '_' chSfx{cI} '.tiff']));
                end
                % --- per islet: same channels with everything outside that islet set to 0 ---
                for kI = 1:nIslets
                    mI = imresize(IsletMask == kI, [size(OverlayStack,1) size(OverlayStack,2)], 'nearest');
                    for cI = 1:nOv
                        imwrite(to16(double(OverlayStack(:,:,cI)).*mI), ...
                            fullfile(output_dir,sprintf('%s_islet%d_%s.tiff',filename,kI,chSfx{cI})));
                    end
                end
                fprintf('Exported %d channel tiff(s) per file/islet (%s) to %s\n',nOv,strjoin(chSfx,', '),output_dir);

                % --- threshold: selected ONCE, written to Excel for every filename_islet# ---
                [thrFile, thrPath] = uigetfile('*.xlsx','Select the EXISTING threshold Excel file (e.g. CaMPARI_Color_Thresholds.xlsx)');
                if ~isequal(thrFile,0)
                    isletKeys = arrayfun(@(kI) sprintf('%s_islet%d',filename,kI),1:nIslets,'UniformOutput',false);
                    nameAns = inputdlg(['Row keys written to the Excel file (separate with commas):'], ...
                        'Filenames',1,{strjoin(isletKeys,', ')});
                    if ~isempty(nameAns) && ~isempty(strtrim(nameAns{1}))
                        keysOut = strtrim(strsplit(nameAns{1},','));
                        keysOut = keysOut(~cellfun(@isempty,keysOut));
                        campariThreshold = selectCampariThreshold(greenIm,redIm, ...
                            keysOut,fullfile(thrPath,thrFile));
                    end
                end
            else
                disp('Channel classification cancelled - skipping the threshold.');
            end
        end
    end
end

