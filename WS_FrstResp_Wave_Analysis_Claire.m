% THIS PROGRAM IMPORTS CALCIUM IMAGING FILES FOR CELL-BY-CELL ANALYSIS OF CALCIUM TRACE TO IDENTIFY FIRST RESPONDER AND WAVE ORIGIN/WAVE END CELLS
%% REAL-TIME USER INPUT REQUIRED
%% CREDIT: VIRA KRAVETS, FEB 2019; WOLFGANG SCHLEICHER, MARCH 2019; JAEANN DWULET, MARCH 2019

%% HOUSEKEEPING
close all
clear all
clc

addpath('/Users/levittcl/Documents/GitHub/UniversalCode');

%% LOADING THE CA IMAGE FILE
disp('Select Calcium Imaging File');
file = '/Volumes/CHL2021/NTEGcampE9.5WT/To Analyze 9 16 24/New Samples/NTE18-e(214.20)_Airyscan Processing.czi'

%R = bfopen('Select Calcium Imaging File'); % Uses bfopen program to open .czi/.lsm image files
R = bfopen (file);

howmanychannel = 1;

% omeMeta = R{1,4};
% voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
% voxelSizeXdouble = voxelSizeX.doubleValue();                                  % The numeric value represented by this object after conversion to type double                             

pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end

images=double(IMG); % converts images to double precision

images = images(:,:,1:howmanychannel:end);
RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable

clear pics R IMG;
clear omeMeta;

disp('Select Folder for Output Directory');
output_dir = uigetdir;
filename = input('Input Filename\n','s');


%% CODE TO MODIFY TIME COURSES
% This block converts the 3rd dimension of array into timecourse for
% later reference
st=1;

%st = 1921;
[xx, yy, zz] = size(images(:,:,:));
ed=zz;

%ed = 710;
T=st:1:ed;

images=images(:,:,st:ed);
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
ImGray = rgb2gray(RGB); %converts image to gray
imagesc(ImGray) %displays in default colormap; easier to see cell borders
disp('Select Background Signal')
bkgrndC=createMask(imfreehand); %allows user to select background area to be removed
close(GrayFig);
BackC=ImGray.*bkgrndC; %multiplies mask by the image
BackC(~logical(BackC))=nan; %takes inverse of masked image and sets it to nan
thresh = nanmax(nanmax(BackC)); %sets threshold based on the max of values after nan's are removed
close(OGFig);

[ImgNon, ~] = RemovingAreaswNoSignal(ImGray,thresh); %Uses function from JD to remove background based on threshold
ImGray(ImgNon)=0; %sets the removed areas to 0
NoSigFig = figure('Name','Background Removed; Draw ROI around Entire Islet');
imagesc(ImGray);
clear OGFig GrayFig;

%% MASKING ISLET
% User draws ROI around islet to plot whole islet timecourse
% User then determines where the start and end points for first responder
% and wave origin analyses
disp('Draw ROI around entire islet');
ROIMask_Start = imfreehand(); % User draws ROI around islet for average intensity
ROIMask_Start = createMask(ROIMask_Start);
sz_FR=ed-st+1;

StartMaskStack = images.*ROIMask_Start;
IsletTC = shiftdim(mean(mean(StartMaskStack)),2);
IsletTCfig = figure('Name','Whole Islet Time Course');
plot(IsletTC); %plots intensity timecourse of entire islet so user can identify first responder and wave origin ranges
datacursormode on % User can click on the plot to see x and y coordinates
stFR = input('Input Starting Frame for First Responder Analysis\n'); %select where first responder analysis should begin (x coordinate value)
edFR = input('Input Ending Frame for First Responder Analysis\n'); % select where first responder analysis should end (x coordinate value)
stWOall = input('Input Starting Frame for WO/WE Analysis\n'); %select where wave origin analysis should begin (x coordinate value)
edWOall = input('Input Ending Frame for WO/WE Analysis \n'); %select where wave origin analysis should end (x coordinate value)

answer = questdlg('Conduct Wave Origin/End Multipeak Analysis', ...
    'Mulitpeak Dialogue', ...
    'Yes','No','No');

switch answer
    case 'Yes'
        disp('Please Indicate Frame Intervals for each Peak')
        WaveInp = 1;
        NumPeak = 3;
        PeakRef = zeros(NumPeak,2);
        for i = 1:NumPeak
            PeakRef(i,1) = input(['Input Starting Frame for Peak ' num2str(i) '\n']);
            PeakRef(i,2) = input(['Input Ending Frame for Peak ' num2str(i) '\n']);
        end
        
    case 'No'
        disp('Multipeak Analysis Will NOT Be Conducted')
        WaveInp = 0;
end

close(IsletTCfig);
close(NoSigFig);
clear StartMaskStack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User Draws ROIs around each cell within the islet
% ROIs are saved in "CellMask" array and called back to throughout analysis
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
        NoSigFig = figure('Name','Draw ROIs Around Cells of Interest');
        imshow(RGB);
        
        k = 1;
        while k > 0
            disp('Draw ROIs Around Cells of Interest')
            ROIMask = imfreehand(); %User draws region around cell
            ROIMask = createMask(ROIMask); %Mask is created from drawn region
            CellMask = CellMask + ROIMask.*numcells; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
            CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
            UInp = input('Select Additonal ROI? 1=Yes, 0=No \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
            if UInp == 1 %User input to determine if another region is drawn
                k = k+1;
                numcells = numcells+1;
            elseif UInp ~= 1
                k = 0;
            end
        end
        
        save([output_dir '/' filename '_' 'Masks' '.mat'],'CellMask'); %Saves CellMask array
        close(NoSigFig);
        clear NoSigFig;
        
        save([output_dir '/' filename '_' 'CellNumber' '.mat'],'numcells'); %Saves numcells array
end

%% CALCULATING/PLOTTING INTENSITY OVER TIME FOR EACH ROI
sz_FR=edFR-stFR+1; % assigns a range of frames to reflect first responder analysis range
sz_WO=edWOall-stWOall+1; % assigns another range of frames to reflect wave origin analysis range
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
    CellTC(:,i) = TC; %Updates CellTC array with mean intensities from each frame of each cell
    PlotLabels{i} = ['Cell' num2str(i)]; %Updates labels with current cell number for legend
end
toc

% Plotting traces for entire time course
TCFig = figure('Name','Average Intensity Over Time');
plot(CellTC);
legend(PlotLabels);
clear images MaskedIMGstack;
saveas(TCFig,[output_dir '/' filename '_' 'TCPlot' '.tif']); %Saves figure of each cell's timecourse
save([output_dir '/' filename 'CaWaveForm.mat'],'CellTC'); % same cawaveform - Cell TC

%% TRANSLATING INPUT PARAMETERS FOR VK ANALYSES
calciumT = CellTC;
CalciumTable = array2table(calciumT, 'VariableNames', PlotLabels);
writetable(CalciumTable, [output_dir '/' filename '_' 'CalciumT' '.xls']);
clear CalciumTable;

%% FIRST/LAST RESPONDER ANALYSIS
IslAvg = mean(calciumT,2);
[RspMx,IndMx]=maxk(calciumT(stFR:edFR,:),1,1);             % returns max values of the Ca intensity for each cell, and corresponding index (time point)
[RspMn,IndMn]=mink(calciumT(stFR:edFR,:),1,1);

NAVG = (IslAvg - min(IslAvg(stFR:edFR,:)))/(max(IslAvg(stFR:edFR,:))-min(IslAvg(stFR:edFR,:)));
[~,NMn] = min(NAVG(stFR:edFR));
[~,NMx] = max(NAVG(stFR:edFR));

calciumN = zeros((edFR-stFR+1),numcells);
for k=1:numcells                                  % itterative index for cells
    currentcell_FR = calciumT(stFR:edFR,k);              % picking column of the calcium table corresponding to current cell, j
    currentcellN_FR =((currentcell_FR-RspMn(k))/(RspMx(k)-RspMn(k)));   % normalizing each cell to be between [0:1]
    calciumN(:,k)= currentcellN_FR;           % writing each normalized cell into array
%     [~, ccNFRMxI] = max(currentcellN_FR);
    currentcellN_FR = currentcellN_FR(1:IndMx(k));
    [~,cHHtime(:,k)] = min(abs(currentcellN_FR-0.5));     % cHHtime - time at which Ca elevation of the k-th cell reaches it's half-height; HHval - not important
end

half_time = cHHtime';

%cHHtime = cHHtime+NMn;
figure(5)
plot(calciumN)                                               % plotting normalized [0:1] Ca timelapses

[FirstVal,FirstInd] = mink(cHHtime,numcells,2);                    % FirstVal - HHElevT of the 1000 (or other number) 1st-responders; FirstInd - cell numbers of the first 1000 (or other number) 1st-responders
FirstInd = FirstInd';                                       % converting cell numbers of the 1st responding cells to column vector
Tresp=cHHtime';                                              % converting time at which Ca elevation of the k-th cell reaches it's half-height to column vector

minTresp = min(Tresp);
maxTresp = max(Tresp);
NormResptest = (Tresp-minTresp)/(maxTresp-minTresp);
NormResp = NormResptest'; 

FirstRespCa = calciumT(:,FirstInd(1,1));          %finding the ca timelapse corresponding to 1st-responder
LastRespCa = calciumT(:,FirstInd(end));         %finding the ca timelapse corresponding to last-responder
RespFig = figure('Name','First Responder, Last Responder, and Islet Avg');
plot(FirstRespCa);
title('Calcium; First/Last Resp. and Islet Average')
hold on
plot(IslAvg);
hold on
plot(LastRespCa)
legend('First Responder','Islet Average','Last Responder');

saveas(RespFig,[output_dir '/' filename '_' 'FR_LR_Avg' '.tif']); %Saves figure of each cell's timecourse
save([output_dir '/' filename '_' 'cHHTimes' '_' num2str(stFR) '_' num2str(edFR) '.mat'],'cHHtime'); %Saves CellMask array

%% VK WAVE ORIGIN/END CELL ANALYSIS

WholeWaveLags = WOWEAnalysis(calciumT,numcells,stWOall,edWOall,output_dir,filename);

if WaveInp == 1
    PeakLags = zeros(numcells,NumPeak);
    for i = 1:NumPeak
        PeakLags(:,i) = WOWEAnalysis(calciumT,numcells,PeakRef(i,1),PeakRef(i,2),output_dir,filename);
    end
end

%% GENERATING HEATMAPS OF RESPONDERS AND WAVE CELLS

Cent = zeros(numcells,2);
for j = 1:numcells
    DummyBW = CellMask; %Pulls in CellMask Array
    DummyBW(find(DummyBW ~= j)) = 0; %Gets rid of all but current mask
    STATS = regionprops(logical(DummyBW),'Centroid'); %Finds centroid of mask
    Cent(j,:) = STATS.Centroid; %Updates array with current centroid coordinates
end

WVCellMask = CellMask; %Calls in CellMask Array
FSCellMask = CellMask; %Calls in CellMask Array

for i = 1:numcells
    WVCellMask(find(WVCellMask == i)) = WholeWaveLags(i)+0.0001; %Assigns normalized lag value to cell mask based on its index
    FSCellMask(find(FSCellMask == i)) = NormResp(i)+0.0001; %Assigns normalized response value to cell mask based on its index
end

MnMxLim = round(numcells/10); % Determines number of 10% of cells

% This section identifies 10% of cells with highest lag and 10% of cells
% with lowest lag each individual peak and the entire wave propagation,
% then compares them to identify "true" wave origin/end cells

for kk = 1:MnMxLim
    [WvRfMin(kk), MnRfInd(kk)] = min(WholeWaveLags); %Indexes minimum of normalized lag values
    [WvRfMax(kk), MxRfInd(kk)] = max(WholeWaveLags); %Indexes maximum of normlalized lag values
    [RespMin(kk), RspMinInd(kk)] = min(NormResp); %Indexes minimum of normalized response values
    [RespMax(kk), RspMaxInd(kk)] = max(NormResp); %Indexes maximum of normalized response values
    WholeWaveLags(MnRfInd(kk)) = NaN; %Removes minimum for next cycle
    WholeWaveLags(MxRfInd(kk)) = NaN; %Removes maximum for next cycle
    NormResp(RspMinInd(kk)) = NaN; %Removes minimum for next cycle
    NormResp(RspMaxInd(kk)) = NaN; %Removes maximum for next cycle
    if kk == MnMxLim
        DupMn = find(WholeWaveLags == WvRfMin(kk));
        DupMn = DupMn';
        MnRfInd = [MnRfInd,DupMn];
        DupMx = find(WholeWaveLags == WvRfMax(kk));
        DupMx = DupMx';
        MxRfInd = [MxRfInd,DupMx];
        DupFrst = find(NormResp == RespMin(kk));
        %DupFrst = DupFrst';
        RspMinInd = [RspMinInd,DupFrst];
        DupLst = find(NormResp == RespMax(kk));
        %DupLst = DupLst';
        RspMaxInd = [RspMaxInd,DupLst];
    end  
end

RespRef = [RspMinInd, RspMaxInd]; %Concatenates max and min of responder cell indices
WaveRef = [MnRfInd, MxRfInd]; %Concatenates max and min of wave cell indices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if WaveInp == 1
    DupIndMn = zeros(NumPeak,numcells);
    DupIndMx = zeros(NumPeak,numcells);
    for i = 1:NumPeak
        for kk = 1:MnMxLim
            [WaveMin(i,kk), WvMinInd(i,kk)] = min(PeakLags(:,i)); %Indexes minimum of normalized lag values
            [WaveMax(i,kk), WvMaxInd(i,kk)] = max(PeakLags(:,i)); %Indexes maximum of normlalized lag values
            PeakLags(WvMinInd(i,kk),i) = NaN; %Removes minimum for next cycle
            PeakLags(WvMaxInd(i,kk),i) = NaN; %Removes maximum for next cycle
            if kk == MnMxLim
                DupInd = find(PeakLags(:,i) == WaveMin(i,kk));
                if ~isempty(DupInd)
                    for jj = 1:size(DupInd)
                        DupIndMn(i,jj) = DupInd(jj,1);
                    end
                end
                DupInd = find(PeakLags(i,:) == WaveMax(i,kk));
                if ~isempty(DupInd)
                    for jj = 1:size(DupInd)
                        DupIndMx(i,jj) = DupInd(jj,1);
                    end
                end
                
            end
        end
    end
    
    WvMinInd = [WvMinInd,DupIndMn];
    WvMaxInd = [WvMaxInd,DupIndMx];
    
    WvMn = zeros(MnMxLim,1); % preallocates WvMn matrix for size
    WvMx = zeros(MnMxLim,1); % preallocates WvMx matrix for size
    for ii = 1:MnMxLim
        [wMmax, wFmax] = mode(WvMaxInd(WvMaxInd>0),'all'); % determines the mode and mode-frequency of 10% of cells with highest lag for each peak
        if wFmax >= NumPeak-1 % reference for mode-frequency threshold
            WvMx(ii,1) = wMmax; % assigns the mode to WvMx matrix if its frequency across multiple peaks meets threshold
        end
        WvMaxInd(find(WvMaxInd == wMmax)) = NaN; % removes mode for next cycle
        [wMmin, wFmin] = mode(WvMinInd(WvMinInd>0),'all'); % determines the mode and mode-frequency of 10% of cells with lowest lag for each peak
        if wFmin >= NumPeak-1 % reference for mode-frequency threshold
            WvMn(ii,1) = wMmin; % assigns the mode to WvMn matrix if its frequency across multiple peaks meets threshold
        end
        WvMinInd(find(WvMinInd == wMmin)) = NaN; % removes values equal to mode for next cycle
    end
    WvMx = [MxRfInd, WvMx']; % concatenates the Max Lag reference matrix with WvMx 
    WvMn = [MnRfInd, WvMn']; % concatenates the Min Lag reference matrix with WvMn
    for j = 1:MnMxLim
        [Twx,Tfx] = mode(WvMx, 'all'); % Takes mode of concatenated matrices to determine which cells agree between WvMx and the reference
        if Tfx > 1
            TruWave(1,j) = Twx; % assigns the mode to a new matrix if it's frequency is greater than one
        end
        WvMx(find(WvMx == Twx)) = NaN; % removes values equal to current mode for next cycle
        [Twn,Tfn] = mode(WvMn,'all'); % Takes mode of concatenated matrices to determine which cells agree between WvMn and the reference
        if Tfn > 1
            TruWave(2,j) = Twn; % assigns the mode to a new matrix if it's frequency is greater than one
        end
        WvMn(find(WvMn == Twn)) = NaN; % removes values equal to current mode for next cycle
    end
    try
        TruWaveRef = sort(nonzeros(TruWave)); % Sorts and removes the zeros from matrix containing Cell ID numbers of "true" wave origin/end cells
    catch
        warning('No True Wave Origin/End Cells Detected.  Assigning "WaveInp" value of 0 to avoid further errors.');
        WaveInp = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
WaveCells = zeros(xx,yy);
RespCells = zeros(xx,yy);

for kk = 1:numel(WaveRef)
    WaveMask = CellMask; %Calls in CellMask Array for wave cell reference
    WaveMask(find(WaveMask ~= WaveRef(kk))) = 0; %Gets rid of all but current mask
    WaveCells = or(WaveMask,WaveCells); %Updates a logical with current mask
    if kk == length(WaveRef)
        WavePerim = bwperim(WaveCells); %Outlines the min and max 10% of wave cells
    end
end

for kk = 1:numel(RespRef)
    RespMask = CellMask; %Calls in CellMask Array for resp cell reference
    RespMask(find(RespMask ~= RespRef(kk))) = 0; %Gets rid of all but current mask
    RespCells = or(RespMask,RespCells); %Updates a logical with current mask
    if kk == length(RespRef)
        RespPerim = bwperim(RespCells); %Outlines the min and max 10% of responder cells
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if WaveInp == 1
    TruWaveCells = zeros(xx,yy);
    for kk = 1:numel(TruWaveRef)
        TruWaveMask = CellMask; %Calls in CellMask Array for wave cell reference
        TruWaveMask(find(TruWaveMask ~= TruWaveRef(kk))) = 0; %Gets rid of all but current mask
        TruWaveCells = or(TruWaveMask,TruWaveCells); %Updates a logical with current mask
        if kk == length(TruWaveRef)
            TruWavePerim = bwperim(TruWaveCells); %Outlines the min and max 10% of wave cells
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wcs = [0,0,1]; %WaveCell colormap start
Wce = [1,0,1]; %WaveCell colormap end 
WaveR = linspace(Wcs(1),Wce(1))';
WaveG = linspace(Wcs(2),Wce(2))';
WaveB = linspace(Wcs(3),Wce(3))';
WaveMap = [WaveR,WaveG,WaveB]; %Generates colormap for WaveCells

Fcs = [0.8,0,0]; %RespCell colormap start
Fce = [0,0.8,0.1]; %RespCell colormap end
FstR = linspace(Fcs(1),Fce(1))';
FstG = linspace(Fcs(2),Fce(2))';
FstB = linspace(Fcs(3),Fce(3))';
FirstMap = [FstR,FstG,FstB]; %Generates colormap for RespCells

pause on
pause(1);
ImgWave = rgb2gray(RGB); %Calls in representative image in grayscale
ImgWave(WavePerim) = 255; %Outlines cells of 10% max and 10% min lag
FinalWave = imoverlayNurin(ImgWave,WVCellMask,[0.0001,1],[0,0.8],WaveMap); %Overlays heatmap of wavecells 
colorbar;

pause(1);

ci = 1;
for c = 1:numcells
    text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w'); %Labels cells in the image with their respective region number
    ci = ci+1;
end

pause(1);

% Scalebar_length = 10/voxelSizeXdouble; % in pixels
% x_location = round(xx.*0.05); % denotes x location in pixels
% y_location = round(yy.*0.95); % denotes y location in pixels
% hold on
% quiver(x_location,y_location,Scalebar_length,0,'ShowArrowHead','off','Color','w','LineWidth',3.0);
% txt = '10um';
% text(x_location,y_location-10,txt,'Color','w');
% 
% set(gcf, 'Position',  [100, 100, 512, 512])

saveas(FinalWave,[output_dir '/' filename '_' 'WaveCells' '.fig']); %saves heatmap image with labels and outlined cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if WaveInp == 1
    pause(1)
    ImgWave = rgb2gray(RGB); %Calls in representative image in grayscale
    ImgWave(TruWavePerim) = 255; %Outlines cells of 10% max and 10% min lag
    FinalWave = imoverlayNurin(ImgWave,WVCellMask,[0.0001,1],[0,0.8],WaveMap); %Overlays heatmap of wavecells
    colorbar;
    
    pause(1);
    
    ci = 1;
    for c = 1:numcells
        text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w'); %Labels cells in the image with their respective region number
        ci = ci+1;
    end
    
    pause(1);
    
%     Scalebar_length = 10/voxelSizeXdouble; % in pixels
%     x_location = round(xx.*0.05); % denotes x location in pixels
%     y_location = round(yy.*0.95); % denotes y location in pixels
%     hold on
%     quiver(x_location,y_location,Scalebar_length,0,'ShowArrowHead','off','Color','w','LineWidth',3.0);
%     txt = '10um';
%     text(x_location,y_location-10,txt,'Color','w');
%     
%     set(gcf, 'Position',  [100, 100, 512, 512])
%     
    saveas(FinalWave,[output_dir '/' filename '_' 'TruWaveCells' '.fig']); %saves heatmap image with labels and outlined cells
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pause(1);

ImgResp = rgb2gray(RGB); %calls in original image in grayscale
ImgResp(RespPerim) = 255; %Outlines cells of 10% max and 10% min response times
FinalResp = imoverlayNurin(ImgResp,FSCellMask,[0.0001,1.0001],[0,1],FirstMap); %overlays heatmap of responder cells
colorbar;

pause(1);

ci = 1;
for c = 1:numcells
    text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w'); %Labels cells in the image with their respective region number
    ci = ci+1;
end

pause(1);

% Scalebar_length = 10/voxelSizeXdouble; % in pixels
% x_location = round(xx.*0.05); % denotes x location in pixels
% y_location = round(yy.*0.95); % denotes y location in pixels
% % hold on
% quiver(x_location,y_location,Scalebar_length,0,'ShowArrowHead','off','Color','w','LineWidth',3.0);
% txt = '10um';
% text(x_location,y_location-10,txt,'Color','w');
% 
% set(gcf, 'Position', [100, 100, 512, 512])

saveas(FinalResp,[output_dir '/' filename '_' 'Responders' '.fig']); %saves heatmap image with labels and outlined cells

pause(1);
%% MEASURING DISTANCE BETWEEN CELLS OF INTEREST
answer = questdlg('Would you like to measure distances between cells?', ...
    'Distance Dialogue', ...
	'Yes','No','No');

switch answer
    case 'Yes'
        disp([answer ': ' 'Measure Distances between Cells']);
        DistFig = figure('Name','Measure Distance Between Cells of Interest');
        imshow(RGB); %Calls in original image in color
        ci = 1;
        for c = 1:j
            text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w'); %Labels cells in image with respective regions
            ci = ci+1;
        end
        n = 1;
        kk = 1;
        while kk > 0
            disp('Input IDs of Cells you wish to Link');
            CellA = num2str(input('ID of Cell A \n')); % User inputs cell number for first cell
            CellB = num2str(input('ID of Cell B \n')); % User inputs cell number for next cell
            CellDist.Label(n,1) = cellstr([CellA 'to' CellB]); % Stores the above user inputs at labels
            disp('Draw Line linking Nuclei of Cells of Interest') % Prompts user to draw line between cells they identified previously
            DistLine = drawline('LineWidth',0.5,'Color','y'); % Allows user to draw lines linking cells
            CellDist.Value(n,1) = pdist(DistLine.Position,'euclidean').*voxelSizeXdouble; % Calculates distance of each line in micrometers
            UInp = input('Measure Another Distance? 1=Yes, 0=No \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
            if kk == UInp
                n = n+1;
            elseif UInp == 0
                kk = 0;
            end
        end
        
        CellDist = [CellDist.Label num2cell(CellDist.Value)]; % Concatenates distance calculation array and label array
        xlswrite([output_dir '\' filename '_' 'Distances'],CellDist); % Exports the concatenated distance and label array as excel sheet.
    case 'No'
        disp('End');
end

%%%%%%%%-----END OF CODE-----%%%%%%%%
