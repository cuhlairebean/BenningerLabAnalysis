
% Modified Code for First responders only
%% HouseKeeping
clear all; close all; clc;
addpath('/Users/levittcl/Documents/GitHub/UniversalCode');
addpath ('/Users/levittcl/Documents/GitHub/UniCode')

%% Loading Calcium Imaging File 

disp('Select Calcium Imaging File');
filename = [''];
R = bfopen(filename); % Uses bfopen program to open .czi/.lsm image files
%R = bfopen('Select Calcium Imaging File'); % Uses bfopen program to open .czi/.lsm image files

% omeMeta = R{1,4};
% voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
% voxelSizeXdouble = voxelSizeX.doubleValue();                                  % The numeric value represented by this object after conversion to type double                             

pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end

images=double(IMG); % converts images to double precision
RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable

clear pics R IMG;
clear omeMeta;

disp('Select Folder for Output Directory');
output_dir = uigetdir;
filename = input('Input Filename\n','s');

%% Modify Timecourse

% This block converts the 3rd dimension of array into timecourse for
% later reference
st = 2715; % This dictates the Starting Frame (modify based on initial time lag)

[xx, yy, zz] = size(images(:,:,:));
ed=zz;
%ed = 1479;
T=st:1:ed;

images=images(:,:,st:ed);

%% Declare Image Properties

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

%% Masking the Islet

% User draws ROI around islet to plot whole islet timecourse
% User then determines where the start and end points for first responder
% analysis only (uncomment for Wave)

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

% For Wave Analysis %
%stWOall = input('Input Starting Frame for WO/WE Analysis\n'); %select where wave origin analysis should begin (x coordinate value)
%edWOall = input('Input Ending Frame for WO/WE Analysis \n'); %select where wave origin analysis should end (x coordinate value)

 disp('Multipeak Analysis Will NOT Be Conducted')
        WaveInp = 0;
 
close(IsletTCfig);
close(NoSigFig);
clear StartMaskStack;

%% Drawing ROIs

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
            %ROIMask = imfreehand(); %User draws region around cell
            ROIMask = drawcircle();
            ROIMask = createMask(ROIMask); %Mask is created from drawn region
            CellMask = CellMask + ROIMask.*numcells; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
            CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
            UInp = input('Select Additonal ROI? 1=Yes, 0=No, 2=Redraw \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
            if UInp == 1 %User input to determine if another region is drawn
                k = k+1;
                numcells = numcells+1;
            elseif UInp == 2
                numcells = numcells;
                k = k;
            elseif UInp ~= 1
                k = 0;
            end
        end
        
        save([output_dir '/' filename '_' 'Masks' '.mat'],'CellMask'); %Saves CellMask array
        close(NoSigFig);
        clear NoSigFig;
        
        save([output_dir '/' filename '_' 'CellNumber' '.mat'],'numcells'); %Saves numcells array
end

%% Calculate and Plot intensity for each ROI

sz_FR=edFR-stFR+1; % assigns a range of frames to reflect first responder analysis range
%sz_WO=edWOall-stWOall+1; % assigns another range of frames to reflect wave origin analysis range
sz_2=ed-st; % includes all frames for exporting of entire trace

PlotHandles = zeros(1,numcells);
PlotLabels = cell(1,numcells);
CellTC = zeros(sz_2,numcells);
TC = zeros(sz_2,1);

tic % How long does it take to do this analysis (tic --> toc)
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

%% calciumT = CellTC;
calciumT = CellTC;
CalciumTable = array2table(calciumT, 'VariableNames', PlotLabels);
writetable(CalciumTable, [output_dir '/' filename '_' 'CalciumT' '.xls']);
clear CalciumTable; % makes it faster to run 

%% First/Last Responder Analysis

IslAvg = mean(calciumT,2);

[RspMx,IndMx]=maxk(calciumT(stFR:edFR,:),1,1); % returns max values of the Ca intensity for each cell, and corresponding index (time point)
[RspMn,IndMn]=mink(calciumT(stFR:edFR,:),1,1);

NAVG = (IslAvg - min(IslAvg(stFR:edFR,:)))/(max(IslAvg(stFR:edFR,:))-min(IslAvg(stFR:edFR,:))); % normalized
[~,NMn] = min(NAVG(stFR:edFR));
[~,NMx] = max(NAVG(stFR:edFR));

%% Calciulating the Half Time %% -- Want to calculate the actual value 
cal_smoothed = smoothdata(calciumT(stFR:edFR,:),1);
currentcellN_FR = (cal_smoothed-min(cal_smoothed))./range(cal_smoothed);
calciumN = zeros((edFR-stFR+1),numcells);

for k=1:numcells                                  % itterative index for cells
    currentcell_FR = calciumT(stFR:edFR,k);              % picking column of the calcium table corresponding to current cell, j
    currentcellN_FR =((currentcell_FR-RspMn(k))/(RspMx(k)-RspMn(k)));   % normalizing each cell to be between [0:1]
    calciumN(:,k)= currentcellN_FR;           % writing each normalized cell into array
%     [~, ccNFRMxI] = max(currentcellN_FR);
    currentcellN_FR = currentcellN_FR(1:IndMx(k));
    [~,cHHtime(:,k)] = min(abs(currentcellN_FR-0.5));   % cHHtime - time at which Ca elevation of the k-th cell reaches it's half-height; HHval - not important
end

% %% what jenny would do :)
% [~,cell_phase] = sort(CHHtime)
% 
% 
% %cHHtime = cHHtime+NMn;
% figure(5)
% plot(calciumN)
% hold on
% plot (cHHtime,val,'o')

%% First Responder

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

%% Creating a heat map of cell responders

Cent = zeros(numcells,2);
for j = 1:numcells
    DummyBW = CellMask; %Pulls in CellMask Array
    DummyBW(find(DummyBW ~= j)) = 0; %Gets rid of all but current mask
    STATS = regionprops(logical(DummyBW),'Centroid'); %Finds centroid of mask
    Cent(j,:) = STATS.Centroid; %Updates array with current centroid coordinates
end
%% 
for i = 1:numcells
FSCellMask = CellMask; %Calls in CellMask Array
end 

MnMxLim = round(numcells/10); % Determines number of 10% of cells

% This section identifies 10% of cells with highest lag and 10% of cells
% with lowest lag each individual peak and the entire wave propagation,
% then compares them to identify "true" wave origin/end cells

for kk = 1:MnMxLim
    %[WvRfMin(kk), MnRfInd(kk)] = min(WholeWaveLags); %Indexes minimum of normalized lag values
    %[WvRfMax(kk), MxRfInd(kk)] = max(WholeWaveLags); %Indexes maximum of normlalized lag values
    [RespMin(kk), RspMinInd(kk)] = min(NormResp); %Indexes minimum of normalized response values
    [RespMax(kk), RspMaxInd(kk)] = max(NormResp); %Indexes maximum of normalized response values
    %WholeWaveLags(MnRfInd(kk)) = NaN; %Removes minimum for next cycle
    %WholeWaveLags(MxRfInd(kk)) = NaN; %Removes maximum for next cycle
    NormResp(RspMinInd(kk)) = NaN; %Removes minimum for next cycle
    NormResp(RspMaxInd(kk)) = NaN; %Removes maximum for next cycle
    if kk == MnMxLim
%         DupMn = find(WholeWaveLags == WvRfMin(kk));
%         DupMn = DupMn';
%         MnRfInd = [MnRfInd,DupMn];
%         DupMx = find(WholeWaveLags == WvRfMax(kk));
%         DupMx = DupMx';
%         MxRfInd = [MxRfInd,DupMx];
        DupFrst = find(NormResp == RespMin(kk));
        %DupFrst = DupFrst';
        RspMinInd = [RspMinInd,DupFrst];
        DupLst = find(NormResp == RespMax(kk));
        %DupLst = DupLst';
        RspMaxInd = [RspMaxInd,DupLst];
    end  
end

RespRef = [RspMinInd, RspMaxInd]; %Concatenates max and min of responder cell indices
%WaveRef = [MnRfInd, MxRfInd]; %Concatenates max and min of wave cell indices

RespCells = zeros(xx,yy);
for kk = 1:numel(RespRef)
    RespMask = CellMask; %Calls in CellMask Array for resp cell reference
    RespMask(find(RespMask ~= RespRef(kk))) = 0; %Gets rid of all but current mask
    RespCells = or(RespMask,RespCells); %Updates a logical with current mask
    if kk == length(RespRef)
        RespPerim = bwperim(RespCells); %Outlines the min and max 10% of responder cells
    end
end

% Color Maps % 
Fcs = [0.8,0,0]; %RespCell colormap start
Fce = [0,0.8,0.1]; %RespCell colormap end
FstR = linspace(Fcs(1),Fce(1))';
FstG = linspace(Fcs(2),Fce(2))';
FstB = linspace(Fcs(3),Fce(3))';
FirstMap = [FstR,FstG,FstB]; %Generates colormap for RespCells


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

saveas(FinalResp,[output_dir '/' filename '_' 'Responders' '.fig']); %saves heatmap image with labels and outlined cells
pause(1);

half = cHHtime';

disp('End');


%%%%%%%%-----END OF CODE-----%%%%%%%%