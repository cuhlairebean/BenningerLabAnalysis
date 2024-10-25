% Calcium Timecourse Two Channel Analysis %

% 10.24.24 Claire H Levitt
% This code is for live calcium imaging with 2 separate channels

clear all; close all; clc

%% HouseKeeping
addpath('/Users/levittcl/Documents/GitHub/UniversalCode');
addpath('/Users/levittcl/Documents/GitHub/Calcium_Analysis_Bucket_CL/Universal_Code/')

%% LOAD CALCIUM IMAGE FILE
% comment out whichever file upload method you aren't using

% 1. SELECT MAUALLY
% disp('Select Calcium Imaging File');
% R = bfopen('Select Calcium Imaging File'); % Uses bfopen program to open .czi/.lsm image files

% 2. INPUT FILE NAME
filename = '/Users/levittcl/Documents/1_Research/7_Example Imaging Files/2_channels_2_calcium/islet3_2mM_11mM_ghrelin_2mM_11mMGLP1.czi';

R = bfopen(filename); % Uses bfopen program to open .czi/.lsm image files

% Number of channels in file - assuming 2
howmanychannel = 2;

pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end

for i=1:length(pics)
    IMG2(:,:,i)=pics{i};
end

images=double(IMG); % converts images to double precision
Ch2_RawImg = images(:,:,2); % FOR SECOND CHANNEL
images = images(:,:,1:howmanychannel:end); %images = first channel
RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable
ref = images(:,:,1:end); % This chooses from a set of frames to create the max intensity map: mostly important if theres movement
MaxIntensityMap = max(ref, [],3);
ref_pic = (mat2gray(MaxIntensityMap));

Ch2_Ref = mat2gray(Ch2_RawImg); % Reference for second channle

clear pics R IMG;
clear omeMeta;

% SELECT LOCATION FOR SAVING FILES
disp('Select Folder for Output Directory');
output_dir = uigetdir;
filename = input('Input Filename  ','s');

%% CODE TO MODIFY TIME COURSES
% This block converts the 3rd dimension of array into timecourse
% change as necessary (st = start, ed = end)

st=1; % STARTING FRAME: Change this per sample if needed
%st = 100;
[xx, yy, zz] = size(images(:,:,:));
%ed = ed; % ENDING FRAME: Change this per sample if needed
ed=zz;

%ed = 1940;
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

disp('Draw ROI around entire islet');
ROIMask_Start = imfreehand(); % User draws ROI around islet for average intensity
ROIMask_Start = createMask(ROIMask_Start);
sz_FR=ed-st+1;

StartMaskStack = images.*ROIMask_Start;
IsletTC = shiftdim(mean(mean(StartMaskStack)),2);
IsletTCfig = figure('Name','Whole Islet Time Course');
plot(IsletTC); %plots intensity timecourse of entire islet

%% INDIVIDUAL CELL MASKS
% User Draws ROIs around each cell within the islet
% ROIs are saved in "CellMask" array and called back to throughout analysis
% Follow the prompts in the Command Window

%figure, imshow(ref_pic)

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
        imshow(ref_pic);
        
        k = 1;
        while k > 0
            disp('Draw ROIs Around Cells of Interest')
            ROIMask = imfreehand(); %User draws region around cell
            ROIMask = createMask(ROIMask); %Mask is created from drawn region
            CellMask = CellMask + ROIMask.*numcells; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
            CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
            UInp = input('Select Additonal ROI? 1=Yes, 0=No, 3 = switch \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
            if UInp == 1 %User input to determine if another region is drawn
                k = k+1;
                numcells = numcells+1;
   
            elseif UInp ~= 1
                k = 0;
            end
        end
        
        save([output_dir '/' filename '_' 'Masks_ch1' '.mat'],'CellMask'); %Saves CellMask array
        close(NoSigFig);
        clear NoSigFig;
        
        save([output_dir '/' filename '_' 'CellNumber_ch1' '.mat'],'numcells'); %Saves numcells array
end

%% EXTRACT TIMECOURSE 1

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
    PlotLabels{i} = [' Cell ' num2str(i)]; %Updates labels with current cell number for legend
end
toc

%All other cells%
Ch1_Ca = CellTC;

% Plotting traces for entire time course
TCFig = figure('Name','Average Intensity Over Time');
plot(Ch1_Ca);
legend(PlotLabels);
clear MaskedIMGstack;
% saveas(TCFig,[output_dir '/' filename '_' 'TCPlot' '.tif']); %Saves figure of each cell's timecourse
% save([output_dir '/CaWaveForm.mat'],'CellTC'); % same cawaveform - Cell TC

% Normalized Timecourse

% View Cells w labels
ImgWave = rgb2gray(RGB); %Calls in representative image in grayscale
Cent = zeros(numcells,2);
for j = 1:numcells
    DummyBW = CellMask; %Pulls in CellMask Array
    DummyBW(find(DummyBW ~= j)) = 0; %Gets rid of all but current mask
    STATS = regionprops(logical(DummyBW),'Centroid'); %Finds centroid of mask
    Cent(j,:) = STATS.Centroid; %Updates array with current centroid coordinates
end
Channel1_Ca = figure, imshow(ref_pic)
title('Channel 1')
ci = 1;
for c = 1:numcells
    text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w'); %Labels cells in the image with their respective region number
    ci = ci+1;
end

saveas(Channel1_Ca,[output_dir '/' filename '_' 'channel 1' '.tif']);

%% CHANNEL 2

clear CellMask
clear numcells

images=double(IMG2); % converts images to double precision
images = images(:,:,2:howmanychannel:end); %images = first channel
RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable
ref = images(:,:,1:end); % This chooses from a set of frames to create the max intensity map: mostly important if theres movement
MaxIntensityMap = max(ref, [],3);
ref_pic = (mat2gray(MaxIntensityMap));

clear IMG2;

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
%% MASKING ISLET CHANNEL 2
% User draws ROI around islet to plot whole islet timecourse

disp('Draw ROI around entire islet');
ROIMask_Start = imfreehand(); % User draws ROI around islet for average intensity
ROIMask_Start = createMask(ROIMask_Start);
sz_FR=ed-st+1;

StartMaskStack = images.*ROIMask_Start;
IsletTC = shiftdim(mean(mean(StartMaskStack)),2);
IsletTCfig = figure('Name','Whole Islet Time Course');
plot(IsletTC); %plots intensity timecourse of entire islet

%% INDIVIDUAL CELL MASKS 2
% User Draws ROIs around each cell within the islet
% ROIs are saved in "CellMask" array and called back to throughout analysis
% Follow the prompts in the Command Window

%figure, imshow(ref_pic)

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
        imshow(ref_pic);
        
        k = 1;
        while k > 0
            disp('Draw ROIs Around Cells of Interest')
            ROIMask = imfreehand(); %User draws region around cell
            ROIMask = createMask(ROIMask); %Mask is created from drawn region
            CellMask = CellMask + ROIMask.*numcells; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
            CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
            UInp = input('Select Additonal ROI? 1=Yes, 0=No,'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
            if UInp == 1 %User input to determine if another region is drawn
                k = k+1;
                numcells = numcells+1;
   
            elseif UInp ~= 1
                k = 0;
            end
        end
        
        save([output_dir '/' filename '_' 'Masks_ch2' '.mat'],'CellMask'); %Saves CellMask array
        close(NoSigFig);
        clear NoSigFig;
        
        save([output_dir '/' filename '_' 'CellNumber_ch2' '.mat'],'numcells'); %Saves numcells array
end

%% EXTRACT TIMECOURSE 2

sz_2=ed-st; % includes all frames for exporting of entire trace

PlotHandles2 = zeros(1,numcells);
PlotLabels2 = cell(1,numcells);
CellTC_2 = zeros(sz_2,numcells);
TC2 = zeros(sz_2,1);

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

    CellTC2(:,i) = TC; %Updates CellTC array with mean intensities from each frame of each cell
    PlotLabels{i} = [' Cell ' num2str(i)]; %Updates labels with current cell number for legend
end
toc

Ch2_TC = CellTC2;

%All other cells%
All_TC = [CellTC Ch2_TC];

% Plotting traces for entire time course
TCFig = figure('Name','Average Intensity Over Time');
plot(Ch2_TC);
legend(PlotLabels);
clear MaskedIMGstack;
% saveas(TCFig,[output_dir '/' filename '_' 'TCPlot' '.tif']); %Saves figure of each cell's timecourse
% save([output_dir '/CaWaveForm.mat'],'CellTC'); % same cawaveform - Cell TC

% Normalized Timecourse

% View Cells w labels
ImgWave = rgb2gray(RGB); %Calls in representative image in grayscale
Cent = zeros(numcells,2);
for j = 1:numcells
    DummyBW = CellMask; %Pulls in CellMask Array
    DummyBW(find(DummyBW ~= j)) = 0; %Gets rid of all but current mask
    STATS = regionprops(logical(DummyBW),'Centroid'); %Finds centroid of mask
    Cent(j,:) = STATS.Centroid; %Updates array with current centroid coordinates
end
Channel2_Ca = figure, imshow(ref_pic)
title('Channel 2')
ci = 1;
for c = 1:numcells
    text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w'); %Labels cells in the image with their respective region number
    ci = ci+1;
end

saveas(Channel2_Ca,[output_dir '/' filename '_' 'channel 2' '.tif']);

%% COMBINE Ch1 and Ch2

%combine channel 1 and channel 2
figure, plot(Ch1_Ca,'Color', [0.3010 0.7450 0.9330]) %channel 1
hold on
plot(Ch2_TC,'Color',[0.6350 0.0780 0.1840]) %channel 2
title ('Ch1 & 2 Combined')
legend

AllChannels = [Ch1_Ca Ch2_TC];
save([output_dir '/mixedCaWaveForm.mat'],'AllChannels'); % same cawaveform - Cell TC