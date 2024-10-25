%% Extracting Single Cell Timecourse (option: from 2 channels)
% Modified from WS_FirstResp_Wave_Analysis
% 08.28.2023 Claire H Levitt

% This code is designed to extract timecourses of fluorescent signals
% For multiple channels (DELTA CELLS EXAMPLE); 
% all 2nd channel properties are labeled "delta cell"
% comment this out if it doens't exist
% LINE 59-60, LINE 204-207

% otherwise: when circling cells, to see next channel 3 = SWITCH channels

% Code is designed to extract the 3rd channel; but if it was not imaged
% together, need to upload colocalization snap separate. 
% LINE 175

% modify start time: LINE 83

%% HOUSEKEEPING
% Add things whatever
close all
clear all
clc

addpath('/Users/levittcl/Documents/GitHub/UniversalCode');
addpath('/Users/levittcl/Documents/GitHub/Calcium_Analysis_Bucket_CL/Universal_Code/')

%% LOAD CALCIUM IMAGE FILE

% comment out whichever file upload method you aren't using

% 1. SELECT MAUALLY
% disp('Select Calcium Imaging File');
% R = bfopen('Select Calcium Imaging File'); % Uses bfopen program to open .czi/.lsm image files

% 2. INPUT FILE NAME
filename = '/Volumes/Claire Hard Drive/2024_07_01 GCaMP Lux Stim/islet1_2mM_11mM_ghrelin_GLP1.czi';

% 3. Snap Filename
id2 = '/Volumes/Claire Hard Drive/2024_07_01 GCaMP Lux Stim/islet1_2mM_11mM_ghrelin_GLP1.czi';

R = bfopen(filename); % Uses bfopen program to open .czi/.lsm image files


% Number of channels in file - important if its multichannel or ref.channel
howmanychannel = 2;

pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end


images=double(IMG); % converts images to double precision

 second_RawImg = images(:,:,3); % FOR SECOND CHANNEL (RED)
 second_ref = mat2gray(second_RawImg);

images = images(:,:,1:howmanychannel:end);
RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable


ref = images(:,:,1:end); % This chooses from a set of frames to create the max intensity map: mostly important if theres movement
MaxIntensityMap = max(ref, [],3);
ref_pic = (mat2gray(MaxIntensityMap));


clear pics R IMG;
clear omeMeta;


%%  SELECT LOCATION FOR SAVING FILES
disp('Select Folder for Output Directory');
output_dir = uigetdir;
filename = input('Input Filename  ','s');

%% CODE TO MODIFY TIME COURSES
% This block converts the 3rd dimension of array into timecourse for
% later reference

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

save([output_dir '/' filename '_' 'backgroundthresh' '.mat'],'thresh'); %Saves threshold value array
close(OGFig);

%%

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
StartMask_NoBackground = StartMaskStack.*ImGray;
IsletTC = shiftdim(mean(mean(StartMask_NoBackground)),2);
IsletTCfig = figure('Name','Whole Islet Time Course');
plot(IsletTC); %plots intensity timecourse of entire islet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For colocalization% - comment this out if you don't have multiple channels
% disp('Select Colocalization File')
% [file2, path2] = uigetfile('*.czi')
% id2 = [path2 file2];
% disp(id2);
% HI_image = bfopen(id2); %Human Islet Image (HI)
% 
% red_cells = HI_image{1,1}{3,1}; 
% red_cells_resize = imresize(red_cells,[xx yy]);
% 
% green_cells = HI_image{1,1}{1,1}; %calcium
% green_cells_resize = imresize(green_cells,[xx yy]);
% 
% ca_cells = HI_image{1,1}{2,1}; %green cells
% ca_cells_resize = imresize(ca_cells,[xx yy]);
% 
% r = imfuse(ca_cells,red_cells); %red cells + calcium
% g = imfuse(ca_cells,green_cells); %green cells + calcium
% 
% b = imfuse(r,g);
% 
% figure, imshow(r)
% title ('RED overlay')
% 
% figure, imshow(g)
% title ('GREEN overlay')

%% Upload SST File
% channel {1,1}:all; non delta,  channel {3,1}: Delta Cells
% load colocalization image; COMMENT THIS OUT IF NO COLOCALIZATION SNAP

 % Manually Upload
% disp('Select Colocalization File')
% [file2, path2] = uigetfile('*.czi')
% id2 = [path2 file2];
% disp(id2);

%second_image = bfopen(id2);
second_cells = second_image{1,1}{3,1}; 
second_ref = imresize(second_cells,[xx yy]);

%% Cell Masking
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
        close all
        NoSigFig = figure('Name','Draw ROIs Around Cells of Interest');
        imshow(ref_pic);
        
        k = 1;
        while k > 0
            disp('Draw ROIs Around Cells of Interest')
            ROIMask = imfreehand(); %User draws region around cell
            ROIMask = createMask(ROIMask); %Mask is created from drawn region
            CellMask = CellMask + ROIMask.*numcells; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
            CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
            UInp = input('Select Additonal ROI? 1=Yes, 0=No, 3 = SWITCH \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
            if UInp == 1 %User input to determine if another region is drawn
                k = k+1;
                numcells = numcells+1;
            elseif UInp == 3
                imshow(delta_ref)
                k = k+1;
                numcells = numcells + 1;
            elseif UInp ~= 1
                k = 0;
            end
        end
        
        save([output_dir '/' filename '_' 'Masks' '.mat'],'CellMask'); %Saves CellMask array
        close(NoSigFig);
        clear NoSigFig;
        
        save([output_dir '/' filename '_' 'CellNumber' '.mat'],'numcells'); %Saves numcells array
        save([output_dir '/' filename '_' 'thresh' '.mat'],'thresh'); %Saves numcells array
end

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
        TCnoZero = TCnoZero(TCnoZero>0);%Accounts for any zeros from preallocation
        TC(ii) = mean(TCnoZero); %Calculates mean intensity of current frame
    end

    CellTC(:,i) = TC; %Updates CellTC array with mean intensities from each frame of each cell
    PlotLabels{i} = [' Cell ' num2str(i)]; %Updates labels with current cell number for legend
end
toc

% Plotting traces for entire time course
TCFig = figure('Name','Average Intensity Over Time');
plot(CellTC);
legend(PlotLabels);
clear MaskedIMGstack;
saveas(TCFig,[output_dir '/' filename '_' 'TCPlot' '.tif']); %Saves figure of each cell's timecourse
save([output_dir '/CaWaveForm.mat'],'CellTC'); % same cawaveform - Cell TC

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
cell_labels = figure, imshow(ref_pic)
title('Whole Islet')
ci = 1;
for c = 1:numcells
    text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w'); %Labels cells in the image with their respective region number
    ci = ci+1;
end

saveas(cell_labels,[output_dir '/' filename '_' 'cell_labels' '.tif']); %Saves figure of cell labels

disp('END! Nice Job :D ')