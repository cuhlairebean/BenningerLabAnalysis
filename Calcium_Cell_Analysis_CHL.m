%% Calcium Duty Cycle and Frequency Analysis %% 
% Claire H Levitt
% 01.19.2023
% Analysis Method from Jenn's code - *note* only allows cell count up to 75

% HouseKeeping
clear all; close all; clc;
addpath('/Users/levittcl/Documents/GitHub/UniversalCode');
addpath ('/Users/levittcl/Documents/GitHub/UniCode');

% where to save the timecourse
savepath = '/Users/levittcl/Documents/Research/Projects/Pseudo Islets/Analysis/2023_05_02 Human/Day 0/WT islet';

%% Load Images %%

% change this depending on how many channels you use
howmanychannel = 1;

disp('Select Imaging File');
R = bfopen('Select Calcium Imaging File'); % Uses bfopen program to open .czi/.lsm image files

pics=R{1};
pics=pics(:,1);

for i = 1:length(pics)
    IMG(:,:,i) = pics{i};
end 

ref_pic = IMG(:,:,1);

images=double(IMG); % converts images to double precision
images = images(:,:,1:howmanychannel:end);
RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable

clear pics R;

disp('Select Folder for Output Directory');
output_dir = uigetdir;
filename = input('Input Filename\n','s');

%% MODIFY TIMECOURSE %%

% This block converts the 3rd dimension of array into timecourse for later reference

st=500; % This dictates the Starting Frame (modify based on initial time lag)

[xx, yy, zz] = size(images(:,:,:));
%ed=zz;   
ed = 700;
T=st:1:ed;

images=images(:,:,st:ed);
RawImg=images(:,:,1);
MaxIntensityMap = max(images, [],3);

imagesRaw=images;

threshhold = 1.2; %threshhold for activity.
% Keep threshhold the same for entire project. good thresholds are between 1.75-2.
% Threshholds might change for different dyes, microscope, etc.

%% DECLARING IMAGE PROPERTIES %%
sx=size(images,1);
sy=size(images,2);
sz=size(images,3);

%% MASKING IMAGES %%

cellmask=[];
bkgcell=[];
fullmask=[];

%%%%%%Selecting "background cell" from peak amp map%%%%%%%%%%
AverageIntensityMap = mean(images, 3);
% figure; imagesc(AverageIntensityMap);
if isempty(bkgcell)
    figure;
    imagesc(AverageIntensityMap)
    title('Average Intensity Map');
    disp('select background cell')
    bkgcell=createMask(imfreehand);
else
    figure
    overlay = imoverlay(mat2gray(RawImg), bkgcell, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay)
    title('Background cell');
end
DataOut.bkgcell=bkgcell;

for i=1:sz
    imagesnosignal(:,:,i)=images(:,:,i).*bkgcell;
end

imagesnosignal(~find(imagesnosignal))=nan;
MeanNoSignalMask = squeeze(nanmean(nanmean(imagesnosignal,1),2));
DataOut.MeanNoSignalMask = MeanNoSignalMask;
MeanNoSignalMask0 = detrend(MeanNoSignalMask);
DataOut.MeanNoSignalMask0 = MeanNoSignalMask0;

averageamp = max(abs(MeanNoSignalMask0));

%% MASKING IMAGES

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
        cellmask = dummyMask{:,:}; % Assigns loaded mask file to CellMask Variable
        disp('Please Select Cell Number File');
        [num,npath] = uigetfile('*.mat','Please select Cell Number file',mpath); % User selects .mat cell number file - only valid for MatLab generated files!!!
        addpath(npath);
        dummyNum = struct2cell(load(num)); % Loads in Cell Number file
        numcells = dummyNum{:,:};
    case 'Draw New Masks'
        disp(answer)
    figure
    imshow(mat2gray(MaxIntensityMap));
    title('Max Intensity Map');
    cellmask=zeros(sx,sy, 75);
    fullmask = zeros(sx,sy);
    disp('select cells')
    done=0;
    i=1;
    while done < 1
            bf=imfreehand();
            useMask=createMask(bf);
            cellmask(:,:,i)=cellmask(:,:,i)+useMask;
            fullmask = fullmask+i*useMask;
            fullmask(find(fullmask>i)) = i;
            i=i+1;
            %%somehow need to deal with overlap
            answer = questdlg('More cells?','Yes','No');
            switch answer
                case 'Yes'
                    done = 0;
                case 'No'
                    done = 1;
            end
    end
 
     cellmasktemp=cellmask(:,:, 1:i-1);
     clear cellmask
     cellmask = cellmasktemp;

    figure
    cmap = hsv(max(max(fullmask)));
    overlay = imoverlay(mat2gray(RawImg), fullmask, 'colormap', cmap, 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay);
    caxis([1 max(max(fullmask))]);
    colormap(cmap)
    colorbar;
    title('All cells - May have overlap');
% else
%     figure
%     cmap = hsv(max(max(fullmask)));
%     overlay = imoverlay(mat2gray(RawImg), fullmask, 'colormap', cmap, 'facealpha', .75, 'ZeroAlpha', 0);
%     imshow(overlay);
%     colormap(cmap)
%     colorbar;
%     caxis([1 max(max(fullmask))]);
%     title('All cells - May have overlap');



end 

DataOut.mask=fullmask;
DataOut.cellmask = cellmask;
 h = fspecial('average', 5);
clear images

 for i=1:sz
    images(:,:,i)=imfilter(imagesRaw(:,:,i),h,'replicate');
 end

%% analyze each cell
for cell = 1:size(cellmask,3)

    indmask = cellmask(:,:,cell);
    for i=1:sz
        imagescell(:,:,i)=images(:,:,i).*indmask;
    end

    imagescell(~find(imagescell))=nan;
    Mean = squeeze(nanmean(nanmean(imagescell,1),2));
    DataOut.timecourse(cell,:) = Mean;
    DataOut.timecourseto0(cell,:) = detrend(Mean);

    % Duty Cycle
    tccell = Mean;
    tclogical = find(detrend(Mean) > threshhold*averageamp);
    DataOut.dutycycle(cell,1) = length(tclogical)/length(tccell);

    % Frequency
    s = Mean-mean(Mean); %use calcium to determine phase
    s = detrend(s);
    t = (0:length(s)-1); %time vector
    Ts = mean(diff(t)); %% Sampling period
    Fs = 1/Ts;  % Sampling Frequency
    Fn = Fs/2;  %Nyquist Frequency
    L  = length(s); %length of signal

    fts = fft(s-mean(s));
    course =s-mean(s);
    f = Fs*(0:(ceil(L/2)))/L;
    P2 = abs(fts/L);
    P1 = P2(1:(ceil(L/2)+1));
    P1(2:end-1) = 2*P1(2:end-1);

    amp_fts = abs(fts);      % Spectrum Amplitude
    phs_fts = angle(fts);
    highestfreq = f(find(P1==max(P1))); %this is the prob
    value = find(f==highestfreq);
    DataOut.freqvec(cell,1) = f(value);
    DataOut.amp(cell,1) = amp_fts(value);
    figure
    plot(t, Mean-mean(Mean), 'r')
    hold on
    plot((t),max(P1)*sin(2*pi*f(value) *(t)), 'b')
    title('Fourier Frequency')

end

figure
plot(DataOut.timecourseto0')
hold on
plot(ones(1,length(Mean))*threshhold*averageamp, 'k', 'LineWidth', 2);
hold on
plot(MeanNoSignalMask0','r', 'LineWidth', 2);
title('Duty Cycle')


figure
plot(DataOut.timecourse')
hold on
plot(ones(1,length(Mean))*threshhold*averageamp, 'k', 'LineWidth', 2);
hold on
plot(MeanNoSignalMask', 'r', 'LineWidth', 2);
title('Real Time Course')

% Save Structure

save([output_dir '\' filename '_' 'data_analysis.mat'], 'DataOut'); %Saves all DataOut info
save([output_dir '\' filename '_' 'fullmask' '.mat'],'fullmask'); %Saves CellMask array
save([output_dir '\' filename '_' 'multiplemasks' '.mat'],'cellmask'); %Saves CellMask array
save([output_dir '\' filename '_' 'CellNum' '.mat'],'cell'); %Saves cell number scalar for running in another code