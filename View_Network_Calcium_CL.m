%% Show network dot stick
% OG Code by Jennifer Briggs
% This script combined network analysis with hub visualization via
% dot-stick model (Claire 5.2023)

close all 
clear all
clc
addpath /Users/levittcl/Documents/GitHub/UniCode
addpath /Users/levittcl/Documents/GitHub/UniversalCode

% THINGS YOU CHANGE
filepath = '/Users/levittcl/Documents/Research/Projects/SST/Calcium Analysis/2023_8_30/islet 1 2mM 11mM/Normalizedwave/' %input directory where CaWaveForm.mat is 
%imagepath = '/Users/levittcl/Documents/Research/Projects/SST/Calcium Analysis/2023_8_30/islet 5 11mM/'%input directory where Imaging.mat is
savename = '/Users/levittcl/Documents/Research/Projects/SST/Calcium Analysis/2023_8_30/islet 1 2mM 11mM/Normalizedwave/'%input where to save the data
Thr = .8   %input correaltion threshold for network analysis
TitleChoice = 'HubAnalysis_5_30_islet3' %Input title choice here
%Input start and end time of video
% starttime = Vidinfo(illy).starttime(lg); 
% endtime = Vidinfo(illy).endtime(lg);

% Load things: 
load([filepath  'islet1_HG_CalciumWave_detrend.mat'])
%load([imagepath 'Hub_Analysis_09_13_2022_cytokine.mat'])
load([filepath 'islet 1 2mM 11mM_Masks.mat'])       %load masks
load([filepath 'islet 1 2mM 11mM_CellNumber.mat']) % cell number

% detrended_ca
CellTC = detrended_ca; 

% 
zstacks = 1;       %how many z stacks
zz = 1;           %where to start on z stack
cachannel = 2;      %where is the calcium channel
howmanychannel = 1;  %how many imaging channels
%ed = 1000;
st = 1;

% Load video
filename = '/Volumes/Claire Hard Drive/DATA/Sst Delta Cells/2023_08_29 sst FRAP/calcium/islet1_2mM_11mM.czi' %file name can direct to '.mat' analysis file or imaging file
R = bfopen([filename]); % Uses bfopen program to open .czi/.lsm image files

pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
pn = length(pics);

for i=1:pn
    IMG(:,:,i)=pics{i};
end

try
    for i=1:length(pics)
        T(i)=R{4}.getPlaneDeltaT(0, i-1).value;
    end
catch
    T=0:0.5:pn*0.5;
end
T = double(T);
T = T(cachannel:howmanychannel:end);
T = T(1:zstacks:end);

% if starttime == -1
%     st=1;
% else
%     st = starttime;
% end
% 
% if endtime == -1
%     ed=length(T);
% else
%     ed=endtime;
% end

T = T(st:end);
images=double(IMG); % converts images to double precision
images = images(:,:,cachannel:howmanychannel:end);
RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable
images = images(:,:,zz:zstacks:end); %change times if you don't want the whole time series
images = images(:,:,st:end-1);

sx=size(images,1);
sy=size(images,2);
sz=length(T);
for i=1:size(images,3)
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]); %applies filter to clean up images
end


ImAv = sum(images,3); %compresses all frames into single array of intensities
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
HSV(:,:,1) = 0.3333;%converts image to green image
RGB2 = hsv2rgb(HSV); %converts to rgb image

%%
for ll = 1:numcells
    [xp,yp] = find(CellMask == ll);
    x(ll) = mean(xp);
    y(ll) = mean(yp);
end

bad = find(isnan(x))

x(bad) = [];
y(bad) = [];
CellTC(:, bad) = [];
numcells = numcells - length(bad)
%% Run the network analysis
% taken from RunNetworkAnalysis script

%load cvs: 
%capath = '/Users/levittcl/Documents/Research/Projects/Pseudo Islets/Cx36 KO/CaWaveForm.mat' %%%% Add .mat path here
%ca = importdata(capath);
ca = CellTC;
st = 1;
%ed = 1674; %change or whatever
ca = ca(st:end,:);

Opts.figs = 1; %Set 1 if you want figures, 0 if not

[degree, Adj, kpercent, histArrayPercShort,histArrayShort, k, pval, Rij, s] = NetworkAnalysis(ca, Thr, Opts);


%% Show figure
fig = figure
imshow(RGB2)%, hold on

% Adj = corr(CellTC);
% Adj = Adj > Thr;
% Adj = Adj - diag(diag(Adj));
 AdjacencyGraph = graph(Adj);
% Conn = sum(Adj);
% 
% Conarray = [0:max(Conn)]; % adjust based on CONTROL
% Conarray2 = Conarray/length(Conarray)*100;
% try
%     hubth = Conarray(Conarray2 >= 59); hubth = min(hubth); %hub threshold
%     Hubby = find(Conn >= hubth);
a = find(kpercent > 60); hubthreshold = (a(1)); %Find degree threshold for hubs
hubthresh = hubthreshold - 1;
Hubs = find(degree>=hubthresh) %Find Hubs: Number of cells > 60% of Islet Links
    Nodec = repmat([.175 .54 .60],numcells,1); %color of nodes
    for lll = 1:length(Hubs) %load nodes 
        Nodec(Hubs(lll),:)=[.98 .122 .157];
    end
% catch
%     Nodec = repmat([.175 .54 .60],numcells,1);
% end

fig2 = figure
p = plot(AdjacencyGraph, 'Xdata',y,'YData',x, 'EdgeColor', 'w', 'NodeColor',Nodec,'MarkerSize',10, 'LineWidth',2)
p.NodeLabel = [];
set(gca, 'YDir','reverse')
title([TitleChoice])

set(gca, 'Visible', 'off')
set(gca, 'xticklabels', [])
set(gca, 'yticklabels', [])
% 
set(fig, 'Position', [100 100 1000 800])
set(fig2, 'Position', [100 100 1000 800])

saveas(fig, [savename '.png'])
saveas(fig2, [savename '.png'])

clearvars x y
close(fig)
close(fig2)

%% plotting HUB Cell timecourse 

hubs_only = ca(:,Hubs(:,:));

grayColor = [.7 .7 .7];

figure, plot (ca, 'Color', grayColor, 'LineWidth', 0.5)
hold on
plot(hubs_only, 'Color', 'b', 'LineWidth', 2)
title('HUB Cell Calcium Trace')

% what happens if we normalize...

ca_norm = normalize(ca);
hubs_only_norm = ca_norm(:,Hubs(:,:));

figure, plot (ca_norm, 'Color', grayColor, 'LineWidth', 0.5)
hold on
plot(hubs_only_norm, 'Color', 'b', 'LineWidth', 2)
title ('Normalized HUB Cell Calcium Trace')

%% Stats
cells = [1:1:width(ca)]'; %array of number of cells
deg = kpercent'; %probability of cells with X amount of links as a percent
cell_degs = histArrayPercShort'; %percent of cells
links = k'; % array of number of links
num_cells = histArrayShort'; % how many cells have X amount of links
coor_coeff = mean(Rij)'; %mean correlation coefficient per cell
abs_corr_coeff = abs(coor_coeff);

cell_stats = table(cells,coor_coeff,abs_corr_coeff,degree);
stats = table(links,num_cells,deg,cell_degs);

