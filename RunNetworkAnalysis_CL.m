% Run Network analysis: 
% Jennifer Briggs 2022
% Modified Claire Levitt 2022
    % add extracting HUB cell timecourses

addpath('/Users/levittcl/Documents/GitHub/Lab_Analysis')

addpath('/Users/levittcl/Documents/GitHub/UniversalCode')

close all
clear all
clc

% SELECT LOCATION FOR SAVING FILES
savepath = ['/Users/levittcl/Documents/1_Research/7_Example Imaging Files/test/'];

%% OLD %%
%%% load MAT file of Calcium Wave: 
     % For extracting CaWaveForm.mat use Extracting_Individual_Cells code
%    % (Function: STD Analysis) 
% 
% capath = ['/Users/levittcl/Documents/Research/Projects/HUB_ANALYSIS/Consistency Analysis/07 01 2022/CaWaveForm.mat'] %%%% Add .mat path here
% ca = importdata(capath);
% %ca = readtable(capath);
% %ca = table2array(ca);
% 
% 
% 
% Threshold = 0.90 %Here you put the correlation threshold to draw an edge
% Opts.figs = 1 %Set 1 if you want figures, 0 if not
% 
% [degree, Adj, kpercent, histArrayPercShort, histArrayShort, k, pval,Rij,s] = NetworkAnalysis(ca, Threshold, Opts)
% 
% %Rij - correlation matrix
% remove diagonal (
% avcorr = mean2(nonzeros(Rij-diag(diag(Rij))));
% 
% 
% a = find(kpercent >= 60); hubthreshold = (a(1)); %Find degree threshold for hubs
% Hubs = find(degree>=hubthreshold) %List each hug <3

%% NEW %%

%load csv: 
capath = ['/Users/levittcl/Documents/1_Research/7_Example Imaging Files/test/CaWaveForm.mat']; %%%% Add .mat path here
Data = importdata(capath);

%% manually select timecourse

IsletTC = Data;
IsletTCfig = figure('Name','Whole Islet Time Course');
plot(IsletTC); %plots intensity timecourse of entire islet so user can identify first responder and wave origin ranges
datacursormode on % User can click on the plot to see x and y coordinates
st = input('Input Starting Frame for First Responder Analysis\n'); %select where first responder analysis should begin (x coordinate value)
ed = input('Input Ending Frame for First Responder Analysis\n'); % select where first responder analysis should end (x coordinate value)

%% preset timecourse %%
%st =1;
  %ed =300;
% 
%ca = Data;
%ca = Data.data(:,:); % for csv only
% 
 ca = IsletTC(st:ed,:);

%To set the threshold, either manually set: 

%1) ----------
Threshold = 0.9 %Here you put the correlation threshold to draw an edge

%Or determine the threshold using the average degree or a threshold that
%gives you a scale free degree distribution

%2) --------------
% Opts.Method = 'Degree'
% Opts.avDeg = 6; %Set the average degree here
% 
% Opts.Method = 'Scale-Free' 
% %%Set the bounds for max and min average degree for scale free: 
% Opts.Max = 50
% Opts.Min = 2
%NOW RUN
%Threshold = findoptRth(ca, Opts)

Opts.figs = 1; %Set 1 if you want figures, 0 if not

[degree, Adj, kpercent, histArrayPercShort,histArrayShort, k, pval, Rij, s, possible] = NetworkAnalysis(ca, Threshold, Opts);


% [N, Adj, kpercent, histArrayPercShort,histArrayShort,k, pval,Rij,s]; pval,Rij,s, k, histArrayShort]

a = find(kpercent > 60); hubthreshold = (a(1)); %Find degree threshold for hubs
hubthresh = hubthreshold - 1;
Hubs = find(degree>=hubthresh) %List each hug <3

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
avg_coorelation = mean(coor_coeff)
% Y = [avg_correlation, '^ average absolute value coorelation coeff'];
% disp(Y)

degree = degree + 1;

cell_stats = table(cells,coor_coeff,degree);
stats = table(links,num_cells,deg,cell_degs);
save([savepath '/correlation_analysis' '.mat'],'cell_stats'); %Saves cell stats (correlation)
save([savepath '/hub_analysis' '.mat'],'stats'); %Saves stats (hub cell info)