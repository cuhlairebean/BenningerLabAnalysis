% DUTY CYCLE, FREQUENCY, and AMPLITUDE
% This script calculates duty cycle as a proportion of time the cell is
% considered " ON " relative to the *entire* timecourse

close all; clear all; clc

% SELECT LOCATION FOR SAVING FILES
savepath = ['/Users/levittcl/Documents/1_Research/1_Research Projects/3_Pseudo Islet Platform/1_Human Pseudoislet Analysis/Donor Data/Donor 27 8_5_25/Untreated Pseudo/Islet1_1mm/'];

% 07 04 2023: ADDED Amplitude and Frequency Analysis
% 1 2026: New Duty Cycle


%% Load Timecourse

% directory where CalciumWaveForm is located
filepath = ['/Users/levittcl/Documents/1_Research/1_Research Projects/3_Pseudo Islet Platform/1_Human Pseudoislet Analysis/Donor Data/Donor 27 8_5_25/Untreated Pseudo/Islet1_1mm/'];

% Load the thing
load([filepath  'CaWaveForm.mat'])

%CellTC = calcium;

% change time frame
CellTC = CellTC(1:end,:);

% disp('Please Select Background Threshold File');
%         [threshhold,mpath] = uigetfile('*.mat','Please select thresh file'); % User selects .mat mask file - only valid for MatLab generated files!!!
%         addpath(mpath);

%% DUTY CYCLE ANALYSIS

for cell = 1:length(CellTC(1,:))

    tccell = CellTC(:,cell);
    
    %Detrend and Normalize to find avgerage amplitude
    det = detrend(tccell);
    %det = smooth(det);
    detrended_ca(:,cell) = det;

   % FLAT TRACE CHECK - make sure cell is truly active

    minStd = 2;   % adjust if needed; minimum standard deviation of detrended timetrace
    if std(det) < minStd % if less than stdev then it is 0
        DataOut.dutycycle(cell,1) = 0;
        continue
    end

    % ----- PEAK DETECTION -----
    [peaks,locs] = findpeaks(det,'MinPeakProminence', std(det)); %only peaks within a standard dev of specified val


    % If no real peaks, cell is OFF
    if isempty(peaks)
        DataOut.dutycycle(cell,1) = 0;
        continue
    end

    % determine what amplitude is considered "on"
    avg_peakamp = mean(peaks);
    on = .3*avg_peakamp; % correct for if avg_amp is slightly biased
    %on = threshhold; % <-- considering avg = 0
%     on = avg_amp;
    
    % find what percent of the time the cell is "on"
    tclogical = det > on;
 
    DataOut.dutycycle(cell,1) = sum(tclogical) / length(det);

    %     % see what this looks like
    figure
    plot(det)
    yline(on,'--','threshold')
    hold on
    plot(CellTC(:,cell))
    title(sprintf('Cell %d | Duty = %.2f', cell, DataOut.dutycycle(cell)))
    % Overlay frequency as text
end 

%save([output_dir '/DutyCycle.mat'],'duty_cycle'); % save duty cycle vals

% %% METHOD 2 - matlab function
% 
% for cell = 1:length(CellTC(1,:))
% 
%     Mean = CellTC(:,cell);
%     
%     %Detrend and Normalize to set avgerage amplitude 
%     tccell = Mean;
% %     det = detrend(tccell);
% %     norm = normalize(det);
%  
%     [d,initCross,finalCross,nextCross,midRef] = dutycycle(tccell);
%     avg_d = mean(d);
% 
%     % take the average to generalize
%     duty_cycle_M2(cell) = avg_d;
% 
%         % see what this looks like
% %     figure, 
% %     plot (norm)
%   
% end 
% 
% duty_cycle_M2(isnan(duty_cycle_M2)) = 0;
% 
% %% Display avg results for both
% 
% M1 = duty_cycle_M1;
% M2 = duty_cycle_M2';
% T = table(M1, M2);

%% Amplitude and Frequency Analysis


for cell = 1:length(CellTC(1,:))
    
    Mean = CellTC(:,cell);
    tccell = Mean;
%     tclogical = find(detrend(Mean) > thresh*avg_amp);
%     DataOut.dutycycle(cell,1) = length(tclogical)/length(tccell);

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
%     figure
%     plot(t, Mean-mean(Mean), 'r')
%     hold on
%     plot((t),max(P1)*sin(2*pi*f(value) *(t)), 'b')
%     title('Fourier Frequency')

end 

%% Correlation Analysis

%active_cells = find(M2);
%Data = CellTC(:,active_cells);
Data = CellTC;

%% manually select timecourse

IsletTC = Data;
IsletTCfig = figure('Name','Whole Islet Time Course');
plot(IsletTC);
ca = CellTC;
Threshold = 0.9
Opts.figs = 1; %Set 1 if you want figures, 0 if not

%% Network Analysis Function

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
save([savepath '/dutycycle' '.mat'], 'DataOut'); %save duty cycle info