% DUTY CYCLE, FREQUENCY, and AMPLITUDE
% This script calculates duty cycle as a proportion of time the cell is
% considered " ON " relative to the *entire* timecourse


% 07 04 2023: ADDED Amplitude and Frequency Analysis

close all; clear all; clc;

%% Load Timecourse

% directory where CalciumWaveForm is located
filepath = ('/Users/ibarevel/Documents/MATLAB/EndoCaImaging/OutputDirectoryEndoCaImaging/');
javaaddpath('/Users/ibarevel/Documents/MATLAB/bioformats_package.jar');
filepath = ('/Users/ibarevel/Documents/MATLAB/EndoCaImaging/OutputDirectoryEndoCaImaging/');

javaaddpath('/Users/ibarevel/Documents/MATLAB/bioformats_package.jar');
load([filepath  'CaWaveForm_Endo3oht1_11-17.mat'])


%% Moving Average Smoothing
neighbors = 3; % Number of neighbors on each side
windowSize = 2 * neighbors + 1; % Total window size

% Apply moving average to each column
SmoothedCellTC = zeros(size(CellTC)); % Preallocate array
for i = 1:size(CellTC, 2)
    SmoothedCellTC(:, i) = movmean(CellTC(:, i), windowSize);
end

% Update CellTC with smoothed data
CellTC = SmoothedCellTC;
% change time frame
CellTC = CellTC(1840:2081,:);

% disp('Please Select Background Threshold File');
%         [threshhold,mpath] = uigetfile('*.mat','Please select thresh file'); % User selects .mat mask file - only valid for MatLab generated files!!!
%         addpath(mpath);

%% Duty Cycle Analysis
% Assuming some level of consistency in duty cycle across timecourse

for cell = 1:length(CellTC(1,:))

    Mean = CellTC(:,cell);
    thresh = 0.3;
    
    %Detrend and Normalize to find avgerage amplitude
    tccell = Mean;
    det = detrend(tccell);
    norm = normalize(det);
    %avg_amp = mean(tccell);
    avg_amp = mean(norm);

    % determine what amplitude is considered "on"
    on = thresh; % correct for if avg_amp is slightly biased
    %on = threshhold; % <-- considering avg = 0
%     on = avg_amp;

%     % see what this looks like
     figure, 
     plot (norm)
     yline(avg_amp,'-','average amplitude')
     yline(on,'--','threshold')
    
    % find what percent of the time the cell is "on"
    logical = find(norm >= on);
    duty_cycle_M1(cell,1) = length(logical)/length(tccell);
end 

%save([output_dir '/DutyCycle.mat'],'duty_cycle'); % save duty cycle vals
%%%
%%%
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
% 
% %% Amplitude and Frequency Analysis
% 
% 
% for cell = 1:length(CellTC(1,:))
% 
%     % Duty Cycle
%     Mean = CellTC(:,cell);
%     tccell = Mean;
%     tclogical = find(detrend(Mean) > thresh*avg_amp);
%     DataOut.dutycycle(cell,1) = length(tclogical)/length(tccell);
% 
%      % Frequency
%     s = Mean-mean(Mean); %use calcium to determine phase
%     s = detrend(s);
%     t = (0:length(s)-1); %time vector
%     Ts = mean(diff(t)); %% Sampling period
%     Fs = 1/Ts;  % Sampling Frequency
%     Fn = Fs/2;  %Nyquist Frequency
%     L  = length(s); %length of signal
% 
%     fts = fft(s-mean(s));
%     course =s-mean(s);
%     f = Fs*(0:(ceil(L/2)))/L;
%     P2 = abs(fts/L);
%     P1 = P2(1:(ceil(L/2)+1));
%     P1(2:end-1) = 2*P1(2:end-1);
% 
%     amp_fts = abs(fts);      % Spectrum Amplitude
%     phs_fts = angle(fts);
%     highestfreq = f(find(P1==max(P1))); %this is the prob
%     value = find(f==highestfreq);
%     DataOut.freqvec(cell,1) = f(value);
%     DataOut.amp(cell,1) = amp_fts(value);
% %     figure
% %     plot(t, Mean-mean(Mean), 'r')
% %     hold on
% %     plot((t),max(P1)*sin(2*pi*f(value) *(t)), 'b')
% %     title('Fourier Frequency')
% 
% end 
%%% 

%%%
%Calculate the normalization of the timeseries




%%%
