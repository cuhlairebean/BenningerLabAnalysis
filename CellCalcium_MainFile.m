% Peak Amplitude Map
% Single Cell Coordination

clear all
close all
clc

addpath /Users/levittcl/Documents/GitHub/Calcium-Analysis
addpath /Users/levittcl/Documents/GitHub/Calcium-Analysis/CalciumAnalysisMaster

data.Location = ['/Volumes/CHL2021/All Microscope Data/2023_11_19 GCaMP/control_11mM_islet3.czi'];
% Correlation Plot; Max Intensity and Activity Map %
samefile = 0; % is calcium with another channel? 1-yes 0-no
cafirst = 1; %is calcium the first channel? 1-yes 0-no, if only 1 channel doesnt matter what it equals
GFP = 0; %GFP Second (boolean yes no)
redobackground = 1; % 1 = yes; 0 = no
threshold = 1.8;  %threshhold for activity. %1.9 OR 1.8
start_frame = 1;
end_frame = 599;
% Keep threshhold the same for entire project. good thresholds are between 1.75-2. 
% Threshholds might change for different dyes, microscope, etc.

pathname = ['/Users/levittcl/Desktop/test/CorrelationAnalysis.mat'];

%[DataOut] = Run_Chr2_lowglucose(data,0,samefile,cafirst); %Run_Chr2_lowglucose_CHL is same
[DataOut] = Run_Nosilentcell(data, 0, samefile, cafirst,start_frame,end_frame,pathname); %(data,useGFPasMaskBoolean, samefile, GFPsecond)
figure; imagesc(DataOut.CoorMap); colorbar; title('Correlation Map')

%save(pathname,"DataOut",'-v7.3')

[SingleCellData] = Run_TimecourseExtraction (data.Location, start_frame, end_frame)