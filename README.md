# BenningerLabAnalysis
Calcium and Islet Analysis
Calcium_Timecourse_1channel_CL.m
•	Use for a single channel (i.e. Fluo4 or GCaMP)
•	This code is designed to extract timecourses of fluorescent signals
•	For multiple channels: skip down to 'Masking Islet > For Colocalization' and change the parameters as needed
•	To use the STD analysis (takes out areas of ROI that might skew signal pixel by pixel) uncomment Line 246sn

Input: czi
Output: CaWaveform  this can be used for network and characteristic analysis

Calcium_Timecourse_2channels_CL.m
•	This code is designed to extract timecourses of fluorescent signals
•	For multiple channels (DELTA CELLS EXAMPLE);  all 2nd channel properties are labeled "delta cell" comment this out if it doens't exist
•	LINE 59-60, LINE 204-207
•	Otherwise: when circling cells, to see next channel 3 = SWITCH channels
•	Code is designed to extract the 3rd channel; but if it was not imaged together, need to upload colocalization snap separate. LINE 175
•	modify start time: LINE 83

Input: czi
Output: CaWaveform  this can be used for network and characteristic analysis

Calcium_Timecourse_2channels_separate_CL.m
•	This code is designed to extract timecourses of fluorescent signals
•	it was designed for 2 channels and **plots them both separately**
•	Second channel labeled 'delta' (change as necessary); this was designed specifically for live endogenous labels with calcium on a single ca channel

Input: czi
Output: mixedCaWaveForm.m  this can be used for network and characteristic analysis; just make sure the file import reflects this name

Calcium_2CaChannels.m
•	This code is designed to extract timecourses of fluorescent signals from 2 separate channels in the same file – both live calcium imaging
•	Ex. GCaMP + Rhod2
•	Both channels merged as 1 waveform in the end

Input: czi
Output: mixedCaWaveForm.m  this can be used for network and characteristic analysis; just make sure the file import reflects this name
excel_to_mat.m
•	This code is for converting excel ROI timetraces to readable matlab array for further analysis

Input: xlsx
Output: CaWaveForm.m

%%%%%%%%% Functional Analysis %%%%%%%%%%%%%%

1.	Area Active

CellCalcium_MainFile.m
•	This is the file that will run a series of functions to analyze calcium activity and coordination relative to islet area
•	Uses: assessing overall coordinated activity of an islet, assessing activity during entrainment

Input: czi
Output: CorrelationAnalysis (Ratio Active), Activity Map

Run_Nosilentcell.m 
•	[FUNCTION] 
•	Output: [DataOut]
•	Input: data, useGFPasMaskBoolean, samefile, GFPsecond, st, ed, pathname
o	*keep in order*
o	All this info is input in the mainfile

2.	Network Analysis – 2nd phase

RunNetworkAnalysis_CL.m
•	Main file
•	This script will be used to assess calcium network characteristics:
o	Correlation Analysis
o	Hub Cell Analysis (can either change based on number of links or threshold)

Input: CaWaveform.m
Output: cell states (cells, correlation coefficient, number of links); stats overall (hub cell properties)

NetworkAnalysis.m
•	[FUNCTION]
•	Input: calcium wave, threshold (for being considered correlated), Opts
•	Output: degree, correlation, percent total links, total percent of links
WS_FrstResp_Wave_Analysis_Claire.m
•	This is a very user-friendly code to assess first responder cells based on half max time and wave ‘leader’ cells at 2nd phase

Input: czi
Output: first responder excel sheet – requires additional analysis to rank cells
![image](https://github.com/user-attachments/assets/06ec6662-fd30-4be2-aac4-9b918ea94052)
