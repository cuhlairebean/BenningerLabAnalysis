%Excel to CaWaveform%
%Claire H Levitt | 10 11 24

clear all; close all; clc

calcium = readtable ('/Users/levittcl/Desktop/Ca_Load.xlsx');
calcium = table2array(calcium);

figure, plot(calcium)

output_dir = '/Users/levittcl/Documents/1_Research/1_Research Projects/3_Pseudo Islet Platform/1_Human Pseudoislet Analysis/Donor Data/Donor 6 _ 12_20_2023/Calcium/Day9_control/islet3_HG/new_9_25_25/';
save([output_dir '/CaWaveForm.mat'],'calcium'); % same cawaveform - Cell TC