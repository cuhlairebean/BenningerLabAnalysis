%Excel to CaWaveform%
%Claire H Levitt | 10 11 24

clear all; close all; clc

calcium = readtable ('/Users/levittcl/Desktop/load_matlab.xlsx');
calcium = table2array(calcium);

figure, plot(calcium)

output_dir = '/Users/levittcl/Documents/1_Research/4_Collaboration/2024_Ca for Dylan/Calcium Anlaysis/Mutant 1/';
save([output_dir '/CaWaveForm.mat'],'calcium'); % same cawaveform - Cell TC