clear all;
close all;
clc;


file_name = input('Please enter Raw Data file name between " " : \n');
size_of_table = size(readtable(file_name));

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = [30, size_of_table(1)];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = "IFRT";
opts.VariableTypes = "double";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
BGT60TR13Crecord = readtable(file_name, opts);

%% Convert to output type
BGT60TR13Crecord = table2array(BGT60TR13Crecord);

%% Clear temporary variables
clear opts

save("BGT60TR13Crecord.mat")