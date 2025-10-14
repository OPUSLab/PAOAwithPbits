clear all
clc
%% ------------------------------------------------------------------------
% Script Name:  Average Optimized Beta Schedules Plot
% Author:       Abdelrahman, Abdelrhaman
% Date:         October 10, 2025
%
% Description:
% This script reads multiple text files containing optimized beta (inverse 
% temperature) schedules from the folder "optimized_beta". Each file 
% corresponds to an optimization instance and includes two annealing 
% schedules (for two layers). The script:
%
%   1. Scans the folder and reads all .txt files.
%   2. Extracts numerical data from each file, reshaping it into a matrix 
%      of size 17×2 (17 layers × 2 schedules).
%   3. Stores all instances in a cell array for averaging.
%   4. Computes the mean schedule over all instances.
%   5. Plots the averaged schedules for the two layers alongside a 
%      reference exponential schedule defined as:
%           optimal_schedule = exp(linspace(log(0.005), log(0.12), 17))
%
% The plot displays:
%   - Light node schedule (blue)
%   - Heavy schedule (orange)
%   - Reference exponential schedule (green)
%
% Axes are labeled with layer indices and inverse temperature values.
%%
% Define the folder path
folder_path = "optimized_beta";
% Get a list of all .mat files in the folder
file_list = dir(fullfile(folder_path, '*.txt'));

% Initialize a matrix to store all the data
num_files = length(file_list);
data_matrix = []; % Each row will correspond to one file
num_rows = 1;
result = cell(length(file_list), 1);
allInstances = cell(num_files,1);
% Loop through each file and load the data
for i = 1:length(file_list)
    % Load the .mat file
    file_path = fullfile(folder_path, file_list(i).name);
        
    fileID = fopen(file_path, 'r');

    % Read all lines of the file
    data = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    lines = data{1};
    

    instance_number = regexp(file_path, '\d?\d*', 'match'); 
    instance_number = str2double(instance_number); % Convert to numeric values
    instance_number = instance_number(end);


    layers = cell(1, 1);  % to store the 3 schedules layers
    idx=1;
    for p = 17:17
        vector = sscanf(lines{idx}, '%f');
        
        % Reshape into a p x 2 matrix (columns: schedule 1 and 2)
        layers{idx} = reshape(vector, [p, 2]);
        idx=  idx+1;
    end
    
    allInstances{instance_number} = layers;  % store the 12-layer cell in the instance cell
end
%% Get the average
optimal_schedule = exp(linspace(log(0.005),log(0.12),17));
basecolor1 = [0 0.4470 0.7410];
basecolor_heavy = [0.8500 0.3250 0.0980];
color1 = [basecolor1, 0.8];
color_heavy = [basecolor_heavy, 0.8];
base_color3 = [0.4660, 0.6740, 0.1880];
colro3 = [0.4660, 0.6740, 0.1880, 0.8];
x = 0;
for i=1:num_files
x=x+allInstances{i}{1};
end
x= x./i;
plot(1:1:17,x(:,2), '-','LineWidth', 2, 'MarkerSize', 10, 'Color',color1);
hold on
scatter(1:1:17,x(:,2), 100, ...
    basecolor1, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)

plot(1:1:17,x(:,1), '-','LineWidth', 2, 'MarkerSize', 10, 'Color',color_heavy);

scatter(1:1:17, x(:,1), 100, ...
    basecolor_heavy, 'd', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)

plot(1:1:17,optimal_schedule, '-','LineWidth', 2, 'MarkerSize', 10, 'Color',colro3);

scatter(1:1:17, optimal_schedule, 100, ...
    base_color3, 's', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)
ax = gca; ax.FontSize=20; ax.FontName= 'arial';
xlim([0.01,17.5])
xticks([1:4:17])
xticklabels([1:4:17]);

xlabel("number of layers (p)")
ylabel("Inverse temp")
%% Extrapolation 

ratio_last_two = mean(x(end-1:end, 2) ./ x(end-1:end, 1));
num_points = 10000;
base_schedule = exp(linspace(log(0.005), log(0.12), num_points));
scaled_schedule = ratio_last_two * base_schedule;

% --- Step 3: Define colors (from PAOA two-schedule figure)
basecolor1       = [0 0.4470 0.7410];       % blue
color1           = [basecolor1, 0.8];
basecolor_heavy  = [0.8500 0.3250 0.0980];   % orange/red
color_heavy      = [basecolor_heavy, 0.8];
base_color3      = [0.4660, 0.6740, 0.1880]; % green
color3           = [base_color3, 0.8];

% --- Step 4: Plot extrapolation
figure('Color', 'w'); hold on; grid on;
plot(1:num_points, base_schedule, '-', 'LineWidth', 2.5, ...
    'Color', color_heavy, 'DisplayName', '\beta^{opt}');
plot(1:num_points, scaled_schedule, '-', 'LineWidth', 2.5, ...
    'Color', color1, 'DisplayName', '\beta^{opt} scaled');

xlabel('Number of layers, p');
ylabel('Inverse temperature, \beta');
title('Extrapolated Schedules up to 10,000 Layers');
legend('Location','northwest', 'Box','off');
ax = gca;
ax.FontSize = 20;
ax.FontName = 'arial';
xlim([1 num_points]);
grid on;