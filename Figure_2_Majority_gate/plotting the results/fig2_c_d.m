clear all
clc
%% Author: Abdelrahman, Abdelrahman-OPUS Lab (Date: May 18, 2025)
% Purpose: plot the histogram for figure 2c,d

%% Plot the evolution of the histograms Analytical (histogram as a function of time)
load("hist_over_time.mat")
color2 = [248/255, 161/255, 0/255]; 



for step=1:3
    figure(step)
    b = bar(hist_over_time(:,step), 0.8);
    b.FaceColor = color2;
    xticks(1:5:32); % Show all numbers from 0 to 31 on x-axis
    xticklabels(string(0:5:31)); % Label each tick with numbers 0 to 31
    ax = gca; ax.FontSize=30; ax.FontName='arial';
    ylim([0, 0.16])
end

%% Plot the histogram for numerical PAOA (histogram as a function of time)

% 1) Read all numbers from the file as a single vector:
data = dlmread('optimized_beta.txt'); 


% 2) Reshape so that each block of 5 values is a row:
beta_optimized_final = reshape(data, 2, []).';  
beta_optimized_final = round(beta_optimized_final, 1);

load("J_weight.mat")

color1 = [169/255, 117/255, 49/255];
time_steps = 3; % N+1 % 4 is also available... 
num_experiments = 1e7;
PDF_evolution_numerical = PDF_time_evolution(num_experiments,...
                                             beta_optimized_final,J, time_steps);
%% plotting
for step=1:3
    figure(step)
    b = bar(PDF_evolution_numerical(:,step), 0.8);
    b.FaceColor = color1;
    
    xticks(1:5:32); % Show all numbers from 0 to 31 on x-axis
    xticklabels(string(0:5:31)); % Label each tick with numbers 0 to 31
    ax = gca; ax.FontSize=30; ax.FontName='arial';
    ylim([0, 0.16])
end