clear all
clc
%% Author: Abdelrahman, Abdelrahman-OPUS Lab (Date: May 18, 2025)
% Purpose: plot the cost over the course of optimization-figure 2b
%% plot the <cost> as a function of iterations
load("cost_history.mat")
load('costHistory_analysis.mat') 
absolute_min = -8*log(1/8); % each desired state has ~100/8 probability  
color1 = [169/255, 117/255, 49/255];
color2 = [248/255, 161/255, 0/255]; 

figure;
subplot(2,1,1);
semilogx(1:length(cost_vs_iterations_analysis), (cost_vs_iterations_analysis-absolute_min),...
    LineWidth=3, Color=color2)
title("Gradient-based optimization")
ax = gca; ax.FontSize=20; ax.FontName = 'Arial';
xlim([0, 300])
ylim([0,23])

subplot(2,1,2);
semilogx(1:length(cost_history), (cost_history-absolute_min),'-',...
    LineWidth=3, Color=color1)
title("Gradient-free optimization")
ax = gca; ax.FontSize=20; ax.FontName = 'Arial';
xlim([0, 300])
ylim([0,23])
