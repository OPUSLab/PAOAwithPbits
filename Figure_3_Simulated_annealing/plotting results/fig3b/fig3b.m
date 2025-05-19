clear all
clc
%% Author: Abdelrahman, Abdelrahman-OPUS Lab (Date: May 18, 2025)
% Purpose: plot the energy histogram for a randomly selected optimization run-figure 3b
%%
load("all_costs_constant_beta_L6_720.mat")
load("all_costs_L6_1e5_720_avg_representative.mat")
for i=1:3
    
% Example data
costs_for_all_layers_PAOA_single = cost_all_720(i,:)';
costs_for_all_layers_SA_single = cost_all_720_SA(i,:)';


% Define bin edges for the histograms
binEdges = linspace(min([costs_for_all_layers_PAOA_single; costs_for_all_layers_SA_single]), ...
                    max([costs_for_all_layers_PAOA_single; costs_for_all_layers_SA_single]), 30);

% Compute histogram counts for each dataset
counts1 = histcounts(costs_for_all_layers_PAOA_single, binEdges);
counts2 = histcounts(costs_for_all_layers_SA_single, binEdges);


% Normalize counts to probabilities
binWidth = binEdges(2) - binEdges(1); % Calculate the width of each bin
prob1 = counts1 / (sum(counts1)); % Normalize to probability
prob2 = counts2 / (sum(counts2) );

% Compute bin centers for plotting
binCenters = binEdges(1:end-1) + diff(binEdges) / 2;


probMatrix = [prob1; prob2]'; % Each column represents one dataset

figure;
bar(binEdges(1:end-1), probMatrix,1, 'grouped'); % Use 'grouped' to display side-by-side bars
hold on
plot(-360*ones(20,1), linspace(0, 1, 20), 'k--', LineWidth=3)
xlabel('3D-Spin Glass System Energy');
ylabel('Occupation Probability');
title_dynamic = sprintf("p=%d", i*5);
title(title_dynamic);
ax = gca; ax.FontSize=20; ax.FontName='arial';
xlim([-362, -340])
ylim([1e-2 1])
yscale(ax,'log')
legend("Optimied $\beta$", "Initial $\beta$", 'interpreter', 'latex');
end
