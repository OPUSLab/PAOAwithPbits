clear all
clc
%% Script Name:  Coupling Strength and Node Weight Analysis in SK Model
%  Author:       Abdelrahman, Abdelrahman – OPUS Lab
%  Date:         May 19, 2025
%
% Description:
% This script analyzes the coupling strength of each node (p-bit) in the
% Sherrington–Kirkpatrick (SK) spin glass model with Lévy-distributed
% interaction bonds. The goal is to compute a "lambda score" that quantifies
% how strongly each node is coupled to all others:
%
%       λ_i = Σ_j |J_ij|
%
% where J_ij are the coupling weights.
%
% Workflow:
%   1. Reads all `.mat` files in the folder "training_instances".
%   2. Each file contains a coupling matrix J for a given SK instance.
%   3. For each instance:
%        - Loads J and converts it to a full matrix.
%        - Computes λ_i for each node as the sum of absolute couplings.
%        - Sorts nodes by coupling strength (indices adjusted for C++ use).
%   4. Stores λ_i values across all instances in `lambda_score_values`.
%   5. Plots the coupling strength distribution:
%        - Top 10 "heavy" nodes (high λ_i) shown in red diamonds.
%        - Remaining nodes shown in blue circles.
%
% Plot details:
%   - y-axis: logarithmic scale (λ_i)
%   - x-axis: sorted node indices
%   - Purpose: visualize the separation between strongly and weakly coupled
%     nodes across multiple random SK instances.
% ------------------------------------------------------------------------
%% calculate lambda scores (weigth of each node)

folder_path = "training_instances";
file_list = dir(fullfile(folder_path, '*.mat'));

num_pbits=  50;
lambda_score_values = zeros(num_pbits, length(file_list));

for j=1:length(file_list)
    file_path = fullfile(folder_path, file_list(j).name);
    instance_number = regexp(file_path, '\d?\d*', 'match');
    instance_number = str2double(instance_number); % Convert to numeric values
    instance_number = instance_number(end);

    J = load(file_path); 
    J=full(J.J);
    for row=1:length(J)
        lambda_score_values(row, instance_number) = sum(abs(J(row, :)));
    end
    lambda_score_values_per_instance = lambda_score_values(:,instance_number);
    [~, sorted_indices] = sort(lambda_score_values_per_instance, 'descend');
    sorted_indices = sorted_indices-1; % for C++ use
   % file_name = sprintf("%d.mat", instance_number);
   % save(file_name, "J")
   % 
end  


%% Plotting the coupling strength
figure;
heavy_nodes=10;
baseColor1 = [0.8500 0.3250 0.0980];
baseColor2 = [0 0.4470 0.7410];
sorted_scores = sort(lambda_score_values, 'descend');
figure(1)
scatter(1:heavy_nodes, sorted_scores(1:heavy_nodes,:), 80, ...
    baseColor1, 'd','filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)
hold on
scatter(heavy_nodes+1:50, sorted_scores(heavy_nodes+1:end,:), 60, ...
    baseColor2, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.2) % average across instances

ylabel('$\lambda_i = \sum_j |J_{ij}|$', 'Interpreter','latex')
xlabel("ith p-bit (sorted)")
set(gca, 'YScale', 'log', 'FontSize', 20)
xlim([0.9,50])

