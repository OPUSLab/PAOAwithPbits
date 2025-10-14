clear all
clc
%% Author: Abdelrahman, Abdelrahman-OPUS lab (Date: 19 May, 2025)
%  Purpose: study the coupling strength in SK model with Levy bonds
%% calculate lambda scores (weigth of each node)

% folder_path = "/Users/abdulrhmanseed/Downloads/Research/p-bits/Levy_SK_Model/Fast_PAOA_2p/training_instances";
folder_path = "/Users/abdulrhmanseed/Downloads/Research/p-bits/Levy_SK_Model/N_500_alpha_0.9";
file_list = dir(fullfile(folder_path, '*.mat'));

num_pbits=  500;
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
   
end  


%% 
figure;
heavy_nodes=100;
baseColor1 = [0.8500 0.3250 0.0980];
baseColor2 = [0 0.4470 0.7410];
sorted_scores = sort(lambda_score_values, 'descend');
figure(1)
scatter(1:heavy_nodes, sorted_scores(1:heavy_nodes,:), 80, ...
    baseColor1, 'd','filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)
hold on
scatter(heavy_nodes+1:500, sorted_scores(heavy_nodes+1:end,:), 60, ...
    baseColor2, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.2) % average across instances

ylabel('$\lambda_i = \sum_j |J_{ij}|$', 'Interpreter','latex')
xlabel("ith p-bit (sorted)")
set(gca, 'YScale', 'log', 'FontSize', 20)
% set(gca, 'XScale', 'log', 'FontSize', 20)

xlim([0.9,500])
%xticks([1, 5:50:500])
%xticklabels([1, 5:50:500]);
