clear all
clc
%% Two Schedules (PAOA)
figure(1)
final_beta = 0.012;
delta_schedule = 0.8;  % difference in two schedules
num_layers = 17;
main_schedule = exp(linspace(log(0.0001),log(final_beta),num_layers)); % mean schedule
up_schedule = main_schedule * (1 + delta_schedule);   % schedule with higher beta
steps = 8:1:num_layers;
basecolor1 = [0 0.4470 0.7410];
color1 = [basecolor1, 0.8];

basecolor_heavy = [0.8500 0.3250 0.0980];
color_heavy = [basecolor_heavy, 0.8];

base_color3 = [0.4660, 0.6740, 0.1880];
colro3 = [0.4660, 0.6740, 0.1880, 0.8];


plot(1:17,main_schedule, '-','LineWidth', 2, 'MarkerSize', 10, 'Color',color_heavy);
hold on
scatter(1:17,main_schedule, 100, ...
    basecolor_heavy,'d', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)


hold on;
plot(1:17,up_schedule, '-','LineWidth', 2, 'MarkerSize', 10, 'Color',color1);
hold on
scatter(1:17,up_schedule, 100, ...
    basecolor1, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)


xlabel('Number of layers, p');
ylabel('Inverse temperature, \beta');
legend({'\beta_2=\beta^{opt}', '\beta_1 = \beta^{opt} (1+\Delta)'},'box','off', Location='northwest');
ax = gca; ax.FontSize=20;
xlim([0.01,17.5])
xticks([1:4:17])
xticklabels([1:4:17]);
hold off;




%%
final_beta = 0.011;
num_pbits = 500;
delta_schedule = 0.005;
num_layers = 17;

main_schedule_heaviest_node = exp(linspace(log(0.0001), log(final_beta), num_layers));
num_heavy_nodes = 1/5 * num_pbits;
beta_schedules_heavy_nodes = zeros(num_heavy_nodes, num_layers);
beta_schedules_heavy_nodes(1, :) = main_schedule_heaviest_node;

for heavy_nodes_idx = 1:99
    beta_schedules_heavy_nodes(heavy_nodes_idx + 1, :) = ...
        main_schedule_heaviest_node * (1 + delta_schedule)^heavy_nodes_idx;
end

light_nodes_schedule = beta_schedules_heavy_nodes(end, :) * (1 + delta_schedule);
num_light_nodes = 4/5 * num_pbits;
light_nodes_schedule_replicated = repmat(light_nodes_schedule, num_light_nodes, 1);

heavy_light_schedules = [beta_schedules_heavy_nodes; light_nodes_schedule_replicated];

% Define colors
basecolor_heavy = [0.8500 0.3250 0.0980];
basecolor1 = [0 0.4470 0.7410];

% Plot
figure;
hold on;

% Plot heavy nodes
for i = 1:num_heavy_nodes
    plot(1:num_layers, beta_schedules_heavy_nodes(i,:), 'd-', ...
        'Color', [basecolor_heavy 0.8], 'LineWidth', 0.5, 'MarkerFaceColor', basecolor_heavy, 'MarkerSize',7);
end

% Plot light node schedule
plot(1:num_layers, light_nodes_schedule, 'o-', ...
    'Color', [basecolor1 0.8], 'LineWidth', 0.5, 'MarkerFaceColor', basecolor1, 'MarkerSize',7);

xlabel('Layer');
ylabel('Beta');
% yscale('log')
title('Beta Schedules for Heavy and Light Nodes');
ax = gca; ax.FontSize=20;
% legend({'Heavy Node', 'Light Node'});
xlim([0.01,17.5])


xticks([1:4:17])
xticklabels([1:4:17]);
hold off;
hold off;
