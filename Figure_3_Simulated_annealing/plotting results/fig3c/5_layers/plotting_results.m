clear all
clc
%% Author: Abdelrahman, Abdelrahman-OPUS Lab (Date: May 18, 2025)
% Purpose: plot the success probability histograms & the corresponding
% schedules -figure 3c (******** 5 layers *********)
%%

folderPath = 'fig3c/5_layers/optimized_params'; % put the absolute path....
filePattern = fullfile(folderPath, '*.mat');
fileList = dir(filePattern);

% Initialize the cell array to store beta values for each file
optimized_beta_all = zeros(100, 5);

% Loop through each file in the directory
for fileIdx = 1:length(fileList)
    % Construct full file path
    filename = fullfile(fileList(fileIdx).folder, fileList(fileIdx).name);
    disp(filename)
    load(filename)   
    fileNumber = regexp(fileList(fileIdx).name,'(?<=run_)\d+', 'match');
    fileNumber = str2double(fileNumber{1}); % Convert to numeric format
   % Append the beta values of the current file to the main cell array
    optimized_beta_all(fileNumber,:) = optimized_beta_schedule{1};
   
    
    
end

% Display results
disp(optimized_beta_all);
 
schedule_5 = optimized_beta_all;
num_sweeps_layer = 720;

%% Processing Success prob.

folderPath = 'success_prob';
filePattern = fullfile(folderPath, '*.mat');
fileList = dir(filePattern);

% Initialize the cell array to store beta values for each file
success_prob = zeros(100, 1);


% Loop through each file in the directory
for fileIdx = 1:length(fileList)
    % Construct full file path
    filename = fullfile(fileList(fileIdx).folder, fileList(fileIdx).name);
    disp(filename)
    load(filename)
    fileNumber = regexp(fileList(fileIdx).name,'(?<=run_)\d+', 'match');
    fileNumber = str2double(fileNumber{1}); % Convert to numeric format
    success_prob(fileNumber) = success_prob_PAOA;
end

% Display results
disp(success_prob);
success_5 = success_prob;
avg_all  = mean(success_5);
%% Analysis on the schedules
n=3; 
sorted_success_prob = sort(success_prob, 'descend');
highest_n_success = sorted_success_prob(1:n);
% best n performing schedules -- Not efficient
schedule_5_best_n = zeros(n, 5);


for i=1:size(highest_n_success,1) % iterate over rows
    schedule_5_best_n(i,:) = schedule_5(highest_n_success(i)==success_prob(:),:); 
end

%% Plot all success probabilities over the optimizer runs-- Histograms

figure(1)
bar(sorted_success_prob); ax=gca; ax.FontSize=20;
xlabel('Optimizer runs');
ylabel('Success Probability');
title("p=5")
legend(arrayfun(@(x) sprintf('<ps> = %f', x), mean(success_prob), 'UniformOutput', false), 'Location', 'northwest');

ax = gca; ax.FontSize= 20; ax.FontName ='Arial';
ylim([0 1])


%% Plot all schedules across the optimizer runs

figure(2)
baseColor = [0, 114, 189]./255;
lowOpacityColor = [baseColor, 0.05]; % RGBA format for transparency (40%)
x_axis  = 1:5;

plot(x_axis, schedule_5_best_n(1,:), '-o','LineWidth', 2, 'Color', 'r',...
    'MarkerFaceColor','r', 'MarkerSize', 8); % Red for individual points
hold on
for run_idx = 1:100
    scatter(x_axis, schedule_5(run_idx, :), 100, ...
[0, 114, 189]./255, 'filled', 'MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.6);

end
hold on
for run_idx=1:100
    plot(x_axis, schedule_5(run_idx, :), '-','LineWidth', 2, 'Color',...
        lowOpacityColor, 'MarkerSize', 8); % Low opacity (0.2)
end

plot(x_axis, schedule_5_best_n(1,:), '-o','LineWidth', 2, 'Color', 'r',...
    'MarkerFaceColor','r', 'MarkerSize', 8); % Red for individual points
    
hold on

% Plot the starting beta schedule
plot(x_axis, 2*ones(5, 1), '-o',...
    'LineWidth', 2, 'Color', [0.96 0.45 0.03],'MarkerFaceColor',[0.96 0.45 0.03],...
    'MarkerSize', 8, 'DisplayName', 'initial \beta (unoptimized)'); % Black for starting values

xlim([0.9 5.1])
xticks(x_axis)
xlabel('number of layers (p)');
ylabel('\beta (Inverse Temperature)');
title("p=5")
ax = gca; ax.FontSize= 20; ax.FontName ='Arial';
ylim([0 6])
