clear all
clc
%%
% Parameters
n = 26; % system size


% Collect files matching pattern
files = dir('out_p*.txt');
num_p = numel(files);

% Initialize container assuming maximum p = 17
max_p = 17;
num_instances = 30;
data_all = NaN(max_p, num_instances);  % Use NaN to avoid uninitialized rows

% Read each file
files = dir('out_p*.txt');
for k = 1:length(files)
    % Extract p value from filename, e.g., out_p13.txt → 13
    name = files(k).name;
    p_match = regexp(name, 'out_p(\d+)\.txt', 'tokens');
    p = str2double(p_match{1});  % Extracted p value

    % Read the file lines
    file_path = fullfile(files(k).folder, name);
    fid = fopen(file_path, 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);

    % Extract energy values
    energy_vals = zeros(1, num_instances);
    for i = 1:num_instances
        line = lines{1}{i};
        match = regexp(line, 'Expected Energy:\s*(-?\d+\.\d+)', 'tokens');
        energy_vals(i) = str2double(match{1});
    end

    % Store the results at the correct row
    data_all(p, :) = energy_vals;
end


% Plotting
figure(2); clf;
hold on; box on;

for p = 1:num_p
    xvals = (1/p) * ones(1, num_instances);
    yvals = data_all(p, :) / n;

    % Plot individual runs
    plot(xvals, yvals, '.', 'MarkerSize', 15, 'Color','blue');
    % Plot mean
    plot(1/p, mean(yvals), 'x', 'MarkerSize', 20, 'Color','r');
end

% Reference energy line
plot(1./(0.5:20), -0.6723 * ones(1, 20), 'g-', 'LineWidth', 2);

% Labeling
title("QAOA (star) vs PAOA (cross)");
xlabel('$1/p$', 'Interpreter','latex');
ylabel('$\langle C/n \rangle \ at\ n = 26$', 'Interpreter','latex');

set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([1/18, 1.03]);
ylim([-0.8, -0.2]);

% Custom ticks
xticks(flip(1 ./ [1:17]));
xticklabels(...
    {"$\frac{1}{17}$","","$\frac{1}{15}$","","", ...
     "$\frac{1}{12}$","","", "$\frac{1}{9}$","", ...
     "","$\frac{1}{6}$","","","$\frac{1}{3}$", ...
     "", "$1$"});
set(gca, 'TickLabelInterpreter', 'latex');

ax = gca;
ax.FontSize = 20;
ax.FontName = 'arial';
