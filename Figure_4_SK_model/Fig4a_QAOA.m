%% ------------------------------------------------------------------------
% Script Name:  QAOA Energy per Spin vs Inverse Layer Depth (1/p)
% Author:       Abdelrahman, Abdelrhaman – OPUS Lab
% Date:         October 10, 2025
%
% Description:
% This script analyzes Quantum Approximate Optimization Algorithm (QAOA) 
% performance by plotting the average energy per spin ⟨C/n⟩ as a function 
% of the inverse number of layers (1/p) for n = 26 spins.
%
% Workflow:
%   1. Load ground-state energies from precomputed data.
%   2. Parse output files named 'out_p#.txt' (e.g., out_p13.txt), 
%      extract energies for each run and depth p.
%   3. Compute and store energies across up to 30 instances for p = 1–17.
%   4. Plot:
%        • Blue dots: individual instance energies (⟨C/n⟩)
%        • Red crosses: mean energy per depth
%        • Green line: ground-state reference
%
% The figure shows how QAOA energy converges toward the ground-state energy 
% as the circuit depth p increases (i.e., as 1/p → 0).
% ------------------------------------------------------------------------

clear all; clc;

%% ------------------------------------------------------------------------
% 1. Load data
% ------------------------------------------------------------------------
load("ground_energy_n_26.mat");
ground_energies = E_average_per_spin;
n = 26;  % system size

%% ------------------------------------------------------------------------
% 2. Parse QAOA output files
% ------------------------------------------------------------------------
folder_path = "QAOA_results";
files = dir(fullfile(folder_path, 'out_p*.txt'));

num_p = numel(files);
max_p = 17;
num_instances = 30;
data_all = NaN(max_p, num_instances);  % NaN ensures clean averaging

for k = 1:num_p
    % Extract p value (e.g., out_p13.txt → p = 13)
    name = files(k).name;
    p_match = regexp(name, 'out_p(\d+)\.txt', 'tokens');
    p = str2double(p_match{1});

    % Read file contents
    file_path = fullfile(files(k).folder, name);
    fid = fopen(file_path, 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);

    % Extract energy per run
    energy_vals = zeros(1, num_instances);
    for i = 1:num_instances
        line = lines{1}{i};
        match = regexp(line, 'Expected Energy:\s*(-?\d+\.\d+)', 'tokens');
        run_number = regexp(line, 'Run = (\d+)\,', 'tokens');
        energy_vals(str2double(run_number{1}{1}) + 1) = str2double(match{1});
    end
    data_all(p, :) = energy_vals;
end

%% ------------------------------------------------------------------------
% 3. Plot results
% ------------------------------------------------------------------------
figure(2); clf;
hold on; box on;

for p = 1:num_p
    xvals = (1/p) * ones(1, num_instances);
    yvals = data_all(p, :) / n;

    % Individual runs (blue)
    plot(xvals, yvals, '.', 'MarkerSize', 12, 'Color', 'blue');

    % Mean across runs (red cross)
    plot(1/p, mean(yvals), 'x', 'MarkerSize', 15, ...
         'Color', 'r', 'LineWidth', 2);
end

% Ground-state reference (green line)
plot(1./(0.5:20), mean(ground_energies) * ones(1, 20), ...
     'g-', 'LineWidth', 2);

%% ------------------------------------------------------------------------
% 4. Axes, labels, and formatting
% ------------------------------------------------------------------------
title("QAOA");
xlabel('$1/p$', 'Interpreter', 'latex');
ylabel('$\langle C/n \rangle \ at\ n = 26$', 'Interpreter', 'latex');

set(gca, 'XScale', 'log');
xlim([1/18, 1.05]);
ylim([-0.8, -0.25]);

% Custom tick labels
xticks(flip(1 ./ [1:17]));
xticklabels(...
    {"$\frac{1}{17}$","","$\frac{1}{15}$","","", ...
     "$\frac{1}{12}$","","", "$\frac{1}{9}$","", ...
     "","$\frac{1}{6}$","","","$\frac{1}{3}$", ...
     "", "$1$"});
set(gca, 'TickLabelInterpreter', 'latex');

% Appearance
ax = gca;
ax.FontSize = 20;
ax.FontName = 'arial';
grid on;
box on;
