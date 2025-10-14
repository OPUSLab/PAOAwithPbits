%% ------------------------------------------------------------------------
% Script Name:  Comparison of Approximation Ratios — PAOA vs QAOA
% Author:       Abdelrahman, Abdelrhaman – OPUS Lab
% Date:         October 10, 2025
%
% Description:
% This script compares the approximation ratio (AR) performance between 
% the Probabilistic Approximate Optimization Algorithm (PAOA) and the 
% Quantum Approximate Optimization Algorithm (QAOA) for n = 26 spins.
%
% Workflow:
%   1. Load:
%        • PAOA energy data per schedule (avg_energy_per_spin_avg_schedule_2p_26_test_instances.mat)
%        • Ground-state energies (ground_energy_n_26.mat)
%   2. Parse QAOA output files: out_p#.txt (e.g., out_p13.txt)
%        • Extract expected energies from 30 runs per depth p.
%   3. Compute the Approximation Ratio (AR):
%        AR = ⟨E⟩ / E_ground
%        for both PAOA and QAOA at each p = 1…17.
%   4. Bootstrap 95% Confidence Intervals (CIs) for each mean AR.
%   5. Plot mean ± CI error bars:
%        • Blue circles — QAOA
%        • Red squares — PAOA
%
% The resulting plot shows how PAOA compares to QAOA across increasing 
% layer depths (p) in terms of approximation ratio quality.
% ------------------------------------------------------------------------

clear all; clc;

%% ------------------------------------------------------------------------
% 1. Load data
% ------------------------------------------------------------------------
load("avg_energy_per_spin_avg_schedule_2p_26_test_instances.mat");
load("ground_energy_n_26.mat");
ground_energies = E_average_per_spin;

n = 26;              % system size
num_instances = 30;  % number of instances

%% ------------------------------------------------------------------------
% 2. Parse QAOA output files
% ------------------------------------------------------------------------
folder_path = "QAOA_results";
files = dir(fullfile(folder_path, 'out_p*.txt'));
num_p = numel(files);
max_p = 17;

data_all = NaN(max_p, num_instances);  % preallocate with NaNs

for k = 1:num_p
    % Extract layer depth (p) from filename
    name = files(k).name;
    p_match = regexp(name, 'out_p(\d+)\.txt', 'tokens');
    p = str2double(p_match{1});

    % Read file contents
    file_path = fullfile(files(k).folder, name);
    fid = fopen(file_path, 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);

    % Extract energy values per run
    energy_vals = zeros(1, num_instances);
    for i = 1:num_instances
        line = lines{1}{i};
        match = regexp(line, 'Expected Energy:\s*(-?\d+\.\d+)', 'tokens');
        run_number = regexp(line, 'Run = (\d+)\,', 'tokens');
        energy_vals(str2double(run_number{1}{1}) + 1) = str2double(match{1});
    end

    % Store at correct p-row
    data_all(p, :) = energy_vals;
end

%% ------------------------------------------------------------------------
% 3. Compute Approximation Ratios and 95% CIs
% ------------------------------------------------------------------------
figure('Color', 'w'); hold on; box on;
alpha = 0.05;  % 95% CI

for p = 1:num_p
    yvals = data_all(p, :) / n;  % normalize QAOA energies
    QAOA_AR = yvals' ./ E_average_per_spin;
    PAOA_AR = avg_energy_per_spin_avg_schedule_2p_26{p}' ./ ground_energies;

    % Bootstrap confidence intervals
    [qCI_lo, qCI_hi] = bootci(1e4, {@mean, QAOA_AR}, 'alpha', alpha, 'type', 'per');
    [pCI_lo, pCI_hi] = bootci(1e4, {@mean, PAOA_AR}, 'alpha', alpha, 'type', 'per');

    % Means
    q_mean = mean(QAOA_AR);
    p_mean = mean(PAOA_AR);

    % Asymmetric error bars
    q_err = [q_mean - qCI_lo; qCI_hi - q_mean];
    p_err = [p_mean - pCI_lo; pCI_hi - p_mean];

    % Plot mean ± CI
    
    errorbar(p, p_mean, p_err(1), p_err(2), 'rs', ...
        'MarkerSize', 10, 'CapSize', 15, 'LineWidth', 1.5);
    errorbar(p, q_mean, q_err(1), q_err(2), 'bo', ...
        'MarkerSize', 10, 'CapSize', 15, 'LineWidth', 1.5);
end

%% ------------------------------------------------------------------------
% 4. Figure Formatting
% ------------------------------------------------------------------------
title("PAOA (red) vs QAOA (blue): 2p", 'FontSize', 20);
xlabel('p', 'FontSize', 22);
ylabel('Approximation Ratio at n=26', 'FontSize', 22);
xlim([0, 17.5]);
ylim([0.4, 1]);

xticks(1:17);
ax = gca;
ax.FontSize = 25;
ax.FontName = 'Arial';
grid on;
box on;
legend({'PAOA', 'QAOA'}, 'Location', 'southeast', 'Box', 'off');
