%% ------------------------------------------------------------------------
% Script Name:  Energy per Spin vs Inverse Layer Depth (1/p)
% Author:       Abdelrahman, Abdelrhaman – OPUS Lab
% Date:         October 10, 2025
%
% Description:
% This script visualizes the relationship between the inverse number of 
% layers (1/p) and the average energy per spin ⟨E⟩/n for the SK model 
% with n = 26 spins. It compares the energies obtained from the averaged 
% 2-schedule PAOA runs against the theoretical ground-state energy.
%
% Workflow:
%   1. Load the simulated energies:
%        • avg_energy_per_spin_avg_schedule_2p_26_test_instances.mat
%        • ground_energy_n_26.mat
%   2. For each layer depth p = 1…17:
%        • Plot individual instance energies (blue points).
%        • Plot their mean (red ‘×’ markers).
%        • Overlay the mean ground-state energy (green line).
%   3. Format the x-axis in terms of 1/p with LaTeX-style tick labels.
%   4. Display ⟨E⟩/n as a function of inverse layer index on a log x-scale.
%
% The plot highlights how energy per spin converges toward the ground-state 
% value as the number of layers increases.
% ------------------------------------------------------------------------

clear all; clc;

%% ------------------------------------------------------------------------
% 1. Load data
% ------------------------------------------------------------------------
load("avg_energy_per_spin_avg_schedule_2p_26_test_instances.mat");
load("ground_energy_n_26.mat");
ground_energies = E_average_per_spin;

%% ------------------------------------------------------------------------
% 2. Plot ⟨E⟩/n vs 1/p for all 17 layers
% ------------------------------------------------------------------------
figure(4); hold on;

for i = 1:17
    % Scatter all instance energies at this depth
    plot(1./i, avg_energy_per_spin_avg_schedule_2p_26{i}, '.', ...
        'MarkerSize', 12, 'Color', 'blue');

    % Mean energy at this depth
    plot(1./i, mean(avg_energy_per_spin_avg_schedule_2p_26{i}), 'x', ...
        'MarkerSize', 15, 'Color', 'r', 'LineWidth', 2);
end

% Ground-state reference (horizontal line)
plot(1./[0.5:20], mean(ground_energies).*ones(20,1), 'g-', 'LineWidth', 2);

xlabel("1/p");
ylabel("<E>/n");
title("<\beta_1(p)>, <\beta_2(p)>");
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ylim([-0.8, -0.25]);
xlim([1/18, 1.05]);

%% ------------------------------------------------------------------------
% 3. Format axis labels and appearance
% -------------------------------------------------------------------------
xticks(flip(1./[1:17]));
xticklabels({"$\frac{1}{17}$","","$\frac{1}{15}$","","","$\frac{1}{12}$",...
              "","","$\frac{1}{9}$","","","$\frac{1}{6}$","","",...
              "$\frac{1}{3}$","","$1$"});
set(gca, 'TickLabelInterpreter', 'latex');

xlabel('$1/p$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\langle E \rangle/n$ at $n=26$', 'Interpreter', 'latex', 'FontSize', 20);

set(gca, 'XScale', 'log');   % Log scale for x-axis
ax.FontName = 'arial';
grid on;
box on;
