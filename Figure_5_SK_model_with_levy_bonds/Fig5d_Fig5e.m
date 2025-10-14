%% ------------------------------------------------------------------------
% Script Name:  Bootstrap Analysis and Comparison of Single vs Double SA
% Author:       Abdelrahman, Abdelrahman – OPUS Lab
% Date:         October 10, 2025
%
% Description:
% This script reproduces Fig. 5d–5e results comparing single-schedule and
% double-schedule simulated annealing (SA) on SK instances.
%
% Workflow:
%   1. Load success probability matrices for 500 instances of single
%      and double SA runs.
%   2. For each sweep count, compute mean ± 95% bootstrap confidence
%      intervals (CI) over 1e5 bootstrap resamples (Fig. 5d).
%   3. Plot the average success probability vs sweep count for both
%      schedules.
%   4. For a selected sweep value (e.g., 5e5), plot pointwise comparisons
%      of success probabilities for each instance, color-coded by whether
%      single or double schedule performed better (Fig. 5e).
%
%   Figures:
%     • Fig. 5d → mean ± 95% CI success probabilities vs sweeps
%     • Fig. 5e → scatter plot of instancewise comparison at fixed sweeps
% ------------------------------------------------------------------------

clear all; clc;

%% ------------------------------------------------------------------------
% 1. Load data
% ------------------------------------------------------------------------
load("p_success_single_SA_putative_ground_500_instances.mat")
load("p_success_double_SA_putative_ground_500_instances_20_80.mat")

P_single_ground = p_success_single_SA_putative_ground_500_instances;
P_double_ground = p_success_double_SA_putative_ground_500_instances_20_80;


%% ------------------------------------------------------------------------
% 2. Compute bootstrap mean and confidence intervals (95%)
% ------------------------------------------------------------------------
xvals = [1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5];  % sweep counts (x-axis)
alpha = 0.05;   % significance level for 95% CI
B = 1e5;        % number of bootstrap resamples

[m_sg, lo_sg, hi_sg] = bootmeanCI(P_single_ground, B, alpha);
[m_dg, lo_dg, hi_dg] = bootmeanCI(P_double_ground, B, alpha);

% Convert to asymmetric error bar lengths
neg_sg = m_sg - lo_sg; pos_sg = hi_sg - m_sg;
neg_dg = m_dg - lo_dg; pos_dg = hi_dg - m_dg;


%% ------------------------------------------------------------------------
% 3. Plot Fig. 5d — Bootstrap confidence intervals
% ------------------------------------------------------------------------
figure('Color','w'); hold on; grid on;
set(gca, 'XScale', 'log');
errorbar(xvals, m_sg, neg_sg, pos_sg, 'o-', ...
    'LineWidth',1.4, 'MarkerSize',6);
errorbar(xvals, m_dg, neg_dg, pos_dg, 's-', ...
    'LineWidth',1.4, 'MarkerSize',6);

xlim([0.5e2, 6e5]); ylim([0,0.25]);
xticks(xvals); 
xticklabels({'10^2','5×10^2','10^3','5×10^3','10^4','5×10^4','10^5','5×10^5'});
xlabel('Sweeps'); ylabel('p_{success}');
title('Putative Ground: mean ± 95% bootstrap CI');
legend('Single schedule','Two schedules (heavy hotter)','Location','southeast');
ax = gca; ax.FontSize = 20; grid off; box on;


%% ------------------------------------------------------------------------
% 4. Plot Fig. 5e — Instancewise comparison at selected sweeps
% ------------------------------------------------------------------------
selected_sweeps = [5e5];
cols_to_plot = find(ismember(xvals, selected_sweeps));

% ---- Colors ----
c_red   = [0.85 0.10 0.10];       % double-SA better (y>x)
c_blue  = [0.00 0.447 0.741];     % single-SA better (y<x)
c_green = [117, 172, 66]./255;    % tie (y≈x)

msz = 42; mkAlpha = 0.85;
tie_tol = 1e-3;                    % tolerance for tie
use_relative = false;              % absolute vs relative gain

[num_instances, ~] = size(P_single_ground);

for c = cols_to_plot
    x = P_single_ground(:, c);   % single schedule
    y = P_double_ground(:, c);   % double schedule
    
    % Remove invalid pairs
    valid = ~(isnan(x) | isnan(y));
    xv = x(valid);  yv = y(valid);
    if isempty(xv), continue; end

    % Masks
    diff = yv - xv;
    mask_tie   = abs(diff) <= tie_tol;
    mask_above = diff >  tie_tol;   % double-SA better
    mask_below = diff < -tie_tol;   % single-SA better

    % ---- Figure ----
    figure('Color','w'); hold on; grid off;
    plot([0, 1], [0, 1], 'k--', 'LineWidth', 1.3);

    % Scatter categories
    hBlue  = scatter(xv(mask_below), yv(mask_below), msz, c_blue,  'filled', 'MarkerFaceAlpha', mkAlpha);
    hRed   = scatter(xv(mask_above), yv(mask_above), msz, c_red,   'filled', 'MarkerFaceAlpha', mkAlpha);
    hGreen = scatter(xv(mask_tie),   yv(mask_tie),   msz, c_green, 'filled', 'MarkerFaceAlpha', mkAlpha);

    % Count annotation
    n_double = sum(mask_above);
    n_single = sum(mask_below);
    n_tie    = sum(mask_tie);
    txt = sprintf('Double better: %d\nSingle better: %d\nTie: %d', ...
        n_double, n_single, n_tie);
    text(0.05, 0.80, txt, 'FontSize', 11, 'Color', [0.1 0.1 0.1]);

    xlabel('Single schedule success probability');
    ylabel('Double schedule success probability');
    title(sprintf('Comparison at sweeps = %g', xvals(c)));
    legend([hRed, hBlue, hGreen], ...
        'y>x (double schedule better)', ...
        'y<x (single schedule better)', ...
        'y≈x (tie)', 'Location','southeast');
    axis square; box on;
    xlim([0 1]); ylim([0 1]);
    ax = gca; ax.FontSize = 20;
end


%% ------------------------------------------------------------------------
% 5. Local function — Bootstrap mean & percentile CI
% ------------------------------------------------------------------------
function [m, lo, hi] = bootmeanCI(M, B, alpha)
    % Computes bootstrap mean and confidence intervals column-wise.
    if nargin < 3, alpha = 0.05; end
    if nargin < 2, B = 5000; end

    [~, K] = size(M);
    m  = nan(1,K); lo = nan(1,K); hi = nan(1,K);
    hasBootci = exist('bootci','file') == 2;

    for k = 1:K
        x = M(:,k); x = x(~isnan(x));
        n = numel(x);
        if n == 0, continue; end

        m(k) = mean(x);
        if n == 1, lo(k) = x; hi(k) = x; continue; end

        if hasBootci
            bootfun = @(z) mean(z,'omitnan');
            ci = bootci(B, {bootfun, x}, 'alpha', alpha);
            lo(k) = ci(1); hi(k) = ci(2);
        else
            bs = zeros(B,1);
            for b = 1:B
                idx = randi(n, n, 1);
                bs(b) = mean(x(idx));
            end
            srt = sort(bs);
            lo(k) = srt(max(1, floor(B*(alpha/2))));
            hi(k) = srt(min(B, ceil(B*(1 - alpha/2))));
        end
    end
end
