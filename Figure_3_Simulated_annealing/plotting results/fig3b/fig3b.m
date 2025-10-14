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
%% Different way of visulaization
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

figure;  hold on

% (1) stem plot for the optimised schedule
h1 = stem(binEdges(1:end-1), prob1, 'filled');          % vertical lines + markers
set(h1, 'LineStyle', '-',  ...                   % solid stems
        'Marker',     'o', ...                   % square markers
        'LineWidth',  1.6, ...
        'MarkerSize', 8, ...
        'Color',      [0 0.4470 0.7410]);        % MATLAB blue

% (2) stem plot for the initial schedule
h2 = stem(binEdges(1:end-1), prob2, 'filled');
set(h2, 'LineStyle', '--', ...                   % dashed stems
        'Marker',     'o', ...                   % diamond markers
        'LineWidth',  1.6, ...
        'MarkerSize', 8, ...
        'Color',      [0.8500 0.3250 0.0980]);   % MATLAB orange

% ----- cosmetic tweaks ------------------------------------------------
% Make sure stems start at a positive base (needed for log‑scale y‑axis)
set([h1 h2], 'BaseValue', 1e-3);   % any small positive number works
% If either prob1 or prob2 contains exact zeros,
% replace them with NaN (so log‑scale doesn’t choke):
prob1(prob1==0) = NaN;
prob2(prob2==0) = NaN;

% Reference vertical line at E = –360
% plot([-360 -360], ylim, 'k-', 'LineWidth', 10);
plot(-360*ones(40,1), linspace(1e-2, 1, 40), 'g--', LineWidth=2.5)

xlabel('3D-Spin Glass System Energy');
ylabel('Occupation Probability');
title(sprintf('p = %d', i*5));
ax = gca;  ax.FontSize = 20;  ax.FontName = 'arial';
xlim([-362 -340]);
ylim([1e-2 1]);
yscale(ax, 'log');

legend({'Optimised $\beta$', 'Initial $\beta$'}, ...
       'Interpreter','latex', 'Location','northwest');
box on;  
ax.XMinorTick = 'on';  % activate minor ticks on x-axis

end




%% Majorization & CDF dominance across p = 5,10,15  (PAOA / optimized β)

p_vals   = [5 10 15];
idxs     = 1:3;                 % i=1->p=5, i=2->p=10, i=3->p=15
nbins    = 8;                  % finer bins → smoother CDF (adjust as you like)

% ----- Common binning over all energies (for fair comparison) -----
allE = reshape(cost_all_720(idxs,:).', [], 1);
Emin = min(allE);  Emax = max(allE);
edges = linspace(Emin, Emax, nbins+1);
centers = edges(1:end-1) + diff(edges)/2;

% ----- Build PMFs and CDFs for each p -----
PMF = zeros(nbins, numel(idxs));
CDF = zeros(nbins, numel(idxs));
for k = 1:numel(idxs)
    Ei       = cost_all_720(idxs(k),:).';
    counts   = histcounts(Ei, edges);
    PMF(:,k) = counts / max(1, sum(counts));
    CDF(:,k) = cumsum(PMF(:,k));
end

% ===== 1) Plot CDFs vs Energy & check stochastic dominance =====
figure('Name','CDF vs Energy'); hold on;
plot(centers, CDF(:,1), '-',  'LineWidth',2);
plot(centers, CDF(:,2), '--', 'LineWidth',2);
plot(centers, CDF(:,3), '-.', 'LineWidth',2);
yl = ylim;  plot([-360 -360], yl, 'k:', 'LineWidth',1.2); ylim(yl);  % ref line (optional)
xlabel('3D-Spin Glass Energy'); ylabel('CDF  P(E \le e)');
title('Energy CDFs  (expect F_5 \le F_{10} \le F_{15})');
legend(compose('p = %d', p_vals), 'Location','southeast');
set(gca,'FontSize',18,'FontName','arial','Box','on');


%% ===== 2) Majorization on PDFs (top-m partial sums) =====
p5s  = sort(PMF(:,1),  'descend');
p10s = sort(PMF(:,2),  'descend');
p15s = sort(PMF(:,3),  'descend');

S5   = cumsum(p5s);
S10  = cumsum(p10s);
S15  = cumsum(p15s);

tol = 1e-12;
maj_15_vs_10_ok = all(S15(1:end-1) - S10(1:end-1) >= -tol);
maj_10_vs_5_ok  = all(S10(1:end-1) - S5(1:end-1)  >= -tol);

fprintf('\nMajorization (descending-PMF partial sums):\n');
fprintf('  P15 majorizes P10? %s\n', yesno(maj_15_vs_10_ok));
fprintf('  P10 majorizes P5?  %s\n', yesno(maj_10_vs_5_ok));

figure('Name','Majorization (Top-m)'); hold on;
plot(1:nbins, S5,  '-',  'LineWidth',2);
plot(1:nbins, S10, '--', 'LineWidth',2);
plot(1:nbins, S15, '-.', 'LineWidth',2);
xlabel('m (from most probable to least probable event)');
ylabel('\Sigma_{j=1}^{m} p_j, m \in \{1, 2, ..., K\} and p_1> p_2> ...>p_K');
title('majorization (higher curve -> more concentrated)');
legend('P_5','P_{10}','P_{15}','Location','southeast');
set(gca,'FontSize',18,'FontName','arial','Box','on');

% ----- small helper -----
function s = yesno(flag), if flag, s='YES'; else, s='NO'; end, end
