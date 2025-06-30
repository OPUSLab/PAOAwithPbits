clc; clearvars; close all;
%% Author: Abdelrahman, Abdelrahman--OPUS lab (Date: 21st June)
% purpose: plot the results for figure 5
%%
% plotting schedules
figure(1)
final_beta = 0.011;
delta_schedule = 0.005;  % difference in two schedules
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
plot(1:17,up_schedule, '--','LineWidth', 2, 'MarkerSize', 10, 'Color',color1);
hold on
scatter(1:17,up_schedule, 100, ...
    basecolor1, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)


xlabel('Number of layers, p');
ylabel('Inverse temperature, \beta');
legend({'\beta_2=\beta^{opt}', '\beta_1 = \beta^{opt} (1+\Delta)'},'box','off', Location='northwest');
ax = gca; ax.FontSize=20; ax.FontName='arial';
xlim([0.01,17.5])
xticks([1:4:17])
xticklabels([1:4:17]);
hold off;
%%
figure(2)

% plot avg E/N for \beta_1 schedule
% for ii = steps
%     load('Levy_1p_SA_bf_0.12_other.mat');
%     kkAll2(ii) = mean(avg_energy_per_spin_avg_schedule_2p_500{ii});
% end

% plot avg E/N for \beta_2 schedule
for ii = 1:length(steps)
    final_beta =0.011;
    fileName = ['Levy_1p_SA_bf_' num2str(final_beta) '_500.mat'];
    load(fileName);

    kkAll(ii) = mean(avg_energy_per_spin_avg_schedule_2p_500{ii});
end

% plot avg E/N for PAOA schedule
for ii = 1:length(steps)
    final_beta = 0.011;
    delta = 0.005;
    fileName = ['Levy_2p_SA_bf_' num2str(final_beta) '_asym_delta' num2str(delta) '_20p_500_spins.mat'];
    load(fileName);

    kkAll3(ii) = mean(avg_energy_per_spin_avg_schedule_2p_500{ii});
end

for ii = 1:length(steps)
    load('Levy_2p_SA_bf_0.012_asym_delta0.8_20p.mat')
    kkAll4(ii) = mean(avg_energy_per_spin_avg_schedule_2p_500{ii});
end


% for ii = 1:length(steps)
%     final_beta = 0.012;
%     delta = 0.8;
%     fileName = ['Levy_2p_SA_bf_' num2str(final_beta) '_asym_delta' num2str(delta) '_20p.mat'];
%     load(fileName);
% 
%     kkAll4(ii) = mean(avg_energy_per_spin_avg_schedule_2p_50{ii});
% end
% 



% plot(steps(8:17),kkAll2(8:17), '--','LineWidth', 2, 'MarkerSize', 10, 'Color',color1);
% hold on
% scatter(steps(8:17),kkAll2(8:17), 100, ...
%     basecolor1, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)


plot(steps,kkAll, '--','LineWidth', 2, 'MarkerSize', 10, 'Color',color_heavy);
hold on
scatter(steps,kkAll, 100, ...
    basecolor_heavy,'d', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)



plot(steps,kkAll3, '-','LineWidth', 2, 'MarkerSize', 10, 'Color',colro3);

scatter(steps,kkAll3, 100, ...
    base_color3,'s', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)


plot(steps,kkAll4, 'b-','LineWidth', 2, 'MarkerSize', 10);

scatter(steps,kkAll4, 100,'bs', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)


ylim([-5344 -5320])
xlabel('Number of layers, p');
ylabel('Average energy per spin, \langleE/N\rangle');
set(gca,'XScale','log');
set(gca,'YScale','log');
ax = gca; ax.FontSize=20;
xlim([7.8,17.5])
xticks([8:1:17])
xticklabels([8:1:17]);

hold off;