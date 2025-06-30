clear all
clc
%%
steps = 1:17;
basecolor1 = [0 0.4470 0.7410];
color1 = [basecolor1, 0.8];

basecolor_heavy = [0.8500 0.3250 0.0980];
color_heavy = [basecolor_heavy, 0.8];

base_color3 = [0.4660, 0.6740, 0.1880];
colro3 = [0.4660, 0.6740, 0.1880, 0.8];



% plot avg E/N for \beta_2 schedule
for ii = steps
    final_beta=0.12;
    fileName = ['Levy_1p_SA_bf_' num2str(final_beta) '.mat'];
    load(fileName);
    kkAll_50(ii) = mean(avg_energy_per_spin_avg_schedule_2p_500{ii});
end

% plot avg E/N for PAOA schedule
for ii = 1:10
    final_beta = 0.12;
    delta = 0.07;
    fileName = ['Levy_2p_SA_bf_' num2str(final_beta) '_asym_delta' num2str(delta) '_20p_50_spins.mat'];
    load(fileName);

    kkAll3_50(ii) = mean(avg_energy_per_spin_avg_schedule_2p_50{ii});
end

for ii = steps
    load('Levy_2p_SA_bf_0.12_asym_delta0.7_20p.mat')
    kkAll4_50(ii) = mean(avg_energy_per_spin_avg_schedule_2p_500{ii});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% n=500 %%%%%%%%%%%%%%%%%%%%%%%%
steps = 8:17;
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


subplot(2, 1, 1)
plot(steps,kkAll_50(8:17), '--','LineWidth', 2, 'MarkerSize', 10, 'Color',color_heavy);
hold on
scatter(steps,kkAll_50(8:17), 100, ...
    basecolor_heavy,'d', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)

plot(steps,kkAll3_50, '-','LineWidth', 2, 'MarkerSize', 10, 'Color',colro3);

scatter(steps,kkAll3_50, 100, ...
    base_color3,'s', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)

plot(steps,kkAll4_50(8:17), 'b-','LineWidth', 2, 'MarkerSize', 10);

scatter(steps,kkAll4_50(8:17), 100,'bs', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)
xlabel('Number of layers, p');
set(gca,'XScale','log');
set(gca,'YScale','log');
ax = gca; ax.FontSize=20;
xlim([7.8,17.5])
xticks([8:1:17])
xticklabels([8:1:17]);
ylim([-839.5 -837.5])


subplot(2, 1, 2)

plot(steps,kkAll, '--','LineWidth', 2, 'MarkerSize', 10, 'Color',color_heavy);
hold on
scatter(steps,kkAll, 100, ...
    basecolor_heavy,'d', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)

plot(steps,kkAll3, '-','LineWidth', 2, 'MarkerSize', 10, 'Color',colro3);

scatter(steps,kkAll3, 100, ...
    base_color3,'s', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)

plot(steps,kkAll4, 'b-','LineWidth', 2, 'MarkerSize', 10);

scatter(steps,kkAll4, 100,'bs', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.5)



xlabel('Number of layers, p');
ylabel('\langleE/N\rangle');
set(gca,'XScale','log');
set(gca,'YScale','log');
ax = gca; ax.FontSize=20;
xlim([7.8,17.5])
xticks([8:1:17])
xticklabels([8:1:17]);
ylim([-5344 -5320])

hold off;