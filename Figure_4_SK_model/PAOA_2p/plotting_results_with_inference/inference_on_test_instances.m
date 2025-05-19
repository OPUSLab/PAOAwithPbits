clear all
clc
%%

% Define the folder path
folder_path = 'replace_with_your_path/Figure_4/PAOA_2p/plotting_results_with_inference/optimized_beta';
 % Get a list of all .mat files in the folder
file_list = dir(fullfile(folder_path, '*.txt'));

% Initialize a matrix to store all the data
num_files = length(file_list);
data_matrix = []; % Each row will correspond to one file
num_rows = 17;
result = cell(length(file_list), 1);
allInstances = cell(num_files,1);
% Loop through each file and load the data
for i = 1:length(file_list)
    % Load the .mat file
    file_path = fullfile(folder_path, file_list(i).name);
        
    fileID = fopen(file_path, 'r');

    % Read all lines of the file
    data = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    lines = data{1};
    

    instance_number = regexp(file_path, '\d?\d*', 'match'); 
    instance_number = str2double(instance_number); % Convert to numeric values
    instance_number = instance_number(end);


    layers = cell(17, 1);  % to store the 12 schedule layers
    
    for p = 1:17
        vector = sscanf(lines{p}, '%f');
        
        % Reshape into a p x 2 matrix (columns: schedule 1 and 2)
        layers{p} = reshape(vector, [p, 2]);
    end
    
    allInstances{instance_number} = layers;  % store the 12-layer cell in the instance cell
end
%% Getting the average schedule (across instances)

avg_schedule = cell(17,1);
for p=1:17
    results = 0;
    for instance=1:num_files
        
        all_results = allInstances{instance};
        results = results + all_results{p};
        
        
    end
    avg_schedule{p} = results./num_files;

end

%% inference on the test instances

folder_path = "/Users/abdulrhmanseed/Downloads/Research/p-bits/github_codes/Figure_4/PAOA_2p/plotting_results_with_inference/test_instances";
file_list = dir(fullfile(folder_path, '*.mat'));



NE=1e6; % number of experiments
num_layers_v = [1:17];

for num_layers_idx = 1:length(num_layers_v)

    selected_2_beta_schedules = avg_schedule{num_layers_idx}; % 
    % ps_generalized = zeros(length(file_list),1);
    num_layers  = num_layers_v(num_layers_idx);

    % test the optimized params for instance 1 on all instances!
    for j=1:length(file_list)
        file_path = fullfile(folder_path, file_list(j).name);
            
        instance_number = regexp(file_path, '\d?\d*', 'match'); 
        instance_number = str2double(instance_number); % Convert to numeric values
        
        J = load(file_path); 
        J = J.J;
        
        if instance_number<=30
            J = J*sqrt(length(J));
        end

       energy_per_spin = MCMC_2p(NE,...
                                selected_2_beta_schedules,J);
       avg_energy_per_spin_avg_schedule(j) = mean(energy_per_spin)

    end  

    avg_energy_per_spin_avg_schedule_2p_26{num_layers_idx} = avg_energy_per_spin_avg_schedule;

end

save('avg_energy_per_spin_avg_schedule_2p_26_test_instances.mat', 'avg_energy_per_spin_avg_schedule_2p_26')
%%
load("avg_energy_per_spin_avg_schedule_2p_26_test_instances.mat")


for i=1:17

figure(4)
plot(1./i, avg_energy_per_spin_avg_schedule_2p_26{i}, 'b.', 'MarkerSize', 15)
hold on
plot(1./i,mean(avg_energy_per_spin_avg_schedule_2p_26{i}), 'rx', 'MarkerSize', 20)
hold on
% plot(1./[1:20],mean(E_average_per_spin).*ones(20,1), 'g-', 'LineWidth', 2)
plot(1./[0.5:20],-0.6723.*ones(20,1), 'g-', 'LineWidth', 2)
% hold on
% plot(1./[0.5:20],-0.763166.*ones(20,1), 'r-', 'LineWidth', 2)
% hold on 
% plot(1./[0.5:20],-2/pi.*ones(20,1), 'k-', 'LineWidth', 2)

xlabel("1/p")
title("<\beta_1(p)>, <\beta_2(p)>")
ylabel("<E>/n")
ax= gca; ax.FontSize=18; ax.FontName='Arial';
ylim([-0.8 -0.2])
xlim([1/18 1.05])

end
xticks(flip([1./[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17]]));
xticklabels({"$\frac{1}{17}$","","$\frac{1}{15}$","","","$\frac{1}{12}$","","","$\frac{1}{9}$","", "","$\frac{1}{6}$","","","$\frac{1}{3}$", "", "$1$"});
% Make sure the interpreter is set to LaTeX so the fractions render properly:
set(gca, 'TickLabelInterpreter','latex')

xlabel('$1 /p$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\langle E \rangle/n$ at $n=26$','Interpreter','latex', 'FontSize', 20);
xscale("log")
yscale("log")


