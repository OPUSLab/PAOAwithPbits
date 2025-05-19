clc
clearvars
close all
%% Author: Abdelrahman, Abdelrahman-OPUS Lab (Date: May 18, 2025)
% Purpose: Finding optimizd beta schedules for 3D spin Glass problems
% (hardware implementation of Fig.3a)
%% Preliminary setting for clk frequencies
ref_clk = 125e6; 
pbit_clk = 15e6;
bram_read_clk = 100e6;
BRAM_depth = 1e6;

%% ************************************************ FPGA interface and addresses/ Memory mappings *************************************
mem = aximanager('Xilinx','interface','pcie');

% optimized C and C++ scripts for bit conversion and cost calculations
mex bit_conversion_mex.c
mex -v COMPFLAGS="$COMPFLAGS /Qopenmp /03" -I"C:path_to_eigen\eigen" compute_avg_energy_optimized.cpp
mex -v COMPFLAGS="$COMPFLAGS /Qopenmp /03" -I"C:path_to_eigen\eigen" compute_cost_vector.cpp


weight_address = hex2dec('003D0920'); % J address
h_address = hex2dec('0040B2A0');
tictoc_counter_limit_address = hex2dec('00000004');

% trigger registers
J_trigger_address = hex2dec('004229A0');
h_trigger_address = hex2dec('004229A4');


% beta_counter registers
num_clk_limit_address = hex2dec('004229B4'); % num of clk cycles




% snapshot registers
take_snapshot_counter_address= hex2dec('004229C0');
write_snapshot_counter_address= hex2dec('004229C4');
num_sweep_address= hex2dec('004229C8');

%% convert from bipolar to binary and reshape based on target adjacency
% J weight using L =3 3D-spin glass problem realization
load L6.mat
J = full(J);


base  =J;
num_replicas = 10;
for i=1:num_replicas-1

J = blkdiag(J,base);

end


h = 0*ones(length(J),1);
clc
save J_10.mat J
save h_10.mat h

load J_10.mat
load h_10.mat
num_pbits = length(J);
num_replicas=1;
% convert from bipolar to binary and reshape based on target adjacency
J_bipolar = J;
h_bipolar = h;
max_num_neighbors = max(degree(graph(J))); % max neighbors in target
target_graph = J;

[~,J_indices] = get_triu_adjacency(target_graph, max_num_neighbors);
[J_binary_reshaped, h_binary] = triu_bip2bin_weight_converter(J_bipolar,...
                                                              h_bipolar,...
                                                              J_indices,...
                                                              max_num_neighbors);



% convert the weights
fixed_a = 4; fixed_b =5; % fixed points as ([s]a.b), a for int, b for fraction
% laod weights
[J_fpga,h_fpga] = fixed_point_weights_j_h(fixed_a, fixed_b, J_binary_reshaped, h_binary);
h_fpga_with_flag = [h_fpga; fi(0, 1, fixed_a+fixed_b+1, fixed_b)];

% load the weights
RESET =1; SHIFT = 1;


% load the weights and beta
writememory(mem,weight_address,J_fpga,'BurstType','Increment'); % get it out of the loop
writememory(mem,h_address,h_fpga_with_flag,'BurstType','Increment');     % get it out of the loo
writememory(mem,J_trigger_address,SHIFT,'BurstType','Increment'); 
writememory(mem,h_trigger_address,SHIFT,'BurstType','Increment'); 


%% ********************************************************* PAOA training ********************************************************* 
% layers/time steps vector
for COBYLA_run=1:100
    layers = 15; % number of layers
    optimized_beta_schedule = cell(length(layers),1);
    putative_ground_state = 360;%207;
    num_words_pbit= ceil(num_pbits*num_replicas/32);
    num_experiments_inside_FPGA = 1e5/10;%floor(BRAM_depth/num_words_pbit); % maximum
    total_number_of_experiments = 1e5;
    num_experiments_outside = total_number_of_experiments/(num_experiments_inside_FPGA*10);
    
    success_prob_PAOA = zeros(length(layers),1);
    
    for i=1:length(layers)
        
        fprintf("The simulation is for %d-layers PAOA.....", layers(i))
        
        % number of betas/layers/time steps
        num_layers = layers(i);
        
        
        % Determine beta living time and loading sim. parameters into FPGA: Time steps/sweeps
        synthesized_beta_length = 61; % the one is for the zero
        num_repetition = (synthesized_beta_length-1)/num_layers; 
        num_sweeps_layer = 720;%number of sweeps/unique beta value 
        num_sweeps_synthezied_beta = num_sweeps_layer/num_repetition; % this has to be integer
        tictoc_counter_limit = num_experiments_inside_FPGA*synthesized_beta_length*num_sweeps_synthezied_beta; % total annealing time
        num_clk_limit = num_sweeps_synthezied_beta; % each synthesized beta living time
        
        % Write and Take counters for snapshot modules
        write_snapshot_counter = ceil(num_words_pbit*pbit_clk/bram_read_clk)+1;
        take_snapshot_counter = 1; % every clk cycle
    
    
    
        % write on FPGA each beta living time and the whole time for all betas 
        % % write on FPGA each beta living time and the whole time for all betas 
        writememory(mem,tictoc_counter_limit_address,uint32(tictoc_counter_limit), 'BurstType','Increment');
        writememory(mem,take_snapshot_counter_address,take_snapshot_counter, 'BurstType','Increment');
        writememory(mem,write_snapshot_counter_address,write_snapshot_counter, 'BurstType','Increment');
        writememory(mem,num_sweep_address,num_experiments_inside_FPGA,'BurstType','Increment'); % experiments inside FPGA
        writememory(mem,num_clk_limit_address,uint32(num_clk_limit), 'BurstType','Increment');
    
        if log2(tictoc_counter_limit) > 32
            error("You exceeded the tictoc counter limit capacity")
        end
    
        % intial guess for cooling schedule
        initial_beta = 2*ones(1,num_layers);%linspace(0.1, 5, num_layers);%ones(1, num_layers);
        
        % Objective function to be minimized
        objective = @(initial_beta) PAOA_circuit_online_annealing_snapshot(num_experiments_outside,...
                                    num_experiments_inside_FPGA, num_replicas, J, initial_beta); %beta is the variable
        options = struct();
        options.maxfun=5000;
        options.radius_init=1;
        options.radius_final=1e-4;
        options.quiet=false;
        
    
        
        %optimization/training 
        [optimized_beta_schedule{i}, fx, exitflag, output] = cobyla(objective, initial_beta, options);
         
        
       
        %inference
        cost_for_all_exp = inference_snapshot(num_experiments_outside,num_experiments_inside_FPGA, 10, J, J(1:216,1:216), optimized_beta_schedule{i});
    
        success_prob_PAOA(i) = length(find(cost_for_all_exp<=-putative_ground_state))/(total_number_of_experiments);
        
        fprintf("The Success Probability for %d layers in percentage is: %.2f\n", layers(i), 100*success_prob_PAOA(i))
        
    
    end
file_name_1 = sprintf('success_prob_L6_1e6_run_%d_15_layers_720', COBYLA_run);
save(file_name_1, "success_prob_PAOA")

file_name = sprintf('optimized_beta_schedule_L6_1e6_run_%d_15_layers_720', COBYLA_run);
save(file_name, "optimized_beta_schedule")

end


%% ***************************************************** Inference PAOA **********************************************************
folderPath = 'L6_1e5_FPGA_optimized_params'; % Replace with the actual path to your files
filePattern = fullfile(folderPath, '*.mat');
fileList = dir(filePattern);

% Initialize the cell array to store beta values for each file
optimized_beta_all = cell(3,100);
% Loop through each file in the directory
for fileIdx = 1:length(fileList)
    % Construct full file path
    filename = fullfile(fileList(fileIdx).folder, fileList(fileIdx).name);
    disp(filename)
    load(filename)    

    % Extract the run number from the filename using regular expressions
    runNumber = regexp(fileList(fileIdx).name, '_run_(\d+)', 'tokens');
    
    if ~isempty(runNumber)
        runNumber = str2double(runNumber{1}{1}); % Convert to a numeric format
        fprintf('Processing run number: %d\n', runNumber);
    else
        warning('Run number not found in filename: %s', fileList(fileIdx).name);
    end

   % Append the beta values of the current file to the main cell array
   for i=1:length(optimized_beta_schedule)
    optimized_beta_all{i,runNumber} = optimized_beta_schedule{i};
   end
    
end



success_prob_PAOA = zeros(length(layers), 1);
costs_for_all_layers = zeros(length(layers), total_number_of_experiments);

% Example 2: Using start and end values
startValue =0.5;     % Initial value
endValue = 5;     % Final value

% Calculate ratio based on start and end values
ratio_calculated = (endValue / startValue)^(1/(layers-1));

% Generate geometric sequence
schedule_end = startValue * ratio_calculated.^(0:layers-1);


for i=1:length(layers)
    optimized_beta_schedule_derived = schedule_end'; % this must be column vector
    %inspace(0.1, 5, 15)';
    %[];
    %optimized_beta_schedule{1};

%optimized_beta_all{i,COBYLA_run};

    fprintf("The simulation is for %d-layers PAOA.....", layers(i))
    
     % number of betas/layers/time steps
    num_layers = layers(i);
    
    
    % Determine beta living time and loading sim. parameters into FPGA: Time steps/sweeps
    synthesized_beta_length = 61; % the one is for the zero
    num_repetition = (synthesized_beta_length-1)/num_layers; 
    num_sweeps_layer = 720;%number of sweeps/unique beta value 
    num_sweeps_synthezied_beta = num_sweeps_layer/num_repetition; % this has to be integer
    tictoc_counter_limit = num_experiments_inside_FPGA*synthesized_beta_length*num_sweeps_synthezied_beta; % total annealing time
    num_clk_limit = num_sweeps_synthezied_beta; % each synthesized beta living time
    
    % Write and Take counters for snapshot modules
    write_snapshot_counter = ceil(num_words_pbit*pbit_clk/bram_read_clk)+1;
    take_snapshot_counter = 1; % every clk cycle



    % write on FPGA each beta living time and the whole time for all betas 
    % % write on FPGA each beta living time and the whole time for all betas 
    writememory(mem,tictoc_counter_limit_address,uint32(tictoc_counter_limit), 'BurstType','Increment');
    writememory(mem,take_snapshot_counter_address,take_snapshot_counter, 'BurstType','Increment');
    writememory(mem,write_snapshot_counter_address,write_snapshot_counter, 'BurstType','Increment');
    writememory(mem,num_sweep_address,num_experiments_inside_FPGA,'BurstType','Increment'); % experiments inside FPGA
    writememory(mem,num_clk_limit_address,uint32(num_clk_limit), 'BurstType','Increment');

    if log2(tictoc_counter_limit) > 32
        error("You exceeded the tictoc counter limit capacity")
    end


    %inference
    cost_for_all_exp = inference_snapshot(num_experiments_outside,num_experiments_inside_FPGA, 10, J, J(1:216,1:216), optimized_beta_schedule_derived);
    costs_for_all_layers(i,:) = cost_for_all_exp;
    success_prob_PAOA(i) = length(find(cost_for_all_exp<=-putative_ground_state))/(total_number_of_experiments);
    

end
disp(success_prob_PAOA)


%% ********************************************************** Flat SA *********************************************************
for run=1:100
    success_prob_SA = zeros(1, length(layers));
    costs_for_all_layers_SA  = zeros(length(layers), total_number_of_experiments);
    for i=1:length(layers)
        
        fprintf("The simulation is for %d-layers Linear SA.....", layers(i))
        
        % number of betas/layers/time steps
        num_layers = layers(i);
        
        
        % Determine beta living time and loading sim. parameters into FPGA: Time steps/sweeps
        synthesized_beta_length = 61; % the one is for the zero
        num_repetition = (synthesized_beta_length-1)/num_layers; 
        num_sweeps_layer = 720;%number of sweeps/unique beta value 
        num_sweeps_synthezied_beta = num_sweeps_layer/num_repetition; % this has to be integer
        tictoc_counter_limit = num_experiments_inside_FPGA*synthesized_beta_length*num_sweeps_synthezied_beta; % total annealing time
        num_clk_limit = num_sweeps_synthezied_beta; % each synthesized beta living time
        
        % Write and Take counters for snapshot modules
        write_snapshot_counter = ceil(num_words_pbit*pbit_clk/bram_read_clk)+1;
        take_snapshot_counter = 1; % every clk cycle
    
    
    
        % write on FPGA each beta living time and the whole time for all betas 
        % % write on FPGA each beta living time and the whole time for all betas 
        writememory(mem,tictoc_counter_limit_address,uint32(tictoc_counter_limit), 'BurstType','Increment');
        writememory(mem,take_snapshot_counter_address,take_snapshot_counter, 'BurstType','Increment');
        writememory(mem,write_snapshot_counter_address,write_snapshot_counter, 'BurstType','Increment');
        writememory(mem,num_sweep_address,num_experiments_inside_FPGA,'BurstType','Increment'); % experiments inside FPGA
        writememory(mem,num_clk_limit_address,uint32(num_clk_limit), 'BurstType','Increment');
    
        if log2(tictoc_counter_limit) > 32
            error("You exceeded the tictoc counter limit capacity")
        end
    
    
    
        beta = 2*ones(num_layers,1);%linspace(2, 2, num_layers);
        
    
        %inference
        %inference
        cost_for_all_exp = inference_snapshot(num_experiments_outside,num_experiments_inside_FPGA, 10, J, J(1:216,1:216), beta);
        %costs_for_all_layers_SA(i,:)=  cost_for_all_exp;
        success_prob_SA(i) = length(find(cost_for_all_exp<=-putative_ground_state))/(total_number_of_experiments);
        disp(success_prob_SA)
    end
    
    disp(success_prob_SA)
    %file_name_3 = sprintf('success_prob_L6_1e5_run_%d_all_2', run);
    %save(file_name_3, "success_prob_SA")
    % disp(mean(costs_for_all_layers_SA, 2))
    file_name_1 = sprintf('all_costs_SA_L6_1e5_run_%d_15_layers', 0);
    save(file_name_1, "cost_for_all_exp")

end





