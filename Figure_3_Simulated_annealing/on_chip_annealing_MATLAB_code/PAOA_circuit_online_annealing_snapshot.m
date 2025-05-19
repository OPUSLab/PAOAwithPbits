function[cost] = PAOA_circuit_online_annealing_snapshot(NE,num_sweeps_inside_FPGA, NR, J, B)

% Description:
% This function solves Eq.1 and Eq.2 using MCMC technique in FPGA.
% It is used in optimization to find a good beta schedule, and can be used
% in inference as well.
%------------------------------------------------------
% Inputs:
% - Input_1: number of experiments (NE)
% - Input_2: J matrix-- problem graph
% - Input_3: beta (B), the cooling schedule
%
%
% -------------------------------------------------------------------
% Outputs:
% - cost/energy: the total energy measured at the end of the cooling
% process
% -----------------------------------------------------
%
% Developer's Name: Abdelrahman, Abdelrahman
% Date: 30th of Mar, 2024


%% FPGA interface and addresses/ Memory mappings
mem = aximanager('Xilinx','interface','pcie');

reset_tictoc_address = hex2dec('00000000');
s_address = hex2dec('00000020');


% trigger registers
weight_trigger_beta_address = hex2dec('004229A8');

% flags 
frozen_flag_address = hex2dec('004229AC');
weight_load_done_flag_address = hex2dec('004229B0');
clamp_flag_address = hex2dec('004229BC');


% beta_counter registers
reset_beta_address = hex2dec('004229B8');
beta_address = hex2dec('00416E20');





%% Extending beta FPGA

synthesized_beta_length = 60;

C = repmat(B', synthesized_beta_length/length(B), 1);
extended_beta = C(:)';
extended_beta = [0 extended_beta];
fixed_a = 4; fixed_b =5; % fixed points as ([s]a.b), a for int, b for fraction

%% PAOA-FPGA Optimization
beta = extended_beta; % : adding zero to ensure independence [0  extended_beta(1:end-1)];%;
beta_fpga = fixed_point_weights_beta(fixed_a, fixed_b, beta);
writememory(mem,beta_address,beta_fpga,'BurstType','Increment');
writememory(mem,weight_trigger_beta_address,1,'BurstType','Increment'); % shifting all weights
% wait until J,h, and beta are loaded
loaded = 0;
while ~loaded
    loaded= readmemory(mem,weight_load_done_flag_address, 1, 'BurstType','Increment');
end

num_words_pbit = ceil(length(J)/32);
total_number_of_words = num_words_pbit*(num_sweeps_inside_FPGA);

s_all = zeros(total_number_of_words, NE);
num_pbits = length(J);

for exper=1:NE

    
    % clear beta counter
    RESET_beta=1;
    writememory(mem,reset_beta_address,RESET_beta, 'BurstType','Increment');
    
    
    
    % % % needs to be removed from the design
    clamp_flag = 0;
    writememory(mem,clamp_flag_address,clamp_flag, 'BurstType','Increment');
    

    % p-bit activation
    ENABLE = 1; % reset the tictoc counter-- start from the begining
    writememory(mem,reset_tictoc_address,ENABLE, 'BurstType','Increment');
    


    % wait for until pbit is frozen PAOA circuit. Frozen will be one once
    % tic_toc counter limit (<= 32 bit) is hit that we have already specified
    
    frozen = 0;
    while ~frozen
        frozen = readmemory(mem,frozen_flag_address, 1, 'BurstType','Increment');
    end
    

    s_all(:,exper) = double(readmemory(mem,s_address,total_number_of_words, 'BurstType','Increment'));
    
   
    
end


% reshaping and conversion
s_all = reshape(s_all, 1, total_number_of_words*NE);
msbfirst=0;
s_bits = bit_conversion_mex(uint32(s_all), msbfirst); % Optimized MEX function
s_reshaped =(reshape(s_bits(:),num_words_pbit*32,num_sweeps_inside_FPGA*NE))';
s = s_reshaped(:,1:num_pbits);
s= double(s);

cost = compute_avg_energy_optimized(s, J)/10; %-0.5*m*J*m';10 is NR here
disp(cost)

return






