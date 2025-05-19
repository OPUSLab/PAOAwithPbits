function[PDF_evolution] = PDF_time_evolution(num_experiments,...
                                             beta_optimized_final,J_cost, time_steps)


h = zeros(length(J_cost),1);
% number of p-bits
NM = size(J_cost,2);

% conversion to decimal
conv = [8 4 2 1];

PDF_evolution = zeros(2^NM, time_steps);


for exper=1:num_experiments
    
    m = sign(rand(1,NM)-0.5)';
    % initial distribution
    state = conv*(m+1)/2;
    PDF_evolution(state+1,1) = PDF_evolution(state+1,1)+1;
    beta_counter = 1;
    for timestep=2:time_steps        
        if time_steps >3
            beta_counter=1;
        end
        % inner loop for updating the whole network
        for i=1:NM
            I = beta_optimized_final(i,beta_counter)*(J_cost(i,:)*m+h(i));
            m(i) = sign(tanh(I)-(-1+2*rand));
        end
        

        state = conv*(m+1)/2;
        PDF_evolution(state+1,timestep) = PDF_evolution(state+1,timestep)+1;
        beta_counter=  beta_counter+1;
        
    end
    

end
PDF_evolution = PDF_evolution/num_experiments;

end