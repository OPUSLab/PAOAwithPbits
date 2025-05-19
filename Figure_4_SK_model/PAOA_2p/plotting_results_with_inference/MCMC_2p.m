function[E] = MCMC_2p(NE, beta, J)


parfor exper=1:NE
    m = ones(length(J),1);
    
    for step=1:size(beta,1)
        % inner loop for updating the whole network
        for i=1:length(J)
            % first schedule
            if i< length(J)/2
                beta_value = beta(step, 1);
            else % second schedule
                beta_value = beta(step, 2);
            end
            I = beta_value*(J(i,:)*m);
            m(i) = sign(tanh(I)-(-1+2*rand));
        end 
    end

    E(exper) =-0.5*m'*J*m/sqrt(length(m));    
end

 E= E./length(J); 

end