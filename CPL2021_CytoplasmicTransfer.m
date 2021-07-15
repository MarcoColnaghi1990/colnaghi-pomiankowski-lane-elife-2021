%% Cytoplasmic Transfer

%%% Parameters & Variables
mut_load_0 = 0.1;                   %%% Initial mutation load
mutrate = 1e-7 ;                    %%% Mutation rate/bp
mu = mutrate *16569;                %%% mtDNA genome-wide mutation rate
n_start=2^19;                       %%% Initial mtDNA copy number
n_iterations = 10;                  %%% Number of iterations
f=.5;                               %%% Fraction of transferred mitochondria per nurse cell
p_m=.25;                            %%% Probability of transfer of mutant mitochondria
p_w=.75;                            %%% Probability of transfer of wildtype mitochondria
PGC_0 = 32;                         %%% Initial number of primordial germ cells
dM = zeros(1,n_iterations);         %%% Mutation load shift between generations

%%% Simulation

for n=1:n_iterations
    
    n_wildtype = round((1-mut_load_0)*n_start);  %%% Initial number of wild-type mtDNA copies
    n_mutants = round(mut_load_0*n_start);       %%% Initial number of mutant mtDNA copies
    S = [n_wildtype; n_mutants];                 %%% State vector
    
    %%% Embryonic development: cell division without mtDNA replication
    t=1;
    
    for j=2:13
        t=t+1;
        S = cell_division(S);
        if t==12
            % Germline segregation
            S=datasample(S,PGC_0,2,'Replace',false);
            t=t+1;
            disp(S)
        end
    end
    
    %%% PGCs proliferation & Nurse cell growth
    for t=15:31
        S = mito_replication2(S,mu);
        S = cell_division(S);
        disp([t,numel(S)/2])
    end

    %%% Cytoplasmic transfer
    S = cytoplasmic_transfer(S,p_m,p_w,f);
    
    %%% Oocyte growth
    for t=32:41
        S = mito_replication2(S,mu);
        disp([t,mean(mean(sum(S)))])
    end
    
    mut_load = S(2,:)./sum(S);
    
    dM(n) = mean(mut_load)-mut_load_0;     %%% Mutation load shift between generations

end


