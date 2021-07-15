%% Mitochondrial replication

function s = mito_replication2(S,mu)
    n_wt=S(1,:);
    newMutants= random('binomial',n_wt,mu);
    s = 2*S + [-newMutants; +newMutants];
end

