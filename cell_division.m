%% Cell division with random segregation of mitochondria

function s = cell_division(S)
    n_wt=S(1,:);                    %%% Number of wildtype
    n_mut=S(2,:);                   %%% Number of mutants
    s=zeros(2,numel(S));   
    M1 = sum(n_mut);                %%% Total number of wildtypes
    W1 = sum(n_wt);                 %%% Total number of mutants
    mutRandom = rand([M1,1]);       %%% Generate list of random numbers...
    wtRandom = rand([W1,1]);        %%% ...for wt and mutants
    x_Mut=1;
    x_Wt=1;
    
    for j=1:length(n_mut)
        if n_mut(j)>0               
                                    %%% If the cell contains mutants...
                                    %%% Segregate wildtype and mutants
                                    %%% independently
                nMutants = sum((mutRandom(x_Mut:x_Mut-1+n_mut(j))>0.5));
                nWildtype = sum((wtRandom(x_Wt:x_Wt-1+n_wt(j))>0.5));
        else
                                    %%% Otherwise, only segregate wt
                nWildtype = binornd(n_wt(j),.5);
                nMutants=0;
        end
                s(1,2*j-1) = nWildtype;
                s(1,2*j) = n_wt(j)-nWildtype;                
                s(2,2*j-1) = nMutants;                
                s(2,2*j) = n_mut(j)-nMutants;     
                x_Mut = x_Mut+n_mut(j);
                x_Wt = x_Wt+n_wt(j);
    end
    s = s(:,not(isnan(sum(s))));
end

