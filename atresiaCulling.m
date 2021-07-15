%% Culling with cell selection during follicular atresia

function s = atresiaCulling(s,k,xi)
    s1 = s(2,:)./sum(s);
    N = round(numel(s1)/8);          %%% sample size
    W = max(1-k*s1.^xi,0);   %%% fitness
    s1 = randsample(1:numel(s1),N,true,W);     %%% culling
    s = s(:,s1);
end