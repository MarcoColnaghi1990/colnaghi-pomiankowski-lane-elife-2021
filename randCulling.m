%% Culling with individual selection

function s = randCulling(s,k,xi)
    s = s(2,:)./sum(s);
    N = round(numel(s)/8);          %%% sample size
    W = max(1-k*s.^xi,0);           %%% fitness
    s = randsample(s,N,true,W);     %%% culling
end