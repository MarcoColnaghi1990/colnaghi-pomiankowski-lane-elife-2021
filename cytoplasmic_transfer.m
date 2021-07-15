%% Simulates cytoplasmic transfer

function s = cytoplasmic_transfer(s,pMt,pWt,proportion)

n_transferred_mito = round(proportion*sum(s));
n_mutants=zeros(1,numel(s)/2);

for i = 1:numel(s)/2
    mito_list = [zeros(1,s(1,i)),ones(1,s(2,i))];
    mito_weight = [repmat(pWt,1,s(1,i)),repmat(pMt,1,s(2,i))];
    transferd = randsample(mito_list,n_transferred_mito(i),true,mito_weight);
    n_mutants(i) = sum(transferd);
end

n_wildtypes = n_transferred_mito - n_mutants;

s=[n_wildtypes;n_mutants];

x = 1:8:numel(s)/2;
new_s = s(:,x);
for i =1:7
    new_s = new_s + s(:,x+i);
end

s = new_s;

end
