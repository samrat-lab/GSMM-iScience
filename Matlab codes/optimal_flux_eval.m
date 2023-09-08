function opt_flux = optimal_flux_eval(model,reaction_exp,sampler)

% Sampling the solution space using gp sampler and evaluating the optimal flux having the maximum
% correlation with input reaction exp
%Input:
%  model - model structure (.mat)
%  reaction_exp - Reaction expression array  (as returned by'mapGeneToRxn.m')
%  sampler type - 'gp' (string identifier)

%Output
%  opt_flx - Optimal flux distribution of input model

initCobraToolbox()


if ~isempty(reaction_exp) 


changeCobraSolver('gurobi6', sampler)
gurobi_setup
storedata = gpSampler(model);

sampled_mat=storedata.points;


%%%FIND OPT

for i =1: size(sampled_mat,2)
    corrmat(i,1)=corr(sampled_mat(:,i),reaction_exp);
    
    
end
   opt_flux=sampled_mat(:,find(corrmat==max(corrmat)));

    
% save storedata storedata

end
return


