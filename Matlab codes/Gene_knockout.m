function corr_mat = Gene_knockout(queried_model,ref_model,ref_flux,non_ess_genes,reaction_exp,pert_rxns,sampler)

% Performing single gene-knockout in the queried model using non-essential genes catalysing reactions perturbed 
% during a transtion, henece finding potential targets for disease reversal.

%%%% Input:
% queried_model - model structure (.mat)

% ref_model - initial model structure (.mat)
% For transition-1: initial model is model corresponding to monolayer
% For transition-2: initial model is model corresponding to 24 hours spheroids (moruloid)

% ref_flux - flux values of reactions in initial model of a transition (.mat)

% non_ess_genes - genes that are non-essential for normal ovarian cell growth (as returned by 'Filtering_non_ess_genes.m')

% reaction_exp - reaction expression array of the ref_model (as returned by 'mapGeneToRxn.m')

% pert_rxns - reactions perturbed during the transition (.mat)

% sampler - 'gp' (string identifier)

%%%% Output:
% corr_mat - correlation (between the flux values of perturbed reactions during a transition) 
%            obtained after gene-knockout

initCobraToolbox()

%% Finding reactions catalysed by non_ess_genelist

rxns={};
idx={};
for i=1:length(non_ess_genes)
    find(string(queried_model.genes)==non_ess_genes(i,1));
    idx{i,1}=find(queried_model.rxnGeneMat(:,ans));
    rxns{i,1}=model.rxns(find(queried_model.rxnGeneMat(:,ans)));    
end

[~,~,ind]=intersect(pert_rxns,string(queried_model.rxns));
[~,~,ind1]=intersect(pert_rxns,string(ref_model.rxns));

init_corr=corr(que_flux(ind),ref_flux(ind1));

%% Performing single gene-knockout

corr_mat=[];
for i=1:length(non_ess_genes)
    model=queried_model;
    knock_rxn=idx{i,1};
    
    model.ub(knock_rxn)=0;
    model.lb(knock_rxn)=0;

    single_gko=optimal_flux_eval(model,reaction_exp,sampler);
   
    corr_mat(i,1)=corr(single_gko(ind,1),ref_flux(ind1,1));
    i
    clear 'model' 'solution'
end

return