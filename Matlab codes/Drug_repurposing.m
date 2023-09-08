function corr_mat = Drug_repurposing(drug_mat, queried_model, ref_model, ref_flux, pert_rxns, reaction_exp,sampler)

% Performing in-silico drug repurposing using drugs and their targets from DrugBank database. 
% Filtered list of drugs targeting reactions in the queried model was used. 

%%%% Input:
% drug_mat - filtered list of drugs (and their targets) targeting reactions along in the queried model (.mat)
% drug_mat = [drug name, drug targets (Gene ID), reactions catalysed by drug targets];

% queried_model - model structure in which drug repurposing has to be performed (.mat)

% ref_model - initial model structure (.mat)

% ref_flux - flux values of reactions in initial model of a transition (.mat)
% For transition-1: initial model is model corresponding to monolayer
% For transition-2: initial model is model corresponding to 24 hours spheroids (moruloid)

% pert_rxns - reactions perturbed during the transition (.mat)

% reaction_exp - reaction expression array of the ref_model (as returned by 'mapGeneToRxn.m')

% sampler - 'gp' (string identifier)

%%%% Output:
% corr_mat - correlation (between the flux values of perturbed reactions during a transition) 
%            obtained after drug repurposing


initCobraToolbox()

[~,~,g1]=intersect(pert_rxns, string(queried_model.rxns));
[~,~,g2]=intersect(pert_rxns, string(ref_model.rxns));

drug_rxns=drug_mat{:,3};

for i=1:length(drug_rxns)
    
    model=queried_model;
    [~,tp,~]=intersect(string(model1.rxns),drug_rxns{i,1});
    
    for j=1:length(tp)
        
        %% for upper bounds
        if model.ub(tp(j))~=0
            model.ub(tp(j))=model.ub(tp(j))-.5*(model.ub(tp(j)));
        else % case of zero upper bound
            model.ub(tp(j))=0;
        end
        
        %% for lower bounds
        if model.lb(tp(j))~=0
            model.lb(tp(j))=model.lb(tp(j))-.5*(model.lb(tp(j)));
        else % case of zero lower bound
            model.lb(tp(j))=0;
        end
    end
    
    flux_value= optimal_flux_eval(model,reaction_exp,sampler);
    corr_mat(i,1)=corr(flux_value(g1,1),ref_flux(g2,1));
end
return
