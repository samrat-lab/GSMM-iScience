function [rxn_mod,corr_opt] = Reaction_modulation(queried_model,que_flux,ref_model,ref_flux,pert_rxns,reaction_exp,sampler)

% Capturing core metabolic reaction modules responsible for driving disease progression

%%%% Input:
% queried_model - model structure (.mat)

% que_flux - flux values of reactions in queried model (.mat)

% ref_model - initial model structure (.mat)

% ref_flux - flux values of reactions in initial model of a transition (.mat)
% For transition-1: initial model is model corresponding to monolayer
% For transition-2: initial model is model corresponding to 24 hours spheroids (moruloid)

% pert_rxns - reactions perturbed during the transition (.mat)

% reaction_exp - reaction expression array of the ref_model (as returned by 'mapGeneToRxn.m')

% sampler - 'gp' (string identifier)

%%%% Output:
% rxn_mod - reaction modules driving disease progression
% corr_opt - correlation values

initCobraToolbox()

[~,~,ind]=intersect(pert_rxns,string(queried_model.rxns));
[~,~,ind1]=intersect(pert_rxns,string(ref_model.rxns));

init_corr=corr(que_flux(ind),ref_flux(ind1));

%% Single reaction modulation

for i=1:length(pert_rxns)
    rand_rxn=i;
    model=queried_model;
    
    if ref_flux(ind1(rand_rxn))<0
        model.ub(ind(rand_rxn))=0;
        model.lb(ind(rand_rxn))=ref_flux(ind1(rand_rxn));
    elseif ref_flux(ind1(rand_rxn))>0
        model.ub(ind(rand_rxn))=ref_flux(ind1(rand_rxn));
        model.lb(ind(rand_rxn))=0;
    end
    
    single_mod=optimal_flux_eval(model,reaction_exp,sampler);
    
    corr_mat(i)=corr(single_mod(ind,1),ref_flux(ind1,1));
    
    clear 'model' 'single_mod'
end

%% Multiple reaction modulation

rxn=[];
rewired_seed_rxn={};
[a1,b1]=sort(corr_mat,'descend');
corr_single_mod=init_corr;
model=queried_model;
seed=b1;
corr_opt=[];
rxn_mod={};

for l=1:length(a1)
    corr_single_mod=a1(l);
    rxn=[];
    
    model=queried_model;
    level_rxn=[];
    rand_rxn=seed(l);
    if ref_flux(ind1(rand_rxn))<0
        model.ub(ind(rand_rxn))=0;
        model.lb(ind(rand_rxn))=ref_flux(ind1(rand_rxn));
    elseif ref_flux(ind1(rand_rxn))>0
        model.ub(ind(rand_rxn))=ref_flux(ind1(rand_rxn));
        model.lb(ind(rand_rxn))=0;
    end
    seed_model=model;
    
    level_rxn=rand_rxn;
    level_corr=1;
    
    %% Leveling
    level=0;
    precorr=corr_single_mod;
    while (level_corr-precorr)>0.01
        precorr=corr_single_mod;
        level=level+1;
        rxn=[];
        
        for i=1:length(pert_rxns)
            
            rand_rxn=b1(i);
            
            premodel_level=model;
            if ref_flux(ind1(rand_rxn))<0
                model.ub(ind(rand_rxn))=0;
                model.lb(ind(rand_rxn))=ref_flux(ind1(rand_rxn));
            elseif ref_flux(ind1(rand_rxn))>0
                model.ub(ind(rand_rxn))=ref_flux(ind1(rand_rxn));
                model.lb(ind(rand_rxn))=0;
            end
            
            postmodel_level=model;
            flux_value=optimal_flux_eval(model,reaction_exp,sampler);
            corr(flux_value(ind,1),ref_flux(ind1,1));
            
            if ans-corr_single_mod>0.001
                corr_single_mod=ans;
                rxn=[rxn;i];
                model=postmodel_level;
            else
                model=premodel_level;
            end
            
            clear 'ans'
        end
        level_rxn=unique([level_rxn; rxn]);
        model=postmodel_level;
        level_corr=corr_single_mod;
    end
    corr_opt=[corr_opt corr_single_mod];
    rxn_mod{l,1}=level_rxn;
    l
end
return


