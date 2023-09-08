function [fastcore_model, flux_fastcore] = Build_FASTCORE(model, met_DATA, num_sample)
 
% Building metabolic models for each morphological stage using fastcore method and flux evaluation

%Input:
% model - Recon_3D (.mat)
% met_DATA - Metabolic gene expression data for each morphological stage
% num_sample - number of expression columns

%Output
% fastcore_model - context-specific model created using fastcore
% flux_fastcore - flux distribution of fastcore models

initCobraToolbox()

%%% Parsed GPR %%%
[parsedGPR,corrRxn] = extractGPRs(model);

model.c(model.c~=0)=0; % removing objective function
[fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, model_consistent] = findFluxConsistentSubset(Recon3D_processed);


for i=1:num_sample 

category={'week1' 'Hrs24' 'mono'};
exp=mapGeneToRxn(model,met_DATA(:,1),double(met_DATA(:,i+1)),parsedGPR,corrRxn);
exp(find(exp==0))=-1;

options.solver = 'fastcore';


%%% Threshold fixing %%%
% Th-1:
[a,b]=sort(exp,'descend');
threshold=a(round(0.1*length(a)));
intersect(find(model_consistent.fluxConsistentRxnBool),find(exp>threshold));
options.core=ans;

% Th-2:
% intersect(find(model_consistent.fluxConsistentRxnBool),find(exp>min(a(a~=0))));
% options.core=ans;

%%% Model Building %%%
fastcore_model = createTissueSpecificModel(model, options);

%%% Flux evaluation %%%
reaction_exp_ref = mapGeneToRxn(fastcore_model,met_DATA(:,1),double(met_DATA(:,i+1)),parsedGPR,corrRxn);
sampler='gp';
flux_fastcore = optimal_flux_eval(fastcore_model,reaction_exp_ref,sampler);

str=['fastcore_new_models_without_biomass_ATP_demand.flux_',category{i},'=','fastcore_model']; 
str1=['fastcore_new_models_without_biomass_ATP_demand.model_',category{i},'=','flux_fastcore'];

eval(str);
eval(str1);

% clear a b model str str1 fastcore_model exp options flux_fastcore reaction_exp_ref

end
return