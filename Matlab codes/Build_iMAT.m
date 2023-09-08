function [iMAT_model, flux_iMAT]=Build_iMAT(model, met_DATA, num_sample)
 
% Building metabolic models for each morphological stage using iMAT method and flux evaluation

%Input:
% model - Recon_3D (.mat)
% met_DATA - Metabolic gene expression data for each morphological stage
% num_sample - number of expression columns

%Output
% iMAT_model - context-specific model created using iMAT
% flux_iMAT - flux distribution of iMAT models

initCobraToolbox()

%%% Parsed GPR %%%
[parsedGPR,corrRxn] = extractGPRs(model);

model.c(model.c~=0)=0; % removing objective function

for i=1:num_sample 
changeCobraSolver ('gurobi', 'MILP');
category={'week1' 'Hrs24' 'mono'};
exp=mapGeneToRxn(model,met_DATA(:,1),double(met_DATA(:,i+1)),parsedGPR,corrRxn);
exp(find(exp==0))=-1;

options.expressionRxns=exp;
options.solver = 'iMAT';

%%% Threshold fixing %%%
% Th-1:
[a,b]=sort(exp,'descend');
options.threshold_lb=min(a(a~=0));
options.threshold_ub=a(round(0.1*length(a)));

% Th-2:
% options.threshold_lb=mean(exp(exp>0))-.5*std(exp(exp>0));
% options.threshold_ub=mean(exp(exp>0))+.5*std(exp(exp>0));

%%% Model Building %%%
iMAT_model = createTissueSpecificModel(model, options);

%%% Flux evaluation %%%
reaction_exp_ref = mapGeneToRxn(iMAT_model,met_DATA(:,1),double(met_DATA(:,i+1)),parsedGPR,corrRxn);
sampler='gp';
flux_iMAT = optimal_flux_eval(iMAT_model,reaction_exp_ref,sampler);

str=['iMAT_new_models_without_biomass_ATP_demand.flux_',category{i},'=','iMAT_model']; 
str1=['iMAT_new_models_without_biomass_ATP_demand.model_',category{i},'=','flux_iMAT'];

eval(str);
eval(str1);

% clear a b model str str1 iMAT_model exp options flux_iMAT reaction_exp_ref

end
return