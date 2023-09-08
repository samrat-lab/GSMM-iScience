function [INIT_model, flux_INIT] = Build_INIT(model, met_DATA, num_sample)
 
% Building metabolic models for each morphological stage using INIT method and flux evaluation

%Input:
% model - Recon_3D (.mat)
% met_DATA - Metabolic gene expression data for each morphological stage
% num_sample - number of expression columns

%Output
% INIT_model - context-specific model created using INIT
% flux_INIT - flux distribution of INIT models

initCobraToolbox()

%%% Parsed GPR %%%
[parsedGPR,corrRxn] = extractGPRs(model);

model.c(model.c~=0)=0; % removing objective function

for i=1:num_sample 
changeCobraSolver ('gurobi', 'MILP');
category={'week1' 'Hrs24' 'mono'};
exp=mapGeneToRxn(model,met_DATA(:,1),double(met_DATA(:,i+1)),parsedGPR,corrRxn);

options.solver = 'INIT';

%%% Threshold fixing %%%
% Th-1:
[a,b]=sort(exp,'descend');
threshold=a(round(0.1*length(a)));

exp1=exp;
exp1(exp>0)=5*log2(exp(exp>0)/threshold);
exp1(exp<=0)=-2;

options.weights=exp1;


% Th-2:
% exp1=exp;
% exp1(exp>min(a(a~=0)))=5*log2(exp(exp>min(a(a~=0)))/min(a(a~=0)));
% exp1(exp==min(a(a~=0)))=-1;
% exp1(exp<=0)=-2;
% 
% options.weights=exp1;

options.tol=1e-6;

%%% Model Building %%%
INIT_model = createTissueSpecificModel(model, options);

%%% Flux evaluation %%%
reaction_exp_ref = mapGeneToRxn(INIT_model,met_DATA(:,1),double(met_DATA(:,i+1)),parsedGPR,corrRxn);
sampler='gp';
flux_INIT = optimal_flux_eval(INIT_model,reaction_exp_ref,sampler);

str=['INIT_new_models_without_biomass_ATP_demand.flux_',category{i},'=','INIT_model']; 
str1=['INIT_new_models_without_biomass_ATP_demand.model_',category{i},'=','flux_INIT'];

eval(str);
eval(str1);

% clear a b model str str1 INIT_model exp options flux_INIT reaction_exp_ref

end
return