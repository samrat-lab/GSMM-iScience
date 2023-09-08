function [met_score, pval] = Metabolic_functionality_score(model, queried_model)

% Evaluating metabolic functionality scores of input model and finding its statistical significance value

%%% Input:
% model - Recon_3D(.mat)
% queried_model - model structure (.mat)
% keep biomass_rxn.xlsx in test directory

%%% Output:
% met_score - Metabolic functionality score of the queried_model
% pval - Significance of the functionality score


%% Metabolic functionality score

initCobraToolbox()

[~,bmMets,~] = xlsread('biomass_rxn.xlsx');
bmMets = setdiff(bmMets,{'atp[c]','adp[c]','pi[c]','h2o[c]','h[c]'});
numTests = numel(bmMets) + 1;
comp = {'[c]','[e]','[g]','[l]','[m]','[n]','[r]','[x]'};
totTest = 0;
for i = 1:numel(bmMets)
    
    cMet = bmMets{i}(1:end-3);
    cTest = 0;
    for j = 1:numel(comp)
        met = [cMet,comp{j}];
        if ~isempty(intersect(queried_model.mets,met))
            
            nm = addReaction(queried_model,['BMS_',met],{met},-0.5,false,0,1000,0,'Biomass Metabolite Sink');
            nm = changeObjective(nm, ['BMS_',met]);
            sol = optimizeCbModel(nm);
            
            if sol.f > 1e-6
                cTest = 1;
                mets{i}=met;
                
                break;
            end
        end
    end
    totTest = cTest + totTest;
    a(i)=sol.f;
end

warning off all
nm = addReaction(model,'BMS_atp(c)',{'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'},[-0.5,-0.5,0.5,0.5,0.5],false,0,1000,0,'Biomass Metabolite Sink');
nm = changeRxnBounds(nm, 'DM_atp_c_',0,'l');
warning on all
nm = changeObjective(nm, 'BMS_atp(c)');
sol = optimizeCbModel(nm);
if sol.f > 1e-6
    totTest = totTest + 1;
end
a(i+1)=sol.f;

met_score=num2str(totTest);

clear 'totTest' 'numTests' 'mets'

%% Functionality score significance

for iter=1:1000
    iter
    a=[];
    
    randperm(length(model.rxnNames));
    tempmod=removeRxns(model, model.rxns(ans(1:(13543-length(queried_model.rxns)))));
    
    [~,bmMets,~] = xlsread('biomass_rxn.xls');
    bmMets = setdiff(bmMets,{'atp[c]','adp[c]','pi[c]','h2o[c]','h[c]'});
    numTests = numel(bmMets) + 1;
    comp = {'[c]','[e]','[g]','[l]','[m]','[n]','[r]','[x]'};
    totTest = 0;
    for i = 1:numel(bmMets)
        
        cMet = bmMets{i}(1:end-3);
        cTest = 0;
        for j = 1:numel(comp)
            met = [cMet,comp{j}];
            if ~isempty(intersect(tempmod.mets,met))
                
                nm = addReaction(tempmod,['BMS_',met],{met},-0.5,false,0,1000,0,'Biomass Metabolite Sink');
                nm = changeObjective(nm, ['BMS_',met]);
                sol = optimizeCbModel(nm);
                
                if sol.f > 1e-6
                    cTest = 1;
                    mets{i}=met;
                    
                    break;
                end
            end
        end
        totTest = cTest + totTest;
        a(i)=sol.f;
    end
    
    warning off all
    nm = addReaction(tempmod,'BMS_atp(c)',{'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'},[-0.5,-0.5,0.5,0.5,0.5],false,0,1000,0,'Biomass Metabolite Sink'); %Make a copy of ATP demand
    nm = changeRxnBounds(nm, 'DM_atp_c_',0,'l');
    warning on all
    nm = changeObjective(nm, 'BMS_atp(c)');
    sol = optimizeCbModel(nm);
    if sol.f > 1e-6
        totTest = totTest + 1;
    end
    a(i+1)=sol.f;
    
    flux_iter(:,iter)=a;
    
    score(iter)=totTest;
end

pval=sum(score>met_score)/1000;

return

