function non_ess_genes = Filtering_non_ess_genes(DATA,model,genelist)

% Building Eflux models corresponding to each normal ovarian patient data from GTEx and evaluating growth rate
% before and after single gene knockout in these models

%%%% Input:
% DATA - matrix representing the expression value of metabolic genes where column denotes samples
% model - Recon3D (.mat)
% genelist - genes catalysing perturbed reactions of each transition

%%%% Output:
% non_ess_genes - genes with growth rate ratio (growth rate after gene knockout/ growth rate before gene knockout)
%                 equal to 1 in atleast on-third normal samples

%% Building reaction expression array for each normal ovarian model

rxn_exp=zeros(length(model.rxns),size(DATA,2));
for m=1:size(DATA,2)
    x=DATA(:,m);
    m
    parfor i=1:length(model.rxns)
        
        if ~isempty(model.grRules{i})
            ans1=strsplit(model.grRules{i},'or');
            A=zeros(length(ans1),1);
            for j=1:length(ans1)
                ans2=strfind(ans1{j},'and');
                if isempty(ans2)
                    ans3=ans1{j}(setdiff(1:length(ans1{j}),union(union(strfind(ans1{j},' '),strfind(ans1{j},')')),strfind(ans1{j},'('))));
                    [a,b]=intersect(model.genes,ans3);
                    A(j,1)=x(b);
                else
                    ans4=strsplit(ans1{j},'and');
                    B=zeros(length(ans4),1);
                    for k=1:length(ans4)
                        ans5=ans4{k}(setdiff(1:length(ans4{k}),union(union(strfind(ans4{k},'('),strfind(ans4{k},')')),strfind(ans4{k},' '))));
                        [a,b]=intersect(model.genes,ans5);
                        B(k,1)=x(b);
                    end
                    if length(find(B))==0
                        A(j,1)=0;
                    else
                        A(j,1)=min(B(find(B)));
                    end
                end
            end
            rxn_exp(i,m)=sum(A);
        end
    end
end

%% Model Building

A=rxn_exp/max(max(rxn_exp));

initCobraToolbox()

for i=1:size(DATA,2)
    i
    model=Recon3D;
    model.ub(find(A(:,i)))=A(find(A(:,i)),i);
    model.ub(find(A(:,i)==0))=1;
    
    model.lb(intersect(find(model.lb),find(A(:,i))))=-A(intersect(find(model.lb),find(A(:,i))),i);
    model.lb(setdiff(find(model.lb),find(A(:,i))))=-1;
    
    modelWT=model;
    sol = optimizeCbModel(modelWT);
    grRateWT(i,1)=sol;
    
    for j=1:length(genelist)
        modelWT=model;
        find(string(modelWT.genes)==genelist(j,1));
        idx=find(modelWT.rxnGeneMat(:,ans));
        rxns=modelWT.rxns(find(modelWT.rxnGeneMat(:,ans)));
        
        modelWT.ub(idx)=0;
        modelWT.lb(idx)=0;
        modelDel=modelWT;
                
        sol1 = optimizeCbModel(modelDel);
        grRateKO(j,i)=sol1;
        
        grRatio(j,i)=grRateKO(j,i)/grRateWT(i,1);
        
        clear modelWT modelDel sol sol1
        
    end
end

%% Filtering non-essential genes

rank=[];
for i=1:size(grRatio,1)
rank(i,1)=sum(grRatio(i,:)==1);
non_ess_genes=genelist(rank>=55);
end

return