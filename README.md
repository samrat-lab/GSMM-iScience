# Targeting metabolic fluxes reverts metastatic transitions in ovarian cancer

## Matlab codes description

**iMAT.m, Build_iMAT.m**: Codes for building iMAT models

**FASTCORE.m, Build_FASTCORE.m**: Codes for building FASTCORE models

**INIT.m, Build_INIT.m**: Codes for building INIT models

**extractGPRs.m**: Maps the GPR rules of the model to a specified format that is used by the model extraction methods

**mapGeneToRxn.m**: Maps gene expression to reaction expression using the GPR rules

**createTissueSpecificModel.m**: Creates draft tissue-specific model from mRNA expression data using different MEMs

**findFluxConsistentSubset.m**: Finds the subset of S that is flux-consistent

**optimal_flux_eval.m**: Evaluates optimal flux using GP sampler

**Metabolic_functionality_score.m**: Evaluates active metabolic functionality score and the statistical significance of the score

**Reaction_modulation.m**: To obtain disease-driving reaction modules

**Filtering_non_ess_genes.m**: Code for filtering non-essential genes

**Gene_knockout.m**: Evaluation of single gene knockout

**Drug_repurposing.m**: Code for performing drug repurposing

**biomass_rxn.xlsx**: Essential metabolites for evaluating metabolic functionality score

## Supplementary data files description

**Data S1**: Proteomics data file containing protein expressions in all biological replicates of each morphological stage

**Data S2**: List of perturbed pathways enriched using protein expression data

**Data S3**: Comparison table of model extraction methods used in the study

**Data S4**: Results of reaction modulation exercise

**Data S5**: Results of in-silico gene knockout

**Data S6**: Results of in-silico drug repurposing

**Data S7**: Pathway analysis results of drugs that target perturbed reactions and are common to both transitions
