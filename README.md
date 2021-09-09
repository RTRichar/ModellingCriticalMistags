# ModellingCriticalMistags
A resource for controlling critical mistag-associated false discoveries in metagenetic data

## Demonstration of obtaining FDR thresholds for a dataset
- Load functions from *CM_Functions.R* into R environment
- Execute step by step analysis provided in *Analysis_Uniform.R* to obtain FDR thresholds for example *Dataset1.csv* and *Dataset2.csv.* 
- For *Dataset1.csv,* column "Dat1_Exp" is the sum of sequences per taxon across all biological samples. Column "Dat1_Cont" is the sum of sequences per taxon across control samples. Columns prefixed with "NL" labels are no-library negative controls and columns prefixed with "PC" labels are positive controls.
- For *Dataset2.csv,* column "Dat2_Exp" is the sum of sequences per taxon across all biological samples. Column "Dat2_Cont" is the sum of sequences per taxon across control samples. Columns prefixed with "Negative" labels are no-library negative controls.
- When using on new datasets, if data show evidence that the distribution of mistags across dual-index bins is not uniform (i.e., hhen the relationship between total sequence abundance and mistag abundance per taxon is weak or there is spatial autocorrelation in sequence abundances across the dual-index matrix), *Analysis_Poisson.R* provides an alternate modelling proceedure to accomodate such data. 
