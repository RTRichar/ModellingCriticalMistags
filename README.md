# ModellingCriticalMistags
A resource for controlling critical mistag-associated false discoveries in metagenetic data

## Demonstration of obtaining FDR thresholds for a dataset
- Load functions from *CM_Functions.R* into R environment
- Execute step by step analysis provided in *Analysis_Uniform.R* to obtain FDR thresholds for example *SimulatedDataset.csv* 
- When using on new datasets, if data show evidence that the distribution of mistags across dual-index bins is not uniform, *Analysis_Poisson.R* provides an alternate modelling proceedure to accomodate such data
  -  A bootstrapped Kolmogorov–Smirnov test can be used to test whether data meet this assumption

## Citation and further information
Richardson, RT. (2022) Controlling critical mistag-associated false discoveries in metagenetic data. *Methods in Ecology and Evolution*, 13: 938-944. https://doi.org/10.1111/2041-210X.13838
