# Alzheimer's Disease Biomarkers - Rare Variant Analyses

------

Analysis scripts and full summary statics data of the manuscript [**"Whole-exome rare variant analysis of Alzheimer’s disease and related biomarker traits"**](https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/alz.12842) by Küçükali _et al._, Alzheimer's and Dementia (2022), https://doi.org/10.1002/alz.12842 .

------

This repository contains the full association results of our study, together with the gene-based rare variant association analysis scripts used to generate ([**PLINKSKAT.R**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/analysis_scripts/PLINKSKAT.R) & [**PLINKMetaSKAT.R**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/analysis_scripts/PLINKMetaSKAT.R)) and annotate & visualize ([**Annotate_PLINKSKAT_Results.R**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/analysis_scripts/Annotate_PLINKSKAT_Results.R) & [**Annotate_PLINKMetaSKAT_Results.R**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/analysis_scripts/Annotate_PLINKMetaSKAT_Results.R)) these results.

More details and descriptions of these are available in respective subdirectories [**association_results**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/association_results) and [**analysis_scripts**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/analysis_scripts).

## Quick Start for Installing Analysis Scripts

You can use the provided [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) environment file [**environment_plinkskat.yml**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/analysis_scripts/environment_plinkskat.yml) to create an environment and load these packages in this environment. You can also clone the repository and add the R scripts into this environment:

```sh
git clone https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses.git
conda env create -f analysis_scripts/environment_plinkskat.yml
conda activate plinkskat
ln -s analysis_scripts/*.R ~/miniconda3/envs/plinkskat/bin/.
```

Or alternatively, using [`mamba`](https://mamba.readthedocs.io/en/latest/installation.html):

```sh
git clone https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses.git
mamba env create -f environment_plinkskat.yml
mamba activate plinkskat
ln -s analysis_scripts/*.R ~/miniconda3/envs/plinkskat/bin/.
```
After this, you can execute the R scripts (PLINKSKAT.R, PLINKMetaSKAT.R, Annotate_PLINKSKAT_Results.R, and Annotate_PLINKMetaSKAT_Results.R) directly from the command line. For more information, see [**analysis_scripts**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/analysis_scripts) subdirectory.
