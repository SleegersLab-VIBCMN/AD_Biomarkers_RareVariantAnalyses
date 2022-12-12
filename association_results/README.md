# AD Biomarkers Rare Variant Analyses - Association Results

Full summary statistics of the gene-based rare variant association studies (RVASs) conducted for the study [**"Whole-exome rare variant analysis of Alzheimer’s disease and related biomarker traits"**](https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/alz.12842) by Küçükali _et al._, Alzheimer's and Dementia (2022), https://doi.org/10.1002/alz.12842 .

RVAS summary statistics for 17 AD-related traits in [**EMIF-AD**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/association_results/EMIF_AD) MBD WES cohort, for 13 traits in [**ADNI**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/association_results/ADNI) WGS cohort, and for the [**meta-analyses**](https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses/tree/main/association_results/Meta_EMIF_AD_ADNI) of 13 traits tested in both cohorts are represented in each subdirectory where protein-altering and LoF model results are provided in separate subdirectories (named **"Protein_altering"** and **"LoF"**).

The name of each TSV file follows the **TRAIT.MODEL.TEST.COHORT.tsv** pattern which indicate:

1) Name of the tested trait (17 AD-related traits)
2) Name of the model (either protein-altering or LoF)
3) Name of the test (either rare SKAT-O or rare MetaSKAT-O tests)
4) Name of the cohort or meta-analysis of cohorts (either EMIF_AD, ADNI, or Meta_EMIF_AD_ADNI)

Summary statistics TSV files have two columns:

1) **gene_name** : Name of the gene tested for this trait and model.
2) **p_value** : p-values of associations (from either SKAT-O or MetaSKAT-O tests)

For more detailed information on these, including the detailed information of top 10 associations in each trait and model ([Supplementary Tables S8-S11](https://alz-journals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Falz.12842&file=alz12842-sup-00002-tables.xlsx)), please see the [**manuscript**](https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/alz.12842).



