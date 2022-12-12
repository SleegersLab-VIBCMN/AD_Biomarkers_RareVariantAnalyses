#!~/miniconda3/bin/Rscript
## Run gene-based (rare variant) association tests as implemented in SKAT package using PLINK format files
## Written by Fahri Kucukali for the analyses conducted for Kucukali et al., Whole‐exome rare‐variant analysis of Alzheimer's disease and related biomarker traits, Alzheimer's & Dementia (2022). http://doi.org/10.1002/alz.12842
## For issues: https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses or e-mail Fahri.Kucukali@uantwerpen.vib.be

suppressWarnings(suppressMessages(library(SKAT)))
suppressWarnings(suppressMessages(library(argparser)))

p <- arg_parser("Run SKAT using PLINK format files, with different modes and options. Publication Version 1.1, December 2022")
p <- add_argument(p, "--prefix", help = "Prefix of your PLINK bed/bim/fam files.", type = "character", nargs = '+')
p <- add_argument(p, "--dummy", help = "Dummy name for your SSD & INFO files. This prevents any overwriting issues if you provide a different name for each analysis, also useful for downstream analyses with PLINK_METASKAT.R as it will be using these SSD & INFO files.", type = "character", nargs = '+')
p <- add_argument(p, "--SetID", help = "A set file (without header) with a 1st column of Gene Names and 2nd column of Variant IDs (assigned to these genes). Gene-based tests will be performed for these sets.", type = "character", nargs = '+')
p <- add_argument(p, "--variable", help = "Select a variable type for association testing: binary or continuous (and if continuous, please specificy --phenotype)", type = "character", nargs = '+')
p <- add_argument(p, "--phenocovar", help = "Phenotype and covariate file with header (including FID and IID columns as first two columns respective to PLINK fam file). It can also contain other phenotypes that can be provided into association testing with --phenotype. All phenotypes and covariates should be numerically coded. Optional.", nargs = '*')
p <- add_argument(p, "--covnames", help = "Space or comma-separated list of covariate names you want to use from your --phenocovar file. As all covariates should be numerically coded, binary covariates can be indicated as e.g. 0 and 1 or 1 and 2, but for more categories please use dummy covariates. If your PLINK fam file contains Sex column, this can be also supplied if you type in Sex. Optional.", type = "character", nargs = '*')
p <- add_argument(p, "--phenotype", help = "Column name of your binary or continuous phenotype variable in the --phenocovar file; if not present, PLINK binary phenotype from the PLINK fam file will be used (where 2=cases and 1=controls encoding is used). In the phenocovar file though, binary phenotypes should be encoded as 1=cases and 0=controls. Optional.", nargs = 1)
p <- add_argument(p, "--out", help = "Prefix of the output results .tsv file and log file. Optional.", type = "character", nargs = 1, default = "output_PLINKSKAT")
p <- add_argument(p, "--method", help = "Option to change test method: default is SKATO, but could be also Burden or SKAT", type = "character", nargs = 1, default = "SKATO")
p <- add_argument(p, "--missing_cutoff", help = "Defines the call rate cutoff for variants, default is 0.15 meaning that keeping all variants present in at least 85% of individuals.", nargs = 1, type = "numeric", default = 0.85)
p <- add_argument(p, "--impute_method", help = "Selects the impute.method to be used in tests, enter either bestguess or fixed.", type = "character", nargs = '+')

arg <- parse_args(p)

# argument checks
if (is.na(arg$prefix) || is.na(arg$dummy) || is.na(arg$SetID) || is.na(arg$variable) || is.na(arg$impute_method)) {
    message("Error: exiting because at least one of the following required arguments is missing: --prefix, --dummy, --SetID, --variable, --impute_method")
    q("no")
}

if ((arg$variable == "continuous") && is.na(arg$phenotype)) {
    message("Error: exiting because --phenotype argument is missing. Please provide it when you use --variable continuous, and include it in your --phenocovar file.")
    q("no")
}

if ((!is.na(arg$phenocovar)) && ((is.na(arg$covnames)) && (is.na(arg$phenotype)))) {
    message("Error: exiting because both of the following missing even though --phenocovar is present: --covnames and/or --phenotype. Please use at least one.")
    q("no")
}

File.Bed <- paste0(arg$prefix, ".bed")
File.Bim <- paste0(arg$prefix, ".bim")
File.Fam <- paste0(arg$prefix, ".fam")
File.SSD <- paste0(arg$dummy, ".ssd")
File.Info <- paste0(arg$dummy, ".info")
File.SetID <- arg$SetID
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFO <- Open_SSD(File.SSD,File.Info)

if (arg$variable == "binary") {
  if (is.na(arg$phenotype)) {
    message("Using default binary phenotype in the PLINK .fam file")
    if(is.na(arg$covnames)){
    FAM <- Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)
    obj <- SKAT_Null_Model(FAM$Phenotype ~ 1, data=FAM, out_type="D")
    } else if (!is.na(arg$covnames)){
    File.Cov <- arg$phenocovar
    covlist = gsub(","," ", arg$covnames)
    covlist_prepared = unlist(strsplit(covlist,split=" "))
    FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=TRUE, cov_header=TRUE)
    f = paste("FAM$Phenotype", paste(names(FAM[covlist_prepared]),collapse="+FAM$"),sep="~FAM$")
    obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="D")}
    output_skato <- SKATBinary.SSD.All(SSD.INFO, obj, method=arg$method, missing_cutoff=arg$missing_cutoff, impute.method=arg$impute_method)
    Close_SSD()
    write.table(output_skato$results, paste0(arg$out, ".tsv"), col.names=TRUE, quote=F, row.names=F, sep="\t")
    } else {
    message("Using indicated binary phenotype in the phenocovar file")
    NewPheno = arg$phenotype
    File.Cov <- arg$phenocovar
    FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=TRUE, cov_header=TRUE)
    if(is.na(arg$covnames)){
    f = paste0("FAM$", NewPheno, " ~ 1")
    obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="D")
    } else if (!is.na(arg$covnames)){
    covlist = gsub(","," ", arg$covnames)
    covlist_prepared = unlist(strsplit(covlist,split=" "))
    f = paste(paste("FAM$", NewPheno, sep=""), paste(names(FAM[covlist_prepared]),collapse="+FAM$"),sep="~FAM$")
    obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="D")}
    output_skato <- SKATBinary.SSD.All(SSD.INFO, obj, method=arg$method, missing_cutoff=arg$missing_cutoff, impute.method=arg$impute_method)
    Close_SSD()
    write.table(output_skato$results, paste0(arg$out, ".tsv"), col.names=TRUE, quote=F, row.names=F, sep="\t")
    }

} else if (arg$variable == "continuous") {
  NewPhenoCont = arg$phenotype
  message("Using indicated continuous phenotype in the phenocovar file")
  File.Cov <- arg$phenocovar
  FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=FALSE, cov_header=TRUE)
  if(is.na(arg$covnames)){
  f = paste0("FAM$", NewPhenoCont, " ~ 1")
  obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="C")
  } else if (!is.na(arg$covnames)){
  covlist = gsub(","," ", arg$covnames)
  covlist_prepared = unlist(strsplit(covlist,split=" "))
  f = paste(paste("FAM$",NewPhenoCont, sep=""), paste(names(FAM[covlist_prepared]),collapse="+FAM$"),sep="~FAM$")
  obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="C")}
  output_skato <- SKAT.SSD.All(SSD.INFO, obj, method=arg$method, missing_cutoff=arg$missing_cutoff, impute.method=arg$impute_method)
  Close_SSD()
  write.table(output_skato$results, paste0(arg$out, ".tsv"), col.names=TRUE, quote=F, row.names=F, sep="\t")
  
}

log_table <- as.data.frame(t(as.data.frame(arg)[,-c(1:3)]))
log_table <- cbind(rownames(log_table),log_table$V1)
colnames(log_table) <- c("argument","value")
write.table(log_table, paste0("log_", arg$out, ".txt"), col.names=T, quote=F, row.names=F, sep="\t")
