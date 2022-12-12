#!~/miniconda3/bin/Rscript
## Run gene-based (rare variant) association meta-analyses as implemented in MetaSKAT package using PLINK format files
## Written by Fahri Kucukali for the analyses conducted for Kucukali et al., Whole‐exome rare‐variant analysis of Alzheimer's disease and related biomarker traits, Alzheimer's & Dementia (2022). http://doi.org/10.1002/alz.12842
## For issues: https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses or e-mail Fahri.Kucukali@uantwerpen.vib.be

suppressWarnings(suppressMessages(library(MetaSKAT)))
suppressWarnings(suppressMessages(library(SKAT)))
suppressWarnings(suppressMessages(library(argparser)))

p <- arg_parser("Run MetaSKAT using PLINK format files, with different modes and options. Publication Version 1.1, December 2022")
p <- add_argument(p, "--studyprefixnames", help = "Space or comma-separated list of study prefixes that are available in PLINK prefix, also will be used to name the intermediate files that were previously assigned with dummy option in PLINKSKAT.R.", type = "character", nargs = "*")
p <- add_argument(p, "--setIDnames", help = "Space or comma-separated list of SetID file names in the same order.", type = "character", nargs = "*")
p <- add_argument(p, "--studygenoOrdose", help = "Space or comma-separated list of study type in terms of genotype or dosage-based, enter either genotype or dosage in the same order of studypredixnames.", type = "character", nargs="*")
p <- add_argument(p, "--variable", help = "Select a variable type for association testing: binary or continuous (and if continuous, please specificy --phenotype).", type = "character", nargs = 1)
p <- add_argument(p, "--phenotypes", help = "Column names of your phenotype of interest variable in the --phenocovars files, should be again in order with list of studies.", type = "character", nargs = "*")
p <- add_argument(p, "--out", help = "Output file name.", type = "character", nargs = 1, default = "output_metaskat")
p <- add_argument(p, "--method", help = "Option to change test method: default is optimal (corresponds to SKATO), but could be also davies or more, see SKAT manual.", type="character", nargs = 1, default = "optimal")
p <- add_argument(p, "--missing_cutoff", help = "Defines the call rate cutoff for variants, default is 0.15 meaning that keeping all variants present in at least 85% of individuals", nargs = 1, type = "numeric", default = 0.15)
p <- add_argument(p, "--phenocovars", help="Phenocovar files with header (including FID and IID respective to PLINK file) in the same order.", type = "character", nargs = '*')
p <- add_argument(p, "--covnames", help="Covariate names you want to use in the same order, it will be separated by number of covariates. Optional.", type = "character", nargs = '*')
p <- add_argument(p, "--covnos", help = "Number of covariates for each study indicated in the same order. Optional.", nargs='*')
p <- add_argument(p, "--rcorr", help = "Only use when you want to specify r.corr value (so not SKAT-O test, change the mode to davies or something else; see SKAT manual), r.corr 0 is full SKAT model and 1 is full Burden models", default = NULL, nargs = '*')
p <- add_argument(p, "--impute_method", help = "Selects the impute.method to be used in SKAT test, enter either bestguess or fixed.", type = "character", nargs = '+')

arg <- parse_args(p)

# argument checks
if (is.na(arg$studyprefixnames) || is.na(arg$setIDnames) || is.na(arg$studygenoOrdose) || is.na(arg$variable) || is.na(arg$phenotypes) || is.na(arg$phenocovars) || is.na(arg$impute_method)) {
    message("Error: exiting because at least one of the following required arguments is missing: --studyprefixnames, --setIDnames, --studygenoOrdose, --variable, --phenotypes , --phenocovars, --impute_method")
    q("no")
}

studyprefixnames_justincase = gsub(","," ", arg$studyprefixnames)

studyprefixnames = unlist(strsplit(studyprefixnames_justincase,split=" "))

nstudyprefixnames <- length(studyprefixnames)

setIDnames_justincase = gsub(","," ", arg$setIDnames)

setIDnames = unlist(strsplit(setIDnames_justincase,split=" "))

studygenoOrdose_justincase = gsub(","," ", arg$studygenoOrdose)

studygenoOrdose = unlist(strsplit(studygenoOrdose_justincase,split=" "))

phenolist_justincase = gsub(","," ", arg$phenotypes)

phenolist = unlist(strsplit(phenolist_justincase,split=" "))

covfile_justincase = gsub(","," ", arg$phenocovars)

covfile_prepared = unlist(strsplit(covfile_justincase,split=" "))

if(!is.na(arg$covnames)){

covlist_justincase = gsub(","," ", arg$covnames)

covlist_prepared = unlist(strsplit(covlist_justincase,split=" "))

covnos_justincase = gsub(","," ", arg$covnos)

covnos_prepared = as.numeric(unlist(strsplit(covnos_justincase,split=" ")))

covnos_prepared_increasing <- cumsum(covnos_prepared)

covs_separated <- ""

for (index in 1:nstudyprefixnames) {

	if (index == 1) {
	first <- 1

	} else if (index > 1) {
	first <- 1 + covnos_prepared_increasing[index - 1]

	}

last <- covnos_prepared_increasing[index]

covs_separated[index] <- paste(covlist_prepared[first:last], collapse=" ")

}

studymerged <- as.data.frame(cbind(studyprefixnames,setIDnames,studygenoOrdose,phenolist,covfile_prepared,covs_separated))

File.MSSD.vec <- rep("", nstudyprefixnames)
File.MInfo.vec <- rep("", nstudyprefixnames)

    for (study.index in 1:nstudyprefixnames) {
     studyprefix <- as.character(studymerged$studyprefixnames)[study.index]
     if(as.character(studymerged$studygenoOrdose)[study.index] == "genotype"){
          File.Bed <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".bed")
          File.Bim <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".bim")
          File.Fam <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".fam")
          File.MSSD <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".mssd")
          File.MInfo <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".minfo")
          File.SetID <- as.character(studymerged$setIDnames)[study.index]
		  File.Cov <- as.character(studymerged$covfile_prepared)[study.index]
          NewPheno <- as.character(studymerged$phenolist)[study.index]
		  covlist_prepared = unlist(strsplit(as.character(studymerged$covs_separated)[study.index],split=" "))

		if (arg$variable == "binary") {
          FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=TRUE, cov_header=TRUE)
          sample_size <- nrow(FAM)
          f = paste(paste("FAM$",NewPheno, sep=""), paste(names(FAM[covlist_prepared]),collapse="+FAM$"),sep="~FAM$")
          obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="D")
          Generate_Meta_Files(obj, File.Bed, File.Bim, File.SetID, File.MSSD, File.MInfo, sample_size, impute.method=arg$impute_method)
          File.MSSD.vec[study.index] <- File.MSSD
          File.MInfo.vec[study.index] <- File.MInfo

          } else if (arg$variable == "continuous") {
          FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=FALSE, cov_header=TRUE)
          sample_size <- nrow(FAM)
          f = paste(paste("FAM$",NewPheno, sep=""), paste(names(FAM[covlist_prepared]),collapse="+FAM$"),sep="~FAM$")
          obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="C")
          Generate_Meta_Files(obj, File.Bed, File.Bim, File.SetID, File.MSSD, File.MInfo, sample_size, impute.method=arg$impute_method)
          File.MSSD.vec[study.index] <- File.MSSD
          File.MInfo.vec[study.index] <- File.MInfo

          }

     } else if (as.character(studymerged$studygenoOrdose)[study.index] == "dosage") {
          File.Fam <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".fam")
          File.MSSD <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".mssd") 
          File.MInfo <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".minfo")
          File.SetID <- as.character(studymerged$setIDnames)[study.index]
          File.Dosage <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".dosage") #this dosage file should be created
		File.Cov <- as.character(studymerged$covfile_prepared)[study.index]
          NewPheno <- as.character(studymerged$phenolist)[study.index]
		  covlist_prepared = unlist(strsplit(as.character(studymerged$covs_separated)[study.index],split=" "))

          if (arg$variable == "binary") {
          FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=TRUE, cov_header=TRUE)
          sample_size <- nrow(FAM)
          f = paste(paste("FAM$",NewPheno, sep=""), paste(names(FAM[covlist_prepared]),collapse="+FAM$"),sep="~FAM$")
          obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="D")
          Generate_Meta_Files_FromDosage(obj, File.Dosage, File.SetID, File.MSSD, File.MInfo, sample_size, impute.method=arg$impute_method)
          File.MSSD.vec[study.index] <- File.MSSD
          File.MInfo.vec[study.index] <- File.MInfo

          } else if (arg$variable == "continuous") {
          FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=FALSE, cov_header=TRUE)
          sample_size <- nrow(FAM)
          f = paste(paste("FAM$",NewPheno, sep=""), paste(names(FAM[covlist_prepared]),collapse="+FAM$"),sep="~FAM$")
          obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="C")
          Generate_Meta_Files_FromDosage(obj, File.Dosage, File.SetID, File.MSSD, File.MInfo, sample_size, impute.method=arg$impute_method)
          File.MSSD.vec[study.index] <- File.MSSD
          File.MInfo.vec[study.index] <- File.MInfo

          } } }

    cohort_info <- Open_MSSD_File_2Read (File.MSSD.vec, File.MInfo.vec)
    output_skato <- MetaSKAT_MSSD_ALL(cohort_info, is.separate = TRUE, r.corr = arg$rcorr, method=arg$method, missing_cutoff=arg$missing_cutoff)
    write.table(output_skato, paste0(arg$out, ".tsv"), col.names=TRUE, quote=F, row.names=F, sep="\t")

} else if (is.na(arg$covnames)) {

studymerged <- as.data.frame(cbind(studyprefixnames,setIDnames,studygenoOrdose,phenolist, covfile_prepared))

File.MSSD.vec <- rep("", nstudyprefixnames)
File.MInfo.vec <- rep("", nstudyprefixnames)

    for (study.index in 1:nstudyprefixnames) {
     studyprefix <- as.character(studymerged$studyprefixnames)[study.index]
     if(as.character(studymerged$studygenoOrdose)[study.index] == "genotype"){
          File.Bed <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".bed")
          File.Bim <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".bim")
          File.Fam <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".fam")
          File.MSSD <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".mssd")
          File.MInfo <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".minfo")
          File.SetID <- as.character(studymerged$setIDnames)[study.index]
		File.Cov <- as.character(studymerged$covfile_prepared)[study.index]
          NewPheno <- as.character(studymerged$phenolist)[study.index]

		if (arg$variable == "binary") {
          FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=TRUE, cov_header=TRUE)
          sample_size <- nrow(FAM)
          f = paste0("FAM$", NewPheno, " ~ 1")
          obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="D")
          Generate_Meta_Files(obj, File.Bed, File.Bim, File.SetID, File.MSSD, File.MInfo, sample_size, impute.method=arg$impute_method)
          File.MSSD.vec[study.index] <- File.MSSD
          File.MInfo.vec[study.index] <- File.MInfo

          } else if (arg$variable == "continuous") {
          FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=FALSE, cov_header=TRUE)
          sample_size <- nrow(FAM)
          f = paste0("FAM$", NewPheno, " ~ 1")
          obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="C")
          Generate_Meta_Files(obj, File.Bed, File.Bim, File.SetID, File.MSSD, File.MInfo, sample_size, impute.method=arg$impute_method)
          File.MSSD.vec[study.index] <- File.MSSD
          File.MInfo.vec[study.index] <- File.MInfo

          }

     } else if (as.character(studymerged$studygenoOrdose)[study.index] == "dosage") {
          File.Fam <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".fam")
          File.MSSD <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".mssd") 
          File.MInfo <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".minfo")
          File.SetID <- as.character(studymerged$setIDnames)[study.index]
          File.Dosage <- paste0(as.character(studymerged$studyprefixnames)[study.index], ".dosage") #this dosage file should be created
		File.Cov <- as.character(studymerged$covfile_prepared)[study.index]
          NewPheno <- as.character(studymerged$phenolist)[study.index]

          if (arg$variable == "binary") {
          FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=TRUE, cov_header=TRUE)
          sample_size <- nrow(FAM)
          f = paste0("FAM$", NewPheno, " ~ 1")
          obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="D")
          Generate_Meta_Files_FromDosage(obj, File.Dosage, File.SetID, File.MSSD, File.MInfo, sample_size, impute.method=arg$impute_method)
          File.MSSD.vec[study.index] <- File.MSSD
          File.MInfo.vec[study.index] <- File.MInfo

          } else if (arg$variable == "continuous") {
          FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=FALSE, cov_header=TRUE)
          sample_size <- nrow(FAM)
          f = paste0("FAM$", NewPheno, " ~ 1")
          obj <- SKAT_Null_Model(as.formula(f), data=FAM, out_type="C")
          Generate_Meta_Files_FromDosage(obj, File.Dosage, File.SetID, File.MSSD, File.MInfo, sample_size, impute.method=arg$impute_method)
          File.MSSD.vec[study.index] <- File.MSSD
          File.MInfo.vec[study.index] <- File.MInfo

          } } }

    cohort_info <- Open_MSSD_File_2Read (File.MSSD.vec, File.MInfo.vec)
    output_skato <- MetaSKAT_MSSD_ALL(cohort_info, is.separate = TRUE, r.corr = arg$rcorr, method=arg$method, missing_cutoff=arg$missing_cutoff)
    write.table(output_skato, paste0(arg$out, ".tsv"), col.names=TRUE, quote=F, row.names=F, sep="\t")

}