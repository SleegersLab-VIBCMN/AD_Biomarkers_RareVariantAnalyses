#!~/miniconda3/bin/Rscript
## Annotate the results of gene-based (rare variant) association meta-analyses as implemented in MetaSKAT package using PLINK format files
## Written by Fahri Kucukali for the analyses conducted for Kucukali et al., Whole‐exome rare‐variant analysis of Alzheimer's disease and related biomarker traits, Alzheimer's & Dementia (2022). http://doi.org/10.1002/alz.12842
## For issues: https://github.com/SleegersLab-VIBCMN/AD_Biomarkers_RareVariantAnalyses or e-mail Fahri.Kucukali@uantwerpen.vib.be

suppressWarnings(suppressMessages(library(MetaSKAT)))
suppressWarnings(suppressMessages(library(SKAT)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(ggbeeswarm)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(qqman)))
suppressWarnings(suppressMessages(library(rmeta)))

p <- arg_parser("Annotate MetaSKAT results, with different modes and options. Publication Version 1.1, December 2022")
p <- add_argument(p, "--resultfile1", help = "Annotated SKAT File 1, from PLINKSKAT.R", type = "character", nargs = '+')
p <- add_argument(p, "--resultfile2", help = "Annotated SKAT File 2, from PLINKSKAT.R", type = "character", nargs = '+')
p <- add_argument(p, "--addNamefile1", help = "Prefix to be added to the columns of File 1.", type = "character", nargs = '+')
p <- add_argument(p, "--addNamefile2", help = "Prefix to be added to the columns of File 2.", type = "character", nargs = '+')
p <- add_argument(p, "--metaresult", help = "MetaSKAT Result file (from PLINKMetaSKAT.R). Only commonly tested genes (tested in both datasets) will be kept in the final annotation file.", type = "character", nargs = '+')
p <- add_argument(p, "--geneAnnotFile", help = "An annotation file for tested genes, with header, first column is gene names as supplied in SetID. Rest of the columns can be any information, but chr (without chr prefix) and start columns are required to generate Manhattan plots.", type="character", nargs = '+')
p <- add_argument(p, "--threshold", help = "Thresholds for plots and tables created for full analysis, or the given gene SetID; also used for labelling the manhattan plots.", type = "numeric", nargs = '+')
p <- add_argument(p, "--testtype", help = "Only for annotating plots, enter whether you tested protein-altering or LoF variants. Optional.", type = "character", nargs = '*', default = "")
p <- add_argument(p, "--phenotype", help = "Only for annotating plots, enter the tested phenotype name. Optional.", type="character", nargs = '*', default = "")

arg <- parse_args(p)

if (is.na(arg$resultfile1) || is.na(arg$resultfile2) || is.na(arg$addNamefile1) || is.na(arg$addNamefile2) || is.na(arg$metaresult) || is.na(arg$geneAnnotFile)) {
    message("Error: exiting because at least one of the following required arguments is missing: --resultfile1, --resultfile2, --addNamefile1, --addNamefile2, --metaresult , --geneAnnotFile")
    q("no")
}

resultfile1 <- fread(arg$resultfile1, header=TRUE)

resultfile2 <- fread(arg$resultfile2, header=TRUE)

metaresult <- fread(arg$metaresult, header=TRUE)

addNamefile1 <- arg$addNamefile1

addNamefile2 <- arg$addNamefile2

geneAnnotFile <- fread(arg$geneAnnotFile, header=TRUE)

colnames(geneAnnotFile)[1] <- "SetID"

sum_resultfile1 <- as.data.table(cbind.data.frame(resultfile1[,c(1:20)], resultfile1[, exome_wide_check]))
sum_resultfile2 <- as.data.table(cbind.data.frame(resultfile2[,c(1:20)], resultfile2[, exome_wide_check]))

colnames(sum_resultfile1)[ncol(sum_resultfile1)] <- "exome_wide_check"
colnames(sum_resultfile2)[ncol(sum_resultfile2)] <- "exome_wide_check"

gene_list <- unique(rbind(sum_resultfile1,sum_resultfile2)[,1])

colnames(sum_resultfile1) <- paste(addNamefile1, colnames(sum_resultfile1), sep = "_")

colnames(sum_resultfile2) <- paste(addNamefile2, colnames(sum_resultfile2), sep = "_")

colnames(sum_resultfile1)[1] <- "SetID"

colnames(sum_resultfile2)[1] <- "SetID"

metaresult_selected <- subset(metaresult, SetID %in% gene_list$SetID)

stats_merged1 <- left_join(metaresult_selected,sum_resultfile1,by="SetID")

stats_merged2 <- left_join(stats_merged1,sum_resultfile2,by="SetID")

annotated_stats_merged2 <- left_join(stats_merged2,geneAnnotFile,by="SetID")

colnames(annotated_stats_merged2)[2] <- "MetaSKATO_P"

noNA_annotated_stats_merged2 <- annotated_stats_merged2[!is.na(annotated_stats_merged2$MetaSKATO_P),]

noNA_annotated_stats_merged2$MetaSKATO_P <- ifelse(noNA_annotated_stats_merged2$MetaSKATO_P >= 1, 1, noNA_annotated_stats_merged2$MetaSKATO_P)

exome_wide_threshold <- 0.05/(nrow(noNA_annotated_stats_merged2))

suggestive_threshold <- 1/(nrow(noNA_annotated_stats_merged2))

noNA_annotated_stats_merged2$METASKATO_exome_wide_check <- ifelse(noNA_annotated_stats_merged2$MetaSKATO_P <= exome_wide_threshold, "exomeWideHit", ifelse(noNA_annotated_stats_merged2$MetaSKATO_P > exome_wide_threshold & noNA_annotated_stats_merged2$MetaSKATO_P <= suggestive_threshold , "suggestiveHit", NA))

final_noNA_annotated_stats_merged2 <- as.data.table(noNA_annotated_stats_merged2 %>% relocate(METASKATO_exome_wide_check, .after=MetaSKATO_P))

commonlyfound <- final_noNA_annotated_stats_merged2[(!is.na(final_noNA_annotated_stats_merged2[[4]]) & !is.na(final_noNA_annotated_stats_merged2[[24]]))]

exome_wide_threshold <- 0.05/(nrow(commonlyfound))

suggestive_threshold <- 1/(nrow(commonlyfound))

commonlyfound$METASKATO_exome_wide_check <- ifelse(commonlyfound$MetaSKATO_P <= exome_wide_threshold, "exomeWideHit", ifelse(commonlyfound$MetaSKATO_P > exome_wide_threshold & commonlyfound$MetaSKATO_P <= suggestive_threshold , "suggestiveHit", NA))

if (colnames(commonlyfound)[18] == paste(addNamefile1, "OR_noOfMutations", sep = "_") & colnames(commonlyfound)[38] == paste(addNamefile2, "OR_noOfMutations", sep = "_")) {

	output_prepareRMETA <- as.data.frame(matrix(0, 1, 15))

	colnames(output_prepareRMETA) <- c("SetID","RE_ORsum","RE_ORsum_L95","RE_ORsum_U95","RE_ORsum_SE", "RE_tau2", "RE_het", "RE_hetP", "FE_ORsum","FE_ORsum_L95","FE_ORsum_U95","FE_ORsum_SE","FE_tau2", "FE_het", "FE_hetP")

	for (i in commonlyfound$SetID){

	table <- commonlyfound[commonlyfound$SetID == i,]

	random <- meta.summaries(c(as.numeric(table[1,21]),as.numeric(table[1,41])), c(as.numeric(table[1,22]),as.numeric(table[1,42])), names=c(addNamefile1,addNamefile2),method="random", logscale=TRUE)

	random$confL <- exp((random$summary)-(1.96*(random$se.summary)))

	random$confU <- exp((random$summary)+(1.96*(random$se.summary)))

	fixed <- meta.summaries(c(as.numeric(table[1,21]),as.numeric(table[1,40])), c(as.numeric(table[1,22]),as.numeric(table[1,42])), names=c(addNamefile1,addNamefile2),method="fixed", logscale=TRUE)

	fixed$confL <- exp((fixed$summary)-(1.96*(fixed$se.summary)))

	fixed$confU <- exp((fixed$summary)+(1.96*(fixed$se.summary)))
	
	SetID <- i

	selectRandomFix <- cbind.data.frame(SetID,round(exp(random$summary), digits = 3),round(random$confL, digits = 3),round(random$confU, digits = 3),round(random$se.summary, digits = 3),random$tau2,random$het[1],random$het[3],round(exp(fixed$summary), digits = 3),round(fixed$confL, digits = 3),round(fixed$confU, digits = 3),round(fixed$se.summary, digits = 3),fixed$tau2,fixed$het[1],fixed$het[3])
	
	colnames(selectRandomFix) <- c("SetID","RE_ORsum","RE_ORsum_L95","RE_ORsum_U95","RE_ORsum_SE", "RE_tau2", "RE_het", "RE_hetP", "FE_ORsum","FE_ORsum_L95","FE_ORsum_U95","FE_ORsum_SE","FE_tau2", "FE_het", "FE_hetP")

	output_prepareRMETA <- rbind.data.frame(output_prepareRMETA, selectRandomFix)

	}

} else if(colnames(commonlyfound)[18] == paste(addNamefile1, "BETA_noOfMutations", sep = "_") & colnames(commonlyfound)[38] == paste(addNamefile2, "BETA_noOfMutations", sep = "_")) {

	output_prepareRMETA <- as.data.frame(matrix(0, 1, 15))

	colnames(output_prepareRMETA) <- c("SetID","RE_BETAsum","RE_BETAsum_L95","RE_BETAsum_U95","RE_BETAsum_SE", "RE_tau2", "RE_het", "RE_hetP", "FE_BETAsum","FE_BETAsum_L95","FE_BETAsum_U95","FE_BETAsum_SE","FE_tau2", "FE_het", "FE_hetP")

	for (i in commonlyfound$SetID){

	table <- commonlyfound[commonlyfound$SetID == i,]

	random <- meta.summaries(c(as.numeric(table[1,21]),as.numeric(table[1,41])), c(as.numeric(table[1,22]),as.numeric(table[1,42])), names=c(addNamefile1,addNamefile2),method="random", logscale=FALSE)

	random$confL <- (random$summary)-(1.96*(random$se.summary))

	random$confU <- (random$summary)+(1.96*(random$se.summary))

	fixed <- meta.summaries(c(as.numeric(table[1,21]),as.numeric(table[1,41])), c(as.numeric(table[1,22]),as.numeric(table[1,42])), names=c(addNamefile1,addNamefile2),method="fixed", logscale=FALSE)

	fixed$confL <- (fixed$summary)-(1.96*(fixed$se.summary))

	fixed$confU <- (fixed$summary)+(1.96*(fixed$se.summary))

	SetID <- i

	selectRandomFix <- cbind.data.frame(SetID,round(random$summary, digits = 3),round(random$confL, digits = 3),round(random$confU, digits = 3),round(random$se.summary, digits = 3),random$tau2,random$het[1],random$het[3],round(fixed$summary, digits = 3),round(fixed$confL, digits = 3),round(fixed$confU, digits = 3),round(fixed$se.summary, digits = 3),fixed$tau2,fixed$het[1],fixed$het[3])
	
	colnames(selectRandomFix) <- c("SetID","RE_BETAsum","RE_BETAsum_L95","RE_BETAsum_U95","RE_BETAsum_SE", "RE_tau2", "RE_het", "RE_hetP", "FE_BETAsum","FE_BETAsum_L95","FE_BETAsum_U95","FE_BETAsum_SE","FE_tau2", "FE_het", "FE_hetP")

	output_prepareRMETA <- rbind.data.frame(output_prepareRMETA, selectRandomFix)
	
	}
	
} else {print("Check your header names, there are not as expected, there might be something wrong.")}

	output_readyRMETA <- output_prepareRMETA[-1,]

part1 <- commonlyfound[,c(1:3)]

part2 <- commonlyfound[,-c(2:3)]

part1_RMETA <- left_join(part1,output_readyRMETA,by="SetID")

commonlyfound2 <- left_join(part1_RMETA,part2, by="SetID")

write.table(commonlyfound2[order(MetaSKATO_P),],paste0("ANNOTATED.",arg$metaresult), col.names=TRUE, quote=F, row.names=F, sep="\t")

	maxP = min(commonlyfound2$MetaSKATO_P)
	b = -log10(maxP)
	ymax = b + 2

	nCHR <- length(unique(commonlyfound2$chr))
	commonlyfound2$StartPoscum <- 0
	s <- 0
	nbp <- c()
	for (i in unique(sort(commonlyfound2$chr))){
	  nbp[i] <- max(commonlyfound2[commonlyfound2$chr == i,]$start)
	  commonlyfound2[commonlyfound2$chr == i,"StartPoscum"] <- commonlyfound2[commonlyfound2$chr == i,"start"] + s
	  s <- s + nbp[i]
	}

	axis.set <- commonlyfound2 %>% 
	  group_by(chr) %>% 
	  summarize(center = (max(StartPoscum) + min(StartPoscum)) / 2)

	is.odd <- function(x) x %% 2 != 0

	commonlyfound2$chrcolor <- ifelse(is.odd(commonlyfound2$chr), 'royalblue3', 'firebrick2')

	commonlyfound2LABEL <- commonlyfound2[commonlyfound2$MetaSKATO_P <= arg$threshold | !is.na(commonlyfound2$METASKATO_exome_wide_check),]

	phenotypename <- arg$phenotype
	
	testtype <- arg$testtype
		
	png(paste0("manhattanMeta.",arg$metaresult,".png"), height=10, width=15, res=600, units="in")
	q1 <- ggplot() + geom_point(data = commonlyfound2, aes(x = StartPoscum, y = -log10(MetaSKATO_P)), color=commonlyfound2$chrcolor) + scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) + scale_y_continuous(limits = c(0,ymax),breaks = seq(0,ymax, by=1)) + geom_hline(yintercept=-log10(exome_wide_threshold), colour = "red", size = 0.5, linetype="dashed") + geom_hline(yintercept=-log10(suggestive_threshold), colour = "blue4", size = 0.25, linetype="dashed") + labs(x="Chromosomal Position",y=expression(-log[10](P))) + geom_text_repel(data=commonlyfound2LABEL, aes(x = StartPoscum, y = -log10(MetaSKATO_P), label = SetID), box.padding = 1, point.padding = 0.4, size = 4.5, segment.color = "black", min.segment.length = unit(0, 'lines'), nudge_y = .25, max.overlaps = getOption("ggrepel.max.overlaps", default = 1000)) + theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(face="bold", size=15), axis.text.y = element_text(face="bold", size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + ggtitle(paste0(phenotypename," ",testtype," MetaSKAT-O Manhattan Plot")) + theme(plot.title = element_text(hjust = 0.5, size=23))
	print(q1)
	dev.off()

	pvalues <- c(t(commonlyfound2[,MetaSKATO_P]))
	qs <- qchisq(pvalues, df=1, lower.tail = TRUE)

	nonaqs <- qs[!is.na(qs)]

	nonapvalues <- pvalues[!is.na(pvalues)]

	lambda <- round(((median(nonaqs))/0.4549364),digits=3)
	png(paste0("qqplotMeta.",arg$metaresult,".png"), width = 6, height = 6, res=600, units="in")
	qq(nonapvalues, main = paste0(phenotypename," ",testtype,"\nMetaSKAT-O QQ Plot\nlambda = ",lambda))
	dev.off()

