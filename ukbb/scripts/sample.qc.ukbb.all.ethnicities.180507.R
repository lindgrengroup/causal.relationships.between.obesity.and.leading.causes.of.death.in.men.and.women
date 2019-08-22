#!/bin/env Rscript
#$ -cwd

#####################################################################
############ Sample QC in UKBB, Jenny Censin, 2018-05-03 ###########
################### Includes all ethnicities #######################

merge.self.reported.ethnicity <- function(data, qc.counts.file) {

	#Check for mismatches between assessments and sort as Europeans or not Europeans
	#Self-reported ethnicity fields: f.21000-0.0; f.21000-1.0; f.21000-2.0 
	# 1	White 				# 1001	British 		# 1002	Irish				# 1003	Any other white background
	# 2	Mixed				# 2001	White/Black Caribbean	# 2002	White and Black African		# 2003	White and Asian			# 2004	Any other mixed background
	# 3	Asian or Asian British		# 3001	Indian			# 3002	Pakistani			# 3003	Bangladeshi			# 3004	Any other Asian background
	# 4	Black or Black British		# 4001	Caribbean		# 4002	African				# 4003	Any other Black background
	# 5	Chinese
	# 6	Other ethnic group
	# -1	Do not know			# -3	Prefer not to answer	#-10 if to exclude
	
	#Recode NA as 0 and make sure is integer vector
	se = c("se.0", "se.1", "se.2")
	data[, se][is.na(data[, se])] = 0
	data[se] <- lapply(data[se], as.integer)	

	#Make vectors of the different codes
	no.value <- c(-1, -3, 0)
	brit = c(1, 1001:1002)
	white = c(1, 1001:1003)
	mixed = c(2, 6, 2001:2004)
	asian = c(3, 5, 3001:3004)
	black = c(4, 4001:4003)
	
	#Start
	data$check.se <- apply(data[se], 1, function(x) {ifelse(any(x == -10), "exclude",
	ifelse(all(x %in% no.value), "no.value",  
	ifelse(any(x %in% white) & any(x %in% c(mixed, asian, black)), "changed.se", 
	ifelse(any(x %in% mixed) & any(x %in% c(white, asian, black)), "changed.se", 
	ifelse(any(x %in% asian) & any(x %in% c(white, mixed, black)), "changed.se", 
	ifelse(any(x %in% black) & any(x %in% c(white, mixed, asian)), "changed.se", 
	ifelse(any(x %in% brit), "brit", 
	ifelse(any(x %in% white), "white", 
	ifelse(any(x %in% mixed), "mixed", 
	ifelse(any(x %in% asian), "asian", 
	ifelse(any(x %in% black), "black", "rest")))))))))))})

	#Recode those that self-report as British but aren't that genetically to white
	data$check.se[data$check.se == "brit" & data$in.white.British.ancestry.subset != 1] <- "recode.white"
	
	#Make the cleaned dataset: 
	cleaned <- subset(data, !(data$check.se %in% c("exclude")))

	#Add to QC counts file: 
        sink(qc.counts.file, append = T)
        cat(paste("**FILTER** Excluded because -10 in self-reported ethnicity: ",
                        length(which(data$check.se == "exclude")), "\n",
		"**FILTER** Has no self-reported ethnicity: ", 
			length(which(data$check.se == "no.value")), "\n",
                "**FILTER** Reported different ethnicities: ",
                        length(which(data$check.se == "changed.se")), "\n",
                "**FILTER** Included, british: ",
                        length(which(data$check.se == "brit")), "\n",
                "**FILTER** Included, white: ", 
			length(which(data$check.se == "white")), "\n", 
		"**FILTER** Included, mismatch between reported white/brit/irish and genetically brit/irish, coded recode.white: ",
                        length(which(data$check.se == "recode.white")), "\n",
		"**FILTER** Included, mixed ethnicities: ", 
			length(which(data$check.se == "mixed")), "\n", 
		"**FILTER** Included, asian: ", 
			length(which(data$check.se == "asian")), "\n", 
		"**FILTER** Included, black: ", 
			length(which(data$check.se == "black")), "\n", 
		"Number of samples which didn't get coded (check if not 0): ", 
			length(which(data$check.se == "rest")), "\n", 
		"Total remaining: ", nrow(cleaned), "\n", sep = ""))
	sink()

        return(cleaned)

}
	
qc.withdraw.ind <- function(data, qc.counts.file) {
		
	#Pathway to UKBB provided list of individuals that have withdrawn consent
	withdraw = read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20181016.csv", 
			header=F)
	
	cleaned = subset(data, !(data$ID %in% withdraw$V1))

        sink(qc.counts.file, append = T)
        cat(paste("**FILTER** Individuals that withdrew consent: ", length(which(withdraw[,1] %in% data$ID)),
                "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
        sink()

	return(cleaned)
}

qc.negative.ids <- function(data, qc.counts.file) {
	
	cleaned = subset(data, data$ID > 0)

        sink(qc.counts.file, append = T)
        cat(paste("**FILTER** Individuals with negative IDs (withdrawn consent): ", length(which(!(data$ID %in% cleaned$ID))), 
		"; REMAINING: ", nrow(cleaned), "\n", sep = ""))
        sink()

        return(cleaned)
}

qc.het.miss <- function(data, qc.counts.file) {

	cleaned = subset(data, !is.na(data$het.missing.outliers) & data$het.missing.outliers != 1)	

	sink(qc.counts.file, append = T)
	cat(paste("**FILTER** Poor heterozygosity or missingness: ", length(which(!(data$ID %in% cleaned$ID))), 
		"; REMAINING: ", nrow(cleaned), "\n", sep = ""))
	sink()

	return(cleaned)
}

qc.excess.related <- function(data, qc.counts.file) {

	cleaned = subset(data, !is.na(data$excess.relatives) & data$excess.relatives != 1)

	sink(qc.counts.file, append = T)
	cat(paste("**FILTER** Excess relatives (>10 3rd degree relatives): ", length(which(!(data$ID %in% cleaned$ID))), 
		"; REMAINING: ", nrow(cleaned) , "\n", sep = ""))
	sink()

	return(cleaned)	
}	

qc.not.in.phasing  <- function(data, qc.counts.file) {

	cleaned = subset(data, !is.na(data$in.Phasing.Input.chr1_22) & data$in.Phasing.Input.chr1_22 != 0)

	sink(qc.counts.file, append = T)
	cat(paste("**FILTER** Not used in autosome phasing: ", length(which(!(data$ID %in% cleaned$ID))), 
		"; REMAINING: ", nrow(cleaned), "\n", sep = ""))
	sink()

	return(cleaned)	
}		
	
qc.sex.chr.aneupl  <- function(data, qc.counts.file) {
	
	cleaned = subset(data, !is.na(data$putative.sex.chromosome.aneuploid) & data$putative.sex.chromosome.aneuploid != 1)
	
	sink(qc.counts.file, append = T)
	cat(paste("**FILTER** Putative sex chr aneuploidy: ", length(which(!(data$ID %in% cleaned$ID))), 
		"; REMAINING: ", nrow(cleaned), "\n", sep = ""))
	sink()

	return(cleaned)	
}	

qc.sex.mismatch <- function(data, qc.counts.file) {

	cleaned = subset(data, !is.na(data$Submitted.Gender) & !is.na(data$Inferred.Gender) & data$Submitted.Gender == data$Inferred.Gender)
 
	sink(qc.counts.file, append=T)
	cat(paste("**FILTER** Sex mismatch ", length(which(!(data$ID %in% cleaned$ID))), 
		"; REMAINING: ", nrow(cleaned), "\n", sep=""))
	sink()

	return(cleaned)	
}	

qc.related <- function(data, qc.counts.file) { 
	
	#Pathway to UKBB list of related individuals
	related <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/ukb1186_rel_s488366.dat", header=T)
	
	#For each pair of related individuals, remove the samples with the highest missingness
	related <- related[related$Kinship > 0.0884 & related$ID1 %in% data$ID & related$ID2 %in% data$ID, ]
	
	related$miss1 = data$sample.qc.missing.rate[match(related$ID1, data$ID)]
	related$miss2 = data$sample.qc.missing.rate[match(related$ID2, data$ID)]
	related$max.miss <- pmax(related$miss1, related$miss2)

	#Make a list
	related$id.remove <- ifelse(is.na(related$miss1) & is.na(related$miss2), related$ID2, 
				ifelse(is.na(related$miss1), related$ID1, 
				ifelse(is.na(related$miss2), related$ID2, 
				ifelse(related$miss1 == related$max.miss, related$ID1, 
				ifelse(related$miss2 == related$max.miss, related$ID2, 
				"problem")))))

	cleaned <- subset(data, !(data$ID %in% related$id.remove))	

	sink(qc.counts.file, append = T)
	cat(paste("**FILTER** Relatedness pairs with problems - look into if not 0: ", 
			length(which(related$id.remove == "problem")), "\n", 
		"**FILTER** Individuals excluded because of relatedness: ", 
			nrow(data[data$ID %in% related$id.remove, ]), "\n",
		"REMAINING NOT RELATED: ", nrow(cleaned), "\n\n", sep = "")) 
	sink()

	return(cleaned)

}

qc.kinship.table <- function(data, qc.counts.file) {

	cleaned <- subset(data, !is.na(data$excluded.from.kinship.inference) & data$excluded.from.kinship.inference == 0)

	sink(qc.counts.file, append = T)
	cat(paste("**FILTER** Excluded from kinship inference: ", length(which(!(data$ID %in% cleaned$ID))), 
		" ; REMAINING: ", nrow(cleaned), "\n", sep = ""))
	sink()
	
	return(cleaned)
}
	
#Make a file with all the details about the sample QC
qc.counts.file = "../ukbb.sample.qc.count.all.ancestries.181023.txt"

#Pathway to main phenotype file from UKBB
pheno = read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
	header = T, sep = ",", na.string = c("NA", "", "."), stringsAsFactors = F)
colnames.pheno <- colnames(pheno)
colnames(pheno) <- gsub("X", "f.", colnames.pheno)
colnames(pheno)[1] <- "f.eid"

sink(qc.counts.file)
cat(paste("SAMPLES IN PHENOTYPE FILE: ", nrow(pheno), "\n", sep = ""))
sink()

#Pathway to main UKBB provided QC file:
qc <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt", header = T, 
	na.string = c("NA", "", "."), stringsAsFactors = F)

sink(qc.counts.file, append = T)
cat(paste("SAMPLES IN QC FILE: ", nrow(qc), "\n\n", sep = ""))
sink()

#Pathway to the fam file corresponding to the QC file provided by the UKBB:
fam <- read.table("/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam", header = F)

#Add IDs to QC file: 
qc$ID <- fam[, 1]

#Add self-reported ethnicity to QC dataframe
qc$se.0 <- pheno$f.21000.0.0[match(qc$ID, pheno$f.eid)]
qc$se.1 <- pheno$f.21000.1.0[match(qc$ID, pheno$f.eid)]
qc$se.2 <- pheno$f.21000.2.0[match(qc$ID, pheno$f.eid)]

#Start running sample QC functions on data
cleaned <- merge.self.reported.ethnicity(qc, qc.counts.file)
cleaned <- qc.withdraw.ind(cleaned, qc.counts.file)
cleaned <- qc.negative.ids(cleaned, qc.counts.file)	
cleaned <- qc.het.miss(cleaned, qc.counts.file)
cleaned <- qc.excess.related(cleaned, qc.counts.file)
cleaned <- qc.not.in.phasing(cleaned, qc.counts.file)
cleaned <- qc.sex.chr.aneupl(cleaned, qc.counts.file)
cleaned <- qc.sex.mismatch(cleaned, qc.counts.file)
cleaned <- qc.related(cleaned, qc.counts.file)
cleaned <- qc.kinship.table(cleaned, qc.counts.file)

#Write file with samples to exclude
sample.qc.excluded.id <- subset(qc, !(qc$ID %in% cleaned$ID))
write.table(sample.qc.excluded.id, "../ukbb.excluded.ids.sample.qc.all.ancestries.180504.txt", quote = F, row.names = F)

#Total excluded
sink(qc.counts.file, append = T)
cat(paste("TOTAL REMOVED STANDARD QC: QC - CLEANED DATASET: ", nrow(qc) - nrow(cleaned), "\n", 
	"POST GENERAL SAMPLE QC CLEANED SAMPLE: ", nrow(cleaned), "\n", sep = ""))
sink()

#Extract columns of interest from phenotype file
pheno.columns <- c("f.eid", "f.34.0.0", "f.48.0.0", "f.49.0.0", "f.54.0.0", "f.3140.0.0", "f.21000.0.0", 
			"f.21000.1.0", "f.21000.2.0", "f.21001.0.0", "f.21003.0.0")
pheno.col.names <- c("ID", "year.of.birth", "wc", "hc", "assessment.centre", "pregnant", "self.ethnicity.0",
			"self.ethnicity.1", "self.ethnicity.2", "bmi", "age.assessment")

#Subset phenotype file to only include samples passing general QC
relevant.pheno <- pheno[ , pheno.columns]
colnames(relevant.pheno) <- pheno.col.names

#Make a data frame with the relevant data
pc <- paste("PC", 1:40, sep = "")
cleaned.columns <- c("ID", "genotyping.array", "Submitted.Gender", pc, "se.0", "se.1", "se.2", "check.se")
relevant.pheno <- merge(relevant.pheno, cleaned[, cleaned.columns], by.x = "ID", by.y = "ID")

#Make columns for whr, WHRadjBMI, age_squared  
relevant.pheno$whr <- as.numeric(relevant.pheno$wc)/as.numeric(relevant.pheno$hc)
relevant.pheno$age.squared <- relevant.pheno$age.assessment^2

#Write table with the relevant phenotypes for all the samples passing QC
write.table(relevant.pheno, "../ukbb.samples.passing.qc.relevant.pheno.all.ancestries.180504.txt", quote = F, row.names = F)

#Write a table with samples were even one in each pair of third degree relatives have been removed
relevant.pheno <- relevant.pheno[relevant.pheno$check.se %in% c("white", "brit", "recode.white"), ]

related <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/ukb1186_rel_s488366.dat", header=T)
related <- related[related$ID1 %in% relevant.pheno$ID & related$ID2 %in% relevant.pheno$ID, ]
relevant.pheno <- relevant.pheno[!(relevant.pheno$ID %in% related$ID1), ]

write.table(relevant.pheno[, c("ID", "ID")], "../ukbb.completely.unrelated.european.samples.190215.txt", 
	quote = F, row.names = F, sep = " ", col.names = F)
