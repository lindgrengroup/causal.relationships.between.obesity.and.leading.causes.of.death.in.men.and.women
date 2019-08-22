#$/bin/env Rscript
#!-cwd

#Read in all SNP-groups
snp_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig",
                "whradjbmi.eur.comb.pulit.sig",
                "bmi.eur.men.pulit", "bmi.eur.men.pulit.sig", "bmi.eur.men.pulit.phet",
                "whr.eur.men.pulit", "whr.eur.men.pulit.sig", "whr.eur.men.pulit.phet",
                "whradjbmi.eur.men.pulit", "whradjbmi.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.phet",
                "bmi.eur.women.pulit", "bmi.eur.women.pulit.sig", "bmi.eur.women.pulit.phet",
                "whr.eur.women.pulit", "whr.eur.women.pulit.sig", "whr.eur.women.pulit.phet",
                "whradjbmi.eur.women.pulit", "whradjbmi.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.phet",
                "bmi.eur.comb.pulit.winner", "whr.eur.comb.pulit.winner",
                "whradjbmi.eur.comb.pulit.winner",
                "bmi.eur.men.pulit.winner", "whr.eur.men.pulit.winner",
                "whradjbmi.eur.men.pulit.winner",
                "bmi.eur.women.pulit.winner", "whr.eur.women.pulit.winner",
                "whradjbmi.eur.women.pulit.winner")

#Read in the the list of SNPs that are tri-allelic and pass QC and are correlated
duplicated <- read.table("/well/lindgren/jc/ukbb/genotype.qc.snplists.180507/non.biallelic.snps.txt", 
		stringsAsFactors = F, header = F)
qc <- read.table("/well/lindgren/jc/ukbb/genotype.qc.snplists.180507/snps_that_pass_all_qc_180508.txt", 
		stringsAsFactors = F, header = F)
qc[qc$V1 == "9:140079779_T_C_C_T", ] <- "rs144926207_C_T"

ld <- read.table("../snps.that.are.correlated.over.0.05.181123.txt", stringsAsFactors = F, header = T)

for (snp_group in snp_groups) {
	input_file <- paste("../index.files.with.trialleles/", snp_group, ".180513.txt", sep = "")
	df <- read.table(input_file, stringsAsFactors = F, header = T)

	print(snp_group)
        print(df[(df$SNP %in% duplicated$V1), ])
	df <- df[!(df$SNP %in% duplicated$V1), ]
	
	print(df[!(df$SNP %in% qc$V1), ])
	df <- df[df$SNP %in% qc$V1, ]

	#Only remove from files if have both the correlated variants. The "remove" SNP is the one with the
	#highest P-value, either combined if in one of the sig/phet/winner SNP-groups, or sex-specific P-values if
	#in one of the index files for a sex. Note, only works for the current set of correlated SNPs
	for (i in 1:nrow(ld)) {
		if (nrow(df[df$SNP %in% ld[i, c("SNP_A", "SNP_B")], ]) >= 2) {
			df <- df[!(df$SNP == ld[i, "remove"]), ]
			print(ld[i, "remove"])
		}
	}
	
	print(nrow(df))

	new_output_file <- paste("../", snp_group, ".180513.txt", sep = "")
	write.table(df, new_output_file, quote = F, row.names = F, sep = " ")

}

#Make unweighted GRSs for the winner SNPs:
winner_groups <- c("bmi.eur.comb.pulit.winner", "whr.eur.comb.pulit.winner",
                "whradjbmi.eur.comb.pulit.winner")

for (winner_group in winner_groups) {
	df <- read.table(paste0("../", winner_group, ".180513.txt"), stringsAsFactors = F, header = T)
	df[, grep("^beta\\.", colnames(df))] <- 1
	write.table(df, paste0("../", winner_group, "_unweighted.180513.txt"), quote = F, row.names = F, sep = " ")
}

sex_het_groups <- c("bmi.pulit", "whr.pulit", "whradjbmi.pulit")
het_or_not_groups <- c("not.sex.het", "sex.het")

for (sex_het_group in sex_het_groups) {
	for (het_or_not_group in het_or_not_groups) {
		input_file <- paste("../sex.het.enrichment/", het_or_not_group, ".snps.", sex_het_group, ".180921.txt", sep = "")
		df <- read.table(input_file, stringsAsFactors = F, header = F)
		
		df$simple <- gsub("_[A-Z]+_[A-Z]+$", "", df$V1)

		print(sex_het_group)
		print(het_or_not_group)
		print(df[df$V1 %in% duplicated$V1, ])
		print(nrow(df))
		df <- df[!(df$V1 %in% duplicated$V1), ]
		print(nrow(df))

		print(df[!(df$V1 %in% qc$V1), ])
        	df <- df[df$V1 %in% qc$V1, ]

		#Since all of these are the combined SNP-groups. Only works for the current set of
		#Correlated SNPs
		if (sex_het_group == "bmi.pulit") {
			print(df[(df$V1 %in% ld[ld$dataset == "pulit" & ld$trait == "bmi" & ld$sex_group == "comb", "remove"]), ])
			df <- df[!(df$V1 %in% ld[ld$dataset == "pulit" & ld$trait == "bmi" & ld$sex_group == "comb", "remove"]), ]		
		} else if (sex_het_group == "whradjbmi.pulit") {
			print(df[(df$V1 %in% ld[ld$dataset == "pulit" & ld$trait == "whradjbmi" & ld$sex_group == "comb", "remove"]), ])
			df <- df[!(df$V1 %in% ld[ld$dataset == "pulit" & ld$trait == "whradjbmi" & ld$sex_group == "comb", "remove"]), ]
		}

		df$simple <- NULL
		print(nrow(df))
		new_output_file <- paste("../sex.het.enrichment/", het_or_not_group, ".snps.", sex_het_group, ".no.tri.180921.txt", sep = "")
		write.table(df, new_output_file, row.names = F, col.names = F, quote = F, sep = " ")
	}
}
