#!/bin/env Rscript
#$-cwd

#This file updates to dbSNP 151
rs_aliases <- read.table("/well/lindgren/jc/resources/rsaliases.181010.txt",
                        stringsAsFactors = F, header = T, sep = " ")

#To clean the fasting glucose summary GWAS:
for (trait in c("FG", "FI")) {

	#According to Reedik at Magic, the chr:pos format SNPs are in build 37
        input_file <- paste("../full.summary.gwas/", trait, "_STAGE_GWAS_SEX.out.gz", sep = "")
        df <- read.table(gzfile(input_file), stringsAsFactors = F, header = T)
	df$chr <- gsub("rs[0-9]+|chr|:[0-9]+|_[0-9]+", "", df$rs_number)
	df$pos <- gsub("rs[0-9]+|chr[0-9]+:|chr[0-9]+_", "", df$rs_number)

	grch <- read.table("/well/lindgren/jc/resources/grch37.all.rsid.txt", stringsAsFactors=F, header =T, sep = " ")

	#Realized that this script might introduce duplicates. Have checked downstreams, and should be unaffected. 
	#No SNPs in the grch file is "" and "" for chr and position, so only SNPs in chr:pos format are updated
	#Since only merge those that have chr and position in the rs_number, the others are NA for the SNP
	df <- merge(df, grch, by = c("chr", "pos"), all.x = T)
	
	df[!is.na(df$snp), "rs_number"] <- df[!is.na(df$snp), "snp"]
	df$chr <- NULL
	df$pos <- NULL

	df <- merge(df, rs_aliases, by.x = "rs_number", by.y = "rsHigh", all.x = T)
	df[!is.na(df$rsCurrent), "rs_number"] <- df[!is.na(df$rsCurrent), "rsCurrent"]

	#Get rs123_A_C instead
	df[, c("alpha_first", "alpha_second")] <- t(apply(df[, c("reference_allele", "other_allele")], 1, sort))
	df$rs_number <- paste(df$rs_number, "_", df$alpha_first, "_", df$alpha_second, sep ="")
	df[, c("alpha_first", "alpha_second")] <- NULL

        df_comb <- df[df$beta != -9, c(1:3, 4:6, 10, 16)]
        df_men <- df[df$male_beta != -9, c(1:3, 18:20, 24, 26)]
        df_women <- df[df$female_beta != -9, c(1:3, 27:29, 33, 35)]
	colnames(df_men) <- gsub("^male_", "", colnames(df_men))
	colnames(df_women) <- gsub("^female_", "", colnames(df_women))

        output_file_comb <- paste("../full.summary.gwas/", trait, "_comb.txt", sep = "")
        output_file_men <- paste("../full.summary.gwas/", trait, "_men.txt", sep = "")
        output_file_women <- paste("../full.summary.gwas/", trait, "_women.txt", sep = "")

        write.table(df_comb, output_file_comb, quote = F, row.names = F, sep = "\t")
        write.table(df_men, output_file_men, quote = F, row.names = F, sep = "\t")
        write.table(df_women, output_file_women, quote = F, row.names = F, sep = "\t")

}

#To get the instruments
get.exp.df <- function(snp_group) {
        input_file <-  paste("../", snp_group, ".180513.txt", sep = "")
        exp_df <- read.table(input_file, stringsAsFactors = F, header = T)
	
	colnames(exp_df) <- gsub("\\..+", "", colnames(exp_df))
	
	colnames(exp_df)[colnames(exp_df) == "snp"] <- "SNP"
	colnames(exp_df)[colnames(exp_df) %in% c("ea", "EA")] <- "exp_ea"
	colnames(exp_df)[colnames(exp_df) %in% c("nea", "NEA")] <- "exp_nea"
	colnames(exp_df)[colnames(exp_df) == "p"] <- "pval"
	colnames(exp_df)[colnames(exp_df) %in% c("n", "nmeta")] <- "samplesize"
	colnames(exp_df)[colnames(exp_df) %in% c("eaf", "frqA1")] <- "exp_eaf"
	colnames(exp_df)[colnames(exp_df) == "beta"] <- "exp_beta"
	colnames(exp_df)[colnames(exp_df) == "se"] <- "exp_se"	

	print(range(exp_df$exp_beta))
        return(exp_df)
}

#To get the MRs: 
get.missing.snps <- function(results_df, snp_group, trait) {
	exp_df <- get.exp.df(snp_group)
	
	print(snp_group)
	sex <- gsub("(bmi|whr|whradjbmi)\\.eur\\.|\\.pulit.+|\\.sig|\\.phet", "", snp_group)

	print(sex)
	
	sex_input_file <- paste("../full.summary.gwas/", trait, "_", sex, ".txt", sep = "")
	out_df <- read.table(sex_input_file, stringsAsFactors = F, header = T, sep = "\t")

	print(nrow(exp_df))
	snp_not_in_file <- exp_df[!(exp_df$SNP %in% out_df$rs_number), "SNP"]
	print(length(snp_not_in_file))
	print(snp_not_in_file)
	results_df <- c(results_df, snp_not_in_file)
	print(results_df)
	return(results_df)

}

snp_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig",
                "whradjbmi.eur.comb.pulit.sig",
                "bmi.eur.men.pulit.sig", "bmi.eur.men.pulit.phet",
                "whr.eur.men.pulit.sig", "whr.eur.men.pulit.phet",
                "whradjbmi.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.phet",
                "bmi.eur.women.pulit.sig", "bmi.eur.women.pulit.phet",
                "whr.eur.women.pulit.sig", "whr.eur.women.pulit.phet",
                "whradjbmi.eur.women.pulit.sig",
                "whradjbmi.eur.women.pulit.phet")

traits <- c("FG", "FI")

results_df <- c()

for (trait in traits) {
	for (snp_group in snp_groups) {
		results_df <- get.missing.snps(results_df, snp_group, trait)
	}
}

results_df <- results_df[!duplicated(results_df)]
results_df<- gsub("_[A-Z]+_[A-Z]+$", "", results_df)
write.table(results_df, "../full.summary.gwas/adiposity.snps.not.in.fg.fi.180921.txt", quote = F, 
			row.names = F, sep = "\t", col.names = F)

