#!/bin/env Rscript
#$-cwd

library(MendelianRandomization)
library(mada)

#To get the instruments
get.exp.df <- function(snp_group) {
	if (length(snp_group[grep("fdr|giukbb|internal", snp_group)]) != 1) {
                input_file <-  paste("../", snp_group, ".180513.txt", sep = "")
        } else if (length(snp_group[grep("fdr|giukbb|internal", snp_group)]) == 1) {
                input_file <- paste("../sensitivity.risk.weights/", snp_group, ".180513.txt", sep = "")
        }

        exp_df <- read.table(input_file, stringsAsFactors = F, header = T)
	
	if (length(snp_group[grep("fdr|pulit", snp_group)]) == 0) {
	        colnames(exp_df) <- gsub("[a-z]+\\.", "", colnames(exp_df))
	} else if (length(snp_group[grep("pulit|fdr|giukbb", snp_group)]) == 1) {
		colnames(exp_df) <- gsub("\\.(.)+", "", colnames(exp_df))
	}
	
	colnames(exp_df)[colnames(exp_df) == "snp"] <- "SNP"
	colnames(exp_df)[colnames(exp_df) %in% c("ea", "EA")] <- "exp_ea"
	colnames(exp_df)[colnames(exp_df) %in% c("nea", "NEA")] <- "exp_nea"
	colnames(exp_df)[colnames(exp_df) %in% c("p", "P")] <- "pval"
	colnames(exp_df)[colnames(exp_df) %in% c("n", "nmeta", "NMISS")] <- "samplesize"
	colnames(exp_df)[colnames(exp_df) %in% c("eaf", "frqA1")] <- "exp_eaf"
	colnames(exp_df)[colnames(exp_df) %in% c("beta", "BETA")] <- "exp_beta"
	colnames(exp_df)[colnames(exp_df) == "se"] <- "exp_se"	

        #exp_df$SNP <- gsub("_[A-Z]+_[A-Z]+$", "", exp_df$SNP)
	print(range(exp_df$exp_beta))
        return(exp_df)
}

#To get the MRs: 
get.mr.results <- function(results_df, snp_group, trait) {
	exp_df <- get.exp.df(snp_group)
	
	print(snp_group)
	sex <- gsub("(bmi|whr|whradjbmi)\\.eur\\.|\\.sig|\\.pulit(.)+|\\.internal|\\.unweighted|\\.pulit|\\.giukbb(.)+|\\.phet|\\.0(.)+", "", snp_group)
	
	n_exp_snps <- nrow(exp_df)
	
	sex_input_file <- paste("../full.summary.gwas/", trait, "_", sex, ".txt", sep = "")
	out_df <- read.table(sex_input_file, stringsAsFactors = F, header = T, sep = "\t")

	#Basically only gain about half of the missing SNPs if using proxies - e.g. for BMI, ~570 SNPs, 
	#538 are available with proxies, 506 without... So, considering the 
	#danger of using proxies suggest not use
	#proxies <- read.table("../full.summary.gwas/proxies.for.fg.fi.180922.txt", stringsAsFactors = F, 
	#		header = T, sep = "\t")
	#proxies[, c("alpha_first", "alpha_second")] <- t(apply(proxies[, c("MAJOR", "MINOR")], 1, sort))
        #proxies$RSID <- paste(proxies$RSID, "_", proxies$alpha_first, "_", proxies$alpha_second, sep ="")
        #proxies[, c("alpha_first", "alpha_second")] <- NULL
	#exp_df$simple_rsid <- gsub("_[A-Z]+_[A-X]+", "", exp_df$SNP)
	#exp_df <- merge(exp_df, proxies, by.x = "simple_rsid", by.y = "QRSID", all.x = T)
	#exp_df[is.na(exp_df$RSID), "RSID"] <- exp_df[is.na(exp_df$RSID), "SNP"]
	#mr_df <- merge(exp_df, out_df, by.x = "RSID", by.y = "rs_number")
	#mr_df <- mr_df[order(mr_df$SNP, -mr_df$R2, abs(mr_df$DIST)), ]
	#mr_df <- mr_df[!(duplicated(mr_df$SNP)), ]

	mr_df <- merge(exp_df, out_df, by.x = "SNP", by.y = "rs_number")

	#Realized that my script to update rsIDs might introduce duplicated SNPs. Double-checked, and no duplicated SNPs
	#are kept all the way here, so should be fine
	print(mr_df[duplicated(mr_df$Pos), ])
	print(mr_df[duplicated(mr_df$pos), ])
	n_exp_snps_in_out <- nrow(mr_df)

	mr_df$flipped_ea <- mr_df$reference_allele
	mr_df$flipped_nea <- mr_df$other_allele
	mr_df$flipped_eaf <- mr_df$eaf
	mr_df$flipped_beta <- mr_df$beta
	
	mr_df[mr_df$reference_allele != mr_df$exp_ea, "flipped_beta"] <- mr_df[mr_df$reference_allele != mr_df$exp_ea, "flipped_beta"] * -1
	mr_df[mr_df$reference_allele != mr_df$exp_ea, "flipped_eaf"] <- 1 - mr_df[mr_df$reference_allele != mr_df$exp_ea, "flipped_eaf"]
	mr_df[mr_df$reference_allele != mr_df$exp_ea, c("flipped_ea", "flipped_nea")] <- mr_df[mr_df$reference_allele != mr_df$exp_ea, c("other_allele", "reference_allele")]

	mr_df <- mr_df[mr_df$exp_ea == mr_df$flipped_ea & mr_df$exp_nea == mr_df$flipped_nea, ]
	n_snps_same_alleles <- nrow(mr_df)
	
        set.seed(42)
	mr_results <- mr_allmethods(mr_input(bx = mr_df$exp_beta, bxse = mr_df$exp_se, 
			by = mr_df$flipped_beta, byse = mr_df$se, exposure = snp_group, outcome = trait), method = "all")$Values
	colnames(mr_results) <- c("Method", "Estimate", "Se", "LCI", "UCI", "P")

	print(mr_results)
	for (i in 1:nrow(mr_results)) {

		if (mr_results[i, "Method"] == "(intercept)") {
			mr_results[i, "Method"] <- paste(mr_results[i-1, "Method"], " ", mr_results[i, "Method"], sep = "")		
		}
		current_row <- data.frame("snp_group" = snp_group, "sex_group" = sex, "trait" = trait,
                "method" = mr_results[i, "Method"], "beta" = mr_results[i, "Estimate"], "se" = mr_results[i, "Se"],
                "beta_lci" = mr_results[i, "LCI"], "beta_uci" = mr_results[i, "UCI"], "p" = mr_results[i, "P"],
		"n_exp_snps" = n_exp_snps, "n_exp_snps_in_out" = n_exp_snps_in_out, "n_snps_same_alleles" = n_snps_same_alleles, 
		stringsAsFactors =F)

		results_df <- merge(results_df, current_row, by = colnames(results_df), all = T)
	}	
	print(results_df)
	return(results_df)
	
}

#Does not include unweigthed since don't have SEs for those
snp_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig",
                "whradjbmi.eur.comb.pulit.sig",
                "bmi.eur.men.pulit", "bmi.eur.men.pulit.sig", "bmi.eur.men.pulit.phet",
                "whr.eur.men.pulit", "whr.eur.men.pulit.sig", "whr.eur.men.pulit.phet",
                "whradjbmi.eur.men.pulit", "whradjbmi.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.phet",
                "bmi.eur.women.pulit", "bmi.eur.women.pulit.sig", "bmi.eur.women.pulit.phet",
                "whr.eur.women.pulit", "whr.eur.women.pulit.sig", "whr.eur.women.pulit.phet",
                "whradjbmi.eur.women.pulit", "whradjbmi.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.phet",
                "bmi.eur.men.0.01.fdr", "bmi.eur.men.0.05.fdr",
                "bmi.eur.men.0.1.fdr", "bmi.eur.women.0.01.fdr",
                "bmi.eur.women.0.05.fdr", "bmi.eur.women.0.1.fdr",
                "whr.eur.men.0.01.fdr", "whr.eur.men.0.05.fdr",
                "whr.eur.men.0.1.fdr", "whr.eur.women.0.01.fdr",
                "whr.eur.women.0.05.fdr", "whr.eur.women.0.1.fdr",
                "whradjbmi.eur.men.0.01.fdr", "whradjbmi.eur.men.0.05.fdr",
                "whradjbmi.eur.men.0.1.fdr", "whradjbmi.eur.women.0.01.fdr",
                "whradjbmi.eur.women.0.05.fdr", "whradjbmi.eur.women.0.1.fdr",
                "bmi.eur.comb.giukbb.sig", "whr.eur.comb.giukbb.sig", "whradjbmi.eur.comb.giukbb.sig",
                "bmi.eur.men.giukbb.sig", "whr.eur.men.giukbb.sig", "whradjbmi.eur.men.giukbb.sig",
                "bmi.eur.women.giukbb.sig", "whr.eur.women.giukbb.sig", "whradjbmi.eur.women.giukbb.sig")

traits <- c("FG", "FI")

results_df <- data.frame("snp_group" = character(), "sex_group" = character(), "trait" = character(), 
		"method" = character(), "beta" = numeric(), "se" = numeric(), 
		"beta_lci" = numeric(), "beta_uci" = numeric(), "p" = numeric(), 
		"n_exp_snps" = integer(), "n_exp_snps_in_out" = integer(), "n_snps_same_alleles" = integer(),
		stringsAsFactors =F)

for (trait in traits) {
	for (snp_group in snp_groups) {
		results_df <- get.mr.results(results_df, snp_group, trait)
	}
}

write.table(results_df, "../results.mr.180730/summary.mr.results.180730.txt", quote = F, 
			row.names = F, sep = "\t")

#Get heterogeneity P-values
women_snp_groups <- snp_groups[grep("\\.women", snp_groups)]
men_snp_groups <- snp_groups[grep("\\.men", snp_groups)]

women_snp_groups <- women_snp_groups[order(women_snp_groups)]
men_snp_groups <- men_snp_groups[order(men_snp_groups)]

#Prepare some columns required for Heterogeneity calculations
results_df$cochran_weights <- 1 / (results_df$se^2)

for (i in 1:length(women_snp_groups)) {
	for (trait in unique(results_df$trait)) {
        	for (method in unique(results_df$method)) {
                	general <- results_df$trait == trait &
                        	results_df$method == method

                        if (nrow(results_df[(general &
                                results_df$snp_group == women_snp_groups[i] &
                                results_df$sex_group == "women") |
                                (general &
                                results_df$snp_group == men_snp_groups[i] &
                                results_df$sex_group == "men"), ]) == 2) {


                                results_df$cochrans_q[(general &
                                	results_df$snp_group == women_snp_groups[i] &
                                        results_df$sex_group == "women") |
                                        (general &
                                        results_df$snp_group == men_snp_groups[i] &
                                        results_df$sex_group == "men")] <- as.numeric(cochran.Q(c(results_df[general &
                                        results_df$snp_group == women_snp_groups[i] &
                                        results_df$sex_group == "women",
                                        "beta"],
                                        results_df[general &
                                        results_df$snp_group == men_snp_groups[i] &
                                        results_df$sex_group == "men",
                                        "beta"]),
                                        c(results_df[general &
                                        results_df$snp_group == women_snp_groups[i] &
                                        results_df$sex_group == "women",
                                        "cochran_weights"],
                                        results_df[general &
                                        results_df$snp_group == men_snp_groups[i] &
                                        results_df$sex_group == "men",
                                        "cochran_weights"]))[1])

                                results_df$cochrans_p[(general &
                                        results_df$snp_group == women_snp_groups[i] &
                                        results_df$sex_group == "women") |
                                        (general &
                                        results_df$snp_group == men_snp_groups[i] &
                                        results_df$sex_group == "men")] <- as.numeric(cochran.Q(c(results_df[general &
                                        results_df$snp_group == women_snp_groups[i] &
                                        results_df$sex_group == "women",
                                        "beta"],
                                        results_df[general &
                                        results_df$snp_group == men_snp_groups[i] &
                                        results_df$sex_group == "men",
                                        "beta"]),
                                        c(results_df[general &
                                        results_df$snp_group == women_snp_groups[i] &
                                        results_df$sex_group == "women",
                                        "cochran_weights"],
                                        results_df[general &
                                        results_df$snp_group == men_snp_groups[i] &
                                        results_df$sex_group == "men",
                                        "cochran_weights"]))[2])

                        }
                }
        }
}

results_df$cochrans_i2 <- (results_df$cochrans_q - 1)/results_df$cochrans_q
results_df$cochrans_i2[results_df$cochrans_i2 < 0] <- 0

write.table(results_df, "../results.mr.180730/summary.mr.results.180730.txt", quote = F,
                        row.names = F, sep = "\t")

