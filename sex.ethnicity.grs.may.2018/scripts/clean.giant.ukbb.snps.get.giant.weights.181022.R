#!/bin/env Rscript
#$-cwd

make.beta.pos <- function(df) {
	for (m in 1:nrow(df)) {
		if (df[m, grep("beta", colnames(df))] < 0) {
			df[m, c("EA", "NEA")] <- df[m, c("NEA", "EA")]
			df[m, "beta"] <- df[m, "beta"] * -1
			df[m, "EAF"] <- 1 - df[m, "EAF"] 
		}
	}
	return(df)
}

loop <- expand.grid(trait = c("bmi", "whr", "whradjbmi"), sex_group = c("comb", "men", "women"),
	snp_group = c("pulit"), eth_group = "eur", stringsAsFactors = F)

results_df <- data.frame("snp_group" = character(), "snps_in_giant_and_ukbb_no_tri" = integer(), 
		"snps_in_giant" = integer(), "n_snps_missing" = integer(), 
		"snps_missing" = character(), "n_added_together" = integer(), stringsAsFactors = F)

for (i in 1:nrow(loop)) {
	input_file <- paste("../", loop[i, "trait"], ".eur.", loop[i, "sex_group"], ".pulit.sig.180513.txt", sep = "") 
	raw_df <- read.table(input_file, stringsAsFactors = F, header = T)

	if (loop[i, "sex_group"] == "comb") {
		sex_group_pulit <- "combined"
	} else if (loop[i, "sex_group"] == "women") {
		sex_group_pulit <- "females" 
	} else if (loop[i, "sex_group"] == "men") {	
		sex_group_pulit <- "males"
	}
 
	giant_input_file <- paste("/well/lindgren/sara/giant-fineMapping/meta-analysis/metal/alignSNPs/giant/giant.2015.", 
			loop[i, "trait"], ".", sex_group_pulit, ".eur.dbsnp151.ukbb.txt.gz", sep = "")
	giant <- read.table(gzfile(giant_input_file), stringsAsFactors = F, header = T, fill = T)
	giant <- giant[!is.na(giant$MarkerName), ]
	#Works because all SNPs are "single alleles"
	giant$SNP <- gsub(":[ATCG]:[ATCG]", "", giant$SNP)
	giant[, c("alpha_first", "alpha_second")] <- t(apply(giant[, c("A1", "A2")], 1, sort))
        giant$SNP <- paste(giant$SNP, "_", giant$alpha_first, "_",
                giant$alpha_second, sep ="")
        giant[, c("alpha_first", "alpha_second")] <- NULL

	df <- merge(raw_df[, c("SNP", "Chr", "Pos")], giant, by = c("SNP"))
	df$MarkerName <- NULL
	colnames(df) <- c("SNP", "Chr", "Pos", "EA", "NEA", "EAF", "beta", "se", "p", "n")
	missing_snps <- raw_df[!(raw_df$SNP %in% df$SNP), "SNP"]

	df <- make.beta.pos(df)
	output_file_sig <- paste("../sensitivity.risk.weights/", loop[i, "trait"], ".", loop[i, "eth_group"], ".", loop[i, "sex_group"], 
			".giukbb.sig.180513.txt", sep = "")		

	write.table(df, output_file_sig, row.names = F, quote = F, sep = " ")
	
	current_row <- data.frame("snp_group" = paste(loop[i, "trait"], ".eur.", loop[i, "sex_group"], ".pulit.sig", sep = ""), 
		"snps_in_giant_and_ukbb_no_tri" = nrow(raw_df),
                "snps_in_giant" = nrow(df), "n_snps_missing" = length(missing_snps),
		"n_added_together" = length(missing_snps) + nrow(df), 
		"snps_missing" = paste(missing_snps, collapse = "; "), stringsAsFactors = F)

	results_df <- rbind(results_df, current_row)
}

write.table(results_df, "../list.of.giant.ukbb.snps.present.in.giant.txt", quote = F, row.names = F, sep = "\t")
