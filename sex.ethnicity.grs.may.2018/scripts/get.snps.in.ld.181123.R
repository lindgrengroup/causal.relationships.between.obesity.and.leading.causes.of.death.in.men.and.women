#!/bin/env Rscript
#$-cwd

traits <- c("bmi", "whr", "whradjbmi")
sex_groups <- c("comb", "men", "women")
datasets <- c("pulit", "pulit_winner")

results_df <- data.frame("trait" = character(), "sex_group" = character(), "dataset" = character(), 
		"CHR_A" = integer(),  "BP_A" = integer(),  "SNP_A" = character(), 
		"CHR_B" = integer(), "BP_B" = integer(), "SNP_B" = character(), "R2" = numeric(), 
		stringsAsFactors = F)

for (trait in traits) {
	for (sex_group in sex_groups) {
		for (dataset in datasets) {
		if (sex_group != "comb" & dataset == "pulit_winner") {
			next()
		}	

		input_file <- paste("../temp.ld.between.snps.181123/", trait, ".", sex_group, ".", dataset, ".snps.ld.181123.ld", sep = "")
		df <- read.table(input_file, stringsAsFactors = F, header = T)
		df$trait <- trait
		df$sex_group <- sex_group
		df$dataset <- dataset

		ld_snps <- df[df$R2 >= 0.05, ]
		results_df <- rbind(results_df, ld_snps)

		if (nrow(ld_snps) >= 1) {
			print(paste(trait, sex_group, dataset, sep = "_"))
			print(ld_snps)
		}
	
		}
	}
}

for (trait in traits) {
	for (sex_group in sex_groups) {
		for (dataset in datasets) {
				comb_input_file <- paste("../original.index.files/", trait, ".eur.comb.pulit.180513.txt", sep = "")
				men_input_file <- paste("../original.index.files/", trait, ".eur.men.pulit.180513.txt", sep = "")
				women_input_file <- paste("../original.index.files/", trait, ".eur.women.pulit.180513.txt", sep = "")

				comb <- read.table(comb_input_file, stringsAsFactors=F, header = T)
				men <- read.table(men_input_file, stringsAsFactors=F, header = T)
				women <- read.table(women_input_file, stringsAsFactors=F, header = T)

				org <- rbind(rbind(comb, men), women)
				org <- org[!duplicated(org$SNP), ]
				org$simple <- gsub(":[ATCG]:[ATCG]", "", org$SNP)
				org[, c("alpha_first", "alpha_second")] <- t(apply(org[, c("A1.combined", "A2.combined")], 1, sort))
				org$SNP <- paste(org$simple, "_", org$alpha_first, "_", org$alpha_second, sep ="")

				pval_column <- ifelse(sex_group == "comb", "pval.combined", 
						ifelse(sex_group == "men", "pval.males", 
						ifelse(sex_group == "women", "pval.females", "MISSING")))

				for (i in 1:nrow(results_df)) {	
					if (results_df[i, c("trait")] == trait & results_df[i, c("sex_group")] == sex_group &
					results_df[i, "dataset"] == dataset) {
					
						#Decide which one to keep based on P-value for that sex group
						results_df[i, "SNP_A_P"] <- org[org$SNP == results_df[i, "SNP_A"], pval_column]
						results_df[i, "SNP_B_P"] <- org[org$SNP == results_df[i, "SNP_B"], pval_column]
						results_df[i, "remove"] <- ifelse(results_df[i, "SNP_A_P"] >= results_df[i, "SNP_B_P"], 
									results_df[i, "SNP_A"], 
									ifelse(results_df[i, "SNP_A_P"] < results_df[i, "SNP_B_P"], 
									results_df[i, "SNP_B"], "MISSING"))
					}
			} 
		}
	}
}				

write.table(results_df, "../snps.that.are.correlated.over.0.05.181123.txt", quote = F, row.names = F, sep = "\t")
