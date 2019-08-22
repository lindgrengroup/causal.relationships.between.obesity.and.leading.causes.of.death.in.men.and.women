#!/bin/env Rscript
#$-cwd

traits <- c("bmi", "whr", "whradjbmi")
sex_groups <- c("comb", "men", "women")

for (trait in traits) {
	for (sex_group in sex_groups) {
	
		input_file_index <- paste("../", trait, ".eur.", sex_group, ".pulit.sig.180513.txt", sep = "")
		df_index <- read.table(input_file_index, stringsAsFactors = F, header = T)

		#Make a separate one for the unweighted weighting
		df_for_unweighted <- df_index
		colnames(df_for_unweighted) <- gsub("\\.combined$|\\.females$|\\.males$", "", colnames(df_for_unweighted))
		df_for_unweighted$TEST <- NA
		df_for_unweighted$STAT <- NA
		df_for_unweighted <- df_for_unweighted[, c("SNP", "EA", "NEA", "Chr", "Pos", "TEST", "nmeta", "beta", "STAT", "pval")]
		colnames(df_for_unweighted) <- c("SNP", "EA", "NEA", "CHR", "BP", "TEST", "NMISS", "BETA", "STAT", "P")

		df_index <- df_index[, c("SNP", "EA", "NEA")]

		input_file_weights <- paste("../temporary.ukbb.weights/", ifelse(trait == "whradjbmi", "res_whr_inv", trait), 
			"_", sex_group, ".assoc.linear", sep = "")
		df_weights <- read.table(input_file_weights, stringsAsFactors = F, header = T)

		df <- merge(df_index, df_weights, by = "SNP")

		#Flip beta if the A1 is not EA, so that beta corresponds to ea
		df[df$A1 != df$EA, "BETA"] <- df[df$A1 != df$EA, "BETA"] * -1
		df$A1 <- NULL

		#For some reason, at least one SNP is NA for BETA
		df <- df[!(is.na(df$BETA)), ]

		#Make sure all are positive
        	for (k in 1:nrow(df)) {
        	        if (df[k, "BETA"] < 0) {
        	                df[k, c("EA", "NEA")] <- df[k, c("NEA", "EA")]
        	                df[k, "BETA"] <- df[k, "BETA"] * -1
        	        }
        	}
	
		#Use "normal" .sig derived from internal weights
		output_file <- paste("../sensitivity.risk.weights/", trait, ".eur.", sex_group, ".internal.normal.180513.txt", sep = "")
		write.table(df, output_file, quote = F, row.names = F)

		#Half all the weights
		df_half <- df
		df_half$BETA <- df_half$BETA / 2
		output_file <- paste("../sensitivity.risk.weights/", trait, ".eur.", sex_group, ".internal.half.180513.txt", sep = "")
		write.table(df_half, output_file, quote = F, row.names = F)

		#Introduce noise
		set.seed(42)
		df_noise <- df
		noise <- rnorm(nrow(df), 0, sd(df$BETA))
		df_noise$BETA <- df_noise$BETA + noise 
		range(df_noise$BETA)

		output_file <- paste("../sensitivity.risk.weights/", trait, ".eur.", sex_group, ".internal.randomnoise.180513.txt", sep = "")
        	write.table(df_noise, output_file, quote = F, row.names = F)

		#Systematic bias
		df_systematic <- df
		df_systematic[, "BETA"] <- ifelse(df_systematic$BETA < median(df_systematic$BETA), df_systematic$BETA /2, df_systematic$BETA *2)
	
        	output_file <- paste("../sensitivity.risk.weights/", trait, ".eur.", sex_group, ".internal.systematicnoise.180513.txt", sep = "")
        	write.table(df_systematic, output_file, quote = F, row.names = F)
	
		#Unweighted scores
		df_unweighted <- df_for_unweighted
		df_unweighted[, "BETA"] <- 1
	
        	output_file <- paste("../sensitivity.risk.weights/", trait, ".eur.", sex_group, ".internal.unweighted.180513.txt", sep = "")
        	write.table(df_unweighted, output_file, quote = F, row.names = F)
	}
}
