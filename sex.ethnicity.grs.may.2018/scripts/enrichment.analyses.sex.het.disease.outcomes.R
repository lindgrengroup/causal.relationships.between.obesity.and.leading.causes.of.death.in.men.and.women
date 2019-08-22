#$/bin/env Rscript
#!-cwd

adiposity_traits <- c("bmi", "whr", "whradjbmi")
traits <- c("t2d_cases_probposs", "cad_cases", "copd_cases", "renal_failure_cases", "ckd_cases")
sex_het_groups <- c("sex_het", "not_sex_het")
datasets <- c("pulit")
eth_groups <- c("all.white")

results_df = data.frame("dataset" = character(), "eth_group" = character(),
	"adiposity_trait" = character(), "trait" = character(),
	"test" = character(), "p" = numeric(), "n_in_input" = integer(), "n_with_p" = integer(), 
	"n_sig_het_adi" = numeric(),
	"n_not_sig_het_adi" = numeric(),
	"n_sig_het_trait" = numeric(),
	"n_sig_not_het_trait" = numeric(),
	stringsAsFactors = F)

for (dataset in datasets) {
	for (eth_group in eth_groups) {
		for (adiposity_trait in adiposity_traits) {
			for (trait in traits) {
	
				#Read in the SNPs that are and are not significantly sexually heterogeneous
				input_file_het <- paste("../sex.het.enrichment/sex.het.snps.", adiposity_trait, 
					".", dataset, ".no.tri.180921.txt", sep = "")
				het_df <- read.table(input_file_het, stringsAsFactors = F, header = F, sep = " ")

				input_file_not_het <- paste("../sex.het.enrichment/not.sex.het.snps.", adiposity_trait,
        				".", dataset, ".no.tri.180921.txt", sep = "")
				not_het_df <- read.table(input_file_not_het, stringsAsFactors = F, header = F, sep = " ")

				n_in_input <- nrow(het_df) + nrow(not_het_df)

				#Get the proportion that are nominally significant in each sex-het category - note
				#that may be that there is beta NA for some...
				trait_input_file <- paste("../sex.het.enrichment/", eth_group, ".comb.", dataset, ".", 
					trait, ".glm.logistic", sep = "")
				trait_df <- read.table(trait_input_file, stringsAsFactors = F, header = F)
				colnames(trait_df) <- c("CHROM", "POS", "ID", "REF", "ALT", "TEST", 
						"OBS_CT", "OR", "SE", "T_STAT", "P")
				trait_df <- trait_df[trait_df$TEST == "ADD" & !is.na(trait_df$P), ]
				n_with_p <- nrow(trait_df[trait_df$ID %in% c(het_df$V1, not_het_df$V1), ])

				trait_df <- trait_df[trait_df$TEST == "ADD" & trait_df$P < 0.05 & !is.na(trait_df$P), ]

				sig_het <- het_df[het_df$V1 %in% trait_df$ID, ]
				sig_not_het <- not_het_df[not_het_df$V1 %in% trait_df$ID, ]
				total_sig <- c(sig_het, sig_not_het)
				
				#Compute the expected count in each cell to decide whether Chi-square or Fischer's 
				#exact test
				raw_prop_sig <- nrow(het_df)/(nrow(not_het_df)+nrow(het_df))

				if (raw_prop_sig * length(total_sig) < 5 | (1-raw_prop_sig)*length(total_sig) < 5) {
					fisher_m <- matrix(c(length(sig_het), length(sig_not_het), nrow(het_df), nrow(not_het_df)), nrow = 2)
					p <- fisher.test(fisher_m, alternative = "two.sided")$p.value
					print(fisher.test(fisher_m, alternative = "two.sided"))
					
					test <- "fisher"				
				} else if (raw_prop_sig * length(total_sig) >= 5 & (1-raw_prop_sig)*length(total_sig) >= 5) {
					p <- chisq.test(x = c(length(sig_het), length(sig_not_het)), p = c(raw_prop_sig, (1-raw_prop_sig)))$p.value
					print(chisq.test(x = c(length(sig_het), length(sig_not_het)), p = c(raw_prop_sig, (1-raw_prop_sig))))
					
					test <- "chi.sq"
				}

				current_row = data.frame("dataset" = dataset, "eth_group" = eth_group,
                                	"adiposity_trait" = adiposity_trait, "trait" = trait,
                                        "test" = test, "p" = p, 
					"n_in_input" = n_in_input, "n_with_p" = n_with_p, 
					"n_sig_het_adi" = nrow(het_df),
                                        "n_not_sig_het_adi" = nrow(not_het_df),
                                        "n_sig_het_trait" = length(sig_het),
                                        "n_sig_not_het_trait" = length(sig_not_het),
                                        stringsAsFactors = F)

				results_df <- rbind(results_df, current_row)
				
				write.table(results_df, "../sex.het.enrichment/results.sex.enrichment.cad.t2d.ukbb.txt", 
					quote = F, row.names = F, sep = "\t")

			}
		}
	}
}

	

				
				
				
	
