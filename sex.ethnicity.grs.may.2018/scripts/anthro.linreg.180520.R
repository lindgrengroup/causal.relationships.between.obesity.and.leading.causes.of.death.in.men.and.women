#!/bin/env Rscript
#$ -cwd

##########################################################
###Sex and ethnicity WHR GRS on anthropometric traits ####
############ Jenny Censin, 2018-05-20 ###################

library(mada)

qc.pregnant <- function(data, qc.counts.file) { 
cleaned <- subset(data, is.na(data$pregnant) | data$pregnant == 0)

sink(qc.counts.file, append = T)
cat(paste("**FILTER** Number of pregnant indidividuals, yes or unsure: ", 
	nrow(subset(data, !(is.na(data$pregnant) | data$pregnant == 0))), 
	" ; REMAINING: ", nrow(cleaned), "\n", sep = ""))
sink()

return(cleaned)
}

qc.below.15.bmi <- function(data, qc.counts.file) {
cleaned <- subset(data, is.na(data$bmi) | data$bmi > 15)

sink(qc.counts.file, append = T)
cat(paste("**FILTER** Number of individuals with BMI below 15: ", 
	nrow(subset(data, !(is.na(data$bmi) | data$bmi> 15))), 
	" ; REMAINING: ", nrow(cleaned), "\n", sep = ""))
sink()

return(cleaned)
} 

linreg <- function(pheno, snp_groups, eth_groups, sex_groups, traits) {

	pheno$dummy_array <- ifelse(pheno$genotyping_array == "UKBB", 0, 1)
	pheno$dummy_sex <- ifelse(pheno$Submitted_Gender == "F", 0, 1)
	results_df <- data.frame("snp_group" = character(), "grs_unit" = character(), "trait_unit" = character(),
				"eth_group" = character(), "sex_group" = character(), "trait" = character(),
                                "grs_r2" = numeric(), "grs_p" = numeric(), "grs_beta" = numeric(),
				"grs_se" = numeric(), "grs_lci" = numeric(), "grs_uci" = numeric(),
				"grs_file" = character(), "n_trait" = integer(),
				"n_complete_cases" = integer(), "f_value" = numeric(),
				"extra_adjustment" = character(), stringsAsFactors = F)

	#Make the covariates
	pcs <- paste("PC", 1:10, sep = "", collapse = " + ")

	for (snp_group in snp_groups) {
	        file <- paste("../bmi.whr.profiles/grs.", snp_group, ".180513.profile", sep = "")
        	scores <- read.table(file, stringsAsFactors = F, header = T)
		df <- merge(x = scores, y = pheno, by.x = "IID", by.y = "ID")		

        	for (eth_group in eth_groups) {
                	if (eth_group == "brit.irish") {
                        	eth_df <- subset(df, df$check_se %in% c("brit"))
                	} else if (eth_group == "all.white") {
                        	eth_df <- subset(df, df$check_se %in% c("brit", "white", "recode.white"))
                	}
              		for (sex_group in sex_groups) {
                        	if (sex_group == "comb") {
                                	sex_df <- eth_df
                        	} else if (sex_group == "women") {
                                	sex_df <- subset(eth_df, eth_df$Submitted_Gender == "F")
                        	} else if (sex_group == "men") {
                                	sex_df <- subset(eth_df, eth_df$Submitted_Gender == "M")
                        	}

				#Make columns with WHR adjusted for BMI - sex is dropped in sex-specific analyses automatic
                        	whrfit <- lm(whr ~ age_assessment + age_squared + dummy_sex + assessment_centre + bmi, 
					data = sex_df, na.action = na.exclude)
                        	sex_df$res_whr <- residuals(whrfit)
                        	sex_df$res_whr_inv <- qnorm((rank(sex_df$res_whr, 
					na.last = "keep") - 0.5) / sum(!is.na(sex_df$res_whr)))
				
				for (grs_unit in c("raw_scoresum")) {
					grs_unit_df <- sex_df
					if (grs_unit == "sd_scoresum") {
						grs_unit_df$SCORESUM <- (grs_unit_df$SCORESUM - mean(grs_unit_df$SCORESUM, 
									na.rm = T))/sd(grs_unit_df$SCORESUM, na.rm = T) 
					}
				
                        		#Make the formula:
                        		for (trait in traits) {
						for (extra_adjustment in c("-", "smoking")) {
							for (trait_unit in c("clin", "sd")) {
								
								#Make RINTed traits - for WHRadjBMI, all use the previous one in the regressions - sex is dropped if only 1 sex
								trait_unit_df <- grs_unit_df
								if (trait != "res_whr_inv") {
									rint.formula <- as.formula(paste(trait,
        	                                                        "~age_assessment+age_squared+dummy_sex+assessment_centre", sep = " "))
									rint_fit <- lm(rint.formula, data = trait_unit_df, na.action = na.exclude)
									trait_unit_df$res_rint <- residuals(rint_fit)
									trait_unit_df$res_rint_inv <- qnorm((rank(trait_unit_df$res_rint,
				        	                                na.last = "keep") - 0.5) / sum(!is.na(trait_unit_df$res_rint)))
								}		

								comb_groups_all <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig", "whradjbmi.eur.comb.pulit.sig", 
										"bmi.eur.comb.giukbb.sig", "whr.eur.comb.giukbb.sig", "whradjbmi.eur.comb.giukbb.sig")
								men_groups_all <- c("bmi.eur.men.pulit", "bmi.eur.men.pulit.sig", "bmi.eur.men.pulit.phet", 
									"whr.eur.men.pulit", "whr.eur.men.pulit.sig", "whr.eur.men.pulit.phet", 
									"whradjbmi.eur.men.pulit", "whradjbmi.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.phet", 
							                "bmi.eur.men.0.01.fdr", "bmi.eur.men.0.05.fdr",
							                "bmi.eur.men.0.1.fdr", 
							                "whr.eur.men.0.01.fdr", "whr.eur.men.0.05.fdr",
							                "whr.eur.men.0.1.fdr", 
							                "whradjbmi.eur.men.0.01.fdr", "whradjbmi.eur.men.0.05.fdr",
							                "whradjbmi.eur.men.0.1.fdr", 
									"bmi.eur.men.giukbb.sig", "whr.eur.men.giukbb.sig", "whradjbmi.eur.men.giukbb.sig")

								women_groups_all <- c("bmi.eur.women.pulit", "bmi.eur.women.pulit.sig", "bmi.eur.women.pulit.phet", 
									"whr.eur.women.pulit", "whr.eur.women.pulit.sig", "whr.eur.women.pulit.phet", 
									"whradjbmi.eur.women.pulit", "whradjbmi.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.phet",
									"bmi.eur.women.0.01.fdr",
							                "bmi.eur.women.0.05.fdr", "bmi.eur.women.0.1.fdr",
									"whr.eur.women.0.01.fdr",
							                "whr.eur.women.0.05.fdr", "whr.eur.women.0.1.fdr",
									"whradjbmi.eur.women.0.01.fdr",
							                "whradjbmi.eur.women.0.05.fdr", "whradjbmi.eur.women.0.1.fdr", 
									"bmi.eur.women.giukbb.sig", "whr.eur.women.giukbb.sig", "whradjbmi.eur.women.giukbb.sig")

								groups_dont_use <- c("bmi.eur.men.pulit", "bmi.eur.men.pulit.phet", "bmi.eur.women.pulit", "bmi.eur.women.pulit.phet", 
									"whr.eur.men.pulit", "whr.eur.men.pulit.phet", "whr.eur.women.pulit", "whr.eur.women.pulit.phet",
									"whradjbmi.eur.men.pulit", "whradjbmi.eur.men.pulit.phet", "whradjbmi.eur.women.pulit", 
									"whradjbmi.eur.women.pulit.phet", 
							                "bmi.eur.men.0.01.fdr", "bmi.eur.men.0.05.fdr",
							                "bmi.eur.men.0.1.fdr", "bmi.eur.women.0.01.fdr",
							                "bmi.eur.women.0.05.fdr", "bmi.eur.women.0.1.fdr",
							                "whr.eur.men.0.01.fdr", "whr.eur.men.0.05.fdr",
							                "whr.eur.men.0.1.fdr", "whr.eur.women.0.01.fdr",
							                "whr.eur.women.0.05.fdr", "whr.eur.women.0.1.fdr",
							                "whradjbmi.eur.men.0.01.fdr", "whradjbmi.eur.men.0.05.fdr",
							                "whradjbmi.eur.men.0.1.fdr", "whradjbmi.eur.women.0.01.fdr",
							                "whradjbmi.eur.women.0.05.fdr", "whradjbmi.eur.women.0.1.fdr")

								if (!((trait %in% c("dbp", "sbp", "wc", "hc") & sex_group == "men" &
									snp_group %in% men_groups_all & !(snp_group %in% groups_dont_use)) |
									(trait %in% c("dbp", "sbp", "wc", "hc") & sex_group == "women" &
                        	                        	        snp_group %in% women_groups_all & !(snp_group %in% groups_dont_use)) |
									(trait %in% c("dbp", "sbp", "wc", "hc") & sex_group == "comb" &
									snp_group %in% comb_groups_all) |
									(trait %in% c("bmi", "whr", "res_whr_inv") & sex_group == "men" &
									snp_group %in% men_groups_all) |
									(trait %in% c("bmi", "whr", "res_whr_inv") & sex_group == "women" &
									snp_group %in% women_groups_all) |
									(trait %in% c("bmi", "whr", "res_whr_inv") & sex_group == "comb" & 
									snp_group %in% comb_groups_all) |
									(trait == "bmi" & snp_group %in% c("bmi.eur.comb.pulit.sig") 
										& sex_group %in% c("men", "women")) |
									(trait == "whr" & snp_group %in% c("whr.eur.comb.pulit.sig") & 
										sex_group %in% c("men", "women")) |
									(trait == "res_whr_inv" & snp_group %in% c("whradjbmi.eur.comb.pulit.sig") & 
										sex_group %in% c("men", "women")))) {
									next()
								} 
										
								if (trait == "res_whr_inv" & extra_adjustment == "-") {
                        		                		linreg.form <- as.formula(paste(trait, 
									"~SCORESUM+dummy_array+", pcs, sep = " "))
                                				} else if (trait_unit == "sd" & extra_adjustment == "-") {
									linreg.form <- as.formula(paste("res_rint_inv~SCORESUM+dummy_array+", pcs, sep = " "))
								} else if (trait_unit == "clin" & extra_adjustment == "-") {
                               		         			linreg.form <- as.formula(paste(trait, 
									"~SCORESUM+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+", 
									pcs, sep = " "))
								} else if (trait == "res_whr_inv" & extra_adjustment == "smoking") {
									linreg.form <- as.formula(paste(trait,
                                                	                "~SCORESUM+smoker_cases+dummy_array+", pcs, sep = " "))
								} else if (trait_unit == "sd" & extra_adjustment == "smoking") {
									linreg.form <- as.formula(paste("res_rint_inv~SCORESUM+smoker_cases+dummy_array+", pcs, sep = " "))
								} else if (trait_unit == "clin" & extra_adjustment == "smoking") {
									linreg.form <- as.formula(paste(trait,
                                                	                "~SCORESUM+smoker_cases+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                                	                pcs, sep = " "))
								}
							
								fit_raw <- lm(formula = linreg.form, data = trait_unit_df, na.action = na.exclude)
                        	                		print(summary(fit_raw))
                                	        		fit <- as.data.frame(summary(fit_raw)$coefficients)
	
								if (trait == "res_whr_inv" & extra_adjustment == "-") {
                                                	                linreg.form2 <- as.formula(paste(trait,
                                                        	        "~dummy_array+", pcs, sep = " "))
								} else if (trait_unit == "sd" & extra_adjustment == "-") {
									linreg.form2 <- as.formula(paste("res_rint_inv~dummy_array+", pcs, sep = " "))
                                                        	} else if (trait_unit == "clin" & extra_adjustment == "-") {
                                                        	        linreg.form2 <- as.formula(paste(trait,
                                                                	"~dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                                                	pcs, sep = " "))
                                                	        } else if (trait == "res_whr_inv" & extra_adjustment == "smoking") {
                                                        	        linreg.form2 <- as.formula(paste(trait,
                                                        	        "~smoker_cases+dummy_array+", pcs, sep = " "))
								} else if (trait_unit == "sd" & extra_adjustment == "smoking") {
									linreg.form2 <- as.formula(paste("res_rint_inv~smoker_cases+dummy_array+", pcs, sep = " "))
                                                        	} else if (trait_unit == "clin" & extra_adjustment == "smoking") {
                                                        	        linreg.form2 <- as.formula(paste(trait,
                                                        	        "~smoker_cases+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                                        	        pcs, sep = " "))
                                                        	}

        	                        		fit2 <- lm(formula = linreg.form2, data = trait_unit_df, na.action = na.exclude)
							iv_r <- summary(fit_raw)$r.squared - summary(fit2)$r.squared 

							current_row <- data.frame("snp_group" = snp_group, "grs_unit" = grs_unit, "trait_unit" = trait_unit,
								"eth_group" = eth_group, "sex_group" = sex_group, 
								"trait" = ifelse(extra_adjustment == "smoking", paste(trait, "_smoking", sep = ""), trait),
								"grs_r2" = round(summary(fit_raw)$adj.r.squared - summary(fit2)$adj.r.squared, 4),
								"grs_p" = fit["SCORESUM", "Pr(>|t|)"],
                                	                	"grs_beta" = fit["SCORESUM", "Estimate"],
	        	                               	        "grs_se" = round(fit["SCORESUM", "Std. Error"], 4),
        	                                        	"grs_lci_beta" = round(fit["SCORESUM", "Estimate"] - (1.96 * fit["SCORESUM", "Std. Error"]), 4),
	                	                                "grs_uci_beta" = round(fit["SCORESUM", "Estimate"] + (1.96 * fit["SCORESUM", "Std. Error"]), 4),
								"grs_file" = file, 
								"n_trait" = sum(!is.na(trait_unit_df[[gsub("_smoking", "", trait)]])), 
								"n_complete_cases" = nobs(fit_raw),
								"f_value" = ( (nobs(fit_raw)-1-1)/1 ) * (iv_r/(1-iv_r)), 
								"extra_adjustment" = extra_adjustment, 
								stringsAsFactors = F)
							print(current_row)
							results_df <- merge(results_df, current_row, by.x = colnames(results_df), 
								by.y = colnames(current_row), all = T)
							print(results_df)
							}
						}
					}
				}
			}
		}
	}
	return(results_df)	
}

#Make a phenotype QC counts file:
pheno.qc.count.file <- "../results.linear.regressions.180514/anthro.ukbb.pheno.qc.sample.count.180513.txt"

#Read in phenotype file:
pheno <- read.table("/well/lindgren/jc/ukbb/ukbb.samples.passing.qc.relevant.pheno.plus.diagnoses.180509.txt",
                 stringsAsFactors = F, header = T)
pheno$age_squared <- pheno$age_assessment * pheno$age_assessment
pheno$assessment_centre <- factor(pheno$assessment_centre)

smoking <- read.table("../ukbb_cases_for_logistic_regression_180620.txt", stringsAsFactors = F, header = T, sep = " ")
smoking <- smoking[, c("ID", "smoker_cases")]
smoking$smoker_cases <- factor(smoking$smoker_cases)

pheno <- merge(pheno, smoking, by.x = "ID", by.y = "ID")

#Make blood pressure columns - mean of baseline
pheno$sbp <- apply(pheno[, c("X4080.0_0", "X4080.0_1")], 1, function(x) {mean(x, na.rm = T)})
pheno$dbp <- apply(pheno[, c("X4079.0_0", "X4079.0_1")], 1, function(x) {mean(x, na.rm = T)})

#Increase the blood pressure if on blood-pressure lowering medications at baseline:
pheno[apply(pheno[, grep("X6177.0_|X6153.0_", names(pheno))], 1, 
	function(x) {any(x == 2, na.rm = T)}), "sbp"] <- pheno[apply(pheno[, grep("X6177.0_|X6153.0_", 
	names(pheno))], 1, function(x) {any(x == 2, na.rm = T)}), "sbp"] + 15

pheno[apply(pheno[, grep("X6177.0_|X6153.0_", names(pheno))], 1,
        function(x) {any(x == 2, na.rm = T)}), "dbp"] <- pheno[apply(pheno[, grep("X6177.0_|X6153.0_",
        names(pheno))], 1, function(x) {any(x == 2, na.rm = T)}), "dbp"] + 10

sink(pheno.qc.count.file, append = F)
cat(paste("Samples in original phenotype file: ", nrow(pheno), "\n", sep = ""))
sink()

#Clean the phenotype file:
pheno <- qc.pregnant(pheno, pheno.qc.count.file)

sink(pheno.qc.count.file, append = T)
cat(paste("Samples in cleaned phenotype file for blood pressure: ", nrow(pheno), "\n", sep = ""))
sink()

#Create the vectors to loop over
snp_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig", 
		"whradjbmi.eur.comb.pulit.sig", 
		"bmi.eur.men.pulit.sig", "whr.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.sig", 
		"bmi.eur.women.pulit.sig", "whr.eur.women.pulit.sig", 
                "whradjbmi.eur.women.pulit.sig", 
		"bmi.eur.comb.giukbb.sig", "whr.eur.comb.giukbb.sig", "whradjbmi.eur.comb.giukbb.sig", 
		"bmi.eur.men.giukbb.sig", "whr.eur.men.giukbb.sig", "whradjbmi.eur.men.giukbb.sig", 
		"bmi.eur.women.giukbb.sig", "whr.eur.women.giukbb.sig", "whradjbmi.eur.women.giukbb.sig")
eth_groups <- c("all.white", "brit.irish")
sex_groups <- c("comb", "women", "men")
traits <- c("sbp", "dbp")

results_df_bp <- linreg(pheno, snp_groups, eth_groups, sex_groups, traits)

#Remove those with a BMI below 15 as well for anthropometric traits
pheno <- qc.below.15.bmi(pheno, pheno.qc.count.file)

sink(pheno.qc.count.file, append = T)
cat(paste("Samples in cleaned phenotype file for anthropometric traits: ", nrow(pheno), "\n", sep = ""))
sink()
	
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
                "bmi.eur.women.giukbb.sig", "whr.eur.women.giukbb.sig", "whradjbmi.eur.women.giukbb.sig", 
	        "bmi.eur.comb.internal.unweighted", "bmi.eur.women.internal.unweighted", "bmi.eur.men.internal.unweighted",
                "whr.eur.comb.internal.unweighted", "whr.eur.women.internal.unweighted", "whr.eur.men.internal.unweighted",
                "whradjbmi.eur.comb.internal.unweighted", "whradjbmi.eur.women.internal.unweighted", "whradjbmi.eur.men.internal.unweighted")

traits <- c("whr", "wc", "hc", "bmi", "res_whr_inv")

results_df_anthro <- linreg(pheno, snp_groups, eth_groups, sex_groups, traits)

results_df <- rbind(results_df_bp, results_df_anthro)

write.table(results_df, "../results.linear.regressions.180514/anthro.results.table.180521.EXTRA.txt",
                row.names = F, sep = "\t", quote = F)

sink("../results.linear.regressions.180514/anthro.results.table.180521.WARNINGS.txt", append = F)
print(warnings())
sink()

comb_snp_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig", 
			"whradjbmi.eur.comb.pulit.sig", "bmi.eur.comb.giukbb.sig", 
			"whr.eur.comb.giukbb.sig", "whradjbmi.eur.comb.giukbb.sig")

women_snp_groups <- snp_groups[grep("\\.women", snp_groups)]
men_snp_groups <- snp_groups[grep("\\.men", snp_groups)]

women_snp_groups <- women_snp_groups[order(women_snp_groups)]
men_snp_groups <- men_snp_groups[order(men_snp_groups)]

#Prepare some columns required for Heterogeneity calculations
results_df$cochran_weights <- 1 / (results_df$grs_se^2)

for (i in 1:length(women_snp_groups)) {
        for (eth_group in unique(results_df$eth_group)) {
                for (trait in unique(results_df$trait)) {
                        for (grs_unit in unique(results_df$grs_unit)) {
				for (trait_unit in unique(results_df$trait_unit)) {
					general <- results_df$eth_group == eth_group &
						results_df$trait == trait &
                                                results_df$grs_unit == grs_unit &
						results_df$trait_unit == trait_unit

					if (i <= length(comb_snp_groups) & 
						nrow(results_df[general & results_df$snp_group == comb_snp_groups[i] &
						results_df$sex_group %in% c("women", "men"), ]) == 2) {
                                	
						results_df$cochrans_q[(general & results_df$snp_group == comb_snp_groups[i] & 
							results_df$sex_group %in% c("women", "men"))] <- as.numeric(cochran.Q(c(results_df[general &
							results_df$snp_group == comb_snp_groups[i] & 
							results_df$sex_group == "women",
                                        		"grs_beta"],
                                        		results_df[general & 
							results_df$snp_group == comb_snp_groups[i] & 
							results_df$sex_group == "men",
                                        		"grs_beta"]),
                                        		c(results_df[general & 
							results_df$snp_group == comb_snp_groups[i] & 
							results_df$sex_group == "women",
                                        		"cochran_weights"],
                                        		results_df[general & 
							results_df$snp_group == comb_snp_groups[i] & 
							results_df$sex_group == "men",
                                        		"cochran_weights"]))[1])

	                                	results_df$cochrans_p[(general &
							results_df$snp_group == comb_snp_groups[i] & 
							results_df$sex_group %in% c("women", "men"))] <- as.numeric(cochran.Q(c(results_df[general & 
							results_df$snp_group == comb_snp_groups[i] & 
							results_df$sex_group == "women",
        	                       	 	        "grs_beta"],
                	                        	results_df[general & 
							results_df$snp_group == comb_snp_groups[i] & 
							results_df$sex_group == "men",
                        	                	"grs_beta"]),
                                	        	c(results_df[general &
							results_df$snp_group == comb_snp_groups[i] & 
							results_df$sex_group == "women",
                                        		"cochran_weights"],
                                        		results_df[general &
							results_df$snp_group == comb_snp_groups[i] & 
							results_df$sex_group == "men",
                                        		"cochran_weights"]))[2])
					
					}
	                 
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
        	        	                        "grs_beta"],
                	        	                results_df[general &
							results_df$snp_group == men_snp_groups[i] & 
							results_df$sex_group == "men",
	                                	        "grs_beta"]),
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
        		                                "grs_beta"],
                		                        results_df[general & 
							results_df$snp_group == men_snp_groups[i] & 
							results_df$sex_group == "men",
                                        		"grs_beta"]),
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
        }
}

results_df$cochrans_i2 <- (results_df$cochrans_q - 1)/results_df$cochrans_q

results_df$cochrans_i2[results_df$cochrans_i2 < 0] <- 0

write.table(results_df, "../results.linear.regressions.180514/anthro.results.table.180521.txt", 
		row.names = F, sep = "\t", quote = F)




