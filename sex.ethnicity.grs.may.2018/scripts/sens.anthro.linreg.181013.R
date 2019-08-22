#!/bin/env Rscript
#$ -cwd

##########################################################
###Sex and ethnicity WHR GRS on anthropometric traits ####
############ Jenny Censin, 2018-05-20 ###################

library(mada)

qc.pregnant <- function(data) { 
cleaned <- subset(data, is.na(data$pregnant) | data$pregnant == 0)
return(cleaned)
}

qc.below.15.bmi <- function(data) {
cleaned <- subset(data, is.na(data$bmi) | data$bmi > 15)
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
                        	whrfit <- lm(whr ~ age_assessment + age_squared + dummy_sex + bmi + assessment_centre, 
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
						for (extra_adjustment in c("-")) {
							for (trait_unit in c("clin", "sd", "normal_sd")) {
								
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
								
								trait_unit_df$normal_sd <- (grs_unit_df[[trait]]- mean(grs_unit_df[[trait]], 
									na.rm = T))/sd(grs_unit_df[[trait]], na.rm = T)

								grs_trait <- gsub("\\.(.)*", "", snp_group)
								if (grs_trait == "whradjbmi") {
									grs_trait <- "res_whr_inv"
								}									
	
								if (grs_trait != trait) {
									next()
								}

        	                                                grs_sex <- gsub("(.)+\\.eur\\.", "", snp_group)
	                                                        grs_sex <- gsub("\\.(.)*", "", grs_sex)
	
        	                                                if (grs_sex != "comb" & grs_sex != sex_group) {
                	                                                next()
                        	                                }

								if (trait == "res_whr_inv" & extra_adjustment == "-") {
                        		                		linreg.form <- as.formula(paste(trait, 
									"~SCORESUM+dummy_array+", pcs, sep = " "))
                                				} else if (trait_unit == "sd" & extra_adjustment == "-") {
									linreg.form <- as.formula(paste("res_rint_inv~SCORESUM+dummy_array+", pcs, sep = " "))
								} else if (trait_unit %in% c("clin") & extra_adjustment == "-") {
                               		         			linreg.form <- as.formula(paste(trait, 
									"~SCORESUM+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+", 
									pcs, sep = " "))
								} else if (trait_unit %in% c("normal_sd") & extra_adjustment == "-") {
                                                                        linreg.form <- as.formula(paste("normal_sd",
                                                                        "~SCORESUM+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
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
                                                        	} else if (trait_unit %in% c("clin") & extra_adjustment == "-") {
                                                        	        linreg.form2 <- as.formula(paste(trait,
                                                                	"~dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                                                	pcs, sep = " "))
                                                        	} else if (trait_unit %in% c("normal_sd") & extra_adjustment == "-") {
                                                                        linreg.form2 <- as.formula(paste("normal_sd",
                                                                        "~dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
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

#Read in phenotype file:
pheno <- read.table("/well/lindgren/jc/ukbb/ukbb.samples.passing.qc.relevant.pheno.plus.diagnoses.180509.txt",
                 stringsAsFactors = F, header = T)
pheno$age_squared <- pheno$age_assessment * pheno$age_assessment
pheno$assessment_centre <- factor(pheno$assessment_centre)

#Clean the phenotype file:
pheno <- qc.pregnant(pheno)
pheno <- qc.below.15.bmi(pheno)

#Create the vectors to loop over
eth_groups <- c("all.white")
sex_groups <- c("comb", "women", "men")
traits <- c("bmi", "whr", "res_whr_inv")
types_internal <- c("normal", "half", "randomnoise", "systematicnoise", "unweighted")
types_fdr <- c("0.01.fdr", "0.05.fdr", "0.1.fdr")
x <- paste(rep(c("bmi", "whr", "whradjbmi"), each = length(sex_groups)), "eur", sex_groups, sep = ".")
internal <- paste(rep(x, each = length(types_internal)), "internal", types_internal, sep = ".")
y <- paste(rep(c("bmi", "whr", "whradjbmi"), each = length(c("men", "women"))), "eur", c("men", "women"), sep = ".")
fdr <- paste(rep(y, each = length(types_fdr)), types_fdr, sep = ".")
snp_groups <- c(internal, fdr)

results_df <- linreg(pheno, snp_groups, eth_groups, sex_groups, traits)

write.table(results_df, "../results.linear.regressions.180514/sens.anthro.results.table.181013.EXTRA.txt",
                row.names = F, sep = "\t", quote = F)

sink("../results.linear.regressions.180514/sens.anthro.results.table.181013.WARNINGS.txt", append = F)
print(warnings())
sink()

comb_snp_groups <- snp_groups[grep("\\.comb\\.", snp_groups)]
women_snp_groups <- snp_groups[grep("\\.women\\.", snp_groups)]
men_snp_groups <- snp_groups[grep("\\.men\\.", snp_groups)]

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

write.table(results_df, "../results.linear.regressions.180514/sens.anthro.results.table.181013.txt", 
		row.names = F, sep = "\t", quote = F)




