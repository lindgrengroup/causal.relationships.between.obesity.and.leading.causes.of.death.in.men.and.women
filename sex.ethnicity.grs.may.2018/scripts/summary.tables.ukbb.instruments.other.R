#!/bin/env Rscript
#$ -cwd

##########################################################
#For describing UKBB characteristics
#########################################################

df <- read.table("../ukbb_cases_for_logistic_regression_180620.txt", 
	stringsAsFactors = F, header = T)

df <- df[df$check_se %in% c("brit", "recode.white", "white"), ]
columns_to_include <- c("ID", "Submitted_Gender", "age_assessment", "check_se",
                        "bmi", "wc", "hc", "whr",
                        "genotyping_array", 
			"t2d_cases_prob", "t2d_cases_probposs", "t1d_cases_prob", 
			"t1d_cases_probposs", "cad_cases", "stroke_cases", 
			"copd_cases", "dementia_cases", "renal_failure_cases", 
			"lungcancer_cases", "colorectal_cancer_cases", 
			"breast_cancer_cases", "female_infertility_cases", 
			"male_infertility_cases", "any_infertility_cases", 
			"nafld_cases", "cld_cases", "aki_cases", "ckd_cases", 
			"haem_stroke_cases", "isch_stroke_cases")
df <- df[, columns_to_include]

#To get blood pressure as well... 
pheno <- read.table("/well/lindgren/jc/ukbb/ukbb.samples.passing.qc.relevant.pheno.plus.diagnoses.180509.txt",
                 stringsAsFactors = F, header = T)

pheno$sbp <- apply(pheno[, c("X4080.0_0", "X4080.0_1")], 1, function(x) {mean(x, na.rm = T)})
pheno$dbp <- apply(pheno[, c("X4079.0_0", "X4079.0_1")], 1, function(x) {mean(x, na.rm = T)})

pheno[apply(pheno[, grep("X6177.0_|X6153.0_", names(pheno))], 1,
        function(x) {any(x == 2, na.rm = T)}), "sbp"] <- pheno[apply(pheno[, grep("X6177.0_|X6153.0_",
        names(pheno))], 1, function(x) {any(x == 2, na.rm = T)}), "sbp"] + 15

pheno[apply(pheno[, grep("X6177.0_|X6153.0_", names(pheno))], 1,
        function(x) {any(x == 2, na.rm = T)}), "dbp"] <- pheno[apply(pheno[, grep("X6177.0_|X6153.0_",
        names(pheno))], 1, function(x) {any(x == 2, na.rm = T)}), "dbp"] + 10
pheno <- pheno[, c("ID", "sbp", "dbp")]

#Merge the two files
df <- merge(df, pheno, by.x = "ID", by.y = "ID")
df$ID <- NULL

results_df <- c("Individuals, N (%)", "British, N (%)", "Age, mean (SD), years", 
		"UK BiLEVE array, N (%)", 
		"Body mass index, mean (SD), kg/m2", 
		"Waist circumference, mean (SD), cm", 
		"Hip circumference, mean (SD), cm", 
		"Waist-hip-ratio, mean (SD)", 
		"Systolic blood pressure, mean (SD), mmHg", 
		"Diastolic blood pressure, mean (SD), mmHg", 
		"Type 2 diabetes cases, N (%)",
		"Coronary artery disease cases, N (%)", 
		"Breast cancer cases, N (%)",
		"Chronic liver disease cases, N (%)",
		"Colorectal cancer cases, N (%)",
		"COPD cases, N (%)", 
		"Dementia cases, N (%)",
		"Infertility cases, N (%)",
		"Lung cancer cases, N (%)",		
		"NAFLD cases, N (%)",
		"Renal failure cases, N (%)",
		"Renal failure, acute, cases, N (%)",
		"Renal failure, chronic, cases, N (%)",
		"Stroke cases, N (%)",
		"Stroke, hemorrhagic, cases, N (%)",
		"Stroke, ischaemic, cases, N (%)",
		"Type 1 diabetes cases, N (%)")

#Can add "combined"
for (sex_group in c("men", "women")) {
	if (sex_group == "men") {
		sex_df <- df[df$Submitted_Gender == "M", ]
	} else if (sex_group == "women") {
		sex_df <- df[df$Submitted_Gender == "F", ]
	} else if (sex_group == "combined") {
		sex_df <- df
	}

	column <- paste(nrow(sex_df), " (", round(nrow(sex_df)/nrow(df)*100, 1), ")", sep = "")
	column <- c(column, paste(nrow(sex_df[sex_df$check_se == "brit", ]), " (", 
			round(nrow(sex_df[sex_df$check_se == "brit", ])/nrow(sex_df) * 100, 1), 
			")", sep = ""))
	column <- c(column, paste(round(mean(sex_df$age_assessment, na.rm = T), 1), " (", 
			round(sd(sex_df$age_assessment, na.rm = T), 1), ")", sep = ""))
	column <- c(column, paste(nrow(sex_df[sex_df$genotyping_array == "UKBL", ]), " (", 
			formatC(nrow(sex_df[sex_df$genotyping_array == "UKBL", ])/nrow(sex_df) * 100, 1, format = "f"),
			")", sep = ""))
	column <- c(column, paste(round(mean(sex_df$bmi, na.rm = T), 1), " (", 
			round(sd(sex_df$bmi, na.rm = T), 1), ")", sep = ""))
	column <- c(column, paste(round(mean(sex_df$wc, na.rm = T), 1), " (", 
			round(sd(sex_df$wc, na.rm = T), 1), ")", sep = ""))
	column <- c(column, paste(round(mean(sex_df$hc, na.rm = T), 1), " (", 
			round(sd(sex_df$hc, na.rm = T), 1), ")", sep = ""))
	column <- c(column, paste(round(mean(sex_df$whr, na.rm = T), 2), " (", 
			round(sd(sex_df$whr, na.rm = T), 2), ")", sep = ""))
	column <- c(column, paste(round(mean(sex_df$sbp, na.rm = T), 1), " (",
                        round(sd(sex_df$sbp, na.rm = T), 1), ")", sep = ""))
	column <- c(column, paste(round(mean(sex_df$dbp, na.rm = T), 1), " (",
                        formatC(sd(sex_df$dbp, na.rm = T), 1, format = "f"), ")", sep = ""))

	traits <- c("t2d_cases_probposs", "cad_cases", "breast_cancer_cases", 
			"cld_cases", "colorectal_cancer_cases", "copd_cases", 
			"dementia_cases", "any_infertility_cases", "lungcancer_cases", 
			"nafld_cases", "renal_failure_cases", "aki_cases", "ckd_cases", 
			"stroke_cases", "haem_stroke_cases", "isch_stroke_cases",
			"t1d_cases_probposs")
	for (trait in traits) {
		column <- c(column, ifelse(nrow(sex_df[sex_df[[trait]] == 1, ]) == 0, paste("-", sep =""), 
			paste(nrow(sex_df[sex_df[[trait]] == 1, ]), " (", 
			formatC(nrow(sex_df[sex_df[[trait]] == 1, ])/nrow(sex_df) * 100, 1, format = "f"), 
			")", sep = "")))
	}

	results_df <- cbind(results_df, column)
}

#Can add combined
results_df <- data.frame(results_df)
colnames(results_df) <- c("Characteristic", "Men", "Women")

write.table(results_df, "/Net/fs1/home/linc4222/summary.characteristics.ukbb.table.txt", sep = "\t", 
		quote = F, row.names = F)

########################################################################
######## SNP GROUP DICTIONARY FOR ALL TABLES ############################
###################################################################

dict_snp_groups_long <- list(bmi.eur.comb.pulit.sig = "BMI,  combined estimates, SNPs = 565",
                         bmi.eur.men.pulit = "BMI, primary SNPs in males, SNPs = 220",
                         bmi.eur.men.pulit.phet = "BMI, male P-heterogeneity Bonferroni, SNPs = 565",
                         bmi.eur.men.pulit.sig = "BMI,  male sex-specific estimates, SNPs = 565",
                         bmi.eur.women.pulit = "BMI, primary SNPs in females, SNPs = 278",
                         bmi.eur.women.pulit.phet = "BMI, female P-heterogeneity Bonferroni, SNPs = 565",
                         bmi.eur.women.pulit.sig = "BMI,  female sex-specific estimates, SNPs = 565",
                         whr.eur.comb.pulit.sig = "WHR,  combined estimates, SNPs = 324",
                         whr.eur.men.pulit = "WHR, primary SNPs in males, SNPs = 77",
                         whr.eur.men.pulit.phet = "WHR, male P-heterogeneity Bonferroni, SNPs = 324",
                         whr.eur.men.pulit.sig = "WHR,  male sex-specific estimates, SNPs = 324",
                         whr.eur.women.pulit = "WHR, primary SNPs in females, SNPs = 203",
                         whr.eur.women.pulit.phet = "WHR, female P-heterogeneity Bonferroni, SNPs = 324",
                         whr.eur.women.pulit.sig = "WHR,  female sex-specific estimates, SNPs = 324",
                         whradjbmi.eur.comb.pulit.sig = "WHRadjBMI,  combined estimates, SNPs = 337",
                         whradjbmi.eur.men.pulit = "WHRadjBMI, primary SNPs in males, SNPs = 90",
                         whradjbmi.eur.men.pulit.phet = "WHRadjBMI, male P-heterogeneity Bonferroni, SNPs = 337",
                         whradjbmi.eur.men.pulit.sig = "WHRadjBMI,  male sex-specific estimates, SNPs = 337",
                         whradjbmi.eur.women.pulit = "WHRadjBMI, primary SNPs in females, SNPs = 264",
                         whradjbmi.eur.women.pulit.phet = "WHRadjBMI, female P-heterogeneity Bonferroni, SNPs = 337",
                         whradjbmi.eur.women.pulit.sig = "WHRadjBMI,  female sex-specific estimates, SNPs = 337",
                         bmi.eur.men.0.01.fdr = "BMI, male P-heterogeneity FDR 1%, SNPs = 565",
                         bmi.eur.men.0.05.fdr = "BMI, male P-heterogeneity FDR 5%, SNPs = 565",
                         bmi.eur.men.0.1.fdr = "BMI, male P-heterogeneity FDR10%, SNPs = 565",
                         bmi.eur.women.0.01.fdr = "BMI, female P-heterogeneity FDR 1%, SNPs = 565",
                         bmi.eur.women.0.05.fdr = "BMI, female P-heterogeneity FDR 5%, SNPs = 565",
                         bmi.eur.women.0.1.fdr = "BMI, female P-heterogeneity FDR10%, SNPs = 565",
                         whr.eur.men.0.01.fdr = "WHR, male P-heterogeneity FDR 1%, SNPs = 324",
                         whr.eur.men.0.05.fdr = "WHR, male P-heterogeneity FDR 5%, SNPs = 324",
                         whr.eur.men.0.1.fdr = "WHR, male P-heterogeneity FDR10%, SNPs = 324",
                         whr.eur.women.0.01.fdr = "WHR, female P-heterogeneity FDR 1%, SNPs = 324",
                         whr.eur.women.0.05.fdr = "WHR, female P-heterogeneity FDR 5%, SNPs = 324",
                         whr.eur.women.0.1.fdr = "WHR, female P-heterogeneity FDR10%, SNPs = 324",
                         whradjbmi.eur.men.0.01.fdr = "WHRadjBMI, male P-heterogeneity FDR 1%, SNPs = 337",
                         whradjbmi.eur.men.0.05.fdr = "WHRadjBMI, male P-heterogeneity FDR 5%, SNPs = 337",
                         whradjbmi.eur.men.0.1.fdr = "WHRadjBMI, male P-heterogeneity FDR10%, SNPs = 337",
                         whradjbmi.eur.women.0.01.fdr = "WHRadjBMI, female P-heterogeneity FDR 1%, SNPs = 337",
                         whradjbmi.eur.women.0.05.fdr = "WHRadjBMI, female P-heterogeneity FDR 5%, SNPs = 337",
                         whradjbmi.eur.women.0.1.fdr = "WHRadjBMI, female P-heterogeneity FDR10%, SNPs = 337",
                         bmi.eur.comb.giukbb.sig = "BMI, combined tgiant weights, SNPs = 478",
                         bmi.eur.women.giukbb.sig = "BMI, female tgiant weights, SNPs = 478",
                         whr.eur.men.giukbb.sig = "WHR, male tgiant weights, SNPs = 253",
                         whradjbmi.eur.comb.giukbb.sig = "WHRadjBMI, combined tgiant weights, SNPs = 275",
                         whradjbmi.eur.women.giukbb.sig = "WHRadjBMI, female tgiant weights, SNPs = 263",
                         bmi.eur.men.giukbb.sig = "BMI, male giant weights, SNPs = 478",
                         whr.eur.comb.giukbb.sig = "WHR, combined tgiant weights, SNPs = 264",
                         whr.eur.women.giukbb.sig = "WHR, female tgiant weights, SNPs = 258",
                         whradjbmi.eur.men.giukbb.sig = "WHRadjBMI, male tgiant weights, SNPs = 248", 
			 bmi.eur.comb.internal.unweighted = "BMI, unweighted, SNPs = 565", 
			 bmi.eur.women.internal.unweighted = "BMI, unweighted, SNPs = 565", 
			 bmi.eur.men.internal.unweighted = "BMI, unweighted, SNPs = 565", 
			 whr.eur.comb.internal.unweighted = "WHR, unweighted, SNPs = 324", 
			 whr.eur.women.internal.unweighted = "WHR, unweighted, SNPs = 324", 
		 	whr.eur.men.internal.unweighted = "WHR, unweighted, SNPs = 324", 
			 whradjbmi.eur.comb.internal.unweighted = "WHRadjBMI, unweighted, SNPs = 337", 
			whradjbmi.eur.women.internal.unweighted = "WHRadjBMI, unweighted, SNPs = 337",
			whradjbmi.eur.men.internal.unweighted = "WHRadjBMI, unweighted, SNPs = 337")

dict_snp_groups_short <- list(bmi.eur.comb.pulit.sig = "BMI, combined, SNPs = 565",
                         bmi.eur.men.pulit.phet = "BMI, men, SNPs = 565",
                         bmi.eur.men.pulit.sig = "BMI, men, SNPs = 565",
                         bmi.eur.women.pulit.phet = "BMI, women, SNPs = 565",
                         bmi.eur.women.pulit.sig = "BMI, women, SNPs = 565",
                         whr.eur.comb.pulit.sig = "WHR, combined, SNPs = 324",
                         whr.eur.men.pulit.phet = "WHR, men, SNPs = 324",
                         whr.eur.men.pulit.sig = "WHR, men, SNPs = 324",
                         whr.eur.women.pulit.phet = "WHR, women, SNPs = 324",
                         whr.eur.women.pulit.sig = "WHR, women, SNPs = 324",
                         whradjbmi.eur.comb.pulit.sig = "WHRadjBMI, combined, SNPs = 337",
                         whradjbmi.eur.men.pulit.phet = "WHRadjBMI, men, SNPs = 337",
                         whradjbmi.eur.men.pulit.sig = "WHRadjBMI, men, SNPs = 337",
                         whradjbmi.eur.women.pulit.phet = "WHRadjBMI, women, SNPs = 337",
                         whradjbmi.eur.women.pulit.sig = "WHRadjBMI, women, SNPs = 337")

######################################################################
#For describing instrument strength: 
#######################################################################

df_raw <- read.table("../results.linear.regressions.180514/anthro.results.table.180521.txt", 
		stringsAsFactors = F, header = T)

pulit_snp_groups <- unique(df_raw$snp_group[grep("pulit|fdr", df_raw$snp_group)])
eth_groups <- c("all.white", "brit.irish")

for (dataset in c("pulit")) {
	df_dataset <- df_raw[df_raw$snp_group %in% pulit_snp_groups, ]	
	for (eth_group in eth_groups) {
		df <- df_dataset
		df <- df[(df$grs_unit == "raw_scoresum" & df$eth_group == eth_group & df$extra_adjustment == "-") &
			((df$snp_group %in% c("bmi.eur.comb.pulit.sig", "bmi.eur.comb.giukbb.sig") &
			df$sex_group %in% c("comb", "men", "women") & 
			df$trait == "bmi") |
			(df$snp_group %in% c("bmi.eur.men.pulit", "bmi.eur.men.pulit.phet", "bmi.eur.men.pulit.sig", 
					"bmi.eur.men.0.01.fdr", "bmi.eur.men.0.05.fdr", "bmi.eur.men.0.1.fdr", "bmi.eur.men.giukbb.sig") &
			df$sex_group == "men" &
			df$trait == "bmi") |
			(df$snp_group %in% c("bmi.eur.women.pulit", "bmi.eur.women.pulit.phet", "bmi.eur.women.pulit.sig", 
					"bmi.eur.women.0.01.fdr", "bmi.eur.women.0.05.fdr", "bmi.eur.women.0.1.fdr", "bmi.eur.women.giukbb.sig") &
			df$sex_group == "women" &
			df$trait == "bmi") |
			(df$snp_group %in% c("whr.eur.comb.pulit.sig", "whr.eur.comb.giukbb.sig") &
			df$sex_group %in% c("comb", "men", "women") &
			df$trait == "whr") |
			(df$snp_group %in% c("whr.eur.men.pulit", "whr.eur.men.pulit.phet", "whr.eur.men.pulit.sig", 
					"whr.eur.men.0.01.fdr", "whr.eur.men.0.05.fdr", "whr.eur.men.0.1.fdr", "whr.eur.men.giukbb.sig") &
			df$sex_group == "men" &
			df$trait == "whr") |
			(df$snp_group %in% c("whr.eur.women.pulit", "whr.eur.women.pulit.phet", "whr.eur.women.pulit.sig", 
					"whr.eur.women.0.01.fdr", "whr.eur.women.0.05.fdr", "whr.eur.women.0.1.fdr", "whr.eur.women.0.1.fdr", "whr.eur.women.giukbb.sig") &
			df$sex_group == "women" &
			df$trait == "whr") |
			(df$snp_group %in% c("whradjbmi.eur.comb.pulit.sig", "whradjbmi.eur.comb.giukbb.sig") &
			df$sex_group %in% c("comb", "men", "women") &
			df$trait == "res_whr_inv") |
			(df$snp_group %in% c("whradjbmi.eur.men.pulit", "whradjbmi.eur.men.pulit.phet", "whradjbmi.eur.men.pulit.sig", 
				"whradjbmi.eur.men.0.01.fdr", "whradjbmi.eur.men.0.05.fdr", 
				"whradjbmi.eur.men.0.1.fdr", "whradjbmi.eur.men.giukbb.sig") &
			df$sex_group == "men" &
			df$trait == "res_whr_inv") |
			(df$snp_group %in% c("whradjbmi.eur.women.pulit", "whradjbmi.eur.women.pulit.phet", "whradjbmi.eur.women.pulit.sig", 
				"whradjbmi.eur.women.0.01.fdr", "whradjbmi.eur.women.0.05.fdr", 
				"whradjbmi.eur.women.0.1.fdr", "whradjbmi.eur.women.giukbb.sig") &
			df$sex_group == "women" &
			df$trait == "res_whr_inv")), ]


		for (i in 1:length(dict_snp_groups_long)) {
			df$instrument <- as.character(replace(df$instrument,
        			df$snp_group == names(dict_snp_groups_long[i]), dict_snp_groups_long[i]))
		}

		df_out <- strsplit(df$instrument, ", SNPs = ")
		df <- data.frame(df, do.call(rbind, df_out))
		df$instrument <- df$X1
		df$n_snps <- df$X2

		df$sex_group <- ifelse(df$sex_group == "comb", "Combined", 
			ifelse(df$sex_group == "men", "Men", 
			"Women"))
		df$trait <- ifelse(df$trait == "bmi", "BMI", 
			ifelse(df$trait == "whr", "WHR", 
			"WHRadjBMI"))

		df$ci <- paste(formatC(df$grs_beta, format = "f", 2), " (",
                                formatC(df$grs_lci, format = "f", 2), ",", formatC(df$grs_uci, format = "f", 2), ")", sep = "")

		#Subset to relevant columns
		df_clin <- df[df$trait_unit == "clin", c("instrument", "n_snps", "sex_group", "trait", "grs_p",
                        "f_value", "grs_r2", "ci", "cochrans_p")]
		df_sd <- df[df$trait_unit == "sd", c("instrument", "n_snps", "sex_group", "trait", "grs_p",
                        "f_value", "grs_r2", "ci", "cochrans_p")]

		df <- merge(df_clin, df_sd, by = c("instrument", "n_snps", "sex_group", "trait"), all = T)

		df$sex_group <- factor(df$sex_group, levels = c("Combined", "Men", "Women"))

		df[, c("grs_p.x", "cochrans_p.x", "grs_p.y", "cochrans_p.y")][df[, c("grs_p.x", 
			"cochrans_p.x", "grs_p.y", "cochrans_p.y")] < 1e-200] <- 0

		df[, c("grs_p.x", "cochrans_p.x", "grs_p.y", "cochrans_p.y")] <- lapply(df[, 
			c("grs_p.x", "cochrans_p.x", "grs_p.y", "cochrans_p.y")], function(x) ifelse(x < 0.001, formatC(x, format = "e", 2),
                                ifelse(x < 0.01, formatC(x, format = "f", 3), formatC(x, format = "f", 2))))

		df[, c("grs_p.x", "cochrans_p.x", "grs_p.y", "cochrans_p.y")][df[, c("grs_p.x", 
			"cochrans_p.x", "grs_p.y", "cochrans_p.y")] == "0.00e+00"] <- "<1e-200"

		df[, c("f_value.x", "f_value.y")] <- lapply(df[, c("f_value.x", "f_value.y")], function(x) round(x, 0))

		df[, c("grs_r2.x", "grs_r2.y")] <- lapply(df[, c("grs_r2.x", "grs_r2.y")], function(x) paste(formatC(x * 100, format = "f", 2), "%", sep = ""))

		#For the linear regressions and WHRadjBMI, SD and clin are the same of course, so put NA
		df[df$trait == "WHRadjBMI", c("grs_p.x", "f_value.x", "grs_r2.x", "ci.x", "cochrans_p.x")] <- "NA"
	
		colnames(df) <- c("GRS", "N SNPs", "Sex-Strata", "Trait",  
			"P, Clin", "F, Clin", "R2, Clin", "Estimate (95% CI), Clin", 
			"Pheterogeneity, Clin", "P, SD", "F, SD", "R2, SD", "Estimate (95% CI), SD",
                        "Pheterogeneity, SD")

		change_column_names <- c("P, Clin", "F, Clin", "R2, Clin", "Estimate (95% CI), Clin",
                        "Pheterogeneity, Clin", "P, SD", "F, SD", "R2, SD", "Estimate (95% CI), SD",
                        "Pheterogeneity, SD")
                eth_group_name <- ifelse(eth_group == "all.white", " - Europeans", " - British")

		colnames(df)[which(names(df) %in% change_column_names)] <- paste(change_column_names, eth_group_name, sep = "")

		df[] <- lapply(df, gsub, pattern = "NA", replacement = "-")
	
		if (eth_group == "all.white") {
                                all_white_df <- df
                } else if (eth_group == "brit.irish") {
                                brit_df <- df
                }

	}
	final_df <- merge(all_white_df, brit_df, by = c("GRS", "N SNPs", "Sex-Strata", "Trait"), all = T)
	final_df <- final_df[order(final_df$Trait, final_df[, "Sex-Strata"], tolower(final_df$GRS)), ]
		
	output_file <- paste("/Net/fs1/home/linc4222/instrument.strength.in.ukbb.table.", dataset, ".180805.txt", sep = "")
	write.table(final_df, output_file, sep = "\t", 
			quote = F, row.names = F, na = "-")

}

########################################################################
#For getting the linear regression results, raw scoresum and both clin
#and SD units in output, with and without adjusting for smoking
########################################################################

df_raw <- read.table("../results.linear.regressions.180514/anthro.results.table.180521.txt",
                stringsAsFactors = F, header = T)

pulit_snp_groups <- unique(df_raw$snp_group[grep("pulit\\.sig", df_raw$snp_group)])

eth_groups <- c("all.white", "brit.irish")
datasets <- c("pulit")

for (extra_adjustment in c("-", "smoking")) {
	for (dataset in datasets) {
	        for (eth_group in eth_groups) {
               		df_dataset <- df_raw[df_raw$snp_group %in% pulit_snp_groups, ]
			df <- df_dataset
			df <- df[(df$grs_unit == "raw_scoresum" & df$eth_group == eth_group & df$extra_adjustment == extra_adjustment &
				df$trait %in% c("sbp", "sbp_smoking", "dbp", "dbp_smoking")) &
		                ((df$snp_group %in% c("bmi.eur.comb.pulit.sig") &
		                df$sex_group %in% c("comb")) |
                		(df$snp_group %in% c("bmi.eur.men.pulit.sig") &
                		df$sex_group == "men") |
                		(df$snp_group %in% c("bmi.eur.women.pulit.sig") &
                		df$sex_group == "women") |
                		(df$snp_group %in% c("whr.eur.comb.pulit.sig") &
                		df$sex_group %in% c("comb")) |
                		(df$snp_group %in% c("whr.eur.men.pulit.sig") &
                		df$sex_group == "men") |
                		(df$snp_group %in% c("whr.eur.women.pulit.sig") &
                		df$sex_group == "women") |
                		(df$snp_group %in% c("whradjbmi.eur.comb.pulit.sig") &
                		df$sex_group %in% c("comb")) |
                		(df$snp_group %in% c("whradjbmi.eur.men.pulit.sig") &
                		df$sex_group == "men") |
                		(df$snp_group %in% c("whradjbmi.eur.women.pulit.sig") &
                		df$sex_group == "women")), ]

			for (i in 1:length(dict_snp_groups_short)) {
        			df$instrument <- as.character(replace(df$instrument,
                			df$snp_group == names(dict_snp_groups_short[i]), dict_snp_groups_short[i]))
			}

			df_out <- strsplit(df$instrument, ", SNPs = ")
			df <- data.frame(df, do.call(rbind, df_out))
			df$instrument <- df$X1
			df$n_snps <- df$X2

			df$sex_group <- ifelse(df$sex_group == "comb", "Combined",
                		        ifelse(df$sex_group == "men", "Men",
                        		"Women"))

			df$trait <- gsub("_smoking", "", df$trait)

			df$unit <- ifelse(df$trait == "bmi", "kg/m2",
                		        ifelse(df$trait == "whr", "Ratio",
                        		ifelse(df$trait == "res_whr_inv", "-",
                        		ifelse(df$trait == "wc", "cm",
                        		ifelse(df$trait == "hc", "cm",
                        		ifelse(df$trait == "sbp", "mmHg",
                        		ifelse(df$trait == "dbp", "mmHg", "MISSING")))))))

			df$trait <- ifelse(df$trait == "bmi", "BMI",
        		                ifelse(df$trait == "whr", "WHR",
                		        ifelse(df$trait == "res_whr_inv", "WHRadjBMI", 
					ifelse(df$trait == "wc", "Waist circumference", 
					ifelse(df$trait == "hc", "Hip circumference", 
					ifelse(df$trait == "sbp", "SBP", 
					ifelse(df$trait == "dbp", "DBP", "MISSING")))))))

			df_clin <- df[df$trait_unit == "clin", ]
			df_sd <- df[df$trait_unit == "sd", ]
			df <- merge(df_clin, df_sd, by = c("snp_group", "grs_unit", "eth_group", "sex_group", 
					"trait", "extra_adjustment", "grs_file", "instrument", "X1", "X2", "n_snps", "unit"), all = T)

			df$ci.x <- paste(formatC(df$grs_beta.x, format = "f", 2), " (", 
				formatC(df$grs_lci.x, format = "f", 2), ",", formatC(df$grs_uci.x, format = "f", 2), ")", sep = "")

			df$ci.y <- paste(formatC(df$grs_beta.y, format = "f", 2), " (",
        		        formatC(df$grs_lci.y, format = "f", 2), ",", formatC(df$grs_uci.y, format = "f", 2), ")", sep = "")

			df <- df[, c("trait", "instrument", "unit", "n_snps", "sex_group",
					"ci.x", "ci.y", "grs_p.x", "grs_p.y", "cochrans_p.x", "cochrans_p.y")]

			df[, c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")][df[, 
				c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")] < 1e-200] <- 0

			df[, c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")] <- lapply(df[, 
				c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")], 
				function(x) ifelse(x < 0.001, formatC(x, format = "e", 2), round(x, 3)))

			df[, c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")][df[,
	        		c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")] == "0.00e+00"] <- "<1E-200"

			colnames(df) <-  c("Outcome", "GRS", "Clinical Unit", "N SNPs", 
				"Sex-Strata", "Estimate, Clinical (95% CI)", "Estimate, SD (95% CI)", 
				"P (Clinical)", "P (SD)", "Pheterogeneity (Clinical)", "Pheterogeneity (SD)")

			df <- df[, c("GRS", "Outcome", "Clinical Unit", "N SNPs",
        		        "Sex-Strata", "Estimate, Clinical (95% CI)", "P (Clinical)", "Pheterogeneity (Clinical)", 
				"Estimate, SD (95% CI)", "P (SD)" , "Pheterogeneity (SD)")]
			
			df[df$Outcome == "WHRadjBMI", c("Estimate, Clinical (95% CI)", "Pheterogeneity (Clinical)", "P (Clinical)")] <- "-"
				
			change_column_names <- c("Estimate, Clinical (95% CI)", "P (Clinical)",
                                "Pheterogeneity (Clinical)", "Estimate, SD (95% CI)", "P (SD)", "Pheterogeneity (SD)")
			eth_group_name <- ifelse(eth_group == "all.white", " - Europeans", " - British")	

			colnames(df)[which(names(df) %in% change_column_names)] <- paste(change_column_names, eth_group_name, sep = "")

			df[] <- lapply(df, gsub, pattern = "NA", replacement = "-")
			
			if (eth_group == "all.white") {
				all_white_df <- df
			} else if (eth_group == "brit.irish") {
				brit_df <- df
			}
		}

		final_df <- merge(all_white_df, brit_df, by = c("GRS", "Outcome", "Clinical Unit", "N SNPs", "Sex-Strata"), all = T)

                final_df$grs_trait_name <- gsub(",.*", "", df$GRS)
                final_df[, "Sex-Strata"] <- factor(final_df[, "Sex-Strata"], levels = c("Combined", "Men", "Women"))
                final_df$Outcome <- factor(final_df$Outcome, levels = c("BMI", "WHR", "WHRadjBMI", "Waist circumference",
                               "Hip circumference", "SBP", "DBP"))
                final_df <- final_df[order(final_df$grs_trait_name, final_df$Outcome, final_df[, "Sex-Strata"]), ]
		final_df$grs_trait_name <- NULL

		adj <- ifelse(extra_adjustment == "-", "not.adj", ifelse(extra_adjustment == "smoking", "adj.smoking", ""))
		output_file <- paste("/Net/fs1/home/linc4222/cleaned.linear.regressions.sex.specific.for.publication.", 
			dataset, ".", adj, ".txt", sep = "")
		write.table(final_df, output_file, row.names = F, quote = F, sep = "\t", na = "-")
		
	}
}

#########################################################################
#Get a table with the logistic regressions for GRSs, with and without adjusting
#for smoking, and the smoking sensititivity analysis
#########################################################################

df_raw <- read.table("../results.logistic.regressions.180514/log.results.table.180627.txt",
                stringsAsFactors = F, header = T)

smoking_names <- unique(df_raw$case_column[grep("_smoking", df_raw$case_column)])

df_raw <- df_raw[df_raw$case_column %in% c("smoker_cases", smoking_names) & 
	!(df_raw$case_column %in% c("t1d_cases_prob_smoking", "t2d_cases_prob_smoking")), ]
eth_groups <- c("all.white", "brit.irish")

for (dataset in c("pulit")) {
	comb_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig",  "whradjbmi.eur.comb.pulit.sig")
        male_groups <- c("bmi.eur.men.pulit.sig", "whr.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.sig")
        female_groups <- c("bmi.eur.women.pulit.sig", "whr.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.sig")

	for (eth_group in eth_groups) {
		df <- df_raw
		#Subset to only keep the actual analyses - no analyses in "wrong" sex even without this extra
		#careful line of code, could just do df$snp_group %in% c(comb_groups, male_groups, female_groups))
		df <- df[(df$grs_unit == "raw_scoresum" & df$eth_group == eth_group) &
        		((df$sex_group == "comb" & df$snp_group %in% comb_groups) |
        		(df$sex_group == "men" & df$snp_group %in% male_groups) |
        		(df$sex_group == "women" & df$snp_group %in% female_groups)), ]

		for (i in 1:length(dict_snp_groups_short)) {
        		df$instrument <- as.character(replace(df$instrument,
        	        	df$snp_group == names(dict_snp_groups_short[i]), dict_snp_groups_short[i]))
		}

		df_out <- strsplit(df$instrument, ", SNPs = ")
		df <- data.frame(df, do.call(rbind, df_out))
		df$instrument <- df$X1
		df$n_snps <- df$X2

		df$sex_group <- ifelse(df$sex_group == "comb", "Combined",
                        ifelse(df$sex_group == "men", "Men",
			ifelse(df$sex_group == "women", "Women", NA)))

		df$ci <- paste(formatC(df$grs_or, format = "f", 2), " (",
                	formatC(df$grs_lci_or, format = "f", 2), ", ", formatC(df$grs_uci_or, format = "f", 2), ")", sep = "")

		df$trait <- df$case_column
		dict_traits <- list(breast_cancer_cases_smoking = "Breast Cancer",
	        	cad_cases_smoking = "CAD",
	        	colorectal_cancer_cases_smoking = "Colorectal Cancer",
	        	copd_cases_smoking = "COPD",
	        	dementia_cases_smoking = "Dementia",
	        	lungcancer_cases_smoking = "Lung Cancer",
	        	renal_failure_cases_smoking = "Renal Failure",
			aki_cases_smoking = "Renal Failure - Acute", 
			ckd_cases_smoking = "Renal Failure - Chronic",
	        	stroke_cases_smoking = "Stroke",
			isch_stroke_cases_smoking = "Stroke - Ischaemic", 
			haem_stroke_cases_smoking = "Stroke - Hemorrhagic",
	        	t2d_cases_probposs_smoking = "Type 2 Diabetes",
	        	t1d_cases_probposs_smoking = "Type 1 Diabetes",
	        	any_infertility_cases_smoking = "Infertility",	
	        	nafld_cases_smoking = "NAFLD",
	        	cld_cases_smoking = "CLD", 
			smoker_cases = "Having smoked/smoker")

		for (i in 1:length(dict_traits)) {
	        	df$trait <- as.character(replace(df$trait,
        	                df$trait == names(dict_traits[i]), dict_traits[i]))
		}

		df$grs_trait <- gsub("\\..*", "", df$snp_group)
		df$order <- ifelse(df$sex_group == "Combined", 1, 
			ifelse(df$sex_group == "Men", 2, 
			ifelse(df$sex_group == "Women", 3, NA)))

		df <- df[order(df$trait, df$grs_trait, df$order), ]

		df <- df[, c("trait", "instrument", "n_snps", "sex_group",
                        "ci", "grs_p", "cochrans_p")]

		df[, c("grs_p", "cochrans_p")][df[, c("grs_p", "cochrans_p")] < 1e-200] <- 0

		df[, c("grs_p", "cochrans_p")] <- lapply(df[, c("grs_p", "cochrans_p")],
        	        function(x) ifelse(x < 0.001, formatC(x, format = "e", 2), 
				ifelse(x < 0.01, formatC(x, format = "f", 3), formatC(x, format = "f", 2))))

		df[, c("grs_p", "cochrans_p")][df[, c("grs_p", "cochrans_p")] == "0.00e+00"] <- "<1e-200"

		colnames(df) <-  c("Outcome", "Sex-Specific Estimates Instrument", "N SNPs",
        	        "Sex-Strata", "OR (95% CI)", "P", "Pheterogeneity")

		change_column_names <- c("OR (95% CI)", "P", "Pheterogeneity")
                eth_group_name <- ifelse(eth_group == "all.white", " - Europeans", " - British")

                colnames(df)[which(names(df) %in% change_column_names)] <- paste(change_column_names, eth_group_name, sep = "")

		if (eth_group == "all.white") {
                                all_white_df <- df
                } else if (eth_group == "brit.irish") {
                                brit_df <- df
                }

	}
	final_df <-  merge(all_white_df, brit_df, by = c("Outcome", "Sex-Specific Estimates Instrument", "N SNPs", 
				"Sex-Strata"), all = T)

	binary_traits <- final_df[!(final_df$Outcome %in% c("Having smoked/smoker")), ]
	binary_traits_file <- paste("/Net/fs1/home/linc4222/associations.for.grs.binary.traits.adj.smoking.", 
		dataset, ".180819.txt", sep = "")
	write.table(binary_traits, binary_traits_file, row.names = F,
                        quote = F, sep = "\t", na = "-")

	smoker <- final_df[(final_df$Outcome %in% c("Having smoked/smoker")), ]
	smoker_file <- paste("/Net/fs1/home/linc4222/associations.for.being.smoker.grs.", dataset, 
		".180820.txt", sep = "")
	write.table(smoker, smoker_file, row.names = F,
                quote = F, sep = "\t", na = "-")
}

###############################################################################
############### MAKE TABLE OF UKBB CASE DEFINITIONS ###########################
###############################################################################

df <- read.table("../ukbb_case_definitions_n_cases_table_180620.txt", 
		stringsAsFactors = F, header = T, sep = "\t")
df <- df[df$diagnose != "smoker", ]

results_df <- data.frame("outcome" = character(), "icd10" = character(), 
			"icd9" = character(), "opcs4" = character(), 
			"ni_noncancer" = character(), 
			"ni_cancer" = character(), 
			"ni_op" = character(),
			"selfreport_doctor" = character(), 
			"description" = character(),
			"add_info" = character(), stringsAsFactors = F)

get.codes <- function(df, trait, diagnose_type) {
	output <- ifelse(length(df[df$diagnose == trait & df$case_from_diagnose_type %in% diagnose_type, "code"]) >=1,
			df[df$diagnose == trait & df$case_from_diagnose_type %in% diagnose_type, "code"], "")
	return(output)	
}		

for (trait in unique(df$diagnose)) {
	current_row <- data.frame("outcome" = trait, 
		"icd10" = get.codes(df, trait, "total_icd10"),
		"icd9" = get.codes(df, trait, "total_icd9"),
		"opcs4" = get.codes(df, trait, "total_opcs4"),
		"ni_noncancer" = get.codes(df, trait, "total_ni_non_cancer"),
		"ni_cancer" = get.codes(df, trait, "total_ni_cancer"),
		"ni_op" = get.codes(df, trait, "total_ni_operation"),
		"selfreport_doctor" = ifelse(get.codes(df, trait, "total_vascular_heart_doctor") != "", 
					get.codes(df, trait, "total_vascular_heart_doctor"), 
					get.codes(df, trait, "total_dvt_copd_doctor")),
		"description" = "",
		"add_info" = "", 
		stringsAsFactors = F) 
	
	results_df <- merge(results_df, current_row, by.x = colnames(results_df),
			by.y = colnames(current_row), all = T)

}

results_df[results_df$outcome == "breast_cancer_women_only", "description"] <- "Codes for breast cancer, including personal history codes"
results_df[results_df$outcome == "exclude_CAD", "description"] <- "Participants were excluded from the CAD-control group if they had the these codes pertaining to heart aneurysm and atherosclerotic cardiovascular disease"
results_df[results_df$outcome == "SOFT_CAD", "description"] <- "Codes for myocardial infarction, percutaneous transluminal coronary angioplasty, coronary artery bypass grafting, chronic ischemic heart disease and angina"
results_df[results_df$outcome == "COPD", "description"] <- "Codes for chronic bronchitis, emphysema, COPD, or complications specified from COPD" 
results_df[results_df$outcome == "chronic_liver_disease", "description"] <- "Codes for fibrosis, sclerosis, cirrhosis of liver and liver failure, including if caused by alcohol, toxic liver disease, biliary cirrhosis, and other causes as well as codes defining complications specified as caused by these"
results_df[results_df$outcome == "colorectal_cancer", "description"] <- "Codes for cancers from caecum to rectum, including appendix"
results_df[results_df$outcome == "dementia", "description"] <- "Codes for dementia in Alzheimer's, vascular dementia, unspecified dementia, senile and presenile dementia"
results_df[results_df$outcome == "any_infertility", "description"] <- "Codes for male or female infertility of different anatomical origins"
results_df[results_df$outcome == "trachea_bronchus_lung_cancer", "description"] <- "Codes for cancers in trachea, bronchi, and lungs, including personal history codes"
results_df[results_df$outcome == "nafld", "description"] <- "Code for fatty liver disease"
results_df[results_df$outcome == "renal_failure", "description"] <- "Codes for both acute, chronic and unspecified renal failure and chronic kidney disease. Also includes renal failure from hypertensive disease and various dialysis procedures in ICD-10 and OPCS-4 codes"
results_df[results_df$outcome == "stroke", "description"] <- "Codes for subarachnoid and intracerebral haemorrhages and cerebral infarctions including cerebral thrombosis and embolism, and unspecified stroke. Does not included transient cerebral ischaemia but includes acute but ill-defined cerebrovascular disease"
results_df[results_df$outcome == "T1D_probposs", "description"] <- "Algorithm sorts participants to likely diabetes status by using information on e.g. self-reported diabetes diagnsosis, age of diagnosis, medications, start of insulin within a year of diagnosis, and other self-reported data at the baseline visit"
results_df[results_df$outcome == "T2D_probposs", "description"] <- "Algorithm sorts participants to likely diabetes status by using information on e.g. self-reported diabetes diagnsosis, age of diagnosis, medications, start of insulin within a year of diagnosis, and other self-reported data at the baseline visit"
results_df[results_df$outcome == "aki", "description"] <- "Codes that are specifically for acute renal failure"
results_df[results_df$outcome == "ckd", "description"] <- "Codes that are specifically for chronic kidney disease"	
results_df[results_df$outcome == "haem_stroke", "description"] <- "Codes that denote hemorrhagic stroke"
results_df[results_df$outcome == "isch_stroke", "description"] <- "Codes that denote ischaemic stroke"
results_df[results_df$outcome == "renal_exclude", "description"] <- "Participants were excluded from the renal failure control groups, including acute renal failure and chronic kidney disease, if they had these codes pertaining to renal failure"
results_df[results_df$outcome == "stroke_exclude", "description"] <- "Participants were excluded from the stroke control groups, including hermorrhagic and ischaemic stroke, if they had these codes pertaining to stroke and transient ischaemic attacks"

results_df[results_df$outcome %in% c("T2D_probposs", "T2D_prob", "T1D_probposs", "T1D_prob"), 
	"add_info"] <- "Outcome definition from Eastwood et al algorithm. Probable and possible cases or probable only. Controls defined as diabetes unlikely. PMID: 27631769"
results_df[results_df$outcome %in% c("SOFT_CAD", "exclude_CAD"), 
	"add_info"] <- "Outcome definition from Nelson et al; SOFT CAD definition including angina. PMID: 28714975"
results_df[results_df$outcome == "breast_cancer_women_only", "add_info"] <- "Run in women only"
results_df[results_df$outcome == "trachea_bronchus_lung_cancer", "add_info"] <- "Includes cancers in trachea and bronchi"
results_df[results_df$outcome == "any_infertility", "add_info"] <- "Sex-specific codes applied in relevant sex only"

for (col in c("icd10", "icd9", "ni_noncancer")) {
	results_df[results_df$outcome == "any_infertility", col] <- paste(results_df[results_df$outcome == 
			"female_infertility", col], results_df[results_df$outcome == "male_infertility", col], 
			sep = ";")
}

results_df <- results_df[!(results_df$outcome %in% c("T2D_prob", "T1D_prob", "female_infertility", "male_infertility")), ]

dict_traits <- list(breast_cancer_women_only = "Breast cancer",
        SOFT_CAD = "CAD",
	exclude_CAD = "CAD - Control group exclusions",
        colorectal_cancer = "Colorectal Cancer",
        COPD = "COPD",
        dementia = "Dementia",
        trachea_bronchus_lung_cancer = "Lung cancer",
        renal_failure = "Renal failure",
        stroke = "Stroke",
        T2D_probposs = "Type 2 diabetes",
	T2D_prob = "Type 2 diabetes",
        T1D_prob = "Type 1 diabetes",
	T1D_probposs = "Type 1 diabetes",
        male_infertility = "Male infertility",
        female_infertility = "Female infertility",
        any_infertility = "Infertility",
        nafld = "NAFLD",
        chronic_liver_disease = "Chronic liver disease", 
	aki = "Renal Failure - acute", 
	ckd = "Renal Failure - chronic", 
	haem_stroke = "Stroke - hemorrhagic", 
	isch_stroke = "Stroke - ischaemic", 
	stroke_exclude = "Stroke - Control group exclusions", 
	renal_exclude = "Renal Failure - Control group exclusions")

for (i in 1:length(dict_traits)) {
        results_df$outcome <- as.character(replace(results_df$outcome,
                        results_df$outcome == names(dict_traits[i]), dict_traits[i]))
}

results_df <- results_df[order(results_df$outcome), ]

results_df[] <- lapply(results_df, gsub, pattern = ";", replacement = "; ")

colnames(results_df) <- c("Outcome", "ICD-10", "ICD-9", "OPCS-4", "NI, Non-cancer", "NI, Cancer", 
		"NI, operation", "Self-reported Diagnosis by Doctor", "Description", "Additional Information")
write.table(results_df, "/Net/fs1/home/linc4222/diagnose.definitions.table.180805.txt", quote = F, sep = "\t", row.names = F)

################################################################################################
#################################### Make the "publication plot" for the continuous MRs
###############################################################################################

df_raw <- read.table("../results.mr.180730/ipd.mr.continuous.results.180815.txt", stringsAsFactors = F, header =T, sep = "\t")

pulit_snp_groups <- unique(df_raw$snp_group[grep("pulit\\.sig|internal\\.unweighted|giukbb", df_raw$snp_group)])

eth_groups <- c("all.white", "brit.irish")

for (dataset in c("pulit")) {
        df_dataset <- df_raw[df_raw$snp_group %in% pulit_snp_groups, ]

        for (eth_group in eth_groups) {
                df <- df_dataset
		df <- df[df$eth_group == eth_group & df$trait %in% c("dbp", "dbp_smoking", "sbp", "sbp_smoking") &
		        df$function_name == "wald", ]

		#Show the results per 1-SD in the exposure
		df_clin <- df[df$exposure_unit == "sd" & df$outcome_unit == "clin", ]
		df_sd <- df[df$exposure_unit == "sd" & df$outcome_unit == "sd", ]

		merge_columns <- c("snp_group", "grs_unit", "eth_group", "sex_group", "trait", "grs_trait", 
				"function_name", "sd_exposure_corresponds_to", "sd_outcome_corresponds_to", 
				"grs_file", "n_complete_cases", "extra_adjustment")
		df <- merge(df_clin, df_sd, by = merge_columns)

		for (i in 1:length(dict_snp_groups_long)) {
        		df$instrument <- as.character(replace(df$instrument,
                		df$snp_group == names(dict_snp_groups_long[i]), dict_snp_groups_long[i]))
		}

		df_out <- strsplit(df$instrument, ", SNPs = ")
		df <- data.frame(df, do.call(rbind, df_out))
		df$instrument <- df$X1
		df$n_snps <- df$X2

		df$sex_group <- ifelse(df$sex_group == "comb", "Combined",
                        ifelse(df$sex_group == "men", "Men",
                        "Women"))

		df$outcome_unit <- ifelse(df$trait %in% c("bmi", "bmi_smoking"), "kg/m2",
                        ifelse(df$trait %in% c("whr", "whr_smoking"), "Ratio",
                        ifelse(df$trait %in% c("res_whr_inv", "res_whr_inv_smoking"), "-",
                        ifelse(df$trait %in% c("wc", "wc_smoking"), "cm",
                        ifelse(df$trait %in% c("hc", "hc_smoking"), "cm",
                        ifelse(df$trait %in% c("sbp", "sbp_smoking"), "mmHg",
                        ifelse(df$trait %in% c("dbp", "dbp_smoking"), "mmHg", "MISSING")))))))

		df$exposure_unit <- "SD"

		df$trait <- ifelse(df$trait == "bmi", "BMI",
                        ifelse(df$trait == "whr", "WHR",
                        ifelse(df$trait == "res_whr_inv", "WHRadjBMI",
                        ifelse(df$trait == "wc", "Waist circumference",
                        ifelse(df$trait == "hc", "Hip circumference",
                        ifelse(df$trait == "sbp", "SBP",
                        ifelse(df$trait == "dbp", "DBP", 
			ifelse(df$trait == "bmi_smoking", "BMI, adjusted smoking", 
			ifelse(df$trait == "whr_smoking", "WHR, adjusted smoking", 
			ifelse(df$trait == "res_whr_inv_smoking", "WHRadjBMI, adjusted smoking",
			ifelse(df$trait == "wc_smoking", "Waist circumference, adjusted smoking", 
			ifelse(df$trait == "hc_smoking", "Hip circumference, adjusted smoking", 
			ifelse(df$trait == "sbp_smoking", "SBP, adjusted smoking", 
			ifelse(df$trait == "dbp_smoking", "DBP, adjusted smoking", "MISSING"))))))))))))))

		df$grs_trait <- ifelse(df$grs_trait == "bmi", "BMI", 
			ifelse(df$grs_trait == "whr", "WHR", 
			ifelse(df$grs_trait == "res_whr_inv", "WHRadjBMI", "MISSING")))

		df$ci.x <- paste(formatC(df$grs_beta.x, format = "f", 2), " (",
                	formatC(df$grs_lci_beta.x, format = "f", 2), ",", formatC(df$grs_uci_beta.x, format = "f", 2), ")", sep = "")

		df$ci.y <- paste(formatC(df$grs_beta.y, format = "f", 2), " (",
                	formatC(df$grs_lci_beta.y, format = "f", 2), ",", formatC(df$grs_uci_beta.y, format = "f", 2), ")", sep = "")

		df$sex_group <- factor(df$sex_group, levels = c("Combined", "Men", "Women"))
		df$trait <- factor(df$trait, levels = c("SBP", "SBP, adjusted smoking", 
			"DBP", "DBP, adjusted smoking"))
		df <- df[order(df$grs_trait, df$trait, df$sex_group), ]

		df[, c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")][df[,
        		c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")] < 1e-200] <- 0

		df[, c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")] <- lapply(df[,
                	c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")],
                	function(x) ifelse(x < 0.001, formatC(x, format = "e", 2), ifelse(x < 0.01, formatC(x, format = "f", 3), 
			formatC(x, format = "f", 2))))

		df[, c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")][df[, c("grs_p.x", "cochrans_p.x", "cochrans_p.y", "grs_p.y")] == "0.00e+00"] <- "<1E-200"

		df <- df[, c("grs_trait", "trait", "instrument", "sex_group", "outcome_unit", "ci.x", "grs_p.x", "cochrans_p.x", "ci.y", "grs_p.y", 
			"cochrans_p.y")] 

		colnames(df) <-  c("Exposure", "Outcome", "GRS", "Sex-strata", "Outcome, Clinical Units", 
			"Estimate, Clinical (95% CI)", "P (Clinical)", "Pheterogeneity (Clinical)", 
			"Estimate, SD (95% CI)", "P (SD)", "Pheterogeneity (SD)")

		df[] <- lapply(df, gsub, pattern = "NA", replacement = "-")
		df[df$Outcome %in% c("WHRadjBMI", "WHRadjBMI, adjusted smoking"), 
			c("Estimate, Clinical (95% CI)", "P (Clinical)", "Pheterogeneity (Clinical)")] <- "-"

		change_column_names <- c("Estimate, Clinical (95% CI)", "P (Clinical)", "Pheterogeneity (Clinical)", 
					"Estimate, SD (95% CI)", "P (SD)", "Pheterogeneity (SD)")
                eth_group_name <- ifelse(eth_group == "all.white", " - Europeans", " - British")

                colnames(df)[which(names(df) %in% change_column_names)] <- paste(change_column_names, eth_group_name, sep = "")

		if (eth_group == "all.white") {
                       all_white_df <- df
                } else if (eth_group == "brit.irish") {
                       brit_df <- df
                }
        }

	final_df <- merge(all_white_df, brit_df, by = c("Exposure", "Outcome", "GRS", "Sex-strata", "Outcome, Clinical Units"), all = T)
	final_df <- final_df[order(final_df$Exposure, final_df$Outcome, final_df[, "Sex-strata"], final_df$GRS), ]

	unadj_outcomes <- c("BMI", "WHR", "WHRadjBMI", "Waist circumference", "Hip circumference",
              "SBP", "DBP")

	unadj_smoking <- final_df[final_df$Outcome %in% unadj_outcomes, ]
	unadj_smoking_output_file <- paste("/Net/fs1/home/linc4222/cleaned.mr.continuous.for.publication.unadj.smoking.", 
			dataset, ".180817.txt", sep = "")
	write.table(unadj_smoking, unadj_smoking_output_file, row.names = F,
                       quote = F, sep = "\t", na = "-")

	adj_smoking <- final_df[!(final_df$Outcome %in% unadj_outcomes), ]
	adj_smoking$Outcome <- gsub(", adjusted smoking", "", adj_smoking$Outcome)
	adj_smoking_output_file <- paste("/Net/fs1/home/linc4222/cleaned.mr.continuous.for.publication.adj.smoking.", 
			dataset, ".180817.txt", sep = "")
	write.table(adj_smoking, adj_smoking_output_file, row.names = F,
                       quote = F, sep = "\t", na = "-")

}

###################################################################################
# Make table with the MRs where we've adjusted for smoking
##################################################################################

df_raw <- read.table("../results.mr.180730/ipd.mr.binary.results.180815.txt", stringsAsFactors = F, header =T, sep = "\t")

pulit_snp_groups <- unique(df_raw$snp_group[grep("pulit|internal|giukbb", df_raw$snp_group)])

eth_groups <- c("all.white", "brit.irish")
datasets <- c("pulit")

for (dataset in datasets) {
        df_dataset <- df_raw[df_raw$snp_group %in% pulit_snp_groups, ]

	for (eth_group in eth_groups) {
		#Subset to relevant rows
		df <- df_dataset
		df <- df[df$exposure_unit == "sd" & df$eth_group == eth_group & df$function_name == "wald" & 
			((df$snp_group %in% unique(df[grep("pulit\\.sig$", df$snp_group), "snp_group"]) &
			!(df$trait %in% c("t1d_cases_prob", "t2d_cases_prob", "t1d_cases_prob_smoking", "t2d_cases_prob_smoking"))) |
			(df$trait == "smoker_cases")), ]

		for (i in 1:length(dict_snp_groups_long)) {
        		df$instrument <- as.character(replace(df$instrument,
                		df$snp_group == names(dict_snp_groups_long[i]), dict_snp_groups_long[i]))
		}

		df_out <- strsplit(df$instrument, ", SNPs = ")
		df <- data.frame(df, do.call(rbind, df_out))
		df$instrument <- df$X1
		df$n_snps <- df$X2

		df$sex_group <- ifelse(df$sex_group == "comb", "Combined",
                        ifelse(df$sex_group == "men", "Men",
                        "Women"))

		dict_traits <- list(breast_cancer_cases_smoking = "Breast Cancer",
        		cad_cases_smoking = "CAD",
			cld_cases_smoking = "CLD",
        		copd_cases_smoking = "COPD",
        		renal_failure_cases_smoking = "Renal Failure",
			aki_cases_smoking = "Renal Failure - Acute", 
			ckd_cases_smoking = "Renal Failure - Chronic",
        		t2d_cases_probposs_smoking = "Type 2 Diabetes",
        		nafld_cases_smoking = "NAFLD",
			stroke_cases_smoking = "Stroke", 
			isch_stroke_cases_smoking = "Stroke - Ischaemic", 
			lungcancer_cases_smoking = "Lung Cancer",
			smoker_cases = "Having smoked/smoker", 
			t1d_cases_probposs_smoking = "Type 1 Diabetes", 
			breast_cancer_cases = "Breast Cancer",
                        cad_cases = "CAD",
                        cld_cases = "CLD",
                        copd_cases = "COPD",
                        renal_failure_cases = "Renal Failure",
                        aki_cases = "Renal Failure - Acute",
                        ckd_cases = "Renal Failure - Chronic",
                        t2d_cases_probposs = "Type 2 Diabetes",
                        nafld_cases = "NAFLD",
                        stroke_cases = "Stroke",
                        isch_stroke_cases = "Stroke - Ischaemic",
                        lungcancer_cases = "Lung Cancer",
                        smoker_cases = "Having smoked/smoker",
                        t1d_cases_probposs = "Type 1 Diabetes")

		for (i in 1:length(dict_traits)) {
        		df$trait <- as.character(replace(df$trait,
        	                df$trait == names(dict_traits[i]), dict_traits[i]))
		}

		df$grs_trait_name <- ifelse(df$grs_trait == "bmi", "BMI",
                        ifelse(df$grs_trait == "whr", "WHR",
                        ifelse(df$grs_trait == "res_whr_inv", "WHRadjBMI", NA)))

		df$unique_combinations <- paste(df$trait, "_", ifelse(df$grs_trait_name == "BMI", "01",
                                ifelse(df$grs_trait_name == "WHR", "02",
                                ifelse(df$grs_trait_name == "WHRadjBMI", "03", "missing"))), ifelse(df$sex_group == "Combined", "01",
                                ifelse(df$sex_group == "Men", "02", ifelse(df$sex_group == "Women", "03", "missing"))),
                                df$sex_group, sep = "")

		df <- df[order(df$unique_combinations), ]

		df$ci <- paste(formatC(df$grs_or, format = "f", 2), " (",
                	formatC(df$grs_lci_or, format = "f", 2), ", ", formatC(df$grs_uci_or, format = "f", 2), ")", sep = "")

		df <- df[, c("grs_trait_name", "trait", "instrument", "sex_group",
                        "ci", "grs_p", "cochrans_p", "extra_adjustment")]

		df[, c("grs_p", "cochrans_p")][df[, c("grs_p", "cochrans_p")] < 1e-200] <- 0

		df[, c("grs_p", "cochrans_p")] <- lapply(df[, c("grs_p", "cochrans_p")],
        	        function(x) ifelse(x < 0.001, formatC(x, format = "e", 2),
                                ifelse(x < 0.01, formatC(x, format = "f", 3), formatC(x, format = "f", 2))))

		df[, c("grs_p", "cochrans_p")][df[, c("grs_p", "cochrans_p")] == "0.00e+00"] <- "<1e-200"

		colnames(df) <-  c("Exposure", "Outcome", "GRS", "Sex-Strata", "OR (95% CI)", "P", "Pheterogeneity", "extra_adjustment")

		change_column_names <- c("OR (95% CI)", "P", "Pheterogeneity")
                eth_group_name <- ifelse(eth_group == "all.white", " - Europeans", " - British")

                colnames(df)[which(names(df) %in% change_column_names)] <- paste(change_column_names, eth_group_name, sep = "")

		if (eth_group == "all.white") {
                        all_white_df <- df
                } else if (eth_group == "brit.irish") {
                        brit_df <- df
                }
	}

	final_df <- merge(all_white_df, brit_df, by = c("Exposure", "Outcome", "GRS", "Sex-Strata", "extra_adjustment"), all = T)

	binary_traits <- final_df[!(final_df$Outcome %in% c("Having smoked/smoker")), ]
	binary_traits_unadj <- binary_traits[binary_traits$extra_adjustment == "-", ]
	binary_traits_adj <- binary_traits[binary_traits$extra_adjustment == "smoking", ]

	binary_traits <- merge(binary_traits_unadj, binary_traits_adj, by = c("Exposure", "Outcome", "GRS", "Sex-Strata"))
	binary_traits[, c("extra_adjustment.x", "extra_adjustment.y")] <- NULL
	colnames(binary_traits) <- c("Exposure", "Outcome", "GRS", "Sex-Strata", 
			"OR (95% CI), unadj - Europeans", "P, unadj - Europeans", "Phet, unadj - Europeans", 
			"OR (95% CI), unadj - British", "P, unadj - British", "Phet, unadj - British", 
			"OR (95% CI), adj - Europeans", "P, adj - Europeans", "Phet, adj - Europeans", 
			"OR (95% CI), adj - British", "P, adj - British", "Phet, adj - British")
 
	binary_traits <- binary_traits[order(binary_traits$Exposure, binary_traits$Outcome, binary_traits[, "Sex-Strata"]), ]
	binary_traits_output_file <- paste("/Net/fs1/home/linc4222/associations.for.mr.binary.traits.adj.smoking.", 
		dataset, ".180819.txt", sep = "")
	write.table(binary_traits, binary_traits_output_file, row.names = F,
                       quote = F, sep = "\t", na = "-")

	smoker <- final_df[(final_df$Outcome %in% c("Having smoked/smoker")), ]
	smoker <- smoker[order(smoker$Exposure, smoker$Outcome, smoker[, "Sex-Strata"], smoker$GRS), ]
	smoker_output_file <- paste("/Net/fs1/home/linc4222/associations.for.being.smoker.mr.", dataset, ".180820.txt", sep = "")
	write.table(smoker, smoker_output_file, row.names = F,
                        quote = F, sep = "\t", na = "-")
	
}

#######################################################
############## FG and FI
#####################################################

df_raw <- read.table("../results.mr.180730/summary.mr.results.180730.txt", 
	stringsAsFactors = F, header = T, sep = "\t")

pulit_snp_groups <- unique(df_raw$snp_group[grep("pulit\\.sig", df_raw$snp_group)])
datasets <- c("pulit")

for (dataset in datasets) {
        df_dataset <- df_raw[df_raw$snp_group %in% pulit_snp_groups, ]

        df <- df_dataset[df_dataset$snp_group %in% 
		df_dataset$snp_group[grep("sig", df_dataset$snp_group)], ]

	methods <- c("IVW", "MR-Egger", "MR-Egger (intercept)", "Weighted median")
	df <- df[df$method %in% methods, ]

	df$sex <- ifelse(df$sex_group == "comb", "Combined", 
			ifelse(df$sex_group == "women", "Women", 
			ifelse(df$sex_group == "men", "Men", "MISSING")))

	df$number_snps <- paste(df$n_snps_same_alleles, " (", df$n_exp_snps, ")", sep = "")
	df$exp_trait <- gsub("\\..+", "", df$snp_group)
	df$exp_trait <- ifelse(df$exp_trait == "bmi", "BMI", 
				ifelse(df$exp_trait == "whr", "WHR", 
				ifelse(df$exp_trait == "whradjbmi", "WHRadjBMI", "MISSING")))

	for (i in 1:length(dict_snp_groups_long)) {
                        df$instrument <- as.character(replace(df$instrument,
                                df$snp_group == names(dict_snp_groups_long[i]), dict_snp_groups_long[i]))
                }
	
	df_out <- strsplit(df$instrument, ", SNPs = ")
        df <- data.frame(df, do.call(rbind, df_out))
        df$instrument <- df$X1
        df$n_snps <- df$X2

        df$outcome_unit <- ifelse(df$trait == "FG", "mmol/L", 
				ifelse(df$trait == "FI", "ln(pmol/L)", "MISSING"))

	df$ci <- paste(formatC(df$beta, format = "f", 3), " (",
                        formatC(df$beta_lci, format = "f", 3), ",", 
			formatC(df$beta_uci, format = "f", 3), ")", sep = "")

	df[, c("p", "cochrans_p")][df[, c("p", "cochrans_p")] < 1e-200] <- 0

	df[, c("p", "cochrans_p")] <- lapply(df[, c("p", "cochrans_p")], 
					function(x) ifelse(x < 0.001, formatC(x, format = "e", 2), 
					ifelse(x < 0.01, formatC(x, format = "f", 3), 
					formatC(x, format = "f", 2))))

        df[, c("p", "cochrans_p")][df[, c("p", "cochrans_p")] == "0.00e+00"] <- "<1E-200"

	df <- df[order(df$exp_trait, df$trait, df$sex_group), ]

	df <- df[, c("exp_trait", "trait", "sex", "number_snps", "outcome_unit",
			"method", "ci", "p", "cochrans_p")]
	colnames(df) <- c("Exposure", "Outcome", "Sex-Strata", "N SNPs available (N SNPs GRS)",
			"Outcome Unit", "Method", "Estimate (95% CI)", "P", "Pheterogeneity")

	output_file <- paste("/Net/fs1/home/linc4222/fg.fi.mr.", dataset, ".180923.txt", sep = "")
	write.table(df, output_file, quote = F, row.names = F, sep = "\t", na = "-")
	
}

######################################################################################
######################### GET TABLE WITH SAME N OF CASES AND CONTROLS
#####################################################################################

df_raw <- read.table("../results.logistic.regressions.180514/log.results.table.same.n.cases.controls.180627.txt",
                stringsAsFactors = F, header = T)
df_raw$log_type <- "same_n"

df_diff <- read.table("../results.logistic.regressions.180514/log.results.table.180627.txt", 
		stringsAsFactors = F, header = T)
df_diff <- df_diff[df_diff$grs_unit == "raw_scoresum", ]
df_diff$grs_unit <- NULL
df_diff$log_type <- "diff_n"

df_raw <- rbind(df_raw, df_diff)

smoking_names <- unique(df_raw$case_column[grep("_smoking", df_raw$case_column)])

df_raw <- df_raw[!(df_raw$case_column %in% c("smoker_cases", smoking_names) |
        df_raw$case_column %in% c("t1d_cases_prob", "t2d_cases_prob")), ]
eth_groups <- c("all.white", "brit.irish")

for (dataset in c("pulit")) {
        comb_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig",  "whradjbmi.eur.comb.pulit.sig")
        male_groups <- c("bmi.eur.men.pulit.sig", "whr.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.sig")
        female_groups <- c("bmi.eur.women.pulit.sig", "whr.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.sig")

        for (eth_group in eth_groups) {
                df <- df_raw

                df <- df[df$eth_group == eth_group & df$snp_group %in% c(male_groups, female_groups) & !(df$case_column == "breast_cancer_cases"), ]

                for (i in 1:length(dict_snp_groups_short)) {
                        df$instrument <- as.character(replace(df$instrument,
                                df$snp_group == names(dict_snp_groups_short[i]), dict_snp_groups_short[i]))
                }

                df_out <- strsplit(df$instrument, ", SNPs = ")
                df <- data.frame(df, do.call(rbind, df_out))
                df$instrument <- df$X1
                df$n_snps <- df$X2

                df$sex_group <- ifelse(df$sex_group == "comb", "Combined",
                        ifelse(df$sex_group == "men", "Men",
                        ifelse(df$sex_group == "women", "Women", NA)))

                df$ci <- paste(formatC(df$grs_or, format = "f", 2), " (",
                        formatC(df$grs_lci_or, format = "f", 2), ", ", formatC(df$grs_uci_or, format = "f", 2), ")", sep = "")

                df$trait <- df$case_column
                dict_traits <- list(breast_cancer_cases = "Breast Cancer",
                        cad_cases = "CAD",
                        colorectal_cancer_cases = "Colorectal Cancer",
                        copd_cases = "COPD",
                        dementia_cases = "Dementia",
                        lungcancer_cases = "Lung Cancer",
                        renal_failure_cases = "Renal Failure",
                        aki_cases = "Renal Failure - Acute",
                        ckd_cases = "Renal Failure - Chronic",
                        stroke_cases = "Stroke",
                        isch_stroke_cases = "Stroke - Ischemic",
                        haem_stroke_cases = "Stroke - Hemorrhagic",
                        t2d_cases_probposs = "Type 2 Diabetes",
                        t1d_cases_probposs = "Type 1 Diabetes",
                        any_infertility_cases = "Infertility",
                        nafld_cases = "NAFLD",
                        cld_cases = "CLD")

                for (i in 1:length(dict_traits)) {
                        df$trait <- as.character(replace(df$trait,
                                df$trait == names(dict_traits[i]), dict_traits[i]))
                }

                df$grs_trait <- gsub("\\..*", "", df$snp_group)
                df$order <- ifelse(df$sex_group == "Combined", 1,
                        ifelse(df$sex_group == "Men", 2,
                        ifelse(df$sex_group == "Women", 3, NA)))

                df <- df[order(df$trait, df$grs_trait, df$order), ]

                df <- df[, c("trait", "instrument", "n_snps", "sex_group",
                        "ci", "grs_p", "cochrans_p", "log_type")]

		df[, c("grs_p", "cochrans_p")][df[, c("grs_p", "cochrans_p")] < 1e-200] <- 0
		
                df[, c("grs_p", "cochrans_p")] <- lapply(df[, c("grs_p", "cochrans_p")],
                        function(x) ifelse(x < 0.001, formatC(x, format = "e", 2),
                                ifelse(x < 0.01, formatC(x, format = "f", 3), formatC(x, format = "f", 2))))

		df[, c("grs_p", "cochrans_p")][df[, c("grs_p", "cochrans_p")] == "0.00e+00"] <- "<1e-200"

                colnames(df) <-  c("Outcome", "Sex-Specific Estimates Instrument", "N SNPs",
                        "Sex-Strata", "OR (95% CI)", "P", "Pheterogeneity", "log_type")

                change_column_names <- c("OR (95% CI)", "P", "Pheterogeneity")
                eth_group_name <- ifelse(eth_group == "all.white", " - Europeans", " - British")

                colnames(df)[which(names(df) %in% change_column_names)] <- paste(change_column_names, eth_group_name, sep = "")

                if (eth_group == "all.white") {
                        all_white_df <- df
                } else if (eth_group == "brit.irish") {
                        brit_df <- df
                }

        }

	final_df <- merge(all_white_df, brit_df, by = c("Outcome", "Sex-Specific Estimates Instrument", "N SNPs", "Sex-Strata", "log_type"), all = T)
	final_df_same <- final_df[final_df$log_type == "same_n", ]
	final_df_diff <- final_df[final_df$log_type == "diff_n", ]

	colnames(final_df_same)[grep("^OR|^P", colnames(final_df_same))] <- paste(colnames(final_df_same)[grep("^OR|^P", colnames(final_df_same))], 
										", same N", sep = "")
	colnames(final_df_diff)[grep("^OR|^P", colnames(final_df_diff))] <- paste(colnames(final_df_diff)[grep("^OR|^P", colnames(final_df_diff))],
                                                                                ", different N", sep = "")

	final_df_out <- merge(final_df_diff, final_df_same, by = c("Outcome", "Sex-Specific Estimates Instrument",
			"N SNPs", "Sex-Strata"))
	final_df_out[, c("log_type.x", "log_type.y")] <- NULL

        output_file <- paste("/Net/fs1/home/linc4222/results.grs.same.n.cases.controls.", dataset, ".180926.txt", sep = "")
        write.table(final_df_out, output_file, row.names = F,
                                quote = F, sep = "\t", na = "-")
	
}

################################################
##### Get MR-Egger intercept
##################################################

#These results show that for the sexually heterogeneous results, the "right"
#sex has higher estimates for MR-Egger, IVW, weighted median
#of weighting strategy
df <- read.table("../results.mr.180730/mr.egger.binary.results.180730.txt", 
		stringsAsFactors = F, header = T, sep = "\t")

df$grs_trait <- gsub("\\.(.)+", "", df$snp_group)
df$no_sex_snp_group <- gsub("women|men", "", df$snp_group)

df <- df[df$method %in% c("IVW", "MR-Egger", "Weighted median") & 
	df$snp_group %in% c("bmi.eur.men.pulit.sig", "bmi.eur.women.pulit.sig", 
	"whr.eur.men.pulit.sig", "whr.eur.women.pulit.sig", 
	"whradjbmi.eur.women.giukbb.sig", "whradjbmi.eur.women.giukbb.sig"), ]
df_men <- df[df$sex_group == "men", ]
df_women <- df[df$sex_group == "women", ]

df <- merge(df_women, df_men, by = c("trait", "method", "cochrans_p", "cochrans_q", "cochrans_i2", "grs_trait", "no_sex_snp_group"))

df_women_highest <- df[(df$trait == "t2d_cases_probposs" & df$grs_trait == "bmi"), ]
df_women_highest$women_higher <- ifelse(df_women_highest$beta.x > df_women_highest$beta.y, T, F)
df_women_highest

df_men_highest <- df[(df$trait == "ckd_cases" & df$grs_trait %in% c("whr", "whradjbmi")) |
	(df$trait == "renal_failure_cases" & df$grs_trait == "whr") |
	(df$trait == "copd_cases" & df$grs_trait == "whr"), ]
df_men_highest$men_higher <- ifelse(df_men_highest$beta.x < df_men_highest$beta.y, T, F)
df_men_highest


#Make table of the MR-Egger intercept test
df <- read.table("../results.mr.180730/mr.egger.binary.results.180730.txt",
                stringsAsFactors = F, header = T, sep = "\t")
df <- df[df$method %in% c("MR-Egger (intercept)", "MR-Egger") & df$snp_group %in% unique(df$snp_group[grep("\\.pulit\\.sig", df$snp_group)]) &
		!(df$trait %in% c("t2d_cases_prob", "t1d_cases_prob", "smoker_cases")), ]

grs_results <- read.table("../results.logistic.regressions.180514/log.results.table.180627.txt", stringsAsFactors = F,
                header = T, sep = "\t")

#Get the significant GRS results and subset so get MR-Egger results for those
comb_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig", "whradjbmi.eur.comb.pulit.sig")
male_groups <- c("bmi.eur.men.pulit.sig", "whr.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.sig")
female_groups <- c("bmi.eur.women.pulit.sig", "whr.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.sig")

adjusted_smoking <- unique(grs_results$case_column[grep("_smoking", grs_results$case_column)])
grs_results <- grs_results[!(grs_results$case_column %in% adjusted_smoking), ]
grs_results <- grs_results[((grs_results$snp_group %in% c(comb_groups) &
                grs_results$sex_group == "comb") |
                (grs_results$snp_group %in% c(male_groups) &
                grs_results$sex_group == "men") |
                (grs_results$snp_group %in% c(female_groups) &
                grs_results$sex_group == "women")) &
                grs_results$grs_unit == "raw_scoresum" &
                grs_results$eth_group %in% c("all.white", "brit.irish"), ]

critical_p <- 0.05/(length(unique(grs_results[!(grs_results$case_column %in% c("t1d_cases_prob", "t2d_cases_prob", "smoker_cases")),
                "case_column"]))*length(unique(gsub("\\..*", "", grs_results$snp_group))))

for (i in 1:length(comb_groups)) {
        for (grs_unit in unique(grs_results$grs_unit)) {
                for (case_column in unique(grs_results$case_column[!(grs_results$case_column %in% c("t1d_cases_prob", "smoker_cases", "t2d_cases_prob"))])) {
                         if (nrow(grs_results[(grs_results$snp_group %in% comb_groups[i]) &
                                  grs_results$grs_unit == grs_unit & grs_results$case_column == case_column & grs_results$grs_p < critical_p, ]) == 0) {

                                grs_results <- grs_results[!(grs_results$snp_group %in% c(comb_groups[i]) &
                                                grs_results$grs_unit == grs_unit & grs_results$case_column == case_column), ]
                        }
                }
        }
}

#Removes the non-significant sex-specific GRSs, except the sensitivity analyses
for (i in 1:length(female_groups)) {
        for (grs_unit in unique(grs_results$grs_unit)) {
                 for (case_column in unique(grs_results$case_column[!(grs_results$case_column %in% c("t1d_cases_prob", "smoker_cases", "t2d_cases_prob"))])) {
                         if (nrow(grs_results[grs_results$snp_group %in% c(female_groups[i], male_groups[i]) &
                                  grs_results$grs_unit == grs_unit & grs_results$case_column == case_column & grs_results$grs_p < critical_p, ]) == 0) {

                                  grs_results <- grs_results[!(grs_results$snp_group %in% c(female_groups[i], male_groups[i]) &
                                                grs_results$grs_unit == grs_unit & grs_results$case_column == case_column), ]

                        }
                }
        }
}

for (snp_group in unique(df$snp_group[grep("\\.pulit\\.sig", df$snp_group)])) {
	for (trait in unique(df$trait)) {
		for (sex_group in unique(df$sex_group)) {
			
			if (nrow(grs_results[grs_results$snp_group == snp_group &
				grs_results$case_column == trait & grs_results$sex_group == sex_group, ]) == 0) {
				
			df <- df[!(df$snp_group == snp_group & df$trait == trait & df$sex_group == sex_group), ]
			}
		}
	}
}

for (i in 1:length(dict_snp_groups_short)) {
       df$instrument <- as.character(replace(df$instrument,
                     df$snp_group == names(dict_snp_groups_short[i]), dict_snp_groups_short[i]))
}

df_out <- strsplit(df$instrument, ", SNPs = ")
df <- data.frame(df, do.call(rbind, df_out))
df$instrument <- df$X1
df$n_snps <- df$X2

df$sex_group <- ifelse(df$sex_group == "comb", "Combined",
                ifelse(df$sex_group == "men", "Men",
                ifelse(df$sex_group == "women", "Women", NA)))

df$ci_beta <- paste(formatC(df$beta, format = "f", 3), " (",
              formatC(df$beta_lci, format = "f", 3), ", ", formatC(df$beta_uci, format = "f", 3), ")", sep = "")

df$ci_or <- paste(formatC(df$or, format = "f", 2), " (",
              formatC(df$or_lci, format = "f", 2), ", ", formatC(df$or_uci, format = "f", 2), ")", sep = "")

dict_traits <- list(breast_cancer_cases = "Breast Cancer",
                        cad_cases = "CAD",
                        colorectal_cancer_cases = "Colorectal Cancer",
                        copd_cases = "COPD",
                        dementia_cases = "Dementia",
                        lungcancer_cases = "Lung Cancer",
                        renal_failure_cases = "Renal Failure",
                        aki_cases = "Renal Failure - Acute",
                        ckd_cases = "Renal Failure - Chronic",
                        stroke_cases = "Stroke",
                        isch_stroke_cases = "Stroke - Ischemic",
                        haem_stroke_cases = "Stroke - Hemorrhagic",
                        t2d_cases_probposs = "Type 2 Diabetes",
                        t1d_cases_probposs = "Type 1 Diabetes",
                        any_infertility_cases = "Infertility",
                        nafld_cases = "NAFLD",
                        cld_cases = "CLD")

for (i in 1:length(dict_traits)) {
          df$trait <- as.character(replace(df$trait,
                       df$trait == names(dict_traits[i]), dict_traits[i]))
}

df$pval <- ifelse(df$p < 0.001, formatC(df$p, format = "e", 2),
                     ifelse(df$p < 0.01, formatC(df$p, format = "f", 3), formatC(df$p, format = "f", 2)))

df$grs_trait_name <- gsub("\\.(.)+", "", df$snp_group)
df$grs_trait_name <- ifelse(df$grs_trait_name == "bmi", "BMI", 
			ifelse(df$grs_trait_name == "whr", "WHR", 
			ifelse(df$grs_trait_name == "whradjbmi", "WHRadjBMI", NA)))

df$number_snps <- paste(df$n_snps_same_alleles, " (", df$n_exp_snps, ")", sep = "")

df_egger <- df[df$method == "MR-Egger", c("trait", "grs_trait_name", "sex_group",
                "number_snps", "ci_or", "pval", "cochrans_p")]
df_intercept <- df[df$method == "MR-Egger (intercept)", c("trait", "grs_trait_name", "sex_group",
                "number_snps", "ci_beta", "pval")]

troll <- merge(df_intercept, df_egger, by = c("trait", "grs_trait_name", "sex_group", 
		"number_snps"))

troll[, c("pval.x", "pval.y", "cochrans_p")] <- lapply(troll[, c("pval.x", "pval.y", "cochrans_p")],
                        function(x) ifelse(x < 0.001, formatC(x, format = "e", 2),
                                ifelse(x < 0.01, formatC(x, format = "f", 3), formatC(x, format = "f", 2))))

colnames(troll) <- c("Outcome", "Exposure", "Sex-strata", "N SNPs used (N SNPs GRS)",  
		"Intercept (95% CI)", "P intercept", "OR (95% CI)", 
		"P estimate", "Pheterogeneity")

write.table(troll, "/Net/fs1/home/linc4222/mr.egger.results.181108.txt", quote = F, sep = "\t", 
		row.names = F, na = "-")
