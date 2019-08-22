#!/bin/env Rscript
#$-cwd


extract_cases <- function(data, codes, field) {	 
	data <- data[apply(data[, grep(field, colnames(data))], 1, function(x) {any(x %in% codes)}), 1]
	return(data)
}

##ICD10 cases
cases_icd10 <- function(data, counts_file, diagnose, codes) {
	cases <- integer()
	for (i in seq_along(codes)) {
		code = codes[i]
		icd10_primary = extract_cases(data, code, "^41202\\.")
		icd10_secondary = extract_cases(data, code, "^41204\\.")
		icd10_cancer_registry = extract_cases(data, code, "^40006\\.")
		icd10_death_primary = extract_cases(data, code, "^40001\\.")
		icd10_death_secondary = extract_cases(data, code, "^40002\\.")
		all_cases_this_code <- unique(c(icd10_primary, icd10_secondary,
                                icd10_cancer_registry, icd10_death_primary,
                                icd10_death_secondary))
		all_diagnose_types_names <- list("icd10_primary", "icd10_secondary", 
				"icd10_cancer_registry", "icd10_death_primary", 
				"icd10_death_secondary", "all_cases_this_code")
		all_diagnose_types <- list(icd10_primary, icd10_secondary, 
				icd10_cancer_registry, icd10_death_primary, 
				icd10_death_secondary, all_cases_this_code) 

		for (j in seq_along(all_diagnose_types)) {
			sink(counts_file, append = T)
			cat(paste(diagnose, "icd10", all_diagnose_types_names[j], code,  
				length(all_diagnose_types[[j]]), "\n", sep = "\t"))
			sink()
		}
			
		cases <- c(cases, all_cases_this_code)
	}
	cases <- unique(cases)
	sink(counts_file, append = T)
	cat(paste(diagnose, "icd10", "total_icd10", paste(codes, collapse = ";"), length(cases), "\n", sep = "\t"))
	sink()
	
	return(cases)
}

##ICD9 cases
cases_icd9 <- function(data, counts_file, diagnose, codes) {
        cases <- integer()
        for (i in seq_along(codes)) {
                code = codes[i]
                icd9_primary = extract_cases(data, code, "^41203\\.")
                icd9_secondary = extract_cases(data, code, "^41205\\.")
                icd9_cancer_registry = extract_cases(data, code, "^40013\\.")
                all_cases_this_code <- unique(c(icd9_primary, icd9_secondary,
                                icd9_cancer_registry))
                all_diagnose_types_names <- list("icd9_primary", "icd9_secondary",
                                "icd9_cancer_registry", "all_cases_this_code")
                all_diagnose_types <- list(icd9_primary, icd9_secondary,
                                icd9_cancer_registry, all_cases_this_code)

                for (j in seq_along(all_diagnose_types)) {
                        sink(counts_file, append = T)
                        cat(paste(diagnose, "icd9", all_diagnose_types_names[j], code,
                                length(all_diagnose_types[[j]]), "\n", sep = "\t"))
                        sink()
                }

                cases <- c(cases, all_cases_this_code)
        }
        cases <- unique(cases)
        sink(counts_file, append = T)
        cat(paste(diagnose, "icd9", "total_icd9", paste(codes, collapse = ";"), length(cases), "\n", sep = "\t"))
        sink()

        return(cases)
}

##Operation cases
cases_opcs4 <- function(data, counts_file, diagnose, codes) {
        cases <- integer()
        for (i in seq_along(codes)) {
                code = codes[i]
                opcs4_main = extract_cases(data, code, "^41200\\.")
                opcs4_secondary = extract_cases(data, code, "^41210\\.")
                all_cases_this_code <- unique(c(opcs4_main, opcs4_secondary))
                all_diagnose_types_names <- list("opcs4_main", "opcs4_secondary",
                                "all_cases_this_code")
                all_diagnose_types <- list(opcs4_main, opcs4_secondary, all_cases_this_code)

                for (j in seq_along(all_diagnose_types)) {
                        sink(counts_file, append = T)
                        cat(paste(diagnose, "opcs4", all_diagnose_types_names[j], code,
                                length(all_diagnose_types[[j]]), "\n", sep = "\t"))
                        sink()
                }

                cases <- c(cases, all_cases_this_code)
        }
        cases <- unique(cases)
        sink(counts_file, append = T)
        cat(paste(diagnose, "opcs4", "total_opcs4", paste(codes, collapse = ";"), length(cases), "\n", sep = "\t"))
        sink()

        return(cases)
}

##Self-reported non_cancer in nurse interview 
cases_ni_non_cancer <- function(data, counts_file, diagnose, codes) {
        cases <- integer()
        for (i in seq_along(codes)) {
                code = codes[i]
                ni_non_cancer = extract_cases(data, code, "^20002\\.")
                all_cases_this_code <- unique(c(ni_non_cancer))
                all_diagnose_types_names <- list("ni_non_cancer",
                                "all_cases_this_code")
                all_diagnose_types <- list(ni_non_cancer, all_cases_this_code)

                for (j in seq_along(all_diagnose_types)) {
                        sink(counts_file, append = T)
                        cat(paste(diagnose, "nurse_interview_non_cancer", all_diagnose_types_names[j], code,
                                length(all_diagnose_types[[j]]), "\n", sep = "\t"))
                        sink()
                }
                cases <- c(cases, all_cases_this_code)
        }
        cases <- unique(cases)
        sink(counts_file, append = T)
        cat(paste(diagnose, "nurse_interview_non_cancer", "total_ni_non_cancer", paste(codes, collapse = ";"), length(cases), "\n", sep = "\t"))
        sink()

        return(cases)
}

##Self-reported cancer in nurse interview
cases_ni_cancer <- function(data, counts_file, diagnose, codes) {
        cases <- integer()
        for (i in seq_along(codes)) {
                code = codes[i]
                ni_cancer = extract_cases(data, code, "^20001\\.")
                all_cases_this_code <- unique(c(ni_cancer))
                all_diagnose_types_names <- list("ni_cancer",
                                "all_cases_this_code")
                all_diagnose_types <- list(ni_cancer, all_cases_this_code)

                for (j in seq_along(all_diagnose_types)) {
                        sink(counts_file, append = T)
                        cat(paste(diagnose, "nurse_interview_cancer", all_diagnose_types_names[j], code,
                                length(all_diagnose_types[[j]]), "\n", sep = "\t"))
                        sink()
                }

                cases <- c(cases, all_cases_this_code)
        }
        cases <- unique(cases)
        sink(counts_file, append = T)
        cat(paste(diagnose, "nurse_interview_cancer", "total_ni_cancer", paste(codes, collapse = ";"), length(cases), "\n", sep = "\t"))
        sink()

        return(cases)
}

##Self-reported operations in nurse interview
cases_ni_operation <- function(data, counts_file, diagnose, codes) {
        cases <- integer()
        for (i in seq_along(codes)) {
                code = codes[i]
                ni_operation = extract_cases(data, code, "^20004\\.")
                all_cases_this_code <- unique(c(ni_operation))
                all_diagnose_types_names <- list("ni_operation",
                                "all_cases_this_code")
                all_diagnose_types <- list(ni_operation, all_cases_this_code)

                for (j in seq_along(all_diagnose_types)) {
                        sink(counts_file, append = T)
                        cat(paste(diagnose, "nurse_interview_operation", all_diagnose_types_names[j], code,
                                length(all_diagnose_types[[j]]), "\n", sep = "\t"))
                        sink()
                }

                cases <- c(cases, all_cases_this_code)
        }
        cases <- unique(cases)
        sink(counts_file, append = T)
        cat(paste(diagnose, "nurse_interview_operation", "total_ni_operation", paste(codes, collapse = ";"), length(cases), "\n", sep = "\t"))
        sink()

        return(cases)
}

##Medications
cases_medications <- function(data, counts_file, diagnose, codes) {
        cases <- integer()
        for (i in seq_along(codes)) {
                code = codes[i]
                medications = extract_cases(data, code, "^20003\\.")
                all_cases_this_code <- unique(c(medications))
                all_diagnose_types_names <- list("medications",
                                "all_cases_this_code")
                all_diagnose_types <- list(medications, all_cases_this_code)

                for (j in seq_along(all_diagnose_types)) {
                        sink(counts_file, append = T)
                        cat(paste(diagnose, "medications", all_diagnose_types_names[j], code,
                                length(all_diagnose_types[[j]]), "\n", sep = "\t"))
                        sink()
                }

                cases <- c(cases, all_cases_this_code)
        }
        cases <- unique(cases)
        sink(counts_file, append = T)
        cat(paste(diagnose, "medications", "total_medications", paste(codes, collapse = ";"), length(cases), "\n", sep = "\t"))
        sink()

        return(cases)
}

##Cases self-reported medications
cases_chol_bp_dm_meds <- function(data, counts_file, diagnose, codes) {
        cases <- integer()
        for (i in seq_along(codes)) {
                code = codes[i]
                chol_bp_dm_meds_men = extract_cases(data, code, "^6177\\.")
		chol_bp_dm_meds_women = extract_cases(data, code, "^6153\\.")
		all_cases_this_code <- unique(c(chol_bp_dm_meds_men, chol_bp_dm_meds_women))
                all_diagnose_types_names <- list("chol_bp_dm_meds_men", "chol_bp_dm_meds_women",
                                "all_cases_this_code")
                all_diagnose_types <- list(chol_bp_dm_meds_men, chol_bp_dm_meds_women, all_cases_this_code)

                for (j in seq_along(all_diagnose_types)) {
                        sink(counts_file, append = T)
                        cat(paste(diagnose, "chol_bp_dm_meds", all_diagnose_types_names[j], code,
                                length(all_diagnose_types[[j]]), "\n", sep = "\t"))
                        sink()
                }

                cases <- c(cases, all_cases_this_code)
        }
        cases <- unique(cases)
        sink(counts_file, append = T)
        cat(paste(diagnose, "chol_bp_dm_meds", "total_chol_bp_dm_meds", paste(codes, collapse = ";"), length(cases), "\n", sep = "\t"))
        sink()

        return(cases)
}

##Cases vascular/heart problems diagnosed by doctor
cases_vascular_heart_doctor <- function(data, counts_file, diagnose, codes) {
        cases <- integer()
        for (i in seq_along(codes)) {
                code = codes[i]
                vascular_heart_doctor = extract_cases(data, code, "^6150\\.")
                all_cases_this_code <- unique(c(vascular_heart_doctor))
                all_diagnose_types_names <- list("vascular_heart_doctor",
                                "all_cases_this_code")
                all_diagnose_types <- list(vascular_heart_doctor, all_cases_this_code)

                for (j in seq_along(all_diagnose_types)) {
                        sink(counts_file, append = T)
                        cat(paste(diagnose, "vascular_heart_doctor", all_diagnose_types_names[j], code,
                                length(all_diagnose_types[[j]]), "\n", sep = "\t"))
                        sink()
                }

                cases <- c(cases, all_cases_this_code)
        }
        cases <- unique(cases)
        sink(counts_file, append = T)
        cat(paste(diagnose, "vascular_heart_doctor", "total_vascular_heart_doctor", paste(codes, collapse = ";"), length(cases), "\n", sep = "\t"))
        sink()

        return(cases)
}

##Cases DVT/COPD, etc diagnosed by doctor
cases_dvt_copd_doctor <- function(data, counts_file, diagnose, codes) {
        cases <- integer()
        for (i in seq_along(codes)) {
                code = codes[i]
                dvt_copd_doctor = extract_cases(data, code, "^6152\\.")
                all_cases_this_code <- unique(c(dvt_copd_doctor))
                all_diagnose_types_names <- list("dvt_copd_doctor",
                                "all_cases_this_code")
                all_diagnose_types <- list(dvt_copd_doctor, all_cases_this_code)

                for (j in seq_along(all_diagnose_types)) {
                        sink(counts_file, append = T)
                        cat(paste(diagnose, "dvt_copd_doctor", all_diagnose_types_names[j], code,
                                length(all_diagnose_types[[j]]), "\n", sep = "\t"))
                        sink()
                }

                cases <- c(cases, all_cases_this_code)
        }
        cases <- unique(cases)
        sink(counts_file, append = T)
        cat(paste(diagnose, "dvt_copd_doctor", "total_dvt_copd_doctor", paste(codes, collapse = ";"), length(cases), "\n", sep = "\t"))
        sink()

        return(cases)
}



#Get all cases
get_case_ids <- function(data, counts_file, diagnose, codes_icd10, 
		codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer, 
		codes_ni_operation, codes_medications, codes_chol_bp_dm_meds, 
		codes_vascular_heart_doctor, codes_dvt_copd_doctor) {
	
	ids_icd10 <- cases_icd10(data, counts_file, diagnose, codes_icd10) 
	ids_icd9 <- cases_icd9(data, counts_file, diagnose, codes_icd9)
	ids_opcs4 <- cases_opcs4(data, counts_file, diagnose, codes_opcs4)
	ids_ni_non_cancer <- cases_ni_non_cancer(data, counts_file, diagnose, codes_ni_non_cancer)
	ids_ni_cancer <- cases_ni_cancer(data, counts_file, diagnose, codes_ni_cancer)
	ids_ni_operation <- cases_ni_operation(data, counts_file, diagnose, codes_ni_operation)
	ids_medications <- cases_medications(data, counts_file, diagnose, codes_medications)
	ids_chol_bp_dm_meds <- cases_chol_bp_dm_meds(data, counts_file, diagnose, codes_chol_bp_dm_meds)
	ids_vascular_heart_doctor <- cases_vascular_heart_doctor(data, counts_file, diagnose, codes_vascular_heart_doctor)
	ids_dvt_copd_doctor <- cases_dvt_copd_doctor(data, counts_file, diagnose, codes_dvt_copd_doctor)

	all_cases <- unique(c(ids_icd10, ids_icd9, ids_opcs4, ids_ni_non_cancer, ids_ni_cancer, 
			ids_ni_operation, ids_medications, ids_chol_bp_dm_meds, ids_vascular_heart_doctor, 
			ids_dvt_copd_doctor))

        sink(counts_file, append = T)
        cat(paste(diagnose, "total_unique_cases", "total_unique_cases", "all_codes", length(all_cases), "\n", sep = "\t"))
        sink()

	return(all_cases)
}

#Prints the number of cases from each code for each diagnosis: 
counts_file <- "../ukbb_case_definitions_n_cases_table_180620.txt"
sink(counts_file, append = F)
cat(paste("diagnose\tdiagnose_coding\tcase_from_diagnose_type\tcode\tn\t\n", sep = ""))
sink()

#File with all the UKBB samples that pass QC, with columns with the data fields as detailed above
data <- read.table("/well/lindgren/jc/ukbb/ukbb.samples.passing.qc.relevant.pheno.plus.diagnoses.180509.txt",
                 stringsAsFactors = F, header = T)
colnames(data) <- gsub("^X", "", colnames(data))

#Subset to only British/Irish, recode white, and white:
data <- data[data$check_se %in% c("brit", "recode.white", "white"), ]

#SOFT CAD definition as per Nelson Deloukas 2017 Nature paper - CHECKED
#INDIVIDUALS TO INCLUDE IN SOFT CAD
diagnose = "SOFT_CAD"
codes_icd10 = c("I20", "I200", "I201", "I208", "I209",
		"I21", "I210", "I211", "I212", "I213", 
		"I214", "I219", "I21X", 
		"I22", "I220", "I221", "I228", "I229", 
		"I23", "I230", "I231", "I232", "I233", 
		"I234", "I235", "I236", "I238", 
		"I24", "I240", "I241", "I248", "I249", 
		"I251", "I252", 
		"I255", "I256", "I258", "I259")
codes_icd9 = c(410, 4109, 411, 4119, 412, 4129,
		413, 4139, 4140, 4148, 4149)
codes_opcs4 = c("K40", "K401", "K402", "K403", "K404", 
		"K408", "K409", 
		"K41", "K411", "K412", "K413", "K414", "K418", 
		"K419", 
		"K42", "K421", "K422", "K423", "K424", "K428", 
		"K429", 
		"K43", "K431", "K432", "K433", "K434", "K438", 
		"K439", 
		"K44", "K441", "K442", "K448", "K449", 
		"K45", "K451", "K452", "K453", "K454", 
		"K455", "K456", "K458", "K459", 
		"K46", "K461", "K462", "K463", "K464", "K465", 
		"K468", "K469", 
		"K49", "K491", "K492", "K493", "K494", "K498", 
		"K499", 
		"K501",  
		"K75", "K751", "K752", "K753", "K754", "K758", 
		"K759")
codes_ni_non_cancer = c(1074, 1075)
codes_ni_cancer = c()
codes_ni_operation = c(1070, 1095, 1523)
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c(1, 2)
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds, 
		codes_vascular_heart_doctor, codes_dvt_copd_doctor)

#INDIVIDUALS TO EXCLUDE FROM CONTROLS IN CAD:

diagnose = "exclude_CAD"
codes_icd10 = c("I250", "I253", "I254")
codes_icd9 = c(4141)
codes_opcs4 = c()
codes_ni_non_cancer = c()
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

exclude_cad_data <- data[!(data$ID %in% all_case_ids), ]

all_exclude_ids <- get_case_ids(exclude_cad_data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)


data$cad_cases <- ifelse(data$ID %in% all_case_ids, 1, ifelse(data$ID %in% all_exclude_ids, 
			2, 0))

#COPD - OURS -CHECKED
diagnose = "COPD"
codes_icd10 = c("J41", "J410", "J411", "J418",
		"J42", "J43", "J431", "J432",
		"J438", "J439", 
		"J44", "J440", "J441", "J448", "J449")
codes_icd9 = c(491, 4910, 4911, 4912, 4918, 4919, 
		492, 4929)
codes_opcs4 = c()
codes_ni_non_cancer = c(1112, 1113, 1472)
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c(6)

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$copd_cases <- ifelse(data$ID %in% all_case_ids, 1, 0)

#Dementia - own - CHECKED
diagnose = "dementia"
codes_icd10 = c("F00", "F000", "F001", "F002", "F009", 
		"F01", "F010", "F011", "F012", "F013", 
		"F018", "F019", 
		"F03", 
		"G30", "G300", "G301", "G308", "G309")
codes_icd9 = c(2900, 2901, 2902, 2903, 2904, 3310)
codes_opcs4 = c()
codes_ni_non_cancer = c(1263)
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$dementia_cases <- ifelse(data$ID %in% all_case_ids, 1, 0)

#Trachea, bronchus, and lung cancer - our definition - CHECKED
diagnose = "trachea_bronchus_lung_cancer"
codes_icd10 = c("C33", 
		"C34", "C340", "C341", "C342", "C343",
		"C348", "C349", 
		"Z851")
codes_icd9 = c(162, 1620, 1622, 1623, 1624, 1625, 1628, 
		1629, "V101")
codes_opcs4 = c()
codes_ni_non_cancer = c()
codes_ni_cancer = c(1001, 1027, 1028, 1080)
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$lungcancer_cases <- ifelse(data$ID %in% all_case_ids, 1, 0)

#Colorectal cancer - own - CHECKED
diagnose = "colorectal_cancer"
codes_icd10 = c("C18", "C180", "C181", "C182", "C183", "C184", "C185",
		"C186", "C187", "C188", "C189", 
		"C19", 
		"C20")
codes_icd9 = c(153, 1530, 1531, 1532, 1533, 1534, 1535, 1536, 1537, 
		1538, 1539, 
		154, 1540, 1541)
codes_opcs4 = c()
codes_ni_non_cancer = c()
codes_ni_cancer = c(1020, 1022, 1023)
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$colorectal_cancer_cases <- ifelse(data$ID %in% all_case_ids, 1, 0)

#Breast cancer - own - CHECKED
diagnose = "breast_cancer_women_only"
codes_icd10 = c("C50", "C500", "C501", "C502", "C503", "C504", 
		"C505", "C506", "C508", "C509", 
		"Z853")
codes_icd9 = c(174, 1740, 1741, 1742, 1743, 1744, 1745, 1746, 
		1748, 1749, "V103")
codes_opcs4 = c()
codes_ni_non_cancer = c()
codes_ni_cancer = c(1002)
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

#Restrict to women only
breastcancer_data <- data[data$Submitted_Gender == "F", ]

all_case_ids <- get_case_ids(breastcancer_data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$breast_cancer_cases <- ifelse(data$Submitted_Gender != "F", 2, ifelse(data$ID %in% all_case_ids, 1, 0))

#Female infertility - own - CHECKED
diagnose = "female_infertility"
codes_icd10 = c("N97", "N970", "N971", "N972", "N973", "N978", "N979")
codes_icd9 = c(628, 6280, 6281, 6282, 6283, 6284, 6288, 6289)
codes_opcs4 = c()
codes_ni_non_cancer = c(1403)
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

female_infertility_data <- data[data$Submitted_Gender == "F", ]
all_case_ids <- get_case_ids(female_infertility_data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$female_infertility_cases <- ifelse(data$Submitted_Gender != "F", 2, ifelse(data$ID %in% all_case_ids, 1, 0))

#Male infertility - own - CHECKED
diagnose = "male_infertility"
codes_icd10 = c("N46")
codes_icd9 = c(606, 6069)
codes_opcs4 = c()
codes_ni_non_cancer = c(1404)
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

male_infertility_data <- data[data$Submitted_Gender == "M", ]
all_case_ids <- get_case_ids(male_infertility_data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$male_infertility_cases <- ifelse(data$Submitted_Gender != "M", 2, ifelse(data$ID %in% all_case_ids, 1, 0))

#General infertility - CHECKED - adds both the male and female cases together
data$any_infertility_cases <- ifelse((data$female_infertility_cases == 1 | data$male_infertility_cases == 1), 1, 0)
sink(counts_file, append = T)
cat(paste("any_infertility", "male_or_female_case", "female_case", "female_case", nrow(data[data$female_infertility_cases == 1, ]), "\n", sep = "\t"))
cat(paste("any_infertility", "male_or_female_case", "male_case", "male_case", nrow(data[data$male_infertility_cases == 1, ]), "\n", sep = "\t"))
cat(paste("any_infertility", "total_unique_cases", "total_unique_cases", "all_codes", nrow(data[data$any_infertility_cases == 1, ]), "\n", sep = "\t"))
sink()

#NAFLD - Own
diagnose = "nafld"
codes_icd10 = c("K760")
codes_icd9 = c()
codes_opcs4 = c()
codes_ni_non_cancer = c()
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$nafld_cases <- ifelse(data$ID %in% all_case_ids, 1, 0)

#CLD - own
diagnose = "chronic_liver_disease"
codes_icd10 = c("K702", "K703", "K704", 
		"K717",
		"K721", 
		"K74", "K740", "K741", "K742", "K743", "K744", 
		"K745", "K746")
codes_icd9 = c(27103, 4562, 571, 5712, 5715, 57150, 57151, 57158, 
		57159, 5716)
codes_opcs4 = c()
codes_ni_non_cancer = c(1604, 1158)
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$cld_cases <- ifelse(data$ID %in% all_case_ids, 1, 0)

#Exclude renal failure/aki/ckd
diagnose = "renal_exclude"
codes_icd10 = c("N17", "N170", "N171", "N172", "N178", "N179",
		"N18", "N180", "N181", "N182", "N183", "N184",
                "N185", "N188", "N189", 
		"N19", "I120", "I131", "I132", 
		"Z992")
codes_icd9 = c(584, 5845, 5846, 5847, 5848, 5849,
		585, 5859,
		586, 5869)
codes_opcs4 = c("L746", "X40", "X401", "X402", "X403", "X404", "X405",
                        "X406", "X407", "X408", "X409",
                        "X41", "X411", "X412", "X418", "X419",
                        "X42", "X421", "X428", "X429")
codes_ni_non_cancer = c(1192, 1193, 1194)
codes_ni_cancer = c()
codes_ni_operation = c(1476, 1580, 1581, 1582)
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

renal_exclude_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

#Kidney disease - own definition  - CHECKED
diagnose = "renal_failure"
codes_icd10 = c("N17", "N170", "N171", "N172", "N178", "N179",
                "N18", "N180", "N181", "N182", "N183", "N184",
                "N185", "N188", "N189",
                "N19",
                "I120", "I131", "I132",
                "Z992")
codes_icd9 = c(584, 5845, 5846, 5847, 5848, 5849,
                585, 5859,
                586, 5869)
codes_opcs4 = c("L746", "X40", "X401", "X402", "X403", "X404", "X405",
                        "X406", "X407", "X408", "X409",
                        "X41", "X411", "X412", "X418", "X419",
                        "X42", "X421", "X428", "X429")
codes_ni_non_cancer = c(1192, 1193, 1194)
codes_ni_cancer = c()
codes_ni_operation = c(1476, 1580, 1581, 1582)
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$renal_failure_cases <- ifelse(data$ID %in% all_case_ids, 1, 
				ifelse(data$ID %in% renal_exclude_ids, 2, 0))

#Include acute kidney disease
diagnose = "aki"
codes_icd10 = c("N17", "N170", "N171", "N172", "N178", "N179")
codes_icd9 = c(584, 5845, 5846, 5847, 5848, 5849)
codes_opcs4 = c()
codes_ni_non_cancer = c()
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$aki_cases <- ifelse(data$ID %in% all_case_ids, 1, 
			ifelse(data$ID %in% renal_exclude_ids, 2, 0))

#Chronic kidney disease, own definition
diagnose = "ckd"
codes_icd10 = c("N18", "N180", "N181", "N182", "N183", "N184",
                "N185", "N188", "N189")
codes_icd9 = c(585, 5859)
codes_opcs4 = c()
codes_ni_non_cancer = c()
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$ckd_cases <- ifelse(data$ID %in% all_case_ids, 1, 
			ifelse(data$ID %in% renal_exclude_ids, 2, 0))

#Stroke exclude ids
diagnose = "stroke_exclude"
codes_icd10 = c("I60", "I600", "I601", "I602", "I603", "I604",
                "I605", "I606", "I607", "I608", "I609",
                "I61", "I610", "I611", "I612", "I613", "I614",
                "I615", "I616", "I618", "I619",
                "I63", "I630", "I631", "I632", "I633", "I634",
                "I635", "I636", "I638", "I639",
                "I64", 
		"G45", "G450", "G451", "G452", "G453", "G454", 
		"G458", "G459")
codes_icd9 = c(430, 4309, 431, 4319, 434, 4340, 4341, 4349,
                436, 4369, 
		435, 4359)
codes_opcs4 = c()
codes_ni_non_cancer = c(1081, 1082, 1086, 1491, 1583)
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c(3)
codes_dvt_copd_doctor = c()

stroke_exclude_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

#Stroke
diagnose = "stroke"
codes_icd10 = c("I60", "I600", "I601", "I602", "I603", "I604",
                "I605", "I606", "I607", "I608", "I609",
                "I61", "I610", "I611", "I612", "I613", "I614",
                "I615", "I616", "I618", "I619",
                "I63", "I630", "I631", "I632", "I633", "I634",
                "I635", "I636", "I638", "I639",
                "I64")
codes_icd9 = c(430, 4309, 431, 4319, 434, 4340, 4341, 4349,
                436, 4369)
codes_opcs4 = c()
codes_ni_non_cancer = c(1081, 1086, 1491, 1583)
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c(3)
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$stroke_cases <- ifelse(data$ID %in% all_case_ids, 1, 
			ifelse(data$ID %in% stroke_exclude_ids, 2, 0))

#Haemorragic stroke cases
diagnose = "haem_stroke"
codes_icd10 = c("I60", "I600", "I601", "I602", "I603", "I604",
                "I605", "I606", "I607", "I608", "I609",
                "I61", "I610", "I611", "I612", "I613", "I614",
                "I615", "I616", "I618", "I619")
codes_icd9 = c(430, 4309, 431, 4319)
codes_opcs4 = c()
codes_ni_non_cancer = c(1086, 1491)
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$haem_stroke_cases <- ifelse(data$ID %in% all_case_ids, 1, 
				ifelse(data$ID %in% stroke_exclude_ids, 2, 0))

#Ischaemic stroke cases
diagnose = "isch_stroke"
codes_icd10 = c("I63", "I630", "I631", "I632", "I633", "I634",
                "I635", "I636", "I638", "I639")
codes_icd9 = c(434, 4340, 4341, 4349)
codes_opcs4 = c()
codes_ni_non_cancer = c(1583)
codes_ni_cancer = c()
codes_ni_operation = c()
codes_medications = c()
codes_chol_bp_dm_meds = c()
codes_vascular_heart_doctor = c()
codes_dvt_copd_doctor = c()

all_case_ids <- get_case_ids(data, counts_file, diagnose, codes_icd10,
                codes_icd9, codes_opcs4, codes_ni_non_cancer, codes_ni_cancer,
                codes_ni_operation, codes_medications, codes_chol_bp_dm_meds,
                codes_vascular_heart_doctor, codes_dvt_copd_doctor)

data$isch_stroke_cases <- ifelse(data$ID %in% all_case_ids, 1, 
				ifelse(data$ID %in% stroke_exclude_ids, 2, 0))

#T2D - probable only, from Eastwood et al script
data$t2d_cases_prob <- ifelse(data$sr_prob_t2_diabetes == "Probable_type_2_diabetes", 1, 
			ifelse(data$sr_unlikely_diabetes == "Diabetes_unlikely", 0, 2))
sink(counts_file, append = T)
cat(paste("T2D_prob", "t2d_stata_script", "t2d_prob_inclusion", "prob_t2d", nrow(data[data$sr_prob_t2_diabetes == "Probable_type_2_diabetes", ]), "\n", sep = "\t"))
cat(paste("T2D_prob", "total_unique_cases", "total_unique_cases", "all_codes", nrow(data[data$sr_prob_t2_diabetes == "Probable_type_2_diabetes", ]), "\n", sep = "\t"))
cat(paste("T2D_prob", "t2d_stata_script", "t2d_prob_exclusion", "not_diabetes_unlikely", nrow(data[data$t2d_cases_prob == 2, ]), "\n", sep = "\t"))
sink()

#T2D - prob and possible, from Anubha's script
data$t2d_cases_probposs <- ifelse((data$sr_prob_t2_diabetes == "Probable_type_2_diabetes") | (data$sr_poss_t2_diabetes == "Possible_type_2_diabetes"), 1,
                        ifelse(data$sr_unlikely_diabetes == "Diabetes_unlikely", 0, 2))
sink(counts_file, append = T)
cat(paste("T2D_probposs", "t2d_stata_script", "t2d_probposs_inclusion", "probable_t2d", nrow(data[data$sr_prob_t2_diabetes == "Probable_type_2_diabetes", ]), "\n", sep = "\t"))
cat(paste("T2D_probposs", "t2d_stata_script", "t2d_probposs_inclusion", "possible_t2d", nrow(data[data$sr_poss_t2_diabetes == "Possible_type_2_diabetes", ]), "\n", sep = "\t"))
cat(paste("T2D_probposs", "total_unique_cases", "total_unique_cases", "all_codes", nrow(data[data$t2d_cases_probposs == 1, ]), "\n", sep = "\t"))
cat(paste("T2D_probposs", "t2d_stata_script", "t2d_prob_exclusion", "not_diabetes_unlikely", nrow(data[data$t2d_cases_probposs == 2, ]), "\n", sep = "\t"))
sink()

#T1D - probable, from Eastwood et als script
data$t1d_cases_prob <- ifelse(data$sr_prob_t1_diabetes == "Probable_type_1_diabetes", 1,
                        ifelse(data$sr_unlikely_diabetes == "Diabetes_unlikely", 0, 2))
sink(counts_file, append = T)
cat(paste("T1D_prob", "t1d_stata_script", "t1d_prob_inclusion", "probable_t1d", nrow(data[data$t1d_cases_prob == 1, ]), "\n", sep = "\t"))
cat(paste("T1D_prob", "total_unique_cases", "total_unique_cases", "all_codes", nrow(data[data$t1d_cases_prob == 1, ]), "\n", sep = "\t"))
cat(paste("T1D_prob", "t1d_stata_script", "t1d_prob_exclusion", "not_diabetes_unlikely", nrow(data[data$t1d_cases_prob == 2, ]), "\n", sep = "\t"))
sink()

#T1D - probable and possible, from Anubha's script
data$t1d_cases_probposs <- ifelse((data$sr_prob_t1_diabetes == "Probable_type_1_diabetes") | (data$sr_poss_t1_diabetes == "Possible_type_1_diabetes"), 1,
                        ifelse(data$sr_unlikely_diabetes == "Diabetes_unlikely", 0, 2))
sink(counts_file, append = T)
cat(paste("T1D_probposs", "t1d_stata_script", "t1d_prob_inclusion", "probable_t1d", nrow(data[data$sr_prob_t1_diabetes == "Probable_type_1_diabetes", ]), "\n", sep = "\t"))
cat(paste("T1D_probposs", "t1d_stata_script", "t1d_prob_inclusion", "possible_t1d", nrow(data[data$sr_poss_t1_diabetes == "Possible_type_1_diabetes", ]), "\n", sep = "\t"))
cat(paste("T1D_probposs", "total_unique_cases", "total_unique_cases", "all_codes", nrow(data[data$t1d_cases_probposs == 1, ]), "\n", sep = "\t"))
cat(paste("T1D_probposs", "t1d_stata_script", "t1d_prob_exclusion", "not_diabetes_unlikely", nrow(data[data$t1d_cases_probposs == 2, ]), "\n", sep = "\t"))
sink()

#Smoking - if at baseline reported to have smoked or be a smoker in any column, even "occasionally"
smoker_cases <- extract_cases(data, c(1,2), "^1239\\.0_0|^1249\\.0_0|^20116\\.0_0")
data$smoker_cases <- ifelse(data$ID %in% smoker_cases, 1, 0)
sink(counts_file, append = T)
cat(paste("smoker", "self_report", "self_report", "1,2", nrow(data[data$smoker_cases == 1, ]), "\n", sep = "\t"))
cat(paste("smoker", "total_unique_cases", "total_unique_cases", "all_codes", nrow(data[data$smoker_cases == 1, ]), "\n", sep = "\t"))
sink()

#Subset dataframe to only contain relevant columns
data <- data[, -grep("^([0-9])+\\.([0-9])+", names(data))]

write.table(data, "../ukbb_cases_for_logistic_regression_180620.txt",
                row.names = F, sep = " ", quote = F)

#To add to this dataset without redoing everything:
#df <- read.table("../ukbb_cases_for_logistic_regression_180620.txt", stringsAsFactors = F, header = T, sep = " ")
#data <- read.table("/well/lindgren/jc/ukbb/ukbb.samples.passing.qc.relevant.pheno.plus.diagnoses.180509.txt",
#                 stringsAsFactors = F, header = T)
#colnames(data) <- gsub("^X", "", colnames(data))
#Subset to only British/Irish, recode white, and white:
#data <- data[data$check_se %in% c("brit", "recode.white", "white"), ]
#data <- data[, grep("^ID$|^([0-9])+\\.([0-9])+", names(data))]
#data <- merge(df, data, by.x ="ID", by.y = "ID")
#
