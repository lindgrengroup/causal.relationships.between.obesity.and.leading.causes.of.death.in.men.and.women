#!/bin/env python
#$ -cwd

######################################################################
######## Script to get the relevant phenotype columns in python ######
###################### Jenny Censin, 2018-04-13 ######################

import pandas as pd

#Make functions for extracting different columns
#Primary ICD10, secondary ICD10, primary death, secondary death
def icd10_columns(pheno):
	df_icd10 = pheno.filter(regex="^41202\-|^41204\-|^40001\-|^40002\-")
	return df_icd10

#Primary ICD9, secondary ICD9
def icd9_columns(pheno):
	df_icd9 = pheno.filter(regex="^41203\-|^41205\-")
	return df_icd9

#Main OPCS-4 codes, secondary OPCS-4 codes
def opcs4_columns(pheno):
	df_opcs4 = pheno.filter(regex="^41200\-|^41210\-")
	return df_opcs4

#Medication codes
def med_columns(pheno):
	df_med = pheno.filter(regex="^20003\-")
	return df_med

#Self-reported nurse-interview non-cancer codes
def selfrep_columns(pheno):
	df_selfrep = pheno.filter(regex="^20002\-")
	return df_selfrep

#Self-reported nurse-interview cancer codes
def selfrep_cancer_columns(pheno):
	df_selfrep_cancer = pheno.filter(regex="^20001\-")
	return df_selfrep_cancer

#Self-reported operation codes
def selfrep_operation_columns(pheno):
	df_selfrep_operation = pheno.filter(regex="^20004\-")
	return df_selfrep_operation

#GDM only
def gdm_columns(pheno):
	df_gdm = pheno.filter(regex="^4041\-")
	return df_gdm

#Age dm diagnosis
def age_dm_columns(pheno):
	df_age_dm = pheno.filter(regex="^2976\-")
	return df_age_dm

#Specific cholesterol/HT/DM medications for men and women, respectively
def chol_bp_dm_med_columns(pheno):
	df_chol_bp_dm_med = pheno.filter(regex="^6177\-|^6153\-")
	return df_chol_bp_dm_med

#Started insulin within one year of diagnosis
def start_insulin_columns(pheno):
	df_start_insulin = pheno.filter(regex="^2986\-")
	return df_start_insulin

#Age non-cancer illness first reported
def age_illness_columns(pheno):
	df_age_illness = pheno.filter(regex="^20009\-")
	return df_age_illness

#Age cancer illness first reported
def age_cancer_columns(pheno):
	df_age_cancer = pheno.filter(regex="^20007\-")
	return df_age_cancer

#Systolic, diastolic blood pressure, automated, and pulse, automated
def bp_pulse_columns(pheno):
	df_bp_pulse = pheno.filter(regex="^4080\-|^4079\-|^102\-")
	return df_bp_pulse

#Source of death
def source_death_columns(pheno):
	df_source_death = pheno.filter(regex="^40018\-")
	return df_source_death

#Cancer ICD-10 codes from cancer registry
def cancer_icd10_registry_columns(pheno):
	df_cancer_icd10_registry = pheno.filter(regex="^40006\-")
	return df_cancer_icd10_registry

#Cancer ICD-9 codes from cancer registry
def cancer_icd9_registry_columns(pheno):
	df_cancer_icd9_registry = pheno.filter(regex="^40013\-")
	return df_cancer_icd9_registry

#Vascular/heart problems diagnosed by doctor
def vascular_heart_doctor_columns(pheno):
	df_vascular_heart_doctor = pheno.filter(regex="^6150\-")
	return df_vascular_heart_doctor

#DVT, COPD, asthma etc problems diagnosed by doctor
def dvt_copd_doctor_columns(pheno):
        df_dvt_copd_doctor = pheno.filter(regex="^6152\-")
        return df_dvt_copd_doctor

#Smoking columns
def smoking_columns(pheno):
	df_smoking = pheno.filter(regex="^1239\-|1249\-|20116\-")
	return df_smoking

#Pathway to the list of samples that have passed QC, made previously
samples = pd.read_csv("/well/lindgren/jc/ukbb/ukbb.samples.passing.qc.relevant.pheno.all.ancestries.180504.txt",
                delimiter = " ", index_col = 'ID')
del samples.index.name

#Pathway to the main UKBB provided phenotype file	 
pheno = pd.read_csv("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv", 
	delimiter = ",", header = 0, index_col = 'eid')
del pheno.index.name

#Pathway to the file made from the Eastwood et al algorithm, with diabetes case/control status:
t2d = pd.read_csv("/well/lindgren/jc/ukbb/T2D_adjudicated_2.txt", 
	delimiter = "\t", header = 0, index_col = "id")
del t2d.index.name
t2d = t2d.replace(to_replace = ", | |,|/", value = "_", regex = True)

df_icd10 = icd10_columns(pheno)
df_icd9 = icd9_columns(pheno)
df_opcs4 = opcs4_columns(pheno)
df_med = med_columns(pheno)
df_selfrep = selfrep_columns(pheno)
df_selfrep_cancer = selfrep_cancer_columns(pheno)
df_selfrep_operation = selfrep_operation_columns(pheno)
df_gdm = gdm_columns(pheno)
df_age_dm = age_dm_columns(pheno)
df_chol_bp_dm_med = chol_bp_dm_med_columns(pheno)
df_start_insulin = start_insulin_columns(pheno)
df_age_illness = age_illness_columns(pheno)
df_age_cancer = age_cancer_columns(pheno)
df_bp_pulse = bp_pulse_columns(pheno)
df_source_death = source_death_columns(pheno)
df_cancer_icd10_registry = cancer_icd10_registry_columns(pheno)
df_cancer_icd9_registry = cancer_icd9_registry_columns(pheno)
df_vascular_heart_doctor = vascular_heart_doctor_columns(pheno)
df_dvt_copd_doctor = dvt_copd_doctor_columns(pheno)
df_smoking = smoking_columns(pheno)

#Merge the datasets
df = pd.concat([samples, df_icd10, df_icd9, df_opcs4, df_med, df_selfrep, 
	df_selfrep_cancer, df_selfrep_operation, df_gdm, df_age_dm, 
	df_chol_bp_dm_med, df_start_insulin, df_age_illness, df_age_cancer, 
	df_bp_pulse, df_source_death, df_cancer_icd10_registry, df_cancer_icd9_registry, 
	df_vascular_heart_doctor, df_dvt_copd_doctor, df_smoking], 
	axis=1, join_axes=[samples.index])

#Substitute dots with underscore to avoid python errors
df.columns = df.columns.str.replace('\.', '_')

df = df.merge(t2d, left_index = True, right_index = True)

#Write to a file that contains all the relevant phenotype information
df.to_csv("../ukbb.samples.passing.qc.relevant.pheno.plus.diagnoses.180509.txt", sep = ' ', 
	index = True, index_label = "ID", na_rep = "NA")



