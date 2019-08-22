#!/bin/env python
#$ -cwd

######################################################################
########## Script to extract columns for T2D STATA script ###########
###################### Jenny Censin, 2018-06-20 ######################

import pandas as pd
import numpy as np

#Pathway to the main phenotype file
df = pd.read_csv("/well/lindgren/UKBIOBANK_DATA_LINDGREN/Phenotype_data/October2017/ukb10844.csv",
        delimiter = ",", header = 0, index_col = 'eid')
del df.index.name
df.columns = df.columns.str.replace('\.', '_')

#Pathway to list of individuals that have withdrawn consent
with open('../ids.withdrawn.consent.180511.w11867_20180503.txt') as f:
    withdrew_consent = f.read().splitlines()
df = df.loc[~df.index.isin(withdrew_consent)]

######
###### MAKE DATAFRAME SMALLER FOR SPEED
######

#Sex columns
def sex_columns(pheno):
	df_sex = pheno.filter(regex="^31\-")
	return df_sex

#Ethnicity columns
def ethnic_columns(pheno):
	df_ethnic = pheno.filter(regex="^21000\-")
	return df_ethnic

#Medication codes
def med_columns(pheno):
        df_med = pheno.filter(regex="^20003\-")
        return df_med

#Self-reported nurse-interview non-cancer codes
def selfrep_columns(pheno):
        df_selfrep = pheno.filter(regex="^20002\-")
        return df_selfrep

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

df_sex = sex_columns(df)
df_ethnic = ethnic_columns(df)
df_med = med_columns(df)
df_selfrep = selfrep_columns(df)
df_gdm = gdm_columns(df)
df_age_dm = age_dm_columns(df)
df_chol_bp_dm_med = chol_bp_dm_med_columns(df)
df_start_insulin = start_insulin_columns(df)
df_age_illness = age_illness_columns(df)

df = pd.concat([df_sex, df_ethnic, df_med, df_selfrep, df_gdm, df_age_dm,
        df_chol_bp_dm_med, df_start_insulin, df_age_illness],
        axis=1, join_axes=[df_sex.index])

df.columns = df.columns.str.replace('-', '_')
df = df.add_prefix("n_")

#Write smaller phenotype file with all the relevant columns for the T2D stata script
df.to_csv("../ukbb_columns_needed_for_t2d_script_180620.txt", sep = ' ',
        index = True, index_label = "ID", na_rep = "NA")




