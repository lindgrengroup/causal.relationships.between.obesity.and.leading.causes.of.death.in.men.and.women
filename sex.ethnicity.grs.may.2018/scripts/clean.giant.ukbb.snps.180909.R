#!/bin/env Rscript
#$-cwd

#Function to remove the secondary signals, based on SNPs that are on the same chromosome, 
#within 1mb up/downstreams of each other, keeping the one with the lowest p-value
#from the p-value column provided
remove.secondary <- function(df, p) {
        df <- df[order(df$Chr, df$Pos), ]
        i <- 1
        if(nrow(df) >=2) {
                while (i < nrow(df)){
			print(i)
                        if ((df[i, "Chr"] == df[i+1, "Chr"]) &
                                (abs(df[i, "Pos"] - df[i+1, "Pos"]) < 1000000) &
					(df[i, p] < df[i+1, p])) {
					print(df[c(i, i+1),])
                                        df <- df[-(i+1), ]
                        } else if ((df[i, "Chr"] == df[i+1, "Chr"]) &
                                (abs(df[i, "Pos"] - df[i+1, "Pos"]) < 1000000) &
                                (df[i, p] >= df[i+1, p])) {
					print(df[c(i, i+1), ])
                                        df <- df[-i, ]
                                if (i > 1) {
                                        i <- i - 1
                                }
                        } else {
                                i <- i + 1
                        }
                }
                return(df)
        }
}

make.index.files.correct <- function(trait, eth_group, sex_group) {
	#Pathway to the original Pulit index files. The files in original.pulit.files.from.github.180911, 
	#eg bmi.eur.comb.pulit.180513.txt, were downloaded from Sara's
	#github page https://github.com/lindgrengroup/fatdistnGWAS on 2018-09-11.
	input_file <- paste("../original.index.files/", trait, ".", eth_group,
                ".", sex_group, ".pulit.180513.txt", sep = "")
        df <- read.table(input_file, header = T,
                stringsAsFactors = F)
	colnames(df)[colnames(df) %in% c("Chr.ref.males", "Pos.ref.males", "A1.combined", "A2.combined")] <- c("Chr", "Pos", "EA", "NEA")

        df$SNP <- gsub(":[ATCG]:[ATCG]", "", df$SNP)

	df[, c("alpha_first", "alpha_second")] <- t(apply(df[, c("EA", "NEA")], 1, sort))
	df$SNP <- paste(df$SNP, "_", df$alpha_first, "_",
                df$alpha_second, sep ="")
        df[, c("alpha_first", "alpha_second")] <- NULL

	return(df)
}

make.beta.pos <- function(df) {
	for (m in 1:nrow(df)) {
		if (df[m, grep("beta\\.", colnames(df))] < 0) {
			df[m, c("EA", "NEA")] <- df[m, c("NEA", "EA")]
			df[m, grep("^beta\\.", colnames(df))] <- df[m, grep("^beta\\.", colnames(df))] * -1
			df[m, grep("^frqA1\\.", colnames(df))] <- 1 - df[m, grep("^frqA1\\.", colnames(df))] 
		}
	}
	return(df)
}

loop <- expand.grid(trait = c("bmi", "whr", "whradjbmi"), sex_group = c("comb", "men", "women"),
                snp_group = c("pulit"), eth_group = "eur", stringsAsFactors = F)

for (k in 1:nrow(loop)) {
	comb <- make.index.files.correct(loop[k, "trait"], loop[k, "eth_group"], "comb")
	men <- make.index.files.correct(loop[k, "trait"], loop[k, "eth_group"], "men")
	women <-  make.index.files.correct(loop[k, "trait"], loop[k, "eth_group"], "women")

        columns <- ifelse(loop[k, "sex_group"] == "comb", "\\.combined$",
                        ifelse(loop[k, "sex_group"] == "women", "\\.females$",
                        ifelse(loop[k, "sex_group"] == "men", "\\.males$", NA)))

	output_file_pulit <- paste("../index.files.with.trialleles/", loop[k, "trait"], ".", loop[k, "eth_group"], ".", loop[k, "sex_group"],
                        ".pulit.180513.txt", sep = "")
	if (loop[k, "sex_group"] == "men") {
		men_pulit <- men[, c(grep("^SNP$|^Chr$|^Pos$|^EA$|^NEA$|^psexdiff$", colnames(men)), grep(columns, colnames(men)))]
		men_pulit <- make.beta.pos(men_pulit)
		write.table(men_pulit, output_file_pulit, row.names = F, quote = F, sep = " ")
	} else if (loop[k, "sex_group"] == "women") {
                women_pulit <- women[, c(grep("^SNP$|^Chr$|^Pos$|^EA$|^NEA$|^psexdiff$", colnames(women)), grep(columns, colnames(women)))]
		women_pulit <- make.beta.pos(women_pulit)
                write.table(women_pulit, output_file_pulit, row.names = F, quote = F, sep = " ")
	} 
	
	output_file_pulit_winner <- paste("../index.files.with.trialleles/", loop[k, "trait"], ".", loop[k, "eth_group"], ".", 
			loop[k, "sex_group"], ".pulit.winner.180513.txt", sep = "")
 	winner_df <- comb[, c(grep("^SNP$|^Chr$|^Pos$|^EA$|^NEA$|^psexdiff$", colnames(comb)), grep(columns, colnames(comb)))]
	winner_df <- make.beta.pos(winner_df)
	write.table(winner_df, output_file_pulit_winner, row.names = F, quote = F, sep = " ")	

	df <- merge(merge(comb, men, by = colnames(comb), all = T), women, by = colnames(women), all = T)
	df <- remove.secondary(df, "pval.combined")

        df_sex <- df[, c(grep("^SNP$|^Chr$|^Pos$|^EA$|^NEA$|^psexdiff$", colnames(df)), grep(columns, colnames(df)))]
	df_sex <- make.beta.pos(df_sex)
	output_file_sig <- paste("../index.files.with.trialleles/", loop[k, "trait"], ".", loop[k, "eth_group"], ".", loop[k, "sex_group"], 
			".pulit.sig.180513.txt", sep = "")		
	write.table(df_sex, output_file_sig, row.names = F, quote = F, sep = " ")
	
	if (loop[k, "sex_group"] %in% c("women", "men")) {
		critical_p <- 0.05/nrow(df)
		df_phet <- df
		
		#Make file with all sexually heterogeneous and not sexually heterogenous SNPs
		sex_het_snps <- df_phet[df_phet$psexdiff < critical_p, "SNP"]
		sex_het_output_file <- paste("../sex.het.enrichment/sex.het.snps.", loop[k, "trait"], ".pulit.180921.txt", sep = "")
		write.table(sex_het_snps, sex_het_output_file, quote = F, row.names = F, sep = "\t", col.names = F)
		not_sex_het_snps <- df_phet[df_phet$psexdiff >= critical_p, "SNP"]
		not_sex_het_output_file <- paste("../sex.het.enrichment/not.sex.het.snps.", loop[k, "trait"], ".pulit.180921.txt", sep = "")
		write.table(not_sex_het_snps, not_sex_het_output_file, quote = F, row.names = F, sep = "\t", col.names = F)		

		phet.columns <- gsub(".combined", ".phet", colnames(df[, grep("\\.combined$", names(df_phet))]))
		df_phet[, phet.columns] <- as.numeric(NA)
		
		for (row in 1:nrow(df_phet)) {
			if (df_phet$psexdiff[row] < critical_p) {
				df_phet[row, phet.columns] <- df_phet[row, grep(columns, names(df_phet))]
			} else if (df_phet$psexdiff[row] >= critical_p) {
				df_phet[row, phet.columns] <- df_phet[row, grep("\\.combined$", names(df_phet))]
			}
		}
		df_phet <- df_phet[, grep("^SNP$|^Chr$|^Pos$|^EA$|^NEA$|^psexdiff$|\\.phet$", names(df_phet))]
		df_phet <- make.beta.pos(df_phet)
		output_file_phet <- paste("../index.files.with.trialleles/", loop[k, "trait"], ".", loop[k, "eth_group"], ".", loop[k, "sex_group"], ".pulit.phet.180513.txt", sep = "")
		write.table(df_phet, output_file_phet, quote = F, row.names = F, sep = " ")
	}
}

#Create a list of SNPs that are in any analysis for extraction
loop <- expand.grid(trait = c("bmi", "whr", "whradjbmi"), sex_group = c("comb", "men", "women"),
	snp_group = c("pulit", "pulit.sig", "pulit.phet", "pulit.winner"), eth_group = "eur", stringsAsFactors = F)
loop <- loop[!(loop$sex_group == "comb" & loop$snp_group %in% c("pulit.phet", "pulit")), ]

results_df <- data.frame("SNP" = character())
for (n in 1:nrow(loop)) {
	input_file <- paste("../index.files.with.trialleles/", loop[n, "trait"], ".", loop[n, "eth_group"],
                ".", loop[n, "sex_group"], ".", loop[n, "snp_group"], ".180513.txt", sep = "")
	print(input_file)
	df <- read.table(input_file, stringsAsFactors = F, header = T)
	df <- data.frame(df[, "SNP"], stringsAsFactors = F)
	colnames(df) <- "SNP"
	results_df <- merge(results_df, df, by = "SNP", all = T)
}

#Read in the list of SNPs that pass all QC -  one SNP rs144926207 is
#called by chr:pos in UKBB: 9:140079779_T_C, hence changed
#when creating all the index SNPs 
results_df$SNP[results_df$SNP == "rs144926207_C_T"] <- "9:140079779_T_C_C_T"

qc <- read.table("/well/lindgren/jc/ukbb/genotype.qc.snplists.180507/snps_that_pass_all_qc_180508.txt", 
	stringsAsFactors = F, header = F, colClasses = "character")
in.qc <- results_df[results_df$SNP %in% qc$V1, ]

in.qc <- gsub("_[ATCG]+_[ATCG]+$", "", in.qc)
write.table(in.qc, "../temporary.bed.pulit.180911/snps.to.extract.txt", sep = "\t", row.names = F, quote = F, col.names = F)

not.in.qc <- results_df[!(results_df$SNP %in% qc$V1), ]
write.table(not.in.qc, "../temporary.bed.pulit.180911/snps.that.are.not.in.qc.list.giukbb.txt", sep = "\t", row.names = F, 
	quote = F, col.names = F)

#Make list of the individuals to include in the calculations
samples <- read.table("/well/lindgren/jc/ukbb/ukbb.samples.passing.qc.relevant.pheno.all.ancestries.180504.txt", 
		stringsAsFactors = F, header = T)
write.table(samples[, "ID"], "../temporary.bed.pulit.180911/samples.to.keep.txt", sep = "\t", row.names = F, quote = F, col.names = F)

