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

duplicated <- read.table("/well/lindgren/jc/ukbb/genotype.qc.snplists.180507/non.biallelic.snps.txt",
                stringsAsFactors = F, header = F)
qc <- read.table("/well/lindgren/jc/ukbb/genotype.qc.snplists.180507/snps_that_pass_all_qc_180508.txt",
                stringsAsFactors = F, header = F)
qc[qc$V1 == "9:140079779_T_C_C_T", ] <- "rs144926207_C_T"

ld <- read.table("../snps.that.are.correlated.over.0.05.181123.txt", stringsAsFactors = F, header = T)

loop <- expand.grid(trait = c("bmi", "whr", "whradjbmi"), sex_group = c("comb", "men", "women"),
                snp_group = c("pulit"), eth_group = "eur", stringsAsFactors = F)

for (k in 1:nrow(loop)) {
	comb <- make.index.files.correct(loop[k, "trait"], loop[k, "eth_group"], "comb")
	men <- make.index.files.correct(loop[k, "trait"], loop[k, "eth_group"], "men")
	women <-  make.index.files.correct(loop[k, "trait"], loop[k, "eth_group"], "women")

        columns <- ifelse(loop[k, "sex_group"] == "comb", "\\.combined$",
                        ifelse(loop[k, "sex_group"] == "women", "\\.females$",
                        ifelse(loop[k, "sex_group"] == "men", "\\.males$", NA)))

	df <- merge(merge(comb, men, by = colnames(comb), all = T), women, by = colnames(women), all = T)
	df <- remove.secondary(df, "pval.combined")

	if (loop[k, "sex_group"] %in% c("women", "men")) {
		for (critical_p in c(0.01, 0.05, 0.1)) {
			df_phet <- df
		
			phet.columns <- gsub(".combined", ".phet", colnames(df[, grep("\\.combined$", names(df_phet))]))
			df_phet[, phet.columns] <- as.numeric(NA)
			df_phet$fdr_p <- p.adjust(df$psexdiff, method = "BH")		

			for (row in 1:nrow(df_phet)) {
				if (df_phet$fdr_p[row] < critical_p) {
					df_phet[row, phet.columns] <- df_phet[row, grep(columns, names(df_phet))]
				} else if (df_phet$fdr_p[row] >= critical_p) {
					df_phet[row, phet.columns] <- df_phet[row, grep("\\.combined$", names(df_phet))]
				}
			}
			df_phet <- df_phet[, grep("^SNP$|^Chr$|^Pos$|^EA$|^NEA$|^psexdiff$|\\.phet$", names(df_phet))]
			df_phet <- make.beta.pos(df_phet)
		
			nrow(df_phet)
			
			#Remove tri-allelic SNPs and SNPs that fail QC and SNPs with long-distance LD:
		        print(df_phet[(df_phet[, colnames(df_phet) %in% c("snp", "SNP")] %in% duplicated$V1), ])
        		df_phet <- df_phet[!(df_phet[, colnames(df_phet) %in% c("snp", "SNP")] %in% duplicated$V1), ]

			nrow(df_phet)

        		print(df_phet[!(df_phet[, colnames(df_phet) %in% c("snp", "SNP")] %in% qc$V1), ])
        		df_phet <- df_phet[df_phet[, colnames(df_phet) %in% c("snp", "SNP")] %in% qc$V1, ]

			nrow(df_phet)

			for (s in 1:nrow(ld)) {
		                if (nrow(df_phet[df_phet[, colnames(df_phet) %in% c("snp", "SNP")] %in% ld[s, c("SNP_A", "SNP_B")], ]) >= 2) {
                		        df_phet <- df_phet[!(df_phet[, colnames(df_phet) %in% c("snp", "SNP")] == ld[s, "remove"]), ]
                        		print(ld[s, "remove"])
                		}
        		}
			print(nrow(df_phet))

			output_file_phet <- paste("../sensitivity.risk.weights/", loop[k, "trait"], ".", loop[k, "eth_group"], 
				".", loop[k, "sex_group"], ".", critical_p, ".fdr.180513.txt", sep = "")
			write.table(df_phet, output_file_phet, quote = F, row.names = F, sep = " ")
		}
	}
}


