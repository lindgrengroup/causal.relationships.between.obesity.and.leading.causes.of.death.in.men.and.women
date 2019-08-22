#/bin/env Rscript
#$-cwd

out_df <- data.frame("V1" = character(), stringsAsFactors = F)

for (chr in 1:22) {
	geno_file <- paste("../genotype.qc.snplists.180507/geno.180507.chr", chr, ".snplist", sep = "")
	hwe_file <- paste("../genotype.qc.snplists.180507/hwe.180507.chr", chr, ".snplist", sep = "")
	maf_file <- paste("../genotype.qc.snplists.180507/maf.180507.chr", chr, ".snplist", sep = "")
	info_file <- paste("../genotype.qc.snplists.180507/info.180507.chr", chr, sep = "")
		
	geno_df <- read.table(geno_file, stringsAsFactors = F, header = F)
	hwe_df <- read.table(hwe_file, stringsAsFactors = F, header = F)
	maf_df <- read.table(maf_file, stringsAsFactors = F, header = F)
	info_df <- read.table(info_file, stringsAsFactors = F, header = F)

	df <- merge(merge(merge(geno_df, hwe_df, by = "V1"), maf_df, by = "V1"), info_df, by = "V1")
	
	out_df <- rbind(out_df, df)
}

out_df <- out_df[!duplicated(out_df$V1), ]
write.table(out_df, "../genotype.qc.snplists.180507/snps_that_pass_all_qc_180508.txt", 
	quote = F, row.names = F, col.names = F, sep = "\t")

