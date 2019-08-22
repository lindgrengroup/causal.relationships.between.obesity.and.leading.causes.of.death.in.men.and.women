#/bin/env Rscript
#$-cwd

#Read in the list of giukbb snps that fail QC
snps_df <- read.table("../temporary.bed.pulit.180911/snps.that.are.not.in.qc.list.giukbb.txt", 
	stringsAsFactors = F, header = F)

#Make a vector of the chromosomes that those SNPs are located on
chromosomes <- c(7, 9, 11, 18)

for (chr in chromosomes) {
	geno_file <- paste("/well/lindgren/jc/ukbb/genotype.qc.snplists.180507/geno.180507.chr", chr, ".snplist", sep = "")
        hwe_file <- paste("/well/lindgren/jc/ukbb/genotype.qc.snplists.180507/hwe.180507.chr", chr, ".snplist", sep = "")
        maf_file <- paste("/well/lindgren/jc/ukbb/genotype.qc.snplists.180507/maf.180507.chr", chr, ".snplist", sep = "")
        info_file <- paste("/well/lindgren/jc/ukbb/genotype.qc.snplists.180507/info.180507.chr", chr, sep = "")

	for (file in c(geno_file, hwe_file, maf_file, info_file)) {
	        df <- read.table(file, stringsAsFactors = F, header = F)
		print(file)
		print(df[df$V1 %in% snps_df$V1, ])
	}
}

#Check what SNP-groups change
snp_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig",
                "whradjbmi.eur.comb.pulit.sig",
                "bmi.eur.men.pulit", "bmi.eur.men.pulit.sig", "bmi.eur.men.pulit.phet",
                "whr.eur.men.pulit", "whr.eur.men.pulit.sig", "whr.eur.men.pulit.phet",
                "whradjbmi.eur.men.pulit", "whradjbmi.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.phet",
                "bmi.eur.women.pulit", "bmi.eur.women.pulit.sig", "bmi.eur.women.pulit.phet",
                "whr.eur.women.pulit", "whr.eur.women.pulit.sig", "whr.eur.women.pulit.phet",
                "whradjbmi.eur.women.pulit", "whradjbmi.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.phet",
                "bmi.eur.comb.pulit.winner", "whr.eur.comb.pulit.winner",
                "whradjbmi.eur.comb.pulit.winner",
                "bmi.eur.men.pulit.winner", "whr.eur.men.pulit.winner",
                "whradjbmi.eur.men.pulit.winner",
                "bmi.eur.women.pulit.winner", "whr.eur.women.pulit.winner",
                "whradjbmi.eur.women.pulit.winner")

for (snp_group in snp_groups) {
        new_file <- paste("../", snp_group, ".180513.txt", sep = "")
        old_file <- paste("../index.files.with.trialleles/", snp_group, ".180513.txt", sep = "")

        new <- read.table(new_file, stringsAsFactors = F, header = T)
        old <- read.table(old_file, stringsAsFactors = F, header = T)

        print(snp_group)
        new$included_new <- T
        old$included_old <- T
        difference <- merge(new, old, all = T)
        if (nrow(difference[is.na(difference$included_new) | is.na(difference$included_old), ]) != 0) {
                print(difference[is.na(difference$included_new) | is.na(difference$included_old), ])
        }
}

