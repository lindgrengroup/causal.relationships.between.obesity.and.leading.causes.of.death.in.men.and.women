#!/bin/env Rscript
#$-cwd

library(ggplot2)
library(gridExtra)
library(lattice)
library(grid)
library(grDevices)
library(ggpubr)

X11(height = 30, width = 40)

####################################################################
########################## Start with GRSs 
####################################################################

for (model in c("grs")) {
	df_raw <- read.table("../results.logistic.regressions.180514/log.results.table.180627.txt",
       	        header = T, stringsAsFactors = F)
	smoking_columns <- unique(df_raw$case_column[grep("_smoking", df_raw$case_column)])
	datasets <- c("pulit", "giukbb")
	units <- "raw_scoresum"

	for (eth_group in c("all.white", "brit.irish")) {
		for (dataset in datasets) {
			for (unit in units) {
			df <- df_raw[!(df_raw$case_column %in% c("t1d_cases_prob", "t2d_cases_prob", "smoker_cases", smoking_columns)), ]
	
			if (dataset == "pulit")  {
        		        comb_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig",  "whradjbmi.eur.comb.pulit.sig")
       		        	male_groups <- c("bmi.eur.men.pulit.sig", "whr.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.sig")
               			female_groups <- c("bmi.eur.women.pulit.sig", "whr.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.sig")
			} else if (dataset == "giukbb") {
				comb_groups <- c("bmi.eur.comb.giukbb.sig", "whr.eur.comb.giukbb.sig",  "whradjbmi.eur.comb.giukbb.sig")
               		        male_groups <- c("bmi.eur.men.giukbb.sig", "whr.eur.men.giukbb.sig", "whradjbmi.eur.men.giukbb.sig")
               		        female_groups <- c("bmi.eur.women.giukbb.sig", "whr.eur.women.giukbb.sig", "whradjbmi.eur.women.giukbb.sig")
			}

			#Subset to only keep the actual analyses
			df <- df[(df$grs_unit == unit & df$eth_group == eth_group) &
        				((df$sex_group == "comb" & df$snp_group %in% comb_groups) |
        				(df$sex_group == "men" & df$snp_group %in% male_groups) |
        				(df$sex_group == "women" & df$snp_group %in% female_groups)), ]
		
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
				haem_stroke_cases = "Stroke - Hemorrhagic",
				isch_stroke_cases = "Stroke - Ischemic",
        			t2d_cases_probposs = "Type 2 Diabetes",
        			t1d_cases_probposs = "Type 1 Diabetes",
        			any_infertility_cases = "Infertility",
        			nafld_cases = "NAFLD",
        			cld_cases = "CLD")

			df$trait <- df$case_column
			df$grs_trait_name <- gsub("\\.(.)+", "", df$snp_group)
			df$sex_group <- ifelse(df$sex_group == "comb", "Combined", 
				ifelse(df$sex_group == "men", "Men", 
				ifelse(df$sex_group == "women", "Women", "missing")))

			for (i in 1:length(dict_traits)) {
				df$trait <- as.character(replace(df$trait,

                	        df$trait == names(dict_traits[i]), dict_traits[i]))
			}

			df$ci <- paste(formatC(df$grs_or, digits = 2, format = "f"), " (",
        		        formatC(df$grs_lci_or, digits = 2, format = "f"), ",",
        	        	formatC(df$grs_uci_or, digits = 2, format = "f"), ")", sep = "")

			df$order <- ifelse(df$trait == "Type 2 Diabetes", "17",
				ifelse(df$trait == "CAD", "16", 
				ifelse(df$trait == "Breast Cancer", "15", 
				ifelse(df$trait == "CLD", "14", 
				ifelse(df$trait == "Colorectal Cancer", "13", 
				ifelse(df$trait == "COPD", "12",
				ifelse(df$trait == "Dementia", "11", 
				ifelse(df$trait == "Infertility", "10", 
				ifelse(df$trait == "Lung Cancer", "09",
				ifelse(df$trait == "NAFLD", "08", 
				ifelse(df$trait == "Renal Failure", "07",
				ifelse(df$trait == "Renal Failure - Acute", "06",
				ifelse(df$trait == "Renal Failure - Chronic", "05", 
				ifelse(df$trait == "Stroke", "04",
				ifelse(df$trait == "Stroke - Hemorrhagic", "03", 
				ifelse(df$trait == "Stroke - Ischemic", "02",  
				ifelse(df$trait == "Type 1 Diabetes", "01", 
				"missing")))))))))))))))))

			df$grs_trait_name <- ifelse(df$grs_trait_name == "bmi", "BMI", 
				ifelse(df$grs_trait_name == "whr", "WHR", 
				ifelse(df$grs_trait_name == "whradjbmi", "WHRadjBMI", NA)))

			df$unique_combinations <- paste(df$order, "_", df$case_column, "_", 
				ifelse(df$sex_group == "Combined", "03", 
				ifelse(df$sex_group == "Men", "02", 
				ifelse(df$sex_group == "Women", "01", "missing"))), 
				df$sex_group, sep = "")

			df$lci_arrow <- ifelse(df$grs_lci_or < 0.5, 0.5, NA)
			df$uci_arrow <- ifelse(df$grs_uci_or > 5, 4, NA)
			title_name <- "Odds ratio (95% CI) per 1 unit higher weighted GRS"
			breaks_number <- c(0.5, 1, 5) 
			ylim_number <- c(0.5, 5)
			df$yend_arrow <- df$uci_arrow+1
			critical_p <- 0.05/(length(unique(df[, "case_column"]))*length(unique(gsub("\\..*", "", df$snp_group))))

			df$sig_or <- ifelse(df$grs_p < critical_p, df$grs_or, NA)
			df$not_sig_or <- ifelse(df$grs_p >= critical_p, df$grs_or, NA)

			critical_het_p <- 0.05/(nrow(df[!is.na(df$cochrans_p), ])/2)
			df$sig_heterogeneity <- ifelse(df$cochrans_p < critical_het_p, df$grs_uci_or, NA)

			print_df <- df[, c("trait", "grs_trait_name", "sex_group", "ci", 
				"grs_p", "cochrans_p", "cochrans_i2", "unique_combinations")]

			print_df$grs_p <- ifelse(print_df$grs_p < 0.001, as.character(formatC(print_df$grs_p,
	                        1, format = "e")), ifelse(print_df$grs_p < 0.01, 
				as.character(formatC(print_df$grs_p, format = "f", 3)),
				as.character(formatC(print_df$grs_p, format = "f", 2))))
			print_df$cochrans_p <- ifelse(print_df$cochrans_p < 0.001, as.character(formatC(print_df$cochrans_p, 
				1, format = "e")), ifelse(print_df$cochrans_p < 0.01, 
				as.character(formatC(print_df$cochrans_p, format = "f", 3)), 
				as.character(formatC(print_df$cochrans_p, format = "f", 2))))
	
			print_df_bmi <- subset(print_df, print_df$grs_trait_name == "BMI")
			colnames(print_df_bmi) <- paste("bmi_", colnames(print_df_bmi), sep = "")
			print_df_whr <- subset(print_df, print_df$grs_trait_name == "WHR")
			colnames(print_df_whr) <- paste("whr_", colnames(print_df_whr), sep ="")
			print_df_whradjbmi <- subset(print_df, print_df$grs_trait_name == "WHRadjBMI")
			colnames(print_df_whradjbmi) <- paste("whradjbmi_", colnames(print_df_whradjbmi), sep ="")
	
			print_df <- merge(print_df_bmi, print_df_whr, by.x = c("bmi_trait", "bmi_sex_group"), 
				by.y = c("whr_trait", "whr_sex_group"))

			print_df <- merge(print_df, print_df_whradjbmi, by.x =  c("bmi_trait", "bmi_sex_group"),
				by.y = c("whradjbmi_trait", "whradjbmi_sex_group"))

			print_df <- print_df[order(print_df$bmi_unique_combinations, decreasing = T), 
				c("bmi_trait", "bmi_sex_group", "bmi_ci", "bmi_grs_p", "bmi_cochrans_p", 
				"whr_ci", "whr_grs_p", "whr_cochrans_p", 
				"whradjbmi_ci", "whradjbmi_grs_p", "whradjbmi_cochrans_p")]
			colnames(print_df) <- c("bold(Outcome)", "bold(Sex-strata)", 
				"BMI OR", "BMI P", "BMI Pheterogeneity", 
				"WHR OR", "WHR P", "WHR Pheterogeneity", 
				"WHRadjBMI OR", "WHRadjBMI P", "WHRadjBMI Pheterogeneity")

			print_df[, c("Outcome", "Sex-strata")] <- print_df[, c("bold(Outcome)", "bold(Sex-strata)")]
			print_df[nrow(print_df)+1, ] <- ""
			print_df[duplicated(print_df[["bold(Outcome)"]]), "bold(Outcome)"] <- ""

			#Make the figure and table
			grob_df_first <- tableGrob(print_df[, 1:2], rows = NULL, 
				theme = ttheme_minimal(base_size = 9, parse = T, base_family = "Times", 
				padding = unit(c(3.5,1.71), "mm"), colhead=list(fg_params=list(hjust=0, x=0.1)),
				core=list(fg_params = list(hjust=0, x=0.1),
				bg_params = list(fill = c(rep("white", 3), 
				rep("gray87", 3), "white", rep("gray87", 3), rep("white", 3),
				rep("gray87", 3), rep("white", 3), rep("gray87", 3), rep("white", 3), 
				rep("gray87", 3), rep("white", 3), rep("gray87", 3), rep("white", 3), 
				rep("gray87", 3), rep("white", 3), rep("gray87", 3), rep("white", 3))))))
	
			plot <- ggplot(df, aes(x=unique_combinations, y=grs_or, ymin=grs_lci_or, ymax=grs_uci_or)) +
			geom_point(aes(color = sex_group), color = "white", shape = 18, size = 3, show.legend =F) +
			geom_vline(xintercept = c(45, 41, 35, 29, 23, 17, 11, 5),
			        colour = "grey", alpha = 0.5, size = 16) +
			geom_hline(yintercept = 1, linetype = "dashed") +
			geom_errorbar(size = 0.3, width=.3) +
			geom_point(aes(y = df$not_sig_or, color = sex_group), shape = 23, fill = "white", 
				size = 1.65, show.legend =F) +
			geom_point(aes(y = df$sig_or, color = sex_group), shape = 18, size = 2.3, show.legend =F) +
			geom_point(aes(y = df$sig_heterogeneity), shape = 8, size = 0.7, 
				position = position_nudge(x = 0, y = 0.12)) +
			geom_segment(aes(yend = lci_arrow, y = lci_arrow,  xend = unique_combinations), size = 0.1,
			        arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
			geom_segment(aes(yend = yend_arrow, y = uci_arrow, xend = unique_combinations), 
				size = 0, arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
			theme_classic() +
			scale_color_brewer(palette = "Dark2") +
			scale_y_continuous(trans = "log",
			        name = title_name,
			        breaks = breaks_number) +
			theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8, family = "Times"), 
				axis.text.y=element_blank(), axis.text = element_text(color = "black", family = "Times"), 
				strip.background = element_blank(), strip.text = element_text(face = "bold", family = "Times", 
				vjust = -0.82)) + 
			coord_flip(ylim = ylim_number, xlim = c(0, 49.6), expand = F) +
			theme(plot.margin=unit(c(0.61, 0.5, 0, 0), "cm")) + 
			facet_grid(. ~grs_trait_name, drop = T) +
			theme(panel.spacing = unit(1, "lines"))

			g <- arrangeGrob(grob_df_first, plot, nrow = 1, widths= c(2.7,3.8))
			output_file <- paste("/Net/fs1/home/linc4222/", model, ".results.", dataset, 
				".", eth_group, ".", unit, ".180819.jpeg", sep = "")
			ggsave(output_file, g, unit = "cm", height = 21.5, 
				width = 14, device = "jpeg")

			print_df[, 1:2] <- print_df[, c("Outcome", "Sex-strata")]
			
			output_file <- paste("/Net/fs1/home/linc4222/", model, ".binary.outcomes.table.", dataset, ".",
                                        eth_group, ".", unit, ".180909.txt", sep = "")
	
			#Same number of cases for BMI, WHR, and WHRadjBMI since GRSs
			n_cases <- df[df$grs_trait_name == "BMI", c("sex_group", "trait", "n_cases")]
			colnames(n_cases) <- c("sex_group", "trait", "N cases")
			merging_df <- print_df[, c(12, 13, 3:11)]
			merging_df$order <- 1:nrow(merging_df)
			n_cases <- merge(merging_df, n_cases, by.x = c("Outcome", "Sex-strata"), 
				by.y = c("trait", "sex_group"), all.x = T)								
			n_cases <- n_cases[order(n_cases$order), ]
			n_cases <- n_cases[, c("Outcome", "Sex-strata", "N cases", "BMI OR", "BMI P", 
				"BMI Pheterogeneity", "WHR OR", "WHR P", "WHR Pheterogeneity", "WHRadjBMI OR", 
				"WHRadjBMI P", "WHRadjBMI Pheterogeneity")]
			write.table(n_cases, output_file,
                                       quote = F, row.names = F, sep = "\t", na = "-")
			}
		}
	}
}

###################################################################################
####################### PLOT COMPARING THE DIFFERENT APPROACHES ####################
####################################################################################

df_raw <- read.table("../results.linear.regressions.180514/anthro.results.table.180521.txt",
                stringsAsFactors = F, header = T)
pulit_snp_groups <- unique(df_raw$snp_group[grep("pulit|fdr", df_raw$snp_group)])

eth_groups <- c("all.white", "brit.irish")

df_dataset <- df_raw[df_raw$snp_group %in% pulit_snp_groups, ]

for (eth_group in eth_groups) {
	 df <- df_dataset
         df <- df[(df$grs_unit == "raw_scoresum" & df$eth_group == eth_group & df$extra_adjustment == "-" &
			df$trait_unit == "sd") &
                        ((df$snp_group %in% c("bmi.eur.comb.pulit.sig") &
                        df$sex_group %in% c("men", "women") &
                        df$trait == "bmi") |
                        (df$snp_group %in% c("bmi.eur.men.pulit.phet", "bmi.eur.men.pulit.sig",
                                        "bmi.eur.men.0.01.fdr", "bmi.eur.men.0.05.fdr", "bmi.eur.men.0.1.fdr") &
                        df$sex_group == "men" &
                        df$trait == "bmi") |
                        (df$snp_group %in% c("bmi.eur.women.pulit.phet", "bmi.eur.women.pulit.sig",
                                        "bmi.eur.women.0.01.fdr", "bmi.eur.women.0.05.fdr", "bmi.eur.women.0.1.fdr") &
                        df$sex_group == "women" &
                        df$trait == "bmi") |
                        (df$snp_group %in% c("whr.eur.comb.pulit.sig") &
                        df$sex_group %in% c("men", "women") &
                        df$trait == "whr") |
                        (df$snp_group %in% c("whr.eur.men.pulit.phet", "whr.eur.men.pulit.sig",
                                        "whr.eur.men.0.01.fdr", "whr.eur.men.0.05.fdr", "whr.eur.men.0.1.fdr") &
                        df$sex_group == "men" &
                        df$trait == "whr") |
                        (df$snp_group %in% c("whr.eur.women.pulit.phet", "whr.eur.women.pulit.sig",
                                        "whr.eur.women.0.01.fdr", "whr.eur.women.0.05.fdr", "whr.eur.women.0.1.fdr", "whr.eur.women.0.1.fdr") &
                        df$sex_group == "women" &
                        df$trait == "whr") |
                        (df$snp_group %in% c("whradjbmi.eur.comb.pulit.sig") &
                        df$sex_group %in% c("men", "women") &
                        df$trait == "res_whr_inv") |
                        (df$snp_group %in% c("whradjbmi.eur.men.pulit.phet", "whradjbmi.eur.men.pulit.sig",
                                "whradjbmi.eur.men.0.01.fdr", "whradjbmi.eur.men.0.05.fdr",
                                "whradjbmi.eur.men.0.1.fdr") &
                        df$sex_group == "men" &
                        df$trait == "res_whr_inv") |
                        (df$snp_group %in% c("whradjbmi.eur.women.pulit.phet", "whradjbmi.eur.women.pulit.sig",
                                "whradjbmi.eur.women.0.01.fdr", "whradjbmi.eur.women.0.05.fdr",
                                "whradjbmi.eur.women.0.1.fdr") &
                        df$sex_group == "women" &
                        df$trait == "res_whr_inv")), ]

		df$sex_group <- ifelse(df$sex_group == "comb", "Combined",
                        ifelse(df$sex_group == "men", "Men",
                        ifelse(df$sex_group == "women", "Women", "missing")))
		
		df$instrument <- gsub("bmi\\.eur\\.|whr\\.eur\\.|whradjbmi\\.eur\\.", "", df$snp_group)
		df$instrument <- gsub("men\\.|women\\.", "", df$instrument)
		df$instrument <- ifelse(df$instrument == "comb.pulit.sig", "Combined weights", 
					ifelse(df$instrument == "pulit", "Sex-specific index SNPs only", 
					ifelse(df$instrument == "pulit.phet", "P-heterogeneity, Bonferroni", 
					ifelse(df$instrument == "0.01.fdr", "P-heterogeneity, FDR 1%", 
					ifelse(df$instrument == "0.05.fdr", "P-heterogeneity, FDR 5%", 
					ifelse(df$instrument == "0.1.fdr", "P-heterogeneity, FDR 10%", 
					ifelse(df$instrument == "pulit.sig", "Sex-specific estimates", "MISSING")))))))
		
		subset_df <- df[, c("instrument", "sex_group", "trait", "grs_r2", "grs_beta", "grs_lci", "grs_uci")]
		
		new_df <- subset_df
		new_df$instrument <- factor(new_df$instrument, levels = c("Combined weights", 
			"P-heterogeneity, Bonferroni", 
			"P-heterogeneity, FDR 1%", "P-heterogeneity, FDR 5%", "P-heterogeneity, FDR 10%", 
			"Sex-specific estimates", "Sex-specific index SNPs only"))

			
		new_df <- new_df[order(new_df$instrument), ]
		new_df$unique_combinations <- paste(gsub(" |,|\\.|-|%", "", new_df$instrument), sep = "_")
		new_df$trait <- ifelse(new_df$trait == "bmi", "BMI", 
					ifelse(new_df$trait == "whr", "WHR", 
					ifelse(new_df$trait == "res_whr_inv", "WHRadjBMI", "MISSING")))

		new_df_women <- new_df[new_df$sex_group == "Women", ]
		new_df_men <- new_df[new_df$sex_group == "Men", ]
		
		new_df <- merge(new_df_women, new_df_men, by = c("instrument", "trait"))		
		new_df$grs_r2.x <- new_df$grs_r2.x *100
		new_df$grs_r2.y <- new_df$grs_r2.y *100

                plot1 <- ggplot(new_df, aes(x=instrument, y=grs_beta.x, ymin=grs_lci.x, ymax=grs_uci.x, group = 1)) +
                geom_errorbar(size = 0.3, width=.1, color = "grey28", position = position_nudge(x=-0.1)) +
                geom_errorbar(aes(ymin = grs_lci.y, ymax = grs_uci.y), size = 0.3, width = .1, color = "grey28",
                        position = position_nudge(x = 0.1)) +
                geom_point(shape = 18, color = "#d95f02",
                        size = 1.65, show.legend =F, position = position_nudge(x=-0.1)) +
                geom_point(aes(y = grs_beta.y), color = "#7570b3", shape = 18,
                        size = 1.65, show.legend = F, position = position_nudge(x = 0.1)) +
                scale_color_brewer(palette = "Dark2") +
                theme(plot.margin=unit(c(1, 2, 0, 1), "cm")) +
                theme(axis.text.x=element_blank()) +
                scale_x_discrete(name = "SNP selection and weighting approach") +
		scale_y_continuous(name = "Estimate (95% CI) in SD-units\nfor respective obesity trait",
                        breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1.0"), limits=c(0, 1.35)) +
                theme(axis.title.x = element_blank(), axis.title.y = element_text(color = "black", family = "Times"), 
                        axis.text = element_text(color = "black", size = 7, family = "Times"), axis.ticks.x = element_blank(), 
                        strip.background = element_blank(), strip.text = element_text(face = "bold", family = "Times", 
                        vjust = -0.5), axis.text.y = element_text(color = "black", family = "Times")) +
                scale_fill_manual(name="Sex",
                        values=c(Women="#d95f02", Men="#7570b3")) +
                facet_grid(. ~trait, drop = T)

		plot2 <- ggplot(new_df, aes(x=instrument, y = grs_beta.x, group = 1)) +
		geom_col(aes(y = grs_r2.x, fill = "Women"), width = 0.2, position = position_nudge(x=-0.1)) +
                geom_col(aes(y = grs_r2.y, fill = "Men"), width = 0.2, position = position_nudge(x=0.1)) +
		theme_gray() +
                scale_color_brewer(palette = "Dark2") +
                theme(plot.margin=unit(c(1, 2, 1, 0.5), "cm")) +
                theme(axis.text.x=element_text(angle=(-30), hjust = 0, vjust=1, family = "Times")) +
                scale_x_discrete(name = "SNP selection and weighting approach") +
                scale_y_continuous(name = "% Trait variance\nexplained", breaks = c(0, 3, 6), labels = c("0", "3.0", "6.0")) +
                theme(axis.title.x = element_text(size = 10, family = "Times"), plot.margin=unit(c(0.1, 2, 1, 1), "cm"),
			axis.title.y = element_text(family = "Times"),
                        axis.text = element_text(color = "black", size = 7, family = "Times", hjust = 0.5),  
			panel.grid.minor.y = element_blank(), axis.text.y = element_text(hjust=1, vjust = -0.5, family = "Times"),
                        strip.background = element_blank(), strip.text = element_blank(), legend.position = "bottom", 
			legend.text = element_text(family = "Times"), legend.title = element_text(family = "Times")) +
                scale_fill_manual(name="Sex",
                        values=c(Women="#d95f02", Men="#7570b3")) +
                facet_grid(. ~trait, drop = T)
		
		g <- arrangeGrob(plot1, plot2, ncol = 1, heights = c(1,1))
                output_file <- paste("/Net/fs1/home/linc4222/comparison.weight.strategies.separate.facets.", eth_group, ".181019.jpeg", sep = "")
                ggsave(output_file, g, unit = "cm", height = 18,
                        width = 15, device = "jpeg")

		output_file <- paste("/Net/fs1/home/linc4222/comparison.weight.strategies.separate.facets.", eth_group, ".181019.pdf", sep = "")
                ggsave(output_file, g, unit = "cm", height = 18,
                        width = 15, device = "pdf")

}

######################################################################
###### Make a sort of heatmap and million deaths ##################
######################################################################

df_raw <- read.table("../results.mr.180730/ipd.mr.binary.results.180815.txt", 
	stringsAsFactors = F, header = T)

snp_groups <- unique(df_raw$snp_group[grep("pulit\\.sig", df_raw$snp_group)])
traits <- unique(df_raw$trait[!(grepl("_smoking|nafld_cases|cld_cases|smoker|infertility|prob$", df_raw$trait))])

df_raw <- df_raw[df_raw$eth_group == "all.white" & df_raw$snp_group %in% snp_groups &
		df_raw$exposure_unit == "sd" & df_raw$outcome_unit == "clin" &
		df_raw$trait %in% traits & df_raw$sex_group %in% c("men", "women"), ]

df_raw$log_p <- ifelse(df_raw$grs_p >= (0.05/51), NA, -log10(df_raw$grs_p))

df_raw <- df_raw[, c("grs_trait", "sex_group", "trait", "grs_or", "log_p")]

df <- data.frame(trait = c(unique(df_raw$trait), "diabetes", "breast_cancer_cases", 
	"colorectal_cancer_cases", "dementia_cases", "haem_stroke_cases"), stringsAsFactors = F)

for (trait in unique(df_raw$grs_trait)) {
	for (sex_group in unique(df_raw$sex_group)) {
		df_subset <- df_raw[df_raw$grs_trait == trait & df_raw$sex_group == sex_group, ]
		df_subset[, c("grs_trait", "sex_group")] <- NULL
		colnames(df_subset)[colnames(df_subset) %in% c("grs_or", "log_p")] <- paste(trait, sex_group, 
			colnames(df_subset)[colnames(df_subset) %in% c("grs_or", "log_p")], sep = "_")
	
		df <- merge(df, df_subset, by = "trait", all = T)
	}
}

#Number of 1000 deaths per sex, men first, then women. 
#Data taken from "2016 Global" sheet of "Global summary estimates" 
#downloadable here: http://www.who.int/healthinfo/global_burden_disease/estimates/en/
#on 2018-11-27. In the dict, first the "nice-looking" disease name, then men deaths, 
#then women deaths (per 1000), then order in the table
dict_traits <- list(cad_cases = list("Coronary artery disease", "4,955", "4,478", 10.8),
             copd_cases = list("Chronic obstructive\npulmonary disease", "1,668", "1,373", 7.2),
             lungcancer_cases = list("Lung cancer", "1,177", 531, 5.4),
             renal_failure_cases = list("Renal failure - chronic\nand acute", 623, 557, 1.8),
             aki_cases = list("Acute:", 6, 6, 0.6),
             ckd_cases = list("Chronic:", 617, 551, 1.2),
	     stroke_cases = list("Stroke - hemorrhagic\nand ischemic", "2,893", "2,887", 9),
             isch_stroke_cases = list("Ischaemic:", "1,338", "1,473", 8.4),
             diabetes = list("Diabetes - type 2 and\ntype 1 diabetes", "737", "862", 3.6),
             t2d_cases_probposs = list("Type 2:", NA, NA, 3.6),
             t1d_cases_probposs = list("Type 1:", NA, NA, 3))

for (i in 1:length(dict_traits)) {
	df$deaths_men <- as.character(replace(df$deaths_men,
                 df$trait == names(dict_traits[i]), dict_traits[[i]][[2]]))
        df$deaths_women <- as.character(replace(df$deaths_women,
                 df$trait == names(dict_traits[i]), dict_traits[[i]][[3]]))
	df$order <- replace(df$order,
                 df$trait == names(dict_traits[i]), dict_traits[[i]][[4]])
	df$trait <- as.character(replace(df$trait,
                 df$trait == names(dict_traits[i]), dict_traits[[i]][[1]]))

}

df$dm_deaths_men[df$trait %in% c("Type 1 diabetes", "Type 2 diabetes")] <- 737
df$dm_deaths_women[df$trait %in% c("Type 1 diabetes", "Type 2 diabetes")] <- 862

subtypes <- c("Acute:", "Chronic:", "Ischaemic:", "Type 2:", "Type 1:")
diseases_investigated <- c("Coronary artery disease:", "Stroke - hemorrhagic\nand ischemic:", 
	"Ischemic stroke:", "Chronic obstructive\npulmonary disease:", "Lung cancer:", "Type 2 diabetes:", 
	"Type 1 diabetes:", "Renal failure - chronic\nand acute: ", "Chronic renal failure:", 
	"Acute renal failure:")
df[df$trait %in% subtypes, c("deaths_men", "deaths_women")] <- NA
df <- df[order(df$order), ]

plot <- ggplot(df, aes(x=order)) +
geom_vline(xintercept = c(1.8, 3, 3.6, 5.4, 7.2, 9, 10.8), colour = "gray87", size = 10) +
geom_point(aes(y=11.5, size=bmi_men_grs_or, fill = bmi_men_log_p), shape = 21) +
geom_point(aes(y=12.5, size=bmi_women_grs_or, fill = bmi_women_log_p), shape = 21) +
geom_point(aes(y=14, size=whr_men_grs_or, fill = whr_men_log_p), shape = 21) +
geom_point(aes(y=15, size=whr_women_grs_or, fill = whr_women_log_p), shape = 21) +
geom_point(aes(y=16.5, size=res_whr_inv_men_grs_or, fill = res_whr_inv_men_log_p), shape = 21) +
geom_point(aes(y=17.5, size=res_whr_inv_women_grs_or, fill = res_whr_inv_women_log_p), shape = 21) +
geom_segment(aes(y =0, yend = (as.integer(gsub(",", "", df$deaths_men))/1000), x = order+0.1, xend = order+0.1), 
	colour = '#d95f02', size = 2.5) +
geom_segment(aes(y =0, yend = (as.integer(gsub(",", "", df$deaths_women))/1000), x = order-0.1, xend = order-0.1), 
	colour = '#7570b3', size = 2.5) +
geom_segment(aes(y=0, yend = 5.2, x = 0.02, xend = 0.02)) +
geom_segment(aes(y=0, yend = 0, x = 0, xend = 11.5)) +
geom_text(aes(y=6.3, x=order, label=deaths_men), size = 3, family = "Times", hjust = 1) +
geom_text(aes(y=7.3, x=order, label=deaths_women), size = 3, family = "Times", hjust = 1) +
geom_segment(aes(y =3.5, yend = 4, x = 0.6, xend = 0.6),
        colour = '#d95f02', size = 2.5) +
geom_segment(aes(y =3.5, yend = 4, x = 0.4, xend = 0.4),
        colour = '#7570b3', size = 2.5) +
geom_text(aes(y=3.5, x= 0.87, label="bold('Sex')"), parse = T, size = 3, family = "Times", hjust = 0) +
geom_text(aes(y=4.1, x= 0.6, label="Men"), size = 2.8, family = "Times", hjust = 0) + # was 2.5
geom_text(aes(y=4.1, x= 0.4, label="Women"), size = 2.8, family = "Times", hjust = 0) +	#was 2.5
annotate("text", y=c(0, 7.9), x = rep(12.2, 2), label = c('bold("A. Number of deaths per disease globally/year,")', 
	'bold("B. Effect of obesity traits on leading mortality causes")'), hjust = 0, size = 3, family = "Times", parse = T) +
annotate("text", y=c(3.72, 12, 14.5, 17), x = rep(11.9, 4), label = c('bold("in 1,000 deaths as estimated by the WHO for 2016")',
	 'bold("BMI")', 'bold("WHR")', 'bold("WHRadjBMI")'), size = 3, family = "Times", parse = T) +
annotate("text", y=c(5.95, 6.85, 9.2, 11.5, 12.5, 14, 15, 16.5, 17.5), x = rep(11.5, 9), 
	label = c("Men", "Women", "Investigated disease", rep(c("Men", "Women"), 3)), size = 2.8, family = "Times") + #was 2.5
annotate("text", y=7.95, x = c(10.8, 9, 8.4, 7.2, 5.4, 3.6, 3, 1.8, 1.2, 0.6), label = diseases_investigated, size = 2.8, family = "Times", #was.25
	hjust = 0) +
scale_fill_gradient(na.value = "white", low = "yellow1", high="red2", name = "-log10 P") +
scale_x_continuous(name = "WHO\nDisease", breaks = df$order, labels=ifelse(df$trait %in% subtypes, "", df$trait)) +
scale_size_continuous(range = c(min(df[, grep("_grs_or", colnames(df))], na.rm = T)*2.5, max(df[, grep("_grs_or", colnames(df))], na.rm = T)*2.5), 
 	breaks = c(0.5, 1, 2, 3, 4), trans = "log", name = "Odds ratio") +
theme_classic() +
scale_y_continuous(name = "1,000 deaths per year", breaks = c(0, 2.5, 5), labels=c("0", "2,500", "5,000")) +
coord_flip(ylim = c(0,18), xlim = c(0, 15), expand = F) +
labs(title = "WHO\ndisease") +
theme(axis.line.y = element_blank(), axis.line.x = element_line(colour = "white"), axis.text = element_text(color = "black", family = "Times"), 
	axis.ticks.y = element_blank(), axis.ticks.x = element_line(colour = "black"), 
	axis.text.y = element_text(face = "bold", family = "Times"), axis.title.y = element_blank(), 
	legend.title=element_text(size=8, family = "Times", face = "bold"), legend.text=element_text(size=8, family = "Times"),
	axis.title.x = element_text(colour = "black", family = "Times", size = 9, hjust = 0.1)) +
theme(plot.title = element_text(face="bold", size = 9, family = "Times", hjust = -0.075, margin=margin(b=-105.5, t = 18))) +
guides(size = guide_legend(order=1))

output_file <- paste("/Net/fs1/home/linc4222/heatmap.of.mr.results.and.million.deaths.181127.jpeg", sep = "")
ggsave(output_file, plot, unit = "cm", height = 16, width = 22.3, device = "jpeg")

output_file <- paste("/Net/fs1/home/linc4222/heatmap.of.mr.results.and.million.deaths.181127.pdf", sep = "")
ggsave(output_file, plot, unit = "cm", height = 16, width = 22.3, device = "pdf")


##############################################################################
############## MAKE PLOTS FOR THESIS AND ARTICLE ###################
#################   MR PLOTS        ######################################
##############################################################################

df_raw <- read.table("../results.mr.180730/ipd.mr.binary.results.180815.txt",
        stringsAsFactors = F, header =T, sep = "\t")

for (dataset in c("pulit", "giukbb", "fdr.0.01", "fdr.0.05", "fdr.0.1", "phet", "index", "unweighted")) {
        pulit_snp_groups <- unique(df_raw$snp_group[grep("pulit\\.sig", df_raw$snp_group)])
        comb_snp_groups <- unique(df_raw$snp_group[grep("comb\\.pulit\\.sig", df_raw$snp_group)])
        giukbb_snp_groups <- unique(df_raw$snp_group[grep("giukbb", df_raw$snp_group)])
        fdr.0.01_snp_groups <- unique(df_raw$snp_group[grep("\\.0\\.01\\.fdr", df_raw$snp_group)])
        fdr.0.05_snp_groups <- unique(df_raw$snp_group[grep("\\.0\\.05\\.fdr", df_raw$snp_group)])
        fdr.0.1_snp_groups <- unique(df_raw$snp_group[grep("\\.0\\.1\\.fdr", df_raw$snp_group)])
        phet_snp_groups <- unique(df_raw$snp_group[grep("pulit.phet", df_raw$snp_group)])
        index_snp_groups <- unique(df_raw$snp_group[grep("pulit$", df_raw$snp_group)])
        unweighted_snp_groups <- unique(df_raw$snp_group[grep("unweighted", df_raw$snp_group)])

        if (dataset == "pulit") {
                df_dataset <- df_raw[df_raw$snp_group %in% pulit_snp_groups, ]
        } else if (dataset == "giukbb") {
                df_dataset <- df_raw[df_raw$snp_group %in% giukbb_snp_groups, ]
        } else if (dataset == "fdr.0.01") {
                df_dataset <- df_raw[df_raw$snp_group %in% c(fdr.0.01_snp_groups, comb_snp_groups), ]
        } else if (dataset == "fdr.0.05") {
                df_dataset <- df_raw[df_raw$snp_group %in% c(fdr.0.05_snp_groups, comb_snp_groups), ]
        } else if (dataset == "fdr.0.1") {
                df_dataset <- df_raw[df_raw$snp_group %in% c(fdr.0.1_snp_groups, comb_snp_groups), ]
        } else if (dataset == "phet") {
                df_dataset <- df_raw[df_raw$snp_group %in% c(phet_snp_groups, comb_snp_groups), ]
        } else if (dataset == "index") {
                df_dataset <- df_raw[df_raw$snp_group %in% c(index_snp_groups, comb_snp_groups), ]
        } else if (dataset == "unweighted") {
                df_dataset <- df_raw[df_raw$snp_group %in% c(unweighted_snp_groups), ]
        }

        for (eth_group in c("all.white", "brit.irish")) {
                for (unit in c("sd")) {
			for (type in c("article", "thesis")) {

                	#Subset to relevant analyses
                	df <- df_dataset
                	df <- df[df$exposure_unit == unit &
                	        df$eth_group == eth_group & !(df$trait %in% c("t1d_cases_prob", "t2d_cases_prob", "smoker_cases")) &
                	        df$function_name == "wald" & df$extra_adjustment == "-", ]

			#The first number is for the article order, the second for the thesis order
        	        dict_traits <- list(breast_cancer_cases = c("Breast cancer", "15", "17"),
                                cad_cases = c("CAD", "16", "16"),
                                colorectal_cancer_cases = c("Colorectal cancer", "13", "14"),
                                copd_cases = c("COPD", "12", "13"),
                                dementia_cases = c("Dementia", "11", "12"),
                                lungcancer_cases = c("Lung cancer", "09", "10"),
                                renal_failure_cases = c("Renal failure", "07", "08"),
                                aki_cases = c("Renal failure - acute", "06", "07"), 
                                ckd_cases = c("Renal failure - chronic", "05", "06"),
                                stroke_cases = c("Stroke", "04", "05"),
                                haem_stroke_cases = c("Stroke - hemorrhagic", "03", "04", "Stroke - haemorrhagic"),
                                isch_stroke_cases = c("Stroke - ischemic", "02", "03", "Stroke - ischaemic"),
                                t2d_cases_probposs = c("Type 2 diabetes", "17", "01"),
                                t1d_cases_probposs = c("Type 1 diabetes", "01", "02"),
                                any_infertility_cases = c("Infertility", "10", "11"),
                                nafld_cases = c("NAFLD", "08", "09"),
                                cld_cases = c("CLD", "14", "15"))

	                for (i in 1:length(dict_traits)) {
				df[df$trait == names(dict_traits[i]), "order_article"] <- dict_traits[names(dict_traits[i])][[1]][2]
				df[df$trait == names(dict_traits[i]), "order_thesis"] <- dict_traits[names(dict_traits[i])][[1]][2]
				df[df$trait == names(dict_traits[i]), "trait_name"] <- dict_traits[names(dict_traits[i])][[1]][1]
			}
		
	                df$grs_trait_name <- gsub("\\.(.)+", "", df$snp_group)
        	        df$sex_group <- ifelse(df$sex_group == "comb", "Combined",
                                ifelse(df$sex_group == "men", "Men",
                                ifelse(df$sex_group == "women", "Women", "missing")))

                	df$ci <- paste(formatC(df$grs_or, digits = 2, format = "f")," (",
                        	formatC(df$grs_lci_or, digits = 2, format = "f"), ",",
                        	formatC(df$grs_uci_or, digits = 2, format = "f"), ")", sep = "")

	                df$grs_trait_name <- ifelse(df$grs_trait_name == "bmi", "BMI",
        	                ifelse(df$grs_trait_name == "whr", "WHR",
                	        ifelse(df$grs_trait_name == "whradjbmi", "WHRadjBMI", NA)))

	                df$unique_combinations <- paste(df[, paste("order_", type, sep = "")], "_", df$trait_name, "_", 
				ifelse(df$grs_trait_name == "BMI", "03",
                                ifelse(df$grs_trait_name == "WHR", "02",
                                ifelse(df$grs_trait_name == "WHRadjBMI", "01", "missing"))), ifelse(df$sex_group == "Combined", "03",
                                ifelse(df$sex_group == "Men", "02", ifelse(df$sex_group == "Women", "01", "missing"))),
                                df$sex_group, sep = "")

       	                df$lci_arrow <- ifelse(df$grs_lci_or < 0.5, 0.5, NA)
               	        df$uci_arrow <- ifelse(df$grs_uci_or > 6.8, 5.8, NA)

        	        critical_p <- 0.05/51
			critical_het_p <- 0.05/48
	                df$original_cochrans_p <- df$cochrans_p
        	        df$original_grs_p <- df$grs_p
                	df$sig_or <- ifelse(df$grs_p < critical_p, df$grs_or, NA)
	                df$not_sig_or <- ifelse(df$grs_p >= critical_p, df$grs_or, NA)		

			save <- df

	                df[, c("grs_p", "cochrans_p")] <- lapply(df[, c("grs_p", "cochrans_p")],
        	                function(x) ifelse(x < 0.01, paste(formatC(x,
                	        1, format = "e")),
	                        formatC(x, 2, format = "f")))

        	        df$print_grs_p <- ifelse(df$original_grs_p < (1*10^-200), "\"<1.0 x\"~10^-200", 
				ifelse(df$original_grs_p < 0.01, gsub("^ ", "", paste("\"", gsub("e-0|e-", 
					" x\"~10^-", formatC(df$grs_p, 1, format = "e")), sep = "")), 
				paste("\"", formatC(df$grs_p, 2, format = "f"), "\"", sep = "")))

			df$print_cochrans_p <- ifelse(is.na(df$original_cochrans_p), "", 
				ifelse(df$original_cochrans_p < (1*10^-200), "\"<1.0 x\"~10^-200",
				ifelse(df$original_cochrans_p < 0.01, gsub("^ ", "",
 				paste("\"", gsub("e-0|e-", " x\"~10^-", formatC(df$cochrans_p, 1, format = "e")), sep = "")), 
				paste("\"", formatC(df$original_cochrans_p, 2, format = "f"), "\"", sep = ""))))

			df <- df[order(df$unique_combinations, decreasing = T), ]
        	        rownames(df) <- NULL

			xlim_number <- c(0.3, nrow(df) +1.5)
			 gray_vline <- nrow(df) - as.integer(rownames(df[df$trait_name %in% c("CAD", "COPD", "NAFLD", "Renal failure - acute",
                                               "Stroke", "Type 1 diabetes"), ])) + 1
			height <- 22.23
			
                        breaks_number <- c(0.5, 1, 5)
                        title_name <- "Odds Ratio (95% CI) per 1-SD higher obesity trait"
			ylim_number <- c(0.015, 100)
			col_placement <- c(0.016, 0.08, 0.22, 16, 40, 85, 95)
			header_placement <- c(0.016, 0.08, 0.22, 1.1, 5.5, 19, 45)
			vjust_number <- ifelse(type == "article", -236, 
						ifelse(type == "thesis", 300, -50))

			if (type == "article") {
				palette_colors <- "Dark2"
				shape_empty <- 23
				shape_fill <- 18
				font <- "Times"
				size_fill <- 2.3
			} else if (type == "thesis") {
				palette_colors <- "Dark2"
				shape_empty <- 23
                                shape_fill <- 18
				font <- "Times"
				size_fill <- 2.3
			}

                        df[duplicated(df[, c("trait", "grs_trait_name")]), "grs_trait_name"] <- ""
			df[duplicated(df$trait_name), "trait_name"] <- ""
			df[df$sex_group == "Men", "print_cochrans_p"] <- ""
			df$cochran_star <- ifelse(!is.na(df$print_cochrans_p) & df$print_cochrans_p != "" & 
				df$original_cochrans_p < critical_het_p, "*", "")		

			plot <- ggplot(df, aes(x=unique_combinations, y=grs_or, ymin=grs_lci_or, ymax=grs_uci_or)) +
                        geom_point(aes(color = sex_group), color = "white", shape = shape_fill, size = 3, show.legend =F) +
                        geom_vline(xintercept = gray_vline,
                                colour = "gray87", size = 4) +
			geom_segment(aes(yend = 1, y = 1, xend = 0, x = nrow(df) + 0.5), linetype = "dashed", size = 0.3) +
			geom_segment(aes(yend = breaks_number[1], y = breaks_number[1], xend = 0, x = nrow(df) + 0.5), size = 0.3) +
			geom_segment(aes(yend = breaks_number[1], y = breaks_number[3], xend = 0.4, x = 0.4), size = 0.3) +
                        geom_errorbar(size = 0.3, width=.3) +
                        geom_point(aes(y = df$not_sig_or, color = sex_group), shape = shape_empty, fill = "white", size = 1.65,
                                show.legend =F) +
                        geom_point(aes(y = df$sig_or, color = sex_group), shape = shape_fill, size = size_fill, show.legend =F) +
                        geom_segment(aes(yend = lci_arrow, y = lci_arrow,  xend = unique_combinations), size = 0.1,
                                arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
                        geom_segment(aes(yend = (uci_arrow+1), y = uci_arrow, xend = unique_combinations), size = 0.1,
                                arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
			geom_text(aes(y = col_placement[1], x = unique_combinations), label = df$trait_name, family = font, size = 2.85, hjust = 0) +
			geom_text(aes(y = col_placement[2], x = unique_combinations), label = df$grs_trait_name, family = font, size = 2.85, hjust = 0) +
			geom_text(aes(y = col_placement[3], x = unique_combinations), label = df$sex_group, family = font, size = 2.85, hjust = 0) +
			geom_text(aes(y = col_placement[4], x = unique_combinations), label = df$ci, family = font, size = 2.85, hjust = 1) +
			geom_text(aes(y = col_placement[5], x = unique_combinations), label = df$print_grs_p, parse = T, family = font, size = 2.85, hjust = 1) +
			geom_text(aes(y = col_placement[6], x = unique_combinations), label = df$print_cochrans_p, parse = T, family = font, size = 2.85, hjust = 1, vjust = -0.5) +
			geom_text(aes(y = col_placement[7], x = unique_combinations), label = df$cochran_star, family = font, size = 3, hjust = 1, vjust = -0.5) +
			annotate("text", y = header_placement, x = rep(nrow(df) +1, length(header_placement)), 
				label = c("bold(Outcome)", "bold(\"Risk factor\")", "bold(Sex)", "bold(Estimate)", "bold(\"OR (95% CI)\")", 
				"bold(P)", "bold(P[het])"), parse = T, size = 2.85, hjust = 0, family = font, fontface = "bold") +
			theme_void() +
			scale_color_brewer(palette = palette_colors) +
                        scale_y_continuous(trans = "log", name = title_name, breaks = breaks_number) +
                        theme(axis.text.y=element_blank(), axis.text = element_text(color = "black", size = 8, family = font)) +
                        coord_flip(ylim =ylim_number, xlim = xlim_number, expand = F) +
                        theme(plot.margin=unit(c(-0.35, 0, 0.4, 0), "cm")) +
			ggtitle(title_name) +	
			theme(plot.title = element_text(family = font, size = 8, vjust = vjust_number, hjust = 0.55, face = "bold")) 

			output_file <- paste("/Net/fs1/home/linc4222/new.mr.results.", type, ".", dataset, ".", eth_group, ".", unit, 
				ifelse(type == "article", ".180816.pdf", ".180816.jpeg"), sep = "")
                        ggsave(output_file, plot,  unit = "cm", height = height, width = 15.3, dpi = 800, 
				device = ifelse(type == "article", "pdf", "jpeg"))

			}                        
                }
        }
}


##############################################################################################
############## MAKE IMAGE WITH SIG. SEX-SPECIFIC ANALYSES ONLY + WINNERS SNPS ##############
#############################################################################################

for (type in c("winner_comb", "winner_sex", "winner_unweighted")) {
	df_raw <- read.table("../results.mr.180730/sens.ipd.mr.binary.results.180815.txt", 
		stringsAsFactors = F, header =T, sep = "\t")
		
	if (type == "winner_comb") {
		df <- df_raw[df_raw$eth_group == "all.white" & df_raw$exposure_unit == "sd" &
			df_raw$snp_group %in% unique(df_raw$snp_group[grep("\\.comb\\.pulit\\.winner$", df_raw$snp_group)]) &
			!is.na(df_raw$cochrans_p) &
			!(df_raw$trait %in% c("smoker_cases", "t2d_cases_prob")), ]
	} else if (type == "winner_sex") {
		df <- df_raw[df_raw$eth_group == "all.white" & df_raw$exposure_unit == "sd" &
                        df_raw$snp_group %in% unique(df_raw$snp_group[grep("men\\.pulit\\.winner$", df_raw$snp_group)]) &
                        !is.na(df_raw$cochrans_p) &
                        !(df_raw$trait %in% c("smoker_cases", "t2d_cases_prob")), ]
	} else if (type == "winner_unweighted") {
                df <- df_raw[df_raw$eth_group == "all.white" & df_raw$exposure_unit == "sd" &
                        df_raw$snp_group %in% unique(df_raw$snp_group[grep("\\.comb\\.pulit\\.winner_unweighted", df_raw$snp_group)]) &
                        !is.na(df_raw$cochrans_p) &
                        !(df_raw$trait %in% c("smoker_cases", "t2d_cases_prob")), ]
	}

	dict_traits <- list(copd_cases = c("COPD", "3"),
		renal_failure_cases = c("Renal failure", "2"),
		ckd_cases = c("Renal failure - chronic", "1"),
		t2d_cases_probposs = c("Type 2 diabetes", "4"))

	for (i in 1:length(dict_traits)) {
		df[df$trait == names(dict_traits[i]), "order_article"] <- dict_traits[names(dict_traits[i])][[1]][2]
		df[df$trait == names(dict_traits[i]), "trait_name"] <- dict_traits[names(dict_traits[i])][[1]][1]
	}	

	df$grs_trait_name <- gsub("\\.(.)+", "", df$snp_group)
	df$sex_group <- ifelse(df$sex_group == "comb", "Combined",
		ifelse(df$sex_group == "men", "Men",
	        ifelse(df$sex_group == "women", "Women", "missing")))

	df$ci <- paste(formatC(df$grs_or, digits = 2, format = "f")," (",
		formatC(df$grs_lci_or, digits = 2, format = "f"), ",",
	        formatC(df$grs_uci_or, digits = 2, format = "f"), ")", sep = "")

	df$grs_trait_name <- ifelse(df$grs_trait_name == "bmi", "BMI",
		ifelse(df$grs_trait_name == "whr", "WHR",
	        ifelse(df$grs_trait_name == "whradjbmi", "WHRadjBMI", NA)))

	df$unique_combinations <- paste(df$order_article, "_", df$trait_name, "_", ifelse(df$grs_trait_name == "BMI", "03",
		ifelse(df$grs_trait_name == "WHR", "02",
		ifelse(df$grs_trait_name == "WHRadjBMI", "01", "missing"))), ifelse(df$sex_group == "Combined", "03",
		ifelse(df$sex_group == "Men", "02", ifelse(df$sex_group == "Women", "01", "missing"))),
		df$sex_group, sep = "")

	critical_p <- 0.05/51
	critical_het_p <- 0.05/48
	df$original_cochrans_p <- df$cochrans_p
	df$original_grs_p <- df$grs_p
	df$sig_or <- ifelse(df$grs_p < critical_p, df$grs_or, NA)
	df$not_sig_or <- ifelse(df$grs_p >= critical_p, df$grs_or, NA)

	df[, c("grs_p", "cochrans_p")] <- lapply(df[, c("grs_p", "cochrans_p")],
		function(x) ifelse(x < 0.01, paste(formatC(x,
	        1, format = "e")),
        	formatC(x, 2, format = "f")))

	df$print_grs_p <- ifelse(df$original_grs_p < (1*10^-200), "\"<1.0 x\"~10^-200",
		ifelse(df$original_grs_p < 0.01, gsub("^ ", "", 
		paste("\"", gsub("e-0|e-", " x\"~10^-", formatC(df$grs_p, 1, format = "e")), sep = "")),
	        paste("\"", formatC(df$grs_p, 2, format = "f"), "\"", sep = "")))

	df$print_cochrans_p <- ifelse(is.na(df$original_cochrans_p), "",
        	ifelse(df$original_cochrans_p < (1*10^-200), "\"<1.0 x\"~10^-200",
	        ifelse(df$original_cochrans_p < 0.01, gsub("^ ", "", 
		paste("\"", gsub("e-0|e-", " x\"~10^-", formatC(df$cochrans_p, 1, format = "e")), sep = "")),
		paste("\"", formatC(df$original_cochrans_p, 2, format = "f"), "\"", sep = ""))))

	df <- df[order(df$unique_combinations, decreasing = T), ]
	rownames(df) <- NULL

	xlim_number <- c(0.3, nrow(df) +1.5)

	breaks_number <- c(0.5, 1, 5)
	title_name <- "Odds Ratio (95% CI) per 1-SD higher obesity trait"
	ylim_number <- c(0.015, 100)
	col_placement <- c(0.016, 0.08, 0.22, 16, 40, 85, 95)
	header_placement <- c(0.016, 0.08, 0.22, 1.1, 6.25, 19.8, 45.1)

	shape_empty <- 23
	shape_fill <- 18
	font <- "Times"
	size_fill <- 2.3

	save <- df
	df[duplicated(df[, c("trait", "grs_trait_name")]), "grs_trait_name"] <- ""
	df[duplicated(df$trait_name), "trait_name"] <- ""
	df[df$sex_group == "Men", "print_cochrans_p"] <- ""

        df$cochran_star <- ifelse(!is.na(df$print_cochrans_p) & df$print_cochrans_p != "" &
                          df$original_cochrans_p < critical_het_p, "*", "")
	
	plot <- ggplot(df, aes(x=unique_combinations, y=grs_or, ymin=grs_lci_or, ymax=grs_uci_or)) +
	geom_point(aes(color = sex_group), color = "white", shape = shape_fill, size = 3, show.legend =F) +
	geom_segment(aes(yend = 1, y = 1, xend = 0, x = nrow(df) + 0.5), linetype = "dashed", size = 0.3) +
	geom_segment(aes(yend = breaks_number[1], y = breaks_number[1], xend = 0, x = nrow(df) + 0.5), size = 0.3) +
	geom_segment(aes(yend = breaks_number[1], y = breaks_number[3], xend = 0.4, x = 0.4), size = 0.3) +
	geom_errorbar(size = 0.3, width=.3) +
	geom_point(aes(y = df$not_sig_or, color = sex_group), shape = shape_empty, fill = "white", size = 1.65, show.legend =F) +
	geom_point(aes(y = df$sig_or, color = sex_group), shape = shape_fill, size = size_fill, show.legend =F) +
	geom_text(aes(y = col_placement[1], x = unique_combinations), label = df$trait_name, family = font, size = 2.85, hjust = 0) +
	geom_text(aes(y = col_placement[2], x = unique_combinations), label = df$grs_trait_name, family = font, size = 2.85, hjust = 0) +
	geom_text(aes(y = col_placement[3], x = unique_combinations), label = df$sex_group, family = font, size = 2.85, hjust = 0) +
	geom_text(aes(y = col_placement[4], x = unique_combinations), label = df$ci, family = font, size = 2.85, hjust = 1) +
	geom_text(aes(y = col_placement[5], x = unique_combinations), label = df$print_grs_p, parse = T, family = font, size = 2.85, hjust = 1) +
	geom_text(aes(y = col_placement[6], x = unique_combinations), label = df$print_cochrans_p, parse = T, family = font, size = 2.85, 
		hjust = 1, vjust = -0.5) +
	geom_text(aes(y = col_placement[6]+8, x = unique_combinations), label = df$cochran_star, family = font, size = 3, hjust = 1, 
		vjust = -0.5) +
	annotate("text", y = header_placement, x = rep(nrow(df) +1, length(header_placement)),
	        label = c("bold(Outcome)", "bold(\"Risk factor\")", "bold(Sex)", "bold(Estimate)", "bold(\"OR (95% CI)\")",
        	"bold(P)", "bold(P[het])"), parse = T, size = 2.85, hjust = 0, family = font, fontface = "bold") +
	theme_void() +
	scale_colour_manual(values=c("#d95f02", "#7570b3")) +
	scale_y_continuous(trans = "log", name = title_name, breaks = breaks_number, ) +
	theme(axis.text.y=element_blank(), axis.text = element_text(color = "black", size = 8, family = font)) +
	coord_flip(ylim =ylim_number, xlim = xlim_number, expand = F) +
	theme(plot.margin=unit(c(-0.35, 0, 0.4, 0), "cm")) +
	ggtitle(title_name) +
	theme(plot.title = element_text(family = font, size = 8, vjust = -50.5, hjust = 0.55, face = "bold"))
	
	output_file <- paste0("/Net/fs1/home/linc4222/grs.pic.sex.het.mr.results.pulit.eur.", type, ".sd.180816.pdf")
	ggsave(output_file, plot, unit = "cm", height = 5, width = 17.3, dpi = 800, device = "pdf")
}

##################################################################
########### PLOT OF THE MRs WITH THE RISK FACTORS ################
##################################################################

#The FG, FI MRs - subset to the pulit.sig, IVW method
summary <- read.table("../results.mr.180730/summary.mr.results.180730.txt",
        stringsAsFactors = F, header = T, sep = "\t")
summary$grs_trait <- gsub("\\.(.)+", "", summary$snp_group)
summary <- summary[summary$snp_group %in% unique(summary$snp_group[grep("pulit.sig", summary$snp_group)]) &
	summary$method == "IVW", c("grs_trait", "sex_group", "trait", "beta", "se", "beta_lci", "beta_uci", "p", "cochrans_p")]

#The SBP, DBP MRs
bp <- read.table("../results.mr.180730/ipd.mr.continuous.results.180815.txt", stringsAsFactors = F, header =T, sep = "\t")
bp <- bp[bp$snp_group %in% unique(bp$snp_group[grep("pulit.sig", bp$snp_group)]) &
	bp$trait %in% c("dbp", "sbp") & bp$exposure_unit == "sd" & bp$outcome_unit == "sd" & bp$eth_group == "all.white", 
	c("grs_trait", "sex_group", "trait", "grs_beta", "grs_se", "grs_lci_beta", "grs_uci_beta", "grs_p", "cochrans_p")]

#Merge summary FG and FI with BP
colnames(bp) <- colnames(summary)
summary_bp <- rbind(summary, bp)

#The smoking MRs - NOTE THAT IT'S NOT BETA, BUT OR!!! Just to make plotting easier with same headings
smok <- read.table("../results.mr.180730/ipd.mr.binary.results.180815.txt", stringsAsFactors = F, header =T, sep = "\t")
smok <- smok[smok$snp_group %in% unique(smok$snp_group[grep("pulit.sig", smok$snp_group)]) &
	smok$trait == "smoker_cases" & smok$exposure_unit == "sd" & smok$eth_group == "all.white", 
	c("grs_trait", "sex_group", "trait", "grs_or", "grs_se", "grs_lci_or", "grs_uci_or", "grs_p", "cochrans_p")]
colnames(smok) <- colnames(summary_bp)

for (risk_factor in c("cont", "smok")) {
	if (risk_factor == "cont") {
		df <- summary_bp
	} else if (risk_factor == "smok") {
		df <- smok
	}

	dict_traits <- list(FG = c("FG", 4), 
			FI = c("FI", 3), 
			dbp = c("DBP", 1), 
			sbp = c("SBP", 2),
			smoker_cases = c("Smoker", 5))

	for (i in 1:length(dict_traits)) {
        	df[df$trait == names(dict_traits[i]), "order"] <- dict_traits[names(dict_traits[i])][[1]][2]
	        df[df$trait == names(dict_traits[i]), "trait_name"] <- dict_traits[names(dict_traits[i])][[1]][1]
	}

	df$grs_trait_name <- ifelse(df$grs_trait == "bmi", "BMI", 
				ifelse(df$grs_trait == "whr", "WHR", 
				ifelse(df$grs_trait == "res_whr_inv", "WHRadjBMI", 
				ifelse(df$grs_trait == "whradjbmi", "WHRadjBMI", NA))))

	df$sex_group_name <- ifelse(df$sex_group == "comb", "Combined",
        	ifelse(df$sex_group == "men", "Men",
        	ifelse(df$sex_group == "women", "Women", "missing")))

	#Note that for smoking is OR
	df$ci <- paste(formatC(df$beta, digits = 2, format = "f")," (",
        	formatC(df$beta_lci, digits = 2, format = "f"), ",",
        	formatC(df$beta_uci, digits = 2, format = "f"), ")", sep = "")

	df$unique_combinations <- paste(df$order, "_", df$trait, "_", ifelse(df$grs_trait_name == "BMI", "03",
        	ifelse(df$grs_trait_name == "WHR", "02",
	        ifelse(df$grs_trait_name == "WHRadjBMI", "01", "missing"))), ifelse(df$sex_group_name == "Combined", "03",
        	ifelse(df$sex_group_name == "Men", "02", ifelse(df$sex_group_name == "Women", "01", "missing"))),
	        df$sex_group, sep = "")

	critical_p <- 0.05/15
	critical_het_p <- 0.05/15
	df$original_cochrans_p <- df$cochrans_p
	df$original_grs_p <- df$p
	df$sig_or <- ifelse(df$p < critical_p, df$beta, NA)
	df$not_sig_or <- ifelse(df$p >= critical_p, df$beta, NA)

	df[, c("p", "cochrans_p")] <- lapply(df[, c("p", "cochrans_p")],
        	function(x) ifelse(x < 0.01, paste(formatC(x,
	        1, format = "e")),
        	formatC(x, 2, format = "f")))

	df$print_grs_p <- ifelse(df$original_grs_p < (1*10^-200), "\"<1.0 x\"~10^-200",
        	ifelse(df$original_grs_p < 0.01, 
		gsub("^ |^  ", "", paste("\"", gsub("e-0|e-", " x\"~10^-", formatC(df$p, 1, format = "e")), sep = "")),
        	gsub(" ", "", paste("\"", formatC(df$p, 2, format = "f"), "\"", sep = ""))))

	df$print_cochrans_p <- ifelse(is.na(df$original_cochrans_p), "",
        	ifelse(df$original_cochrans_p < (1*10^-200), "\"<1.0 x\"~10^-200",
	        ifelse(df$original_cochrans_p < 0.01, gsub("^ ", "",
        	paste("\"", gsub("e-0|e-", " x\"~10^-", formatC(df$cochrans_p, 1, format = "e")), sep = "")),
	        paste("\"", formatC(df$original_cochrans_p, 2, format = "f"), "\"", sep = ""))))

	df$cochrans_p_star <- ifelse(df$original_cochrans_p < critical_het_p, "*", "")
	df <- df[order(df$unique_combinations, decreasing = T), ]
	rownames(df) <- NULL

        save <- df
        df[duplicated(df[, c("trait", "grs_trait_name")]), "grs_trait_name"] <- ""
        df[duplicated(df$trait_name), "trait_name"] <- ""
        df[df$sex_group == "men", c("print_cochrans_p", "cochrans_p_star")] <- ""

        if (risk_factor == "cont") {
		dashed_line_place <- 0
		axis_ends <- c(-0.02, 0.35)
		axis_xend <- 0
		xlim_number  <- c(-0.5, nrow(df) + 1.5)
		ylim_number <- c(-0.5, 0.9)
		breaks_number <- c(0, 0.1, 0.2, 0.3)
		estimate_name <- "bold(\"Beta (95% CI)\")"
		col_placement <- c(-0.5, -0.35, -0.15, 0.55, 0.7, 0.82)
		header_placement <- c(-0.5, -0.35, -0.15, 0.17, 0.405, 0.595, 0.77)
		trans <- "identity"
		height <- 15
		vline <- c(10:18, 28:36)
		estimate_label <- "Beta (95% CI) per 1-SD higher obesity trait"
		estimate_vjust <- -148
		estimate_hjust <- 0.47
        } else if (risk_factor == "smok") {
		dashed_line_place <- 1
		axis_ends <- c(0.95, 1.6)
		axis_xend <- 0.4
		xlim_number <- c(0.3, nrow(df) +3)
		ylim_number <- c(0.4, 5)
		breaks_number <- c(1, 1.5)
		estimate_name <- "bold(\"Odds Ratio (95% CI)\")"
		col_placement <- c(0.42, 0.58, 0.78, 2.2, 3.3, 4.5)
		header_placement <- c(0.42, 0.58, 0.78, 1.225, 1.67, 2.7, 3.7)
		trans <- "log"
		height <- 5
		vline <- 30
		estimate_label <- "Odds Ratio (95% CI) per 1-SD higher obesity trait"
		estimate_vjust <- -41.5
		estimate_hjust <- 0.44
        }

	plot <- ggplot(df, aes(x=unique_combinations, y=beta, ymin=beta_lci, ymax=beta_uci)) +
	geom_point(aes(color = sex_group), color = "white", shape = 18, size = 3, show.legend =F) +
	geom_vline(xintercept = vline, colour = "gray87", size = 5.1) +
	geom_segment(aes(yend = dashed_line_place, y = dashed_line_place, xend = 0, x = nrow(df) + 0.5), linetype = "dashed", size = 0.3) +
	geom_segment(aes(yend = axis_ends[1], y = axis_ends[1], xend = axis_xend, x = nrow(df) + 0.5), size = 0.3) +
	geom_segment(aes(yend = axis_ends[1], y = axis_ends[2], xend = axis_xend, x = axis_xend), size = 0.3) +
	annotate("segment", y = breaks_number, yend = breaks_number, x = rep(0, length(breaks_number)), 
		xend = rep(-0.3, length(breaks_number)), size = 0.3) +
	geom_errorbar(size = 0.3, width=.3) +
	geom_point(aes(y = df$not_sig_or, color = sex_group), shape = 23, fill = "white", size = 1.65, show.legend =F) +
	geom_point(aes(y = df$sig_or, color = sex_group), shape = 18, size = 2.3, show.legend =F) +
	geom_text(aes(y = col_placement[1], x = unique_combinations), label = df$trait_name, family = "Times", size = 2.85, hjust = 0) +
	geom_text(aes(y = col_placement[2], x = unique_combinations), label = df$grs_trait_name, family = "Times", size = 2.85, hjust = 0) +
	geom_text(aes(y = col_placement[3], x = unique_combinations), label = df$sex_group_name, family = "Times", size = 2.85, hjust = 0) +
	geom_text(aes(y = col_placement[4], x = unique_combinations), label = df$ci, family = "Times", size = 2.85, hjust = 1) +
	geom_text(aes(y = col_placement[5], x = unique_combinations), label = df$print_grs_p, parse = T, 
		family = "Times", size = 2.85, hjust = 1) +
	geom_text(aes(y = col_placement[6], x = unique_combinations), label = df$print_cochrans_p, parse = T, 
		family = "Times", size = 2.85, hjust = 1, vjust = -0.5) +
	geom_text(aes(y = col_placement[6] + 0.01, x = unique_combinations), label = df$cochrans_p_star, family = "Times", 
		size = 3.5, hjust = 0, vjust = -0.5) +
	annotate("text", y = header_placement, x = rep(nrow(df) +1, length(header_placement)),
	        label = c("bold(Outcome)", "bold(\"Risk factor\")", "bold(Sex)", "bold(Estimate)", estimate_name, 
	        "bold(P)", "bold(P[het])"), parse = T, size = 2.85, hjust = c(rep(0, 3), 0.5, rep(0, 3)), family = "Times", fontface = "bold") +
	ggtitle(estimate_label)+
	theme_void() +
	scale_color_brewer(palette = "Dark2") +
	scale_y_continuous(trans = trans, breaks = breaks_number) +
	theme(axis.text.y=element_blank(), axis.text = element_text(color = "black", size = 8, family = "Times")) +
	coord_flip(ylim = ylim_number, xlim = xlim_number, expand = F) +
	theme(plot.title = element_text(family = "Times", size = 8, vjust = estimate_vjust, hjust = estimate_hjust, face = "bold"), 
		plot.margin=unit(c(0, 0, 1, 0), "cm")) 

	output_file <- paste("/Net/fs1/home/linc4222/mr.", risk_factor, ".mr.risk.estimates.all.white.sd.pdf", sep = "")
	ggsave(output_file, plot, unit = "cm", height = height, width = 17.3, dpi = 800, device = "pdf")


}

