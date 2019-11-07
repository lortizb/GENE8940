suppressMessages({
 library("sleuth")
})

#set input and output dirs
datapath <- "/work/gene8940/lo50936/HMW06/kallisto_results"  # you need to modify this line
resultdir <- "/work/gene8940/lo50936/HMW06/kallisto_results"   # you need to modify this line
setwd(resultdir)

#create a sample to condition metadata description
sample_id <- dir(file.path(datapath))
kallisto_dirs <- file.path(datapath, sample_id)
print(kallisto_dirs)
sample <- c("SRR5344681", "SRR5344682", "SRR5344683", "SRR5344684")
condition <- c("WT", "WT", "dFNR", "dFNR")
samples_to_conditions <- data.frame(sample,condition)
samples_to_conditions <- dplyr::mutate(samples_to_conditions, path = kallisto_dirs)
print(samples_to_conditions)

#run sleuth on the data
sleuth_object <- sleuth_prep(samples_to_conditions, extra_bootstrap_summary = TRUE) # read data into sleuth_object
sleuth_object <- sleuth_fit(sleuth_object, ~condition, 'full') # estimate parameters for the full linear model that includes the conditions as factors
sleuth_object <- sleuth_fit(sleuth_object, ~1, 'reduced') # estimate parameters for the reduced linear model that assumes equal transcript abundances in both conditions
sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full') # perform likelihood ratio test to identify transcripts whose fit is significantly better under full model relative to reduced model
models(sleuth_object)

#summarize the sleuth results and view 20 most significant DE transcripts
sleuth_table <- sleuth_results(sleuth_object, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 2)

#plot an example DE transcript result
differentially_expressed_gene1="lcl|NC_000913.3_cds_NP_414965.1_430"     # you need to modify this line
differentially_expressed_gene2="lcl|NC_000913.3_cds_NP_414964.1_429" # you need to modify this line

p1 <- plot_bootstrap(sleuth_object, differentially_expressed_gene1, units = "est_counts", color_by = "condition")
p2 <- plot_bootstrap(sleuth_object, differentially_expressed_gene2, units = "est_counts", color_by = "condition")

#Print out the plots created above and store in a single PDF file
pdf(file="SleuthResults.pdf")
print(p1)
print(p2)
dev.off()

#Quit R
quit(save="no")
