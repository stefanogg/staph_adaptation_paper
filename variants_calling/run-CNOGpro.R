# Script to run CNOGpro
args <- commandArgs(trailingOnly = TRUE)
iso <- args[1]
if (length(args) > 1) {
	window <- as.numeric(args[2])
} else {
	window <- 1000
}
# prob <- as.numeric(args[3])
prob <- 1e-4


# Run CNOGpro
cnv<- CNOGpro(hitsfile=paste0(iso,".hits"), gbkfile="reference/ref.gbk", windowlength=window)
cnv_norm <- normalizeGC(cnv)
cnv_hmm <- runHMM(cnv_norm, changeprob = prob)
saveRDS(cnv_hmm, file = paste0(iso, ".cnv.Rda"))
#cnv_bootstrap <- runBootstrap(cnv_norm)
#saveRDS(cnv_bootstrap, file = paste0(iso, .cnv.Rda")

# Extract HMM table and save
#cnv_hmm <- cnv_bootstrap
df_cnv <- cnv_hmm$HMMtable
write.table(df_cnv, file = paste0(iso, ".cnv.hmm.tab"))

# Filter regions with > 1 copy and save in bed format. This will raise an error and terminte the script if there are no cnv
df_cnv_clean <- subset(df_cnv, State > 1)
df_cnv_clean["chrom"] <- cnv_hmm$accession
df_cnv_clean <- df_cnv_clean[c("chrom","Startpos","Endpos","State")]
write.table(df_cnv_clean, file = paste0(iso, ".cnv.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

# Extract gene table and save
df_cnv_genes <- cnv_hmm$genes
write.table(df_cnv_genes, file = paste0(iso, ".cnv.genes.tab"), sep = "\t", quote = F, row.names = F, col.names = T)

# Filter regions with > 1 copy and save in bed format. This will raise an error and terminate the script if there are no cnv
df_cnv_genes_clean <- subset(df_cnv_genes, CN_HMM > 1)
df_cnv_genes_clean["chrom"] <- cnv_hmm$accession
df_cnv_genes_clean <- df_cnv_genes_clean[c("chrom","Left","Right","Locus","Strand","Type","Length","CN_HMM")]
write.table(df_cnv_genes_clean, file = paste0(iso, ".cnv.genes.bed"), sep = "\t", quote = F, row.names = F, col.names = F)






