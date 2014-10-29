# code to combine p values from numbers and spacing of DNMs

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis"
DATA_DIR = file.path(CODE_DIR, "data")

fishersMethod <- function(x) {
    # function to combine p values
    pchisq(-2 * sum(log(x)), df=2 * length(x), lower=FALSE)
}


# read in p values from clustering analysis, only for genes with >1 mutation
spacing = read.delim(file.path(DATA_DIR, "de_novo_distance_simulations.increased_set_all.geometric_mean.txt"), header=TRUE)

# read in p values from mupit/daly analyses
probs = read.delim(file.path(DATA_DIR, "Mup-it_Daly_020514_META_undiagnosed.output.txt"), header=TRUE)

# merge the datasets
merged = merge(probs, spacing, by.x=1, by.y=1, all.x=TRUE)

# calculate a combined p-value for each gene
p_values = merged[, c("func_dist_probability", "daly.p.DNM.func")]
combined.p.dist = apply(p_values, 1, fishersMethod)

# adjust the P values by false-discovery rate
num.tests = 18500
combined.fdr.dist = p.adjust(combined.p.dist, method="BH", n=num.tests)

# write the results to a file
output = cbind(merged, combined.p.dist, combined.fdr.dist)
write.table(output, file=file.path(DATA_DIR, "Mup-it_Daly_020514_META_undiagnosed.output.combined.txt"), row.names=F, quote=F, sep="\t")

