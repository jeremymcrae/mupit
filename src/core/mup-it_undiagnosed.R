# program to calculate the significance of seeing N DNMs of a specific
# combination of functional types in a particular gene in M trios

# RATIONALE: use gene coding sequence to predict rate of DNMs in coding sequence
# for each gene of different functional classes, then estimate the probability
# of seeing the observed combination of different functional classes of DNMs
# assuming number of DNMs in each class is Poisson distributed

# initial implementation: use genome-wide mutation rate and scale by length of
# coding sequence, use genome-wide average of functional consequences of coding
# mutations from Kryukov et al 2007


# POTENTIAL FUTURE IMPROVEMENTS (highest priority first):
#  - adapt for chrX, lower mutation rate, and number of transmissions required.
#        Assume male/female of proband = 1:1, unless known [DONE]
#  - incorporate FDR estimation [DONE]
#  - convert to analytical approach from permutation approach [DONE]
#  - use actual number of exons to predict essential_splice_site mutations,
#  - use base composition of coding sequence to predict gene-specific mutation
#        rate for each class of mutation, rather than just scaling the
#        genome-average.
#  - use estimate of de novo mutation discovery power in a gene to better
#        estimate gene-specific mutation rate
#  - get CDS length for all genes, not just this subset, 455/478 in test data.
#        use longest transcript if >1 [DONE]
#  - calculate coding sequence length according to intersection of exome
#        targeted regions and the union of all transcripts for a gene.
#  - account for incomplete sensitivity for DNMs, especially indels
#  - look at clustering of de novos in genes with recurrent mutations, within
#        protein space


# format: read in table of genes and the numbers of families with SNV DNMs
# within different functiona classes and indel DNMS in different functional
# classes, output same table with added column of probability of seeing that
# combination of DNMs.

# input data (validated SNV DNMs in TSV format, HGNC_ID, NUM_LOF, NUM_MISSENSE,
# NUM_LOF_INDEL, NUM_MISSENSE_INDEL):

# note: current test input file is numbers of mutations, not number of families,
# some families have >1 mutation in the same gene

library(Cairo)
library(reshape)

######### SETTINGS ################

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis"
DATA_DIR = file.path(CODE_DIR, "data")
SRC_DIR = file.path(CODE_DIR, "src")
DE_NOVO_DIR = file.path(DATA_DIR, "de_novo_datasets")
source(file.path(SRC_DIR, "core", "mutation_rates_daly.R"))
source(file.path(SRC_DIR, "core", "mutation_rates_length.R"))
source(file.path(SRC_DIR, "core", "open_de_novo_datasets.R"))

CQ.LOF = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant")
CQ.MISSENSE = c("missense_variant", "initiator_codon_variant", "stop_lost", "inframe_deletion", "inframe_insertion")

get_ddd_diagnosed <- function() {
    # find diagnosed probands in the DDD study, to exclude them from our data
    # 
    # Returns:
    #     A list containing vectors with DECIPHER IDs, and sex of the diagnosed
    #     probands
    
    # read in samples that have been diagnosed, so as to remove from our data
    diagnoses = read.delim(file.path(DATA_DIR, "Diagnoses_1133.txt"), header=TRUE)
    diagnosed.index = which(rowSums(diagnoses[, c(4:14)]) > 0)
    
    diagnosed = list()
    diagnosed$id = diagnoses$decipher_id[diagnosed.index]
    diagnosed$sex = diagnoses$Sex[diagnosed.index]
    
    return(diagnosed)
}

get_de_novo_counts <- function(de_novos, lof_cq, missense_cq) {
    # tallies the mutation types observed for each gene
    # 
    # Args:
    #     de_novos: data frame listing all the de novo mutations, with columns
    #         for HGNC symbol, consequence type (VEP style predictions), and a
    #         column indicating SNV, or indel.
    #     lof_cq: vector of types for loss-of-function consequences
    #     missense_cq: vector of types for missense consequences
    # 
    # Returns:
    #     data frame with tally of de novo mutations for each of the mutation 
    #     types
    
    lof_regex = paste(lof_cq, collapse = "|")
    missense_regex = paste(missense_cq, collapse = "|")
    all_regex = paste(c(lof_cq, missense_cq), collapse = "|")
    
    # count the number of de novos for each type/consequence combination
    counts = cast(de_novos, HGNC ~ TYPE + CQ, value = "STUDY", length)
    indels = counts[, grep("INDEL", names(counts))]
    snvs = counts[, grep("SNV", names(counts))]
    
    # and sum the de novo counts in each type/consequence category
    de_novo_counts = data.frame(HGNC = counts$HGNC)
    de_novo_counts$lof.snvs = rowSums(data.frame(snvs[, grep(lof_regex, names(snvs))]))
    de_novo_counts$missense.snvs = rowSums(data.frame(snvs[, grep(missense_regex, names(snvs))]))
    de_novo_counts$lof.indels = rowSums(data.frame(indels[, grep(lof_regex, names(indels))]))
    de_novo_counts$missense.indels = rowSums(data.frame(indels[, grep(missense_regex, names(indels))]))
    
    return(de_novo_counts)
}

get_p_values <- function(rates, de_novos, counts, num.tests) {
    # tests whether genes are enriched with de novo mutations
    # 
    # Args:
    #     rates: gene mutation rates per consequence type
    #     de_novos: data frame containing all the observed de novos for all the 
    #         genes
    #     counts: data frame with tally of de novo mutations for each of the
    #         mutation types.
    #     num.tests: number of tests performed (used for multiple correction).
    # 
    # Returns:
    #     data frame with gene info, mutation rates and P values from testing
    #     for enrichment.
    
    # set-up vectors to store gene-specific information only for observed genes
    observed = data.frame(matrix(NA, nrow = nrow(counts), ncol = 8))
    names(observed) = c("chr", "coord", "snv.missense.rate", 
        "snv.lof.rate", "indel.missense.rate", "indel.lof.rate", "p.func", 
        "p.lof")

    # loop for each observed gene, test for functional variants and lof 
    # variants, using the gene mutation rates
    for (i in 1:nrow(counts)) {
        gene = as.character(counts$HGNC[i])
        
        # continue to next gene if mutation rates not available for the gene
        if (!(gene %in% rates$HGNC)) { next }
        
        gene.index = which(rates$HGNC == gene)
        
        # get the mutation rates for the gene
        observed[i, c("snv.missense.rate", "snv.lof.rate", "indel.missense.rate", "indel.lof.rate")] = rates[gene.index, c("snv.missense.rate", "snv.lof.rate", "indel.missense.rate", "indel.lof.rate")]
        
        # figure out the chromosome and nucleotide position of the gene
        data.index = which(de_novos$HGNC == gene)[1]
        observed$chr[i] = de_novos$CHROM[data.index]
        observed$coord[i] = de_novos$POS[data.index]
        
        # count the observed de novos in each functional category
        lof_count = sum(counts[i, c("lof.snvs", "lof.indels")])
        missense_count = sum(counts[i, c("missense.snvs", "missense.indels")])
        func_count = lof_count + missense_count
        
        # get the mutation rates for each funcational category
        lof_rate = sum(observed[i, c("snv.lof.rate", "indel.lof.rate")])
        missense_rate = sum(observed[i, c("snv.missense.rate", "indel.missense.rate")]) 
        func_rate = lof_rate + missense_rate
        
        # calculate the probability of observing said de novos, given the 
        # gene mutation rates
        observed$p.lof[i] = dpois(lof_count, lambda=lof_rate)
        observed$p.func[i] = dpois(func_count, lambda=func_rate)
        
        if (i %% 100 == 0) {
            print(paste(i, " out of ", nrow(counts), " genes", sep = ""))
        }
    }
    
    # correct the P values for multiple testing by false discovery rate
    observed$fdr.lof = p.adjust(observed$p.lof, method="BH", n=num.tests)
    observed$fdr.func = p.adjust(observed$p.func, method="BH", n=num.tests)
    
    return(observed)
}

label_genes <- function(enriched, p_values, num.tests) {
    # make plot labels for genes with fdr > threshold
    # 
    # Args:
    #     enriched: data frame containing HGNC symbols
    #     p_values: vector of p-values, sorted as per enriched data frame
    #     num.tests: number of tests performed (used for multiple correction).
    
    fdr = p.adjust(p_values, method="BH", n=num.tests)
    fdr.thresh = 0.05
    label.thresh = max(p_values[which(fdr < fdr.thresh)])
    thresh.index = which(p_values < label.thresh)
    num.thresh = length(thresh.index)
    
    # sometimes we don't have any genes with P values more significant than the 
    # FDR threshold, simply return, rather than trying to plot (also since when
    # num.thresh is zero, we get some zero-length errors)
    if (num.thresh == 0) {
        return()
    }
    
    for (i in 1:num.thresh) {
        text(x=thresh.index[i], -log10(p_values[thresh.index[i]]), labels=enriched$HGNC[thresh.index[i]], pos=3, cex=0.5)
    }
}

plot_graphs <- function(enriched, num.tests) {
    # make Manhattan plots for LOF and Func variants separately
    # 
    # Args:
    #     enriched: data frame containing columns for chr, coord (position), and
    #         p values from testing for enrichment with loss-of-function and 
    #         functional consequence variants. Both of these have been tested 
    #         with two sets of mutation rate data, one derived from length 
    #         based rates, and the other from mutation rates provided by Mark 
    #         Daly.
    #     num.tests: number of tests performed (used by Bonferroni and FDR 
    #         correction).
    
    plot_values <- function(length_p_vals, daly_p_vals, title, color_index) {
        # plot the results from de novos
        
        # set up the plot, starting with the length based P values
        plot(-log10(length_p_vals), col=color_index, pch=19, cex=0.75, 
            ylab="-log10(p)", xaxt="n", main=title, xlab="genome position",  
            ylim=c(0, max(-log10(length_p_vals), na.rm=TRUE) + 1))
        
        # plot line showing where Bonferroni correction would be
        abline(h = -log10(0.05 / num.tests), col="red", lty=2)
        
        # add the results from using the alternative mutation rates
        points(-log10(daly_p_vals), col=color_index, pch=1, cex=0.75)
        legend("topleft", legend=c("Length-based rates", "Daly group rates"), 
            pch=c(19, 1), col="darkblue", cex=0.75)
        # label genes which are significant following FDR correction
        label_genes(enriched, length_p_vals, num.tests)
    }
    
    Cairo(file="temp.pdf", type="pdf", width=20, height=20, units = "cm")
    
    enriched = na.omit(enriched)
    
    # set up alternating colors for successive chromosomes
    enriched$chr[enriched$chr == "X"] = "23"
    enriched = enriched[order(as.numeric(as.character(enriched$chr)), as.numeric(as.character(enriched$coord))), ]
    color_index = rep("lightblue3", nrow(enriched))
    odd.chr = c("1", "3", "5", "7", "9", "11", "13", "15", "17", "19", "21", "23")
    color_index[enriched$chr %in% odd.chr] = "darkblue"
    
    # plot the results from loss-of-function de novos
    plot_values(enriched$p.lof, enriched$daly.p.lof, "Loss-of-Function DNMs", color_index)
    plot_values(enriched$p.func, enriched$daly.p.func, "Functional DNMs", color_index)
    
    dev.off()
}

analyse_gene_enrichment <- function(de_novos, num.trios.male, num.trios.female) {
    # run the analysis of whether de novo mutations are enriched in genes
    # 
    # Args:
    #     de_novos: data frame containing all the observed de novos for all the 
    #         genes
    #     num.trios.male: number of trios with male offspring in the dataset
    #     num.trios.female: number of trios with female offspring in the dataset
    
    # tally the de novos by consequence and variant type
    de_novo_counts = get_de_novo_counts(de_novos, CQ.LOF, CQ.MISSENSE)
    
    # get the length-based, and Daly-based mutation rates for each gene
    cds_rates = get_length_based_rates(num.trios.male, num.trios.female)
    daly_rates = get_mutation_rates(num.trios.male, num.trios.female)
    
    # calculate p values for each gene using the different mutation rates
    num.tests = 18500
    p_vals_length = get_p_values(cds_rates, de_novos, de_novo_counts, num.tests)
    p_vals_daly = get_p_values(daly_rates, de_novos, de_novo_counts, num.tests)
    
    # write out results table
    enriched = cbind(de_novo_counts, p_vals_length, p_vals_daly)
    names(enriched) = c("HGNC", "LOF.snvs", "NS.snvs", "LOF.indels",
        "NS.indels", "chr", "coord", "snv.missense.rate",
        "snv.lof.rate", "indel.missense.rate", "indel.lof.rate", "p.func",
        "p.lof", "fdr.lof", "fdr.func", "chr", "coord",
        "daly.snv.missense.rate", "daly.snv.lof.rate",
        "daly.indel.missense.rate", "daly.indel.lof.rate", "daly.p.func",
        "daly.p.lof", "daly.fdr.lof", "daly.fdr.func")
    
    # write.table(enriched, file=file.path(DATA_DIR, "Mup-it_Daly_100414_META_undiagnosed.output.txt"), row.names=F, quote=F, sep="\t")
    
    plot_graphs(enriched, num.tests)
}


main <- function() {
    # here's an example of how to use the functions in this script
    diagnosed = get_ddd_diagnosed()
    num = get_trio_counts(diagnosed)
    num.trios.male = num$male
    num.trios.female = num$female
    
    # open the de novos, and 
    de_novos = open_datasets(diagnosed)
    
    analyse_gene_enrichment(de_novos, num.trios.male, num.trios.female)
}



