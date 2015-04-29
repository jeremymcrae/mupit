# estimate probability of seeing NN genes with >1 functional mutations given
# number of mutations seen
# number of silent DNMs in DDD families, 1 gene with >1 mutation
# can adapt to get number expected with 3, 4, 5, 6, ... mutations to work out
# PPV of each number

library(Cairo)
library(mupit)
library(plyr)

DE_NOVOS_PATH = "/nfs/users/nfs_j/jm33/apps/mupit/data-raw/de_novo_datasets/de_novos.ddd_4k.ddd_only.txt"
RATES_PATH = "/nfs/users/nfs_j/jm33/apps/de_novo_clustering/results/de_novo_gene_rates.ddd_4k.meta-analysis.txt"

#' defines the cohort sizes, used to get the overall population size
#'
#' @param diagnosed list of sex and ID for probands diagnosed in the DDD study
#' @param meta true/false for whether to include meta-analysis populations
#'
#' @return list with total counts of trios with male and female offspring
get_trio_counts <- function(diagnosed=NULL, meta=FALSE) {
    
    # number of trios studied in our data
    male = 2408 # male probands
    female = 1887 # female probands
    
    # remove diagnosed patients, if maximising power
    male = male - sum(diagnosed$sex %in% c("Male", "male", "M", "m"))
    female = female - sum(diagnosed$sex %in% c("Female", "female", "F", "f"))
    
    if (meta) {
        male = male + sum(cohorts$unique_male)
        female = female + sum(cohorts$unique_female)
    }
    
    return(list(male=male, female=female))
}

#' load the mutation rates matched to the DDD de novos
get_rates_dataset <- function(rates_path) {
    rates = read.table(rates_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # convert from my column names to those used when estimating the gene
    # mutation rates given the cohort size
    rates$hgnc = rates$transcript_id
    rates$mis = rates$missense_rate
    rates$non = rates$nonsense_rate
    rates$splice_site = rates$splice_lof_rate
    rates$syn = rates$synonymous_rate
    rates$frameshift = rates$frameshift_rate
    
    rates = rates[, c("hgnc", "chrom", "length", "mis", "non", "splice_site", "syn", "frameshift")]
    
    return(rates)
}

#' load in the DDD de novos
#'
#' We load the de novos from the source file, so that we can futher restrict the
#' de novos to a high confidence set, since if we use the full set of denovos, we
#' overestimate the excess of recurrently mutated genes.
load_ddd_de_novos <- function(path, diagnosed=NULL) {
    de_novos = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # make sure we have the necessary columns for the gene counting
    de_novos$person_id = de_novos$person_stable_id
    de_novos$hgnc = de_novos$symbol
    de_novos$start_pos = de_novos$pos
    de_novos$type = "snv"
    de_novos$type[de_novos$var_type == "DENOVO-INDEL"] = "indel"
    
    # restrict the de novos to the highest quality set
    de_novos = de_novos[de_novos$pp_dnm > 0.9, ]
    
    # exclude the diagnosed probands
    if (!is.null(diagnosed)) {
        de_novos = de_novos[!de_novos$person_id %in% diagnosed$id, ]
    }
    
    return(de_novos)
}

#' simulate the counts of recurrent mutations using the mutation rates
#'
#' @param mutation_rates rates of mutation for each of the genes in the genome
#' @param num.sims number of simulations to perform
#' @param num.DNM.func number of genes observed with at least one functional
#'     mutation
#' @param max.obs highest number of recurrences in a gene that we expect
#'
#' @return matrix with the tallies of how many genes had one de novo, two de
#'     novos, three de novos etc at each simulation
simulate_recurrent_mutations <- function(mutation_rates, num.sims, num.DNM.func, max.obs) {
    
    # get the cumulative distribution of mutation rates for the genes
    summed.rates = cumsum(mutation_rates)
    
    # scale the cumulative distrubution between 0 and 1, so that we can randomly
    # select genes using runif, which defaults to between 0 and 1
    scaled.summed.rates = as.double(summed.rates/max(summed.rates))
    
    # randomly simulate de novos in genes, to count recurrently mutated genes
    simulated = matrix(nrow=num.sims, ncol=max.obs)
    for (i in 1:num.sims) {
        # get as many random numbers as ther are genes with functional mutations
        x = runif(num.DNM.func)
        # select genes, weighted by their mutation rate
        genes_selected = findInterval(x, scaled.summed.rates)
        
        # tally the number of genes by how many mutations they have, then store
        # these values from the current simulation
        recurrent_count = table(table(genes_selected))
        simulated[i, as.numeric(names(recurrent_count))] = recurrent_count
    }
    
    return(simulated)
}

# plot simulations of expected recurrence versus observed
plot_recurrent_distribution <- function(recurrent, func_n, observed_n) {
    Cairo(file="recurrent_coding_genes.pdf", type="pdf", height=15, width=15, units="cm")
    title = paste("Recurrent genes by chance (", func_n, " functional de novos)", sep="")
    hist(recurrent,
        breaks=seq(-0.5, max(observed_n, max(recurrent)) + 5, 5),
        col="green",
        main=title, xlab="Number of genes with recurrent DNMs", las=1)
    abline(v=observed_n, col="red")
    
    dev.off()
}

main <- function() {
    # could calculate mean analytically using mutation rate scaled to fit
    # observed number of genes
    # sum(ppois(1, gene_func_rate*num.DNM.func/sum(gene_func_rate), lower.tail=F))
    # but how to calculate variance without simulation?
    
    # define the diagnosed probands
    diagnosed_path = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx"
    diagnosed = get_likely_diagnosed(diagnosed_path)
    # diagnosed = NULL
    
    de_novos = load_ddd_de_novos(DE_NOVOS_PATH, diagnosed)
    trios = get_trio_counts(diagnosed)
    counts = get_de_novo_counts(de_novos)
    counts$func = apply(counts[, c("lof_indel", "lof_snv", "missense_indel", "missense_snv")], 1, sum)
    
    # set number of functional DNMs and number of genes with >1 functional DNMs
    func_n = sum(counts$func > 0) # de novos from 4k trios cover 3774 genes
    observed_n = sum(counts$func > 1) # de novos from 4k trios are recurrent in 919 genes
    num.sims = 10000 # number of simulations to perform
    max.obs = as.integer(func_n/60) # max times a gene might be expected to be mutated
    
    # ddd_rates = get_rates_dataset(RATES_PATH)
    # rates = get_gene_based_mutation_rates(trios, ddd_rates)
    rates = get_gene_based_mutation_rates(trios)
    gene_func_rate = rates$snv.missense.rate + rates$snv.lof.rate +
        rates$indel.missense.rate + rates$indel.lof.rate
    
    store = simulate_recurrent_mutations(gene_func_rate, num.sims, func_n, max.obs)
    recurrent = rowSums(store[, 2:max.obs], na.rm=TRUE)
    
    if(any(!is.na(store[, max.obs]))) {
        stop("warning maximum observations reached, extend max obs")
    }
    
    plot_recurrent_distribution(recurrent, func_n, observed_n)
    
    # get the 95% confidence interval for the number of excess genes
    error = qnorm(0.975) * sd(recurrent)
    high = observed_n - (mean(recurrent) - error)
    low = observed_n - (mean(recurrent) + error)
    print(c(low, high))
}

main()
