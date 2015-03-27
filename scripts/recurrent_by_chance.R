# estimate probability of seeing NN genes with >1 functional mutations given 
# number of mutations seen
# number of silent DNMs in DDD families, 1 gene with >1 mutation
# can adapt to get number expected with 3, 4, 5, 6, ... mutations to work out 
# PPV of each number

library(mupit)
library(plyr)

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

main <- function(){
    # could calculate mean analytically using mutation rate scaled to fit 
    # observed number of genes
    # sum(ppois(1, gene.func.rate*num.DNM.func/sum(gene.func.rate), lower.tail=F))
    # but how to calculate variance without simulation?
    
    de_novos = get_ddd_de_novos()
    counts = get_de_novo_counts(de_novos)
    counts$func = apply(counts[, c("lof_indel", "lof_snv", "missense_indel", "missense_snv")], 1, sum)
    
    # set number of functional DNMs and number of genes with >1 functional DNMs
    num.DNM.func = sum(counts$func > 0) # de novos from 4k trios cover 3774 genes
    num.recurr.genes = sum(counts$func > 1) # de novos from 4k trios are recurrent in 919 genes
    num.sims = 10000 # number of simulations to perform
    max.obs = as.integer(num.DNM.func/60) # max times a gene might be expected to be mutated
    
    trios = list()
    trios$male = 2408
    trios$female = 1887
    
    rates = get_mutation_rates(trios)
    gene.func.rate = rates$snv.missense.rate + rates$snv.lof.rate + 
        rates$indel.missense.rate + rates$indel.lof.rate
    
    store = simulate_recurrent_mutations(gene.func.rate, num.sims, num.DNM.func, max.obs)
    recurrent = rowSums(store[, 2:max.obs], na.rm=TRUE)
    
    if(any(!is.na(store[, max.obs]))) {
        stop("warning maximum observations reached, extend max obs")
    }
    
    # plot simulations versus observed
    #observed = number of genes with recurrent DNMs
    title = paste("Recurrent genes by chance (Functional DNMs=", num.DNM.func, ")", sep="")
    hist(recurrent, breaks=seq(-0.5, num.recurr.genes + 5, 1), col="green", 
        main=title, xlab="Number of genes with recurrent DNMs")
    abline(v=num.recurr.genes, col="red")
    
    # write out table for plotting boxplot of the same data
    write.table(recurrent, file="num_recurrent_simulation.txt", quote=FALSE, row.names=FALSE, sep="\t")
}

main()


