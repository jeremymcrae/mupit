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
source(file.path(SRC_DIR, "mutation_rates_daly.R"))

CQ.LOF = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant")
CQ.MISSENSE = c("missense_variant", "initiator_codon_variant", "stop_lost", "inframe_deletion", "inframe_insertion")

get_ddd_diagnosed <- function() {
    # find diagnosed probands in the DDD study, to exclude them from our data
    # 
    # Returns:
    #     A list containing vectors with DECIPHER IDs, and sex of the diagnosed
    #     probands
    
    # read in samples that have been diagnosed, so as to remove from our data
    diagnoses = read.delim(file.path(DATA_DIR, "Diagnoses_1133.txt"), header=T)
    diagnosed.index = which(rowSums(diagnoses[, c(4:14)]) > 0)
    
    diagnosed = list()
    diagnosed$id = diagnoses$decipher_id[diagnosed.index]
    diagnosed$sex = diagnoses$Sex[diagnosed.index]
    
    return(diagnosed)
}

get_trio_counts <- function(diagnosed) {
    # defines the cohort sizes, used to get the overall population size
    # 
    # Args:
    #     diagnosed: list of sex and ID for probands diagnosed in the DDD study
    # 
    # Returns:
    #     list with total counts of trios with male and female offspring
    
    # number of trios studied in our data
    male.ddd = 582 # trios with male offspring
    female.ddd = 548 # trios with female offspring
    
    # remove diagnosed patients, if maximising power
    male.ddd = male.ddd - length(which(diagnosed$sex == "Male"))
    female.ddd = female.ddd - length(which(diagnosed$sex == "Female"))
    
    # number of trios studied in deligt data
    male.deligt = 47
    female.deligt = 53
    
    # number of trios studied in autism data
    male.autism = 764
    female.autism = 183
    
    # number of trios studied in rauch data
    male.rauch = 32
    female.rauch = 19
    
    # number of trios studied in fromer data
    male.fromer = 317
    female.fromer = 306
    
    # number of trios studied in epi4k data
    male.epi4k = 156
    female.epi4k = 108
    
    # number of trios studied in zaidi data
    male.zaidi = 220
    female.zaidi = 142
    
    # sum up males and females across studies
    male = male.ddd + male.deligt + male.autism + male.rauch + male.fromer + male.epi4k + male.zaidi
    female = female.ddd + female.deligt + female.autism + female.rauch + female.fromer + female.epi4k + female.zaidi
    
    return(list(male = male, female = female))
}

open_datasets <- function(diagnosed) {
    # combine datasets listing de novo mutations into a single data frame
    # 
    # Args:
    #     diagnosed: list of IDs and sex for probands with diagnoses in the DDD
    # 
    # Returns:
    #     data frame containing HGNC, chrom, position, consequence, SNV or INDEL
    #     type, and study ID.
    
    ########### DDD dataset ###############
    # read in DNMs, genes, CQ and type from file
    our.data = read.delim(file.path(DE_NOVO_DIR, "DNG_Variants_20Feb2014_NonRed_Clean_NoTwins_NoClusters.txt"), header=T, colClasses = "character")
    # remove diagnosed patients, if maximising power
    our.data = our.data[-which(our.data$DECIPHER_ID %in% diagnosed$id), ]
    TYPE.index = which(our.data$snp_or_indel == "DENOVO-SNP")
    our.data$TYPE = "INDEL"
    our.data$TYPE[TYPE.index] = "SNV"
    
    # read in other datasets and calculate numbers of LoF and NS, SNVs and indels
    ########### ID datasets ###############
    rauch = read.delim(file.path(DE_NOVO_DIR, "rauch_v2.txt"), header=T, colClasses = "character")
    deligt = read.delim(file.path(DE_NOVO_DIR, "deligt_v2.txt"), header=T, colClasses = "character")
    
    ########### Epi4K dataset ############
    epi4k = read.delim(file.path(DE_NOVO_DIR, "epi4k_v2.txt"), header=T, colClasses = "character")
    TYPE.index = which(epi4k$Type == "snv")
    epi4k$TYPE = "INDEL"
    epi4k$TYPE[TYPE.index] = "SNV"
    
    ########### Autism dataset ###############
    autism = read.delim(file.path(DE_NOVO_DIR, "autism_v3_PJ.txt"), header=T, colClasses = "character")
    # select only de novos in probands
    autism = autism[which(autism$pheno == "Pro"), ]
    TYPE.index = which(abs(nchar(autism$ref.1) - nchar(autism$var)) == 0)
    autism$TYPE = "INDEL"
    autism$TYPE[TYPE.index] = "SNV"
    
    ########### Schizophrenia dataset ###############
    fromer = read.delim(file.path(DE_NOVO_DIR, "fromer_v2.txt"), header=T, colClasses = "character")
    TYPE.index = which(abs(nchar(fromer$Reference.allele) - nchar(fromer$Alternate.allele)) == 0)
    fromer$TYPE = "INDEL"
    fromer$TYPE[TYPE.index] = "SNV"
    
    ########### CHD dataset ###############
    # could only include syndromic DNMs
    zaidi = read.delim(file.path(DE_NOVO_DIR, "zaidi_VEP.txt"), header=T, colClasses = "character")
    # remove DNMs in controls
    zaidi = zaidi[-which(zaidi$Primary_Cardiac_Class == "Control"), ]
    TYPE.index = which(abs(nchar(zaidi$ref) - nchar(zaidi$alt)) != 0)
    TYPE.index = sort(unique(c(TYPE.index, which(zaidi$ref == "-"), which(zaidi$alt == "-"))))
    zaidi$TYPE = "SNV"
    zaidi$TYPE[TYPE.index] = "INDEL"
    
    # get complete lists of genes
    merged.HGNC = c(our.data$curated_HGNC, deligt$INFO.HGNC, rauch$INFO.HGNC, autism$INFO.HGNC, fromer$INFO.HGNC, epi4k$INFO.HGNC, zaidi$INFO.HGNC)
    merged.CQ = c(our.data$curated_CQ, deligt$INFO.CQ, rauch$INFO.CQ, autism$INFO.CQ, fromer$INFO.CQ, epi4k$INFO.CQ, zaidi$INFO.CQ)
    merged.POS = c(our.data$pos, deligt$POS, rauch$POS, autism$pos, fromer$pos, epi4k$pos, zaidi$pos)
    merged.CHROM = c(our.data$chr, deligt$CHROM, rauch$CHROM, autism$CHROM, fromer$chrom, epi4k$chrom, zaidi$chrom)
    merged.TYPE = c(our.data$TYPE, deligt$TYPE, rauch$TYPE, autism$TYPE, fromer$TYPE, epi4k$TYPE, zaidi$TYPE)
    merged.STUDY = c(rep("DDD", nrow(our.data)), rep("deligt", nrow(deligt)), rep("rauch", nrow(rauch)), rep("autism", nrow(autism)), rep("fromer", nrow(fromer)), rep("epi4k", nrow(epi4k)), rep("zaidi", nrow(zaidi)))
    data = data.frame(merged.HGNC, merged.CQ, merged.POS, merged.CHROM, merged.TYPE, merged.STUDY)
    names(data) = c("HGNC", "CQ", "POS", "CHROM", "TYPE", "STUDY")
    
    # write.table(data, file="/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/DNMs_280114/Meta_analysis_other_DNM_studies/Meta_DNMs_4Jeremy_100314.txt",  quote=F, row.names=F, sep="\t")
    
    return(data)
}

get_length_based_rates <- function(num.trios.male, num.trios.female) {
    # calculate mutation rates based on the gene length
    # 
    # Args:
    #     num.trios.male: number of trios with male offspring in the dataset
    #     num.trios.female: number of trios with female offspring in the dataset
    # 
    # Returns:
    #     list containing vectors of mutation rates
    
    # read-in length of coding sequence of each gene, from Ensembl biomart
    gene.info = read.delim(file.path(DATA_DIR, "CDS_LENGTH_B37_chr.txt"), header=T)
    cds.length = gene.info$CDS_LENGTH
    gene.info$gene = gene.info$ID
    chrX_index = which(gene.info$chr == "X")
    
    auto.transmissions = 2 * (num.trios.male + num.trios.female)
    female.transmissions = num.trios.male + num.trios.female
    male.transmissions = num.trios.female
    
    # get scaling factors using the alpha from the most recent SFHS (Scottish
    # Family Health Study) phased de novo data.
    alpha = 3.4
    male.chrx.scaling = 2 / (1 + (1 / alpha))
    female.chrx.scaling = 2 / (1 + alpha)
    
    # set mutation rates for snvs and indels
    snv.mut.rate = 1.5E-8 # higher than genome-wide mutation rate, due to higher GC ...
    # TODO: check that the following indel mutation rate is correct, since the
    # TODO: comment claims it should be 10% of the SNV mutation rate, but it
    # TODO: currently looks  like 30%
    indel.mut.rate = 0.53E-9 # ~10% of genome-wide SNV mutation rate, no reason to think higher in exome
    
    male.snv.mut.rate = snv.mut.rate * male.chrx.scaling
    female.snv.mut.rate = snv.mut.rate * female.chrx.scaling
    
    male.indel.mut.rate = indel.mut.rate * male.chrx.scaling
    female.indel.mut.rate = indel.mut.rate * female.chrx.scaling
    
    # specify proportion of coding mutations of different types
    snv.prop.lof = 0.0485 # from Daly
    snv.prop.missense = 0.6597 # from Daly
    
    indel.prop.lof = 0.9 # non-3n, from size distribution in neutral sequence
    indel.prop.missense = 0.1
    
    # calculate rates of missense and lof mutations, multiply by 2 for two 
    # transmissions per child and number of trios
    cds_rates = list()
    cds_rates$snv.missense.rate = cds.length * snv.mut.rate * snv.prop.missense * auto.transmissions
    cds_rates$snv.lof.rate = cds.length * snv.mut.rate * snv.prop.lof * auto.transmissions
    cds_rates$indel.missense.rate = cds.length * indel.mut.rate * indel.prop.missense * auto.transmissions
    cds_rates$indel.lof.rate = cds.length * indel.mut.rate * indel.prop.lof * auto.transmissions
    
    # correct non-PAR chrX genes for  fewer transmissions and lower rate 
    # (dependent on alpha) currently doing for all chrX genes, not just non-PAR 
    # genes
    x_snv_scaling = (male.transmissions * male.snv.mut.rate + female.transmissions * female.snv.mut.rate)
    x_indel_scaling = (male.transmissions * male.indel.mut.rate + female.transmissions * female.indel.mut.rate)
    cds_rates$snv.missense.rate[chrX_index] = cds.length[chrX_index] * snv.prop.missense * x_snv_scaling
    cds_rates$snv.lof.rate[chrX_index] = cds.length[chrX_index] * snv.prop.lof * x_snv_scaling
    cds_rates$indel.missense.rate[chrX_index] = cds.length[chrX_index] * indel.prop.missense * x_indel_scaling
    cds_rates$indel.lof.rate[chrX_index] = cds.length[chrX_index] * indel.prop.lof * x_indel_scaling
    
    # checked chrX rate is now slower with plot(gene.snv.missense.rate, cds.length)
    
    values = list(cds_rates = cds_rates, gene.info = gene.info)
    
    return(values)
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
    de_novo_counts$lof.snvs = rowSums(snvs[, grep(lof_regex, names(snvs))])
    de_novo_counts$missense.snvs = rowSums(snvs[, grep(missense_regex, names(snvs))])
    de_novo_counts$lof.indels = rowSums(indels[, grep(lof_regex, names(indels))])
    de_novo_counts$missense.indels = rowSums(indels[, grep(missense_regex, names(indels))])
    
    return(de_novo_counts)
}

get_p_values <- function(rates, de_novos, de_novo_counts, gene_data, num.tests) {
    # tests whether genes are enriched with de novo mutations
    # 
    # Args:
    #     rates: gene mutation rates per consequence type
    #     de_novos: data frame containing all the observed for all the genes
    #     de_novo_counts: data frame with tally of de novo mutations for each 
    #         of the mutation types.
    #     gene_data: data frame with gene info, sorted as per vectors in rates
    #     num.tests: number of tests performed (used for multiple correction).
    # 
    # Returns:
    #     data frame with gene info, mutation rates and P values from testing
    #     for enrichment.
    
    # set-up vectors to store gene-specific information only for observed genes
    observed = data.frame(matrix(NA, nrow = nrow(de_novo_counts), ncol = 8))
    names(observed) = c("chr", "coord", "snv.missense.rate", 
        "snv.lof.rate", "indel.missense.rate", "indel.lof.rate", "p.func", 
        "p.lof")

    # loop for each observed gene, test for functional variants and lof 
    # variants, using the gene mutation rates
    for (i in 1:nrow(de_novo_counts)) {
        gene = as.character(de_novo_counts$HGNC[i])
        gene_info = gene_data$gene == gene
        # catch if info not available for that gene
        if (any(gene_info)) {
            gene.index = which(gene_info)
            
            # get the mutation rates for the gene
            observed$snv.missense.rate[i] = rates$snv.missense.rate[gene.index]
            observed$snv.lof.rate[i] = rates$snv.lof.rate[gene.index]
            observed$indel.missense.rate[i] = rates$indel.missense.rate[gene.index]
            observed$indel.lof.rate[i] = rates$indel.lof.rate[gene.index]
            
            # figure out the chromosome and nucleotide position of the gene
            data.index = which(de_novos$HGNC == de_novo_counts$HGNC[i])[1]
            observed$chr[i] = de_novos$CHROM[data.index]
            observed$coord[i] = de_novos$POS[data.index]
            
            # count the observed de novos in each functional category
            lof_count = de_novo_counts$lof.snvs[i] + de_novo_counts$lof.indels[i]
            func_count = de_novo_counts$lof.snvs[i] + de_novo_counts$lof.indels[i] + de_novo_counts$missense.snvs[i] + de_novo_counts$missense.indels[i]
            
            # calculate the probability of observing said de novos, given the 
            # gene mutation rates
            observed$p.lof[i] = dpois(lof_count, lambda=observed$snv.lof.rate[i] + observed$indel.lof.rate[i])
            observed$p.func[i] = dpois(func_count, lambda=observed$snv.lof.rate[i] + observed$indel.lof.rate[i] + observed$snv.missense.rate[i] + observed$indel.missense.rate[i])
        }
        if (i %% 100 == 0) {
            print(paste(i, " out of ", nrow(de_novo_counts), " genes", sep = ""))
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
            ylim=c(0, max(-log10(length_p_vals), na.rm=T) + 1))
        
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

main <- function() {
    diagnosed = get_ddd_diagnosed()
    num = get_trio_counts(diagnosed)
    num.trios.male = num$male
    num.trios.female = num$female
    num.trios = num.trios.male + num.trios.female
    
    # open the de novos, and tally them by consequence and variant type
    de_novos = open_datasets(diagnosed)
    de_novo_counts = get_de_novo_counts(de_novos, CQ.LOF, CQ.MISSENSE)
    
    # get the length-based mutation rates for each gene
    length_values = get_length_based_rates(num.trios.male, num.trios.female)
    cds_rates = length_values$cds_rates
    gene.info = length_values$gene.info
    cds.length = length_values$cds.length
    
    # calculate the Daly mutation rates
    rates = get_mutation_rates(num.trios.male, num.trios.female)
    
    # calculate p values for each gene using the different mutation rates
    num.tests = 18500
    p_vals_length = get_p_values(cds_rates, de_novos, de_novo_counts, gene.info, num.tests)
    p_vals_daly = get_p_values(rates, de_novos, de_novo_counts, rates$daly, num.tests)
    
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


main()



