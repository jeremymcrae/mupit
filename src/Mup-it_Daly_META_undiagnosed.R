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

######### SETTINGS ################

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis"
DATA_DIR = file.path(CODE_DIR, "data")
SRC_DIR = file.path(CODE_DIR, "src")
DE_NOVO_DIR = file.path(DATA_DIR, "de_novo_datasets")
source(file.path(SRC_DIR, "mutation_rates_daly.R"))

CQ.LOF = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant")
CQ.NS = c( "missense_variant", "initiator_codon_variant", "stop_lost", "inframe_deletion", "inframe_insertion")

get_ddd_diagnosed <- function() {
    # read in samples that have been diagnosed, so as to remove from our data
    diagnoses = read.delim(file.path(DATA_DIR, "Diagnoses_1133.txt"), header=T)
    diagnosed.index = which(rowSums(diagnoses[, c(4:14)]) > 0)
    
    diagnosed = list()
    diagnosed$id = diagnoses$decipher_id[diagnosed.index]
    diagnosed$sex = diagnoses$Sex[diagnosed.index]
    
    return(diagnosed)
}

get_observed_values <- function(diagnosed) {
    # number of trios studied in our data
    male.ddd = 582 # trios with male offspring
    female.ddd = 548 # trios with female offspring
    
    # number of trios studied in deligt data
    male.deligt = 47 # trios with male offspring
    female.deligt = 53 # trios with female offspring
    
    # number of trios studied in autism data
    male.autism = 764 # trios with male offspring
    female.autism = 183 # trios with female offspring
    
    # number of trios studied in rauch data
    male.rauch = 32 # trios with male offspring
    female.rauch = 19 # trios with female offspring
    
    # number of trios studied in fromer data
    male.fromer = 317 # trios with male offspring
    female.fromer = 306 # trios with female offspring
    
    # number of trios studied in epi4k data
    male.epi4k = 156 # trios with male offspring
    female.epi4k = 108 # trios with female offspring
    
    # number of trios studied in zaidi data
    male.zaidi = 220 # trios with male offspring
    female.zaidi = 142 # trios with female offspring
    
    # remove diagnosed patients, if maximising power
    male.ddd = male.ddd - length(which(diagnosed$sex == "Male"))
    female.ddd = female.ddd - length(which(diagnosed$sex == "Female"))
    
    # sum up males and females across studies
    male = male.ddd +  male.deligt + male.autism + male.rauch + male.fromer + male.epi4k + male.zaidi
    female = female.ddd +  female.deligt + female.autism + female.rauch + female.fromer + female.epi4k + female.zaidi
    
    return(list(num.trios.male = male, num.trios.female = female))
}

open_datasets <- function(diagnosed) {
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
    ######### Calculate MUPit mutation rates ###########
    
    # read-in length of coding sequence of each gene, from Ensembl biomart
    gene.info = read.delim(file.path(DATA_DIR, "CDS_LENGTH_B37_chr.txt"), header=T)
    cds.length = gene.info$CDS_LENGTH
    gene.info$gene = gene.info$ID
    gene.index.chrx = which(gene.info$chr == "X")
    
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
    
    # calculate rates of missense and lof mutations, multiply by 2 for 2 transmissions/child and number of trios
    cds_rates = list()
    cds_rates$snv.missense.rate = cds.length * snv.mut.rate * snv.prop.missense * auto.transmissions
    cds_rates$snv.lof.rate = cds.length * snv.mut.rate * snv.prop.lof * auto.transmissions
    cds_rates$indel.missense.rate = cds.length * indel.mut.rate * indel.prop.missense * auto.transmissions
    cds_rates$indel.lof.rate = cds.length * indel.mut.rate * indel.prop.lof * auto.transmissions
    
    # correct non-PAR chrX genes for  fewer transmissions and lower rate (dependent on alpha)
    # currently doing for all chrX genes, not just non-PAR genes
    cds_rates$snv.missense.rate[gene.index.chrx] = cds.length[gene.index.chrx] * snv.prop.missense * (male.transmissions * male.snv.mut.rate + female.transmissions * female.snv.mut.rate)
    cds_rates$snv.lof.rate[gene.index.chrx] = cds.length[gene.index.chrx] * snv.prop.lof * (male.transmissions * male.snv.mut.rate + female.transmissions * female.snv.mut.rate)
    cds_rates$indel.missense.rate[gene.index.chrx] = cds.length[gene.index.chrx] * indel.prop.missense * (male.transmissions * male.indel.mut.rate + female.transmissions * female.indel.mut.rate)
    cds_rates$indel.lof.rate[gene.index.chrx] = cds.length[gene.index.chrx] * indel.prop.lof * (male.transmissions * male.indel.mut.rate + female.transmissions * female.indel.mut.rate)
    
    # checked chrX rate is now slower with plot(gene.snv.missense.rate, cds.length)
    
    values = list(cds_rates = cds_rates, gene.info = gene.info)
    
    return(values)
}

get_de_novo_counts <- function(data) {
    num.variants = nrow(data)

    LOF.variants = rep(0, num.variants)
    NS.variants = rep(0, num.variants)

    LOF.variants[which(data$CQ %in% CQ.LOF)] = 1
    NS.variants[which(data$CQ %in% CQ.NS)] = 1

    LOF.counts = table(data$HGNC, data$TYPE, LOF.variants)
    NS.counts = table(data$HGNC, data$TYPE, NS.variants)

    LOF.snvs = LOF.counts[, 2, 2]
    LOF.indels = LOF.counts[, 1, 2]
    NS.snvs = NS.counts[, 2, 2]
    NS.indels = NS.counts[, 1, 2]
    tot.NS.LOF = rowSums(cbind(LOF.snvs, NS.snvs, LOF.indels, NS.indels))

    # TODO: how are data and input.data different?
    input.data = data.frame(cbind(names(LOF.snvs), LOF.snvs, NS.snvs, LOF.indels, NS.indels, tot.NS.LOF), stringsAsFactors=F)
    
    return(input.data)
}

get_p_values <- function(rates, data, input.data, genes, num.tests) {
    
    # set-up vectors to store gene-specific information only for observed genes
    observed = data.frame(matrix(NA, nrow = nrow(input.data), ncol = 8))
    names(observed) = c("chr", "coord", "snv.missense.rate", 
        "snv.lof.rate", "indel.missense.rate", "indel.lof.rate", "p.DNM.func", 
        "p.DNM.lof")

    # loop for each observed gene, test for functional variants and lof 
    # variants, MUPit mutation rates
    for (i in 1:nrow(input.data)) {
        gene = input.data$V1[i]
        gene_info = genes$gene == gene
        # catch if info not available for that gene
        if (any(gene_info)) {
            gene.index = which(gene_info)
            
            observed$snv.missense.rate[i] = rates$snv.missense.rate[gene.index]
            observed$snv.lof.rate[i] = rates$snv.lof.rate[gene.index]
            observed$indel.missense.rate[i] = rates$indel.missense.rate[gene.index]
            observed$indel.lof.rate[i] = rates$indel.lof.rate[gene.index]
            
            data.index = which(as.character(data$HGNC) == as.character(input.data$V1[i]))[1]
            observed$chr[i] = as.character(data$CHROM[data.index])
            observed$coord[i] = as.character(data$POS[data.index])
            
            lof = as.numeric(input.data$LOF.snvs[i]) + as.numeric(input.data$LOF.indels[i])
            func = as.numeric(input.data$LOF.snvs[i]) + as.numeric(input.data$LOF.indels[i]) + as.numeric(input.data$NS.snvs[i]) + as.numeric(input.data$NS.indels[i])
            
            observed$p.DNM.lof[i] = dpois(lof, lambda=observed$snv.lof.rate[i] + observed$indel.lof.rate[i])
            observed$p.DNM.func[i] = dpois(func, lambda=observed$snv.lof.rate[i] + observed$indel.lof.rate[i] + observed$snv.missense.rate[i] + observed$indel.missense.rate[i])
        }
        if (i %% 100 == 0) {
            print(paste("Checked ", i, " out of ", nrow(input.data), " genes", sep = ""))
        }
    }
    
    # correct the P values for multiple testing by false discovery rate
    observed$fdr.lof = p.adjust(observed$p.DNM.lof, method="BH", n=num.tests)
    observed$fdr.func = p.adjust(observed$p.DNM.func, method="BH", n=num.tests)
    
    return(observed)
}

label_genes <- function(data, p_values, num.tests) {
    
    # label genes with fdr > threshold
    fdr = p.adjust(p_values, method="BH", n=num.tests)
    fdr.thresh = 0.05
    label.thresh = max(p_values[which(fdr < fdr.thresh)])
    thresh.index = which(p_values < label.thresh)
    data.thresh = data[thresh.index,]
    num.thresh = length(thresh.index)
    
    for (i in 1:num.thresh) {
        text(x=thresh.index[i], -log10(p_values[thresh.index[i]]), labels=data$V1[thresh.index[i]], pos=3, cex=0.5)
    }
}

plot_graphs <- function(data, num.tests) {
    
    # make Manhattan plots for LOF and Func variants separately
    Cairo(file="temp.pdf", type="pdf", width=20, height=20, units = "cm")
    
    data$chr[data$chr == "X"] = "23"
    data = data[order(as.numeric(as.character(data$chr)), as.numeric(as.character(data$coord))), ]
    col.chr.index = rep("lightblue3", nrow(data))
    odd.chr = c("1", "3", "5", "7", "9", "11", "13", "15", "17", "19", "21", "23")
    col.chr.index[which(data$chr %in% odd.chr)] = "darkblue"
    
    # CHECK COLOUR LABELS HERE
    plot(-log10(data$p.DNM.lof), col=col.chr.index, pch=19, cex=0.75, ylab="-log10(p)", xaxt="n", main="Loss-of-Function DNMs", xlab="genome position",  ylim=c(0, max(-log10(data$p.DNM.lof), na.rm=T) + 1))
    
    abline(h = -log10(0.05 / num.tests), col="red", lty=2) # Bonferroni correction
    points(-log10(data$daly.p.DNM.lof), col=col.chr.index, pch=1, cex=0.75)
    legend("topleft", legend=c("Length-based rates", "Daly group rates"), pch=c(19,1), col="darkblue", cex=0.75)
    label_genes(data, data$p.DNM.lof, num.tests)
    
    plot(-log10(data$p.DNM.func), col=col.chr.index, pch=19, cex=0.75, ylab="-log10(p)", xaxt="n", main="Functional DNMs", xlab="genome position", ylim=c(0, max(-log10(data$p.DNM.func), na.rm=T) + 1))
    
    abline(h=-log10(0.05 / num.tests), col="red", lty=2)
    points(-log10(data$daly.p.DNM.func), col=col.chr.index, pch=1, cex=0.75)
    legend("topleft", legend=c("Length-based rates", "Daly group rates"), pch=c(19,1), col="darkblue", cex=0.75)
    label_genes(data, data$p.DNM.func, num.tests)
    
    dev.off()
}

main <- function() {
    diagnosed = get_ddd_diagnosed()
    nums = get_observed_values(diagnosed)
    num.trios.male = nums$num.trios.male
    num.trios.female = nums$num.trios.female
    num.trios = num.trios.male + num.trios.female
    
    data = open_datasets(diagnosed)
    input.data = get_de_novo_counts(data)
    
    # get the length-based mutation rates for each gene
    length_values = get_length_based_rates(num.trios.male, num.trios.female)
    cds_rates = length_values$cds_rates
    gene.info = length_values$gene.info
    cds.length = length_values$cds.length
    
    # calculate the Daly mutation rates
    rates = get_mutation_rates(num.trios.male, num.trios.female)
    
    # calculate FDR using what looks like the number of genes (minus some untested genes?)
    num.tests = 18500
    p_vals_length = get_p_values(cds_rates, data, input.data, gene.info, num.tests)
    p_vals_daly = get_p_values(rates, data, input.data, rates$daly, num.tests)
    
    # write out results table
    input.data.plus = cbind(input.data, p_vals_length, p_vals_daly)
    names(input.data.plus) = c("V1", "LOF.snvs", "NS.snvs", "LOF.indels", 
        "NS.indels", "tot.NS.LOF", "chr", "coord",
        "snv.missense.rate", "snv.lof.rate", "indel.missense.rate", 
        "indel.lof.rate", "p.DNM.func", "p.DNM.lof", "fdr.lof", "fdr.func", 
     "chr", "coord", "daly.snv.missense.rate", "daly.snv.lof.rate", 
        "daly.indel.missense.rate", "daly.indel.lof.rate", "daly.p.DNM.func", "daly.p.DNM.lof", 
        "daly.fdr.lof", "daly.fdr.func")
    
    write.table(input.data.plus, file=file.path(DATA_DIR, "Mup-it_Daly_100414_META_undiagnosed.output.txt"), row.names=F, quote=F, sep="\t")
    plot_graphs(input.data.plus, num.tests)
}


main()



