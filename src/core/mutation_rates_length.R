# Calculate mutation rates based on the length of transcripts for a gene
# 
# This has been discontinued in favor of using rates derived from de novo 
# mutation rates, as provided by Mark Daly


CODE_DIR = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis"
DATA_DIR = file.path(CODE_DIR, "data")
SRC_DIR = file.path(CODE_DIR, "src")

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
    chrX = which(gene.info$chr == "X")
    
    auto.transmissions = 2 * (num.trios.male + num.trios.female)
    female.transmissions = num.trios.male + num.trios.female
    male.transmissions = num.trios.female
    
    # get scaling factors using the alpha from the most recent SFHS (Scottish
    # Family Health Study) phased de novo data.
    # NOTE: I'm note sure what alpha represents, or why the scaling factor is
    # generated like it is.
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
    rates = list()
    rates$snv.missense.rate = cds.length * snv.mut.rate * snv.prop.missense * auto.transmissions
    rates$snv.lof.rate = cds.length * snv.mut.rate * snv.prop.lof * auto.transmissions
    rates$indel.missense.rate = cds.length * indel.mut.rate * indel.prop.missense * auto.transmissions
    rates$indel.lof.rate = cds.length * indel.mut.rate * indel.prop.lof * auto.transmissions
    
    # correct non-PAR chrX genes for fewer transmissions and lower rate 
    # (dependent on alpha) currently doing for all chrX genes, not just non-PAR 
    # genes
    x_snv_scaling = (male.transmissions * male.snv.mut.rate + female.transmissions * female.snv.mut.rate)
    x_indel_scaling = (male.transmissions * male.indel.mut.rate + female.transmissions * female.indel.mut.rate)
    rates$snv.missense.rate[chrX] = cds.length[chrX] * snv.prop.missense * x_snv_scaling
    rates$snv.lof.rate[chrX] = cds.length[chrX] * snv.prop.lof * x_snv_scaling
    rates$indel.missense.rate[chrX] = cds.length[chrX] * indel.prop.missense * x_indel_scaling
    rates$indel.lof.rate[chrX] = cds.length[chrX] * indel.prop.lof * x_indel_scaling
    
    # checked chrX rate is now slower with plot(gene.snv.missense.rate, cds.length)
    
    values = list(cds_rates = rates, gene.info = gene.info)
    
    return(values)
}
