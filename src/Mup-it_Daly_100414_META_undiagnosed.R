
# program to calculate the significance of seeing N DNMs of a specific combination of functional types in a particular gene in M trios

# RATIONALE: use gene coding sequence to predict rate of DNMs in coding sequence for each gene of different functional classes, then estimate the probability of seeing the observed combination of different functional classes of DNMs assuming number of DNMs in each class is Poisson distributed

# initial implementation: use genome-wide mutation rate and scale by length of coding sequence, use genome-wide average of functional consequences of coding mutations from Kryukov et al 2007


# POTENTIAL FUTURE IMPROVEMENTS (highest priority first): 
#	adapt for chrX, lower mutation rate, and number of transmissions required. Assume male/female of proband = 1:1, unless known [DONE]
#	incorporate FDR estimation [DONE]
#	convert to analytical approach from permutation approach [DONE]
#	use actual number of exons to predict essential_splice_site mutations, 
#	use base composition of coding sequence to predict gene-specific mutation rate for each class of mutation, rather than just scaling the genome-average
#	use estimate of de novo mutation discovery power in a gene to better estimate gene-specific mutation rate
#	get CDS length for all genes, not just this subset, 455/478 in test data.use longest transcript if >1 [DONE]
#	calculate coding sequence length according to intersection of exome targeted regions and the union of all transcripts for a gene
# account for incomplete sensitivity for DNMs, especially indels
# look at clustering of de novos in genes with recurrent mutations, within protein space


# format: read in table of genes and the numbers of families with SNV DNMs within different functiona classes and indel DNMS in different functional classes, output same table with added column of probability of seeing that combination of DNMs.

#input data (validated SNV DNMs in TSV format, HGNC_ID, NUM_LOF, NUM_MISSENSE, NUM_LOF_INDEL, NUM_MISSENSE_INDEL):

# note: current test input file is numbers of mutations, not number of families, some families have >1 mutation in the same gene

######### SETTINGS ################


CQ.LOF<-c("stop_gained", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant")

CQ.NS<-c( "missense_variant", "initiator_codon_variant", "stop_lost", "inframe_deletion", "inframe_insertion")

#number of trios studied in our data
num.trios.male.ddd<-582 # trios with male offspring
num.trios.female.ddd<-548 # trios with female offspring

#number of trios studied in deligt data
num.trios.male.deligt<-47 # trios with male offspring
num.trios.female.deligt<-53 # trios with female offspring

#number of trios studied in autism data
num.trios.male.autism<-764 # trios with male offspring
num.trios.female.autism<-183 # trios with female offspring

#number of trios studied in rauch data
num.trios.male.rauch <-32 # trios with male offspring
num.trios.female.rauch <-19 # trios with female offspring

#number of trios studied in fromer data
num.trios.male.fromer<-317 # trios with male offspring
num.trios.female.fromer<-306 # trios with female offspring

#number of trios studied in epi4k data
num.trios.male.epi4k<-156 # trios with male offspring
num.trios.female.epi4k<-108 # trios with female offspring

#number of trios studied in zaidi data
num.trios.male.zaidi<-220 # trios with male offspring
num.trios.female.zaidi<-142 # trios with female offspring



# read in sample that have been diagnosed to remove from our data

diagnoses<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/DNMs_280114/Diagnoses_1133.txt", header=T)

diagnosed.index<-which(rowSums(diagnoses[,c(4:14)])>0)

diagnosed.id<-diagnoses$decipher_id[diagnosed.index]

diagnosed.sex<-diagnoses$Sex[diagnosed.index]

# remove diagnosed patients, if maximising power
num.trios.male.ddd<-num.trios.male.ddd - length(which(diagnosed.sex=="Male"))

num.trios.female.ddd<-num.trios.female.ddd - length(which(diagnosed.sex=="Female"))


# sum up males and females across studies

#WITH ZAIDI
num.trios.male<-num.trios.male.ddd +  num.trios.male.deligt + num.trios.male.autism + num.trios.male.rauch + num.trios.male.fromer + num.trios.male.epi4k + num.trios.male.zaidi

num.trios.female<-num.trios.female.ddd +  num.trios.female.deligt + num.trios.female.autism + num.trios.female.rauch + num.trios.female.fromer + num.trios.female.epi4k + num.trios.female.zaidi

#WITHOUT ZAIDI
#num.trios.male<-num.trios.male.ddd +  num.trios.male.deligt + num.trios.male.autism + num.trios.male.rauch + num.trios.male.fromer + num.trios.male.epi4k

#num.trios.female<-num.trios.female.ddd +  num.trios.female.deligt + num.trios.female.autism + num.trios.female.rauch + num.trios.female.fromer + num.trios.female.epi4k


num.trios<-num.trios.male + num.trios.female # to give total of 1130 (1133 minus 1 twin of 3 identical twins)

auto.transmissions<-2*num.trios

female.transmissions<-num.trios # for chrX analyses

male.transmissions<-num.trios.female # for chrX analyses


#specify ratio of male and female mutation rate to estimate sex-specific rate to allow rate estimation for chrX

alpha=3.4

male.chrx.scaling<-2/(1+1/alpha)

female.chrx.scaling<-2/(1+alpha)


#read-in length of coding sequence of each gene, from Ensembl biomart

gene.info<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/CDS_LENGTH_B37_chr.txt", header=T)

cds.length<-gene.info$CDS_LENGTH

gene.index.chrx<-which(gene.info$chr=="X")



##################################

# read in DNMs, genes, CQ and type from file

our.data<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/DNMs_280114/DNG_Variants_20Feb2014_NonRed_Clean_NoTwins_NoClusters.txt", header=T)

# remove diagnosed patients, if maximising power
our.data<-our.data[-which(our.data$DECIPHER_ID %in% diagnosed.id),]

TYPE.index<-which(our.data$snp_or_indel=="DENOVO-SNP")

TYPE<-rep("INDEL", length(our.data[,1]))

TYPE[TYPE.index]<-"SNV"

our.data<-cbind(our.data, TYPE)



# read in other datasets and calculate numbers of LoF and NS, SNVs and indels

########### ID datasets ###############

rauch<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/DNMs_280114/Meta_analysis_other_DNM_studies/rauch_v2.txt", header=T)

deligt<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/DNMs_280114/Meta_analysis_other_DNM_studies/deligt_v2.txt", header=T)


########### Epi4K dataset ############

epi4k<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/DNMs_280114/Meta_analysis_other_DNM_studies/epi4k_v2.txt", header=T)

TYPE.index<-which(epi4k$Type=="snv")

TYPE<-rep("INDEL", length(epi4k[,1]))

TYPE[TYPE.index]<-"SNV"

epi4k<-cbind(epi4k, TYPE)


########### Autism dataset ###############

autism<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/DNMs_280114/Meta_analysis_other_DNM_studies/autism_v3_PJ.txt", header=T)

# select only de novos in probands
autism<-autism[which(autism$pheno=="Pro"),]

TYPE.index<-which(abs(nchar(as.character(autism$ref.1))-nchar(as.character(autism$var)))==0)

TYPE<-rep("INDEL", length(autism[,1]))

TYPE[TYPE.index]<-"SNV"

autism<-cbind(autism, TYPE)


########### Schizophrenia dataset ###############

fromer<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/DNMs_280114/Meta_analysis_other_DNM_studies/fromer_v2.txt", header=T)

TYPE.index<-which(abs(nchar(as.character(fromer$Reference.allele))-nchar(as.character(fromer$Alternate.allele)))==0)

TYPE<-rep("INDEL", length(fromer[,1]))

TYPE[TYPE.index]<-"SNV"

fromer<-cbind(fromer, TYPE)


########### CHD dataset ###############

# could only include syndromic DNMs

zaidi<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/DNMs_280114/Meta_analysis_other_DNM_studies/zaidi_VEP.txt", header=T)

TYPE<-rep("SNV", length(zaidi[,1]))

TYPE.index<-which(abs(nchar(as.character(zaidi$ref))-nchar(as.character(zaidi$alt)))!=0)

TYPE.index<-sort(unique(c(TYPE.index, which(zaidi$ref=="-"), which(zaidi$alt=="-"))))

TYPE[TYPE.index]<-"INDEL"

zaidi<-cbind(zaidi, TYPE)


# remove DNMs in controls

zaidi<-zaidi[-which(zaidi$Primary_Cardiac_Class=="Control"),]



# get complete lists of genes
#WITH ZAIDI

merged.HGNC<-c(as.character(our.data$curated_HGNC), as.character(deligt$INFO.HGNC), as.character(rauch$INFO.HGNC), as.character(autism$INFO.HGNC), as.character(fromer$INFO.HGNC), as.character(epi4k$INFO.HGNC), as.character(zaidi$INFO.HGNC))

merged.CQ<-c(as.character(our.data$curated_CQ), as.character(deligt$INFO.CQ), as.character(rauch$INFO.CQ), as.character(autism$INFO.CQ), as.character(fromer$INFO.CQ), as.character(epi4k$INFO.CQ), as.character(zaidi$INFO.CQ))

merged.POS<-c(as.character(our.data$pos), as.character(deligt$POS), as.character(rauch$POS), as.character(autism$pos), as.character(fromer$pos), as.character(epi4k$pos), as.character(zaidi$pos))

merged.CHROM<-c(as.character(our.data$chr), as.character(deligt$CHROM), as.character(rauch$CHROM), as.character(autism$CHROM), as.character(fromer$chrom), as.character(epi4k$chrom), as.character(zaidi$chrom))

merged.TYPE<-c(as.character(our.data$TYPE), as.character(deligt$TYPE), as.character(rauch$TYPE), as.character(autism$TYPE), as.character(fromer$TYPE), as.character(epi4k$TYPE), as.character(zaidi$TYPE))

merged.STUDY<-c(rep("DDD", length(our.data[,1])), rep("deligt", length(deligt[,1])), rep("rauch", length(rauch[,1])), rep("autism", length(autism[,1])), rep("fromer", length(fromer[,1])), rep("epi4k", length(epi4k[,1])), rep("zaidi", length(zaidi[,1])))

raw.data<-data.frame(merged.HGNC, merged.CQ, merged.POS, merged.CHROM, merged.TYPE, merged.STUDY) 

names(raw.data)<-c("HGNC", "CQ", "POS", "CHROM", "TYPE", "STUDY")

#WITHOUT ZAIDI

#merged.HGNC<-c(as.character(our.data$curated_HGNC), as.character(deligt$INFO.HGNC), as.character(rauch$INFO.HGNC), as.character(autism$INFO.HGNC), as.character(fromer$INFO.HGNC), as.character(epi4k$INFO.HGNC))

#merged.CQ<-c(as.character(our.data$curated_CQ), as.character(deligt$INFO.CQ), as.character(rauch$INFO.CQ), as.character(autism$INFO.CQ), as.character(fromer$INFO.CQ), as.character(epi4k$INFO.CQ))

#merged.POS<-c(as.character(our.data$pos), as.character(deligt$POS), as.character(rauch$POS), as.character(autism$POS), as.character(fromer$POS), as.character(epi4k$POS))

#merged.CHROM<-c(as.character(our.data$chr), as.character(deligt$CHROM), as.character(rauch$CHROM), as.character(autism$CHROM), as.character(fromer$CHROM), as.character(epi4k$CHROM))

#merged.TYPE<-c(as.character(our.data$TYPE), as.character(deligt$TYPE), as.character(rauch$TYPE), as.character(autism$TYPE), as.character(fromer$TYPE), as.character(epi4k$TYPE))

#raw.data<-data.frame(merged.HGNC, merged.CQ, merged.POS, merged.CHROM, merged.TYPE, merged.STUDY) 

#names(raw.data)<-c("HGNC", "CQ", "POS", "CHROM", "TYPE", "STUDY")


#write.table(raw.data, file="/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/DNMs_280114/Meta_analysis_other_DNM_studies/Meta_DNMs_4Jeremy_100314.txt",  quote=F, row.names=F, sep="\t")

############


num.variants<-length(raw.data[,1])

LOF.variants<-rep(0,num.variants)
NS.variants<-rep(0,num.variants)

LOF.variants[which(raw.data$CQ %in% CQ.LOF)]<-1
NS.variants[which(raw.data$CQ %in% CQ.NS)]<-1

LOF.counts<-table(as.character(raw.data$HGNC), raw.data$TYPE, LOF.variants)
NS.counts<-table(as.character(raw.data$HGNC), raw.data$TYPE, NS.variants)

LOF.snvs<-LOF.counts[,2,2]

LOF.indels<-LOF.counts[,1,2]

NS.snvs<-NS.counts[,2,2]

NS.indels<-NS.counts[,1,2]

tot.NS.LOF<-rowSums(cbind(LOF.snvs, NS.snvs, LOF.indels, NS.indels))

input.data<-data.frame(cbind(dimnames(data.frame(LOF.snvs))[[1]] , LOF.snvs, NS.snvs, LOF.indels, NS.indels, tot.NS.LOF), stringsAsFactors=F)

num.genes<-length(input.data[,1])





######### Calculate MUPit mutation rates ###########

#set mutation rates for snvs and indels

snv.mut.rate<-1.5E-8 # higher than genome-wide mutation rate, due to higher GC ...

indel.mut.rate<-0.53E-9 # ~10% of genome-wide SNV mutation rate, no reason to think higher in exome

male.snv.mut.rate<-snv.mut.rate*male.chrx.scaling
female.snv.mut.rate<-snv.mut.rate*female.chrx.scaling

male.indel.mut.rate<-indel.mut.rate*male.chrx.scaling
female.indel.mut.rate<-indel.mut.rate*female.chrx.scaling

# specify proportion of coding mutations of different types

snv.prop.lof<-0.0485 # from Daly
snv.prop.missense<-0.6597 # from Daly

indel.prop.lof<-0.9 # non-3n, from size distribution in neutral sequence
indel.prop.missense<-0.1


#calculate rates of missense and lof mutations, multiply by 2 for 2 transmissions/child and number of trios

gene.snv.missense.rate<-cds.length*snv.mut.rate*snv.prop.missense*auto.transmissions

gene.snv.lof.rate<-cds.length*snv.mut.rate*snv.prop.lof*auto.transmissions

gene.indel.missense.rate<-cds.length*indel.mut.rate*indel.prop.missense*auto.transmissions

gene.indel.lof.rate<-cds.length*indel.mut.rate*indel.prop.lof*auto.transmissions


#correct non-PAR chrX genes for  fewer transmissions and lower rate (dependent on alpha)
# currently doing for all chrX genes, not just non-PAR genes


gene.snv.missense.rate[gene.index.chrx]<-cds.length[gene.index.chrx]*snv.prop.missense*(male.transmissions*male.snv.mut.rate + female.transmissions*female.snv.mut.rate)

gene.snv.lof.rate[gene.index.chrx]<-cds.length[gene.index.chrx]*snv.prop.lof*(male.transmissions*male.snv.mut.rate + female.transmissions*female.snv.mut.rate)

gene.indel.missense.rate[gene.index.chrx]<-cds.length[gene.index.chrx]*indel.prop.missense*(male.transmissions*male.indel.mut.rate + female.transmissions*female.indel.mut.rate)

gene.indel.lof.rate[gene.index.chrx]<-cds.length[gene.index.chrx]*indel.prop.lof*(male.transmissions*male.indel.mut.rate + female.transmissions*female.indel.mut.rate)

# checked chrX rate is now slower with plot(gene.snv.missense.rate, cds.length)


##### Calculate Daly mutation rates #######

daly<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/Compare_Daly_MUPIT/fixed_mut_prob_fs_adjdepdiv.txt", header=T)

daly<-merge(daly, gene.info, by.x=2, by.y=4, all.x=T) #add chromosome annotation

daly.gene.snv.missense.rate<-(10^daly$mis+10^daly$rdt)*auto.transmissions

daly.gene.snv.lof.rate<-(10^daly$non+10^daly$css)*auto.transmissions

daly.gene.indel.missense.rate<-(10^daly$frameshift)/9*auto.transmissions

daly.gene.indel.lof.rate<-(10^daly$frameshift)*auto.transmissions


# could scale to take account of longer transcript

#adapt indel rates to take account of lower rate estimate from validated de novos

valid.nonsense<-102
valid.frameshift<-95
daly.scaling<-1.25 # ratio of frameshift to nonsense

daly.gene.indel.missense.rate<-daly.gene.indel.missense.rate/1.25*valid.frameshift/valid.nonsense

daly.gene.indel.lof.rate<-daly.gene.indel.lof.rate/1.25*valid.frameshift/valid.nonsense


# catch cases where there is no rdt or css mutation rate, resulting in an NA for composite rates

daly.gene.snv.missense.rate[is.na(daly.gene.snv.missense.rate)]<-(10^daly$mis[is.na(daly.gene.snv.missense.rate)])*auto.transmissions

daly.gene.snv.lof.rate[is.na(daly.gene.snv.lof.rate)]<-(10^daly$non[is.na(daly.gene.snv.lof.rate)])*auto.transmissions


#correct non-PAR chrX genes for  fewer transmissions and lower rate (dependent on alpha)


daly.gene.index.chrx<-which(daly$chr=="X")


daly.gene.snv.missense.rate[daly.gene.index.chrx]<-daly.gene.snv.missense.rate[daly.gene.index.chrx]/auto.transmissions*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

#daly.gene.snv.missense.rate[daly.gene.index.chrx]<-(10^daly$mis[daly.gene.index.chrx]+10^daly$rdt[daly.gene.index.chrx])*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

daly.gene.snv.lof.rate[daly.gene.index.chrx]<-daly.gene.snv.lof.rate[daly.gene.index.chrx]/auto.transmissions*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

#daly.gene.snv.lof.rate[daly.gene.index.chrx]<-(10^daly$non[daly.gene.index.chrx]+10^daly$css[daly.gene.index.chrx])*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

daly.gene.indel.missense.rate[daly.gene.index.chrx]<-daly.gene.indel.missense.rate[daly.gene.index.chrx]/auto.transmissions*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

#daly.gene.indel.missense.rate[daly.gene.index.chrx]<-(10^daly$frameshift[daly.gene.index.chrx])/9*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

daly.gene.indel.lof.rate[daly.gene.index.chrx]<-daly.gene.indel.lof.rate[daly.gene.index.chrx]/auto.transmissions*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

#daly.gene.indel.lof.rate[daly.gene.index.chrx]<-(10^daly$frameshift[daly.gene.index.chrx])*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)



# catch cases where there is no rdt or css mutation rate, resulting in an NA for composite rates

#if (is.na(daly.gene.snv.missense.rate[daly.gene.index.chrx])


#set-up vectors to store gene-specific information only for observed genes

observed.cds.length<-rep(0,num.genes)
observed.chr<-rep(0,num.genes)
observed.coord<-rep(0,num.genes)

observed.snv.missense.rate<-rep(0,num.genes)
observed.snv.lof.rate<-rep(0,num.genes)
observed.indel.missense.rate<-rep(0,num.genes)
observed.indel.lof.rate<-rep(0,num.genes)

p.DNM.func<-rep(0,num.genes)
p.DNM.lof<-rep(0,num.genes)

daly.observed.snv.missense.rate<-rep(0,num.genes)
daly.observed.snv.lof.rate<-rep(0,num.genes)
daly.observed.indel.missense.rate<-rep(0,num.genes)
daly.observed.indel.lof.rate<-rep(0,num.genes)

daly.p.DNM.func<-rep(0,num.genes)
daly.p.DNM.lof<-rep(0,num.genes)


#loop for each observed gene, test for functional variants and lof variants, MUPit mutation rates

for (i in seq(1,num.genes)) {
  
  #catch if info not available for that gene
  
  if (length(which(as.character(gene.info$ID) == as.character(input.data$V1[i])))==0) {
    
    observed.cds.length[i]<-NA
    observed.chr[i]<-NA
    observed.coord[i]<-NA
    
    observed.snv.missense.rate[i]<-NA
    observed.snv.lof.rate[i]<-NA
    observed.indel.missense.rate[i]<-NA
    observed.indel.lof.rate[i]<-NA
    
    p.DNM.func[i]<-NA
    p.DNM.lof[i]<-NA
    
  } else {
    
    gene.index<-which(as.character(gene.info$ID) == as.character(input.data$V1[i]))
    
    
    observed.cds.length[i]<-cds.length[gene.index]
    
    observed.snv.missense.rate[i]<-gene.snv.missense.rate[gene.index]
    observed.snv.lof.rate[i]<-gene.snv.lof.rate[gene.index]
    observed.indel.missense.rate[i]<-gene.indel.missense.rate[gene.index]
    observed.indel.lof.rate[i]<-gene.indel.lof.rate[gene.index]
    
    raw.data.index<-which(as.character(raw.data$HGNC) == as.character(input.data$V1[i]))[1]
    
    observed.chr[i]<-as.character(raw.data$CHROM[raw.data.index])
    observed.coord[i]<-as.character(raw.data$POS[raw.data.index])
    
    observed.lof<-as.numeric(input.data$LOF.snvs[i])+as.numeric(input.data$LOF.indels[i])
    
    # observed.lof<-as.numeric(input.data$LOF.snvs[i])
    
    observed.func<-as.numeric(input.data$LOF.snvs[i]) + as.numeric(input.data$LOF.indels[i]) + as.numeric(input.data$NS.snvs[i]) + as.numeric(input.data$NS.indels[i])
    
    # observed.func<-as.numeric(input.data$LOF.snvs[i]) + as.numeric(input.data$NS.snvs[i])		
    
    p.DNM.lof[i]<-dpois(observed.lof, lambda=observed.snv.lof.rate[i] + observed.indel.lof.rate[i])
    
    # 	p.DNM.lof[i]<-dpois(observed.lof, lambda=observed.snv.lof.rate[i])
    
    p.DNM.func[i]<-dpois(observed.func, lambda=observed.snv.lof.rate[i] + observed.indel.lof.rate[i] + observed.snv.missense.rate[i] + observed.indel.missense.rate[i])
    
    #	p.DNM.func[i]<-dpois(observed.func, lambda=observed.snv.lof.rate[i] + observed.snv.missense.rate[i])
    
  }
  print(i)
}


#loop for each observed gene, test for functional variants and lof variants, daly mutation rates

for (i in seq(1,num.genes)) {
  
  #catch if info not available for that gene
  
  if (length(which(as.character(daly$gene) == as.character(input.data$V1[i])))==0) {
   
    
    daly.observed.snv.missense.rate[i]<-NA
    daly.observed.snv.lof.rate[i]<-NA
    daly.observed.indel.missense.rate[i]<-NA
    daly.observed.indel.lof.rate[i]<-NA
    
    daly.p.DNM.func[i]<-NA
    daly.p.DNM.lof[i]<-NA
    
  } else {
    
    gene.index<-which(as.character(daly$gene) == as.character(input.data$V1[i]))
        
    daly.observed.snv.missense.rate[i]<-daly.gene.snv.missense.rate[gene.index]
    daly.observed.snv.lof.rate[i]<-daly.gene.snv.lof.rate[gene.index]
    daly.observed.indel.missense.rate[i]<-daly.gene.indel.missense.rate[gene.index]
    daly.observed.indel.lof.rate[i]<-daly.gene.indel.lof.rate[gene.index]
    
    observed.lof<-as.numeric(input.data$LOF.snvs[i])+as.numeric(input.data$LOF.indels[i]) 
           
    observed.func<-as.numeric(input.data$LOF.snvs[i]) + as.numeric(input.data$LOF.indels[i]) + as.numeric(input.data$NS.snvs[i]) + as.numeric(input.data$NS.indels[i])
    
        
    daly.p.DNM.lof[i]<-dpois(observed.lof, lambda=daly.observed.snv.lof.rate[i] + daly.observed.indel.lof.rate[i])
        
    daly.p.DNM.func[i]<-dpois(observed.func, lambda=daly.observed.snv.lof.rate[i] + daly.observed.indel.lof.rate[i] + daly.observed.snv.missense.rate[i] + daly.observed.indel.missense.rate[i])
    
    #	p.DNM.func[i]<-dpois(observed.func, lambda=observed.snv.lof.rate[i] + observed.snv.missense.rate[i])
    
  }
  print(i)
}





# calculate FDR

num.tests=18500

fdr.lof<-p.adjust(p.DNM.lof, method="BH", n=num.tests)

fdr.func<-p.adjust(p.DNM.func, method="BH", n=num.tests)

daly.fdr.lof<-p.adjust(daly.p.DNM.lof, method="BH", n=num.tests)

daly.fdr.func<-p.adjust(daly.p.DNM.func, method="BH", n=num.tests)


#write out results table


input.data.plus<-cbind(input.data, observed.chr, observed.coord, observed.cds.length, observed.snv.missense.rate, observed.snv.lof.rate, observed.indel.missense.rate, observed.indel.lof.rate, p.DNM.lof, p.DNM.func, fdr.lof, fdr.func, daly.observed.snv.missense.rate, daly.observed.snv.lof.rate, daly.observed.indel.missense.rate, daly.observed.indel.lof.rate, daly.p.DNM.lof, daly.p.DNM.func, daly.fdr.lof, daly.fdr.func )

# input.data.plus<-cbind(input.data, observed.chr, observed.coord, observed.cds.length, observed.snv.missense.rate, observed.snv.lof.rate, p.DNM.lof, p.DNM.func)

write.table(input.data.plus, file="/Volumes/DDD_meh/Analysis/Exome/For_1133_trio_ms/Meta_DD_undiagnosed/Mup-it_Daly_100414_META_undiagnosed.output.txt", row.names=F, quote=F, sep="\t")

#found.gene.index<-which(gene.info$HGNC %in% input.data$HGNC)


# make Manhattan plots for LOF and Func variants separately

quartz(width=15, height=7)

order.index<-order(as.numeric(as.character(input.data.plus$observed.chr))*1e09+as.numeric(as.character(input.data.plus$observed.coord)))

ordered.data<-input.data.plus[order.index,]

col.chr.index<-as.character(ordered.data$observed.chr)

odd.chr<-c("1","3","5","7","9","11","13","15","17","19","21", "X")

even.chr<-c("2", "4", "6", "8", "10", "12", "14", "16", "18", "20", "22")

col.chr.index[which(col.chr.index %in% odd.chr)]<-"darkblue"

col.chr.index[which(col.chr.index %in% even.chr)]<-"lightblue3"

# CHECK COLOUR LABELS HERE

plot(-log10(ordered.data$p.DNM.lof), col=col.chr.index, pch=19, cex=0.75, ylab="-log10(p)", xaxt="n", main="Loss-of-Function DNMs", xlab="genome position",  ylim=c(0,max(-log10(ordered.data$p.DNM.lof), na.rm=T)+1))

abline(h=-log10(0.05/num.tests), col="red", lty=2) # Bonferroni correction

points(-log10(ordered.data$daly.p.DNM.lof), col=col.chr.index, pch=1, cex=0.75)

legend("topleft", legend=c("Length-based rates", "Daly group rates"), pch=c(19,1), col="darkblue", cex=0.75)

# label genes with fdr > threshold

fdr<-p.adjust(ordered.data$p.DNM.lof, method="BH", n=num.tests)

fdr.thresh<-0.05

label.thresh<-max(ordered.data$p.DNM.lof[which(fdr<fdr.thresh)])

thresh.index<-which(ordered.data$p.DNM.lof<label.thresh)

ordered.data.thresh<-ordered.data[thresh.index,]

num.thresh<-length(thresh.index)

for (i in seq(1,num.thresh)) {
  
  text(x=thresh.index[i], -log10(ordered.data$p.DNM.lof[thresh.index[i]]), labels=ordered.data$V1[thresh.index[i]], pos=3, cex=0.5)
  
}

quartz(width=15, height=7)


plot(-log10(ordered.data$p.DNM.func), col=col.chr.index, pch=19, cex=0.75, ylab="-log10(p)", xaxt="n", main="Functional DNMs", xlab="genome position", ylim=c(0,max(-log10(ordered.data$p.DNM.func), na.rm=T)+1))

abline(h=-log10(0.05/num.tests), col="red", lty=2)

points(-log10(ordered.data$daly.p.DNM.func), col=col.chr.index, pch=1, cex=0.75)

legend("topleft", legend=c("Length-based rates", "Daly group rates"), pch=c(19,1), col="darkblue", cex=0.75)

# label genes with fdr > threshold

fdr<-p.adjust(ordered.data$p.DNM.func, method="BH", n=num.tests)

fdr.thresh<-0.05

label.thresh<-max(ordered.data$p.DNM.func[which(fdr<fdr.thresh)])

thresh.index<-which(ordered.data$p.DNM.func<label.thresh)

ordered.data.thresh<-ordered.data[thresh.index,]

num.thresh<-length(thresh.index)

for (i in seq(1,num.thresh)) {
  
  text(x=thresh.index[i], -log10(ordered.data$p.DNM.func[thresh.index[i]]), labels=ordered.data$V1[thresh.index[i]], pos=3, cex=0.5)
  
}


# fishers method for combining p values

# fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)


