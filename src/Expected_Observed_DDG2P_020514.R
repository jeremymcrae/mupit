# R code to estimate expected and observed numbers of mutations in DDG2P genes

observed.ddg2p.ar.lof<-6
observed.ddg2p.ar.missense<-43
observed.ddg2p.ar.silent<-12
observed.ddg2p.ar.inframe<-0
observed.ddg2p.ar.frameshift<-0

observed.ddg2p.nonar.lof<-45
observed.ddg2p.nonar.missense<-97
observed.ddg2p.nonar.silent<-9
observed.ddg2p.nonar.inframe<-1
observed.ddg2p.nonar.frameshift<-38

observed.nonddg2p.lof<-86
observed.nonddg2p.missense<-724
observed.nonddg2p.silent<-239
observed.nonddg2p.inframe<-9
observed.nonddg2p.frameshift<-57

# read in DDG2P, only Confirmed and Probable, split into non-AR and AR

ddg2p<-read.delim("~/Desktop/DDD/DDG2P/DDG2P_Confirmed_Probable_20131107.txt", header=T)

non.AR.CQ<-c("Both", "Hemizygous", "Monoallelic", "X-linked dominant")

chrx.CQ<-c("Hemizygous", "X-linked dominant")

non.AR.auto.CQ<-c("Both", "Monoallelic")

AR.CQ<-c("Biallelic")

ddg2p.non.AR<-ddg2p[which(ddg2p$mode %in% non.AR.CQ),] #select only non AR genes

ddg2p.AR<-ddg2p[which(ddg2p$mode %in% AR.CQ),] #selection only AR genes

ddg2p.non.AR.uniq<-unique(ddg2p.non.AR$gene) # get non-redundant list of HGNC names for non.AR 

ddg2p.AR.uniq<-unique(ddg2p.AR$gene) # get non-redundant list of HGNC names for AR

ddg2p.non.AR.chrx.uniq<-unique(ddg2p$gene[which(ddg2p$mode %in% chrx.CQ)])

ddg2p.non.AR.auto.uniq<-unique(ddg2p$gene[which(ddg2p$mode %in% non.AR.auto.CQ)])



gene.info<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/CDS_LENGTH_B37_chr.txt", header=T)


#number of trios studied in our data
num.trios.male.ddd<-582 # trios with male offspring
num.trios.female.ddd<-548 # trios with female offspring

num.trios.male<-num.trios.male.ddd
num.trios.female<-num.trios.female.ddd

num.trios<-num.trios.male + num.trios.female # to give total of 1130 (1133 minus 1 twin of 3 identical twins)

auto.transmissions<-2*num.trios

female.transmissions<-num.trios # for chrX analyses

male.transmissions<-num.trios.female # for chrX analyses

#specify ratio of male and female mutation rate to estimate sex-specific rate to allow rate estimation for chrX

alpha=3.4 # from most recent SFHS phased de novo data

male.chrx.scaling<-2/(1+1/alpha)

female.chrx.scaling<-2/(1+alpha)



##### Calculate Daly mutation rates #######

daly<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/Compare_Daly_MUPIT/fixed_mut_prob_fs_adjdepdiv.txt", header=T)

daly<-merge(daly, gene.info, by.x=2, by.y=4, all.x=T) #add chromosome annotation

daly.gene.snv.silent.rate<-(10^daly$syn)*auto.transmissions

daly.gene.snv.missense.rate<-(10^daly$mis+10^daly$rdt)*auto.transmissions

daly.gene.snv.lof.rate<-(10^daly$non+10^daly$css)*auto.transmissions

daly.gene.indel.missense.rate<-(10^daly$frameshift)/9*auto.transmissions

daly.gene.indel.lof.rate<-(10^daly$frameshift)*auto.transmissions

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

daly.gene.snv.silent.rate[daly.gene.index.chrx]<-daly.gene.snv.silent.rate[daly.gene.index.chrx]/auto.transmissions*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

daly.gene.snv.missense.rate[daly.gene.index.chrx]<-daly.gene.snv.missense.rate[daly.gene.index.chrx]/auto.transmissions*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

daly.gene.snv.lof.rate[daly.gene.index.chrx]<-daly.gene.snv.lof.rate[daly.gene.index.chrx]/auto.transmissions*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

daly.gene.indel.missense.rate[daly.gene.index.chrx]<-daly.gene.indel.missense.rate[daly.gene.index.chrx]/auto.transmissions*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)

daly.gene.indel.lof.rate[daly.gene.index.chrx]<-daly.gene.indel.lof.rate[daly.gene.index.chrx]/auto.transmissions*(male.transmissions*male.chrx.scaling + female.transmissions*female.chrx.scaling)


#identify which genes in DDG2P

daly.non.AR.auto.index<-which(daly$gene %in% ddg2p.non.AR.auto.uniq)

daly.non.AR.chrx.index<-which(daly$gene %in% ddg2p.non.AR.chrx.uniq)

daly.non.AR.index<-unique(c(daly.non.AR.auto.index, daly.non.AR.chrx.index))

daly.AR.index<-which(daly$gene %in% ddg2p.AR.uniq)

daly.ddg2p.index<-unique(c(daly.non.AR.index, daly.AR.index))



# calculated expected numbers

expected.ddg2p.ar.lof<-sum(daly.gene.snv.lof.rate[daly.AR.index])
expected.ddg2p.ar.missense<-sum(daly.gene.snv.missense.rate[daly.AR.index])
expected.ddg2p.ar.silent<-sum(daly.gene.snv.silent.rate[daly.AR.index])
expected.ddg2p.ar.inframe<-sum(daly.gene.indel.missense.rate[daly.AR.index])
expected.ddg2p.ar.frameshift<-sum(daly.gene.indel.lof.rate[daly.AR.index])

expected.ddg2p.nonar.lof<-sum(daly.gene.snv.lof.rate[daly.non.AR.index])
expected.ddg2p.nonar.missense<-sum(daly.gene.snv.missense.rate[daly.non.AR.index])
expected.ddg2p.nonar.silent<-sum(daly.gene.snv.silent.rate[daly.non.AR.index])
expected.ddg2p.nonar.inframe<-sum(daly.gene.indel.missense.rate[daly.non.AR.index])
expected.ddg2p.nonar.frameshift<-sum(daly.gene.indel.lof.rate[daly.non.AR.index])

expected.nonddg2p.lof<-sum(daly.gene.snv.lof.rate[-daly.ddg2p.index])
expected.nonddg2p.missense<-sum(daly.gene.snv.missense.rate[-daly.ddg2p.index])
expected.nonddg2p.silent<-sum(daly.gene.snv.silent.rate[-daly.ddg2p.index])
expected.nonddg2p.inframe<-sum(daly.gene.indel.missense.rate[-daly.ddg2p.index])
expected.nonddg2p.frameshift<-sum(daly.gene.indel.lof.rate[-daly.ddg2p.index])

observed.ddg2p.ar.lof.p<-ppois(observed.ddg2p.ar.lof-1, expected.ddg2p.ar.lof, F)
observed.ddg2p.ar.missense.p<-ppois(observed.ddg2p.ar.missense-1, expected.ddg2p.ar.missense, F)
observed.ddg2p.ar.silent.p<-ppois(observed.ddg2p.ar.silent-1, expected.ddg2p.ar.silent, F)
observed.ddg2p.ar.inframe.p<-ppois(observed.ddg2p.ar.inframe-1, expected.ddg2p.ar.inframe, F)
observed.ddg2p.ar.frameshift.p<-ppois(observed.ddg2p.ar.frameshift-1, expected.ddg2p.ar.frameshift, F)

observed.ddg2p.nonar.lof.p<-ppois(observed.ddg2p.nonar.lof-1, expected.ddg2p.nonar.lof, F)
observed.ddg2p.nonar.missense.p<-ppois(observed.ddg2p.nonar.missense-1, expected.ddg2p.nonar.missense, F)
observed.ddg2p.nonar.silent.p<-ppois(observed.ddg2p.nonar.silent-1, expected.ddg2p.nonar.silent, F)
observed.ddg2p.nonar.inframe.p<-ppois(observed.ddg2p.nonar.inframe-1, expected.ddg2p.nonar.inframe, F)
observed.ddg2p.nonar.frameshift.p<-ppois(observed.ddg2p.nonar.frameshift-1, expected.ddg2p.nonar.frameshift, F)

observed.nonddg2p.lof.p<-ppois(observed.nonddg2p.lof-1, expected.nonddg2p.lof, F)
observed.nonddg2p.missense.p<-ppois(observed.nonddg2p.missense-1, expected.nonddg2p.missense, F)
observed.nonddg2p.silent.p<-ppois(observed.nonddg2p.silent-1, expected.nonddg2p.silent, F)
observed.nonddg2p.inframe.p<-ppois(observed.nonddg2p.inframe-1, expected.nonddg2p.inframe, F)
observed.nonddg2p.frameshift.p<-ppois(observed.nonddg2p.frameshift-1, expected.nonddg2p.frameshift, F)


rbind(expected.ddg2p.nonar.lof, expected.ddg2p.nonar.missense, expected.ddg2p.nonar.silent, expected.ddg2p.nonar.inframe, expected.ddg2p.nonar.frameshift)

rbind(expected.ddg2p.ar.lof, expected.ddg2p.ar.missense, expected.ddg2p.ar.silent, expected.ddg2p.ar.inframe, expected.ddg2p.ar.frameshift)

rbind(expected.nonddg2p.lof, expected.nonddg2p.missense, expected.nonddg2p.silent, expected.nonddg2p.inframe, expected.nonddg2p.frameshift)

rbind(observed.ddg2p.nonar.lof.p, observed.ddg2p.nonar.missense.p, observed.ddg2p.nonar.silent.p, observed.ddg2p.nonar.inframe.p, observed.ddg2p.nonar.frameshift.p)

rbind(observed.ddg2p.ar.lof.p, observed.ddg2p.ar.missense.p, observed.ddg2p.ar.silent.p, observed.ddg2p.ar.inframe.p, observed.ddg2p.ar.frameshift.p)

rbind(observed.nonddg2p.lof.p, observed.nonddg2p.missense.p, observed.nonddg2p.silent.p, observed.nonddg2p.inframe.p, observed.nonddg2p.frameshift.p)





