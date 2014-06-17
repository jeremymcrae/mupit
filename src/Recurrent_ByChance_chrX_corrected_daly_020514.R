# estimate probability of seeing NN genes with >1 functional mutations given number of mutations seen

# number of silent DNMs in DDD families, 1 gene with >1 mutation

# can adapt to get number expected with 3, 4, 5, 6, ... mutations to work out PPV of each number

#ddg2p<-read.delim("~/Desktop/DDD/DDG2P/DDG2P_ConfProbPossBoth_130521_trim.txt", header=T)


# could calculate mean analytically using mutation rate scaled to fit observed number of genes
# sum(ppois(1, gene.func.rate*num.DNM.func/sum(gene.func.rate), lower.tail=F))
# but how to calculate variance without simulation?


gene.info<-read.delim("/Volumes/DDD_meh/Analysis/Exome/Recurrent_DNM_signif/CDS_LENGTH_B37_chr.txt", header=T)

cds.length<-gene.info$CDS_LENGTH


#set number of functional DNMs and number of genes with >1 functional DNMs

num.DNM.func<-1106

num.recurr.genes<-96

num.sims<-10000 # number of simulations to perform


# calculate numbers of transmissions for autosomes and chrX

num.trios.male<-582 # trios with male offspring

num.trios.female<-548 # trios with female offspring

num.trios<-num.trios.male + num.trios.female # to give total of 1130 (1133 minus 1 twin of 3 identical twins)

auto.transmissions<-2*num.trios

female.transmissions<-num.trios # for chrX analyses

male.transmissions<-num.trios.female # for chrX analyses

alpha=3.4

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


gene.func.rate<-daly.gene.snv.missense.rate + daly.gene.snv.lof.rate + daly.gene.indel.missense.rate + daly.gene.indel.lof.rate

num.genes<-length(gene.func.rate)

summed.rates<-rep(0,num.genes)

for (i in seq(0,num.genes)) summed.rates[i]<-sum(gene.func.rate[1:i])

scaled.summed.rates<-summed.rates/max(summed.rates)

#index<-round(runif(num.DNM.func)*max(round.summed.rates))

#find.index<-findInterval(index, round.summed.rates)

#table(table(find.index))

# table(table(findInterval(round(runif(num.DNM.func)*sum(cds.length[,3])), cds.length$sum)))

max.obs<-15 # max times a gene might be expected to be mutated

store<-matrix(nrow=num.sims, ncol=max.obs)

for (i in seq(1,num.sims)) store[i,]<-hist(table(findInterval(runif(num.DNM.func),scaled.summed.rates)), 0:max.obs, plot=FALSE)$counts

count.means<-colMeans(store)

recurrent<-rowSums(store[,2:max.obs])

if(count.means[max.obs]>0) print("warning maximum observations reached, extend max obs")


#plot simulations versus observed

#observed = number of genes with recurrent DNMs


hist(recurrent, breaks=seq(-0.5,num.recurr.genes+5,1), col="green", main="Recurrent genes by chance (Functional DNMs=1106)", xlab="Number of genes with recurrent DNMs")

abline(v=num.recurr.genes, col="red")


# write out table for plotting boxplot of the same data

write.table(recurrent, file="/Volumes/DDD_meh/Analysis/Exome/For_1133_trio_ms/Recurrent_by_chance/Num_recurr_sim_1106.txt", quote=F, row.names=F, sep="\t")

