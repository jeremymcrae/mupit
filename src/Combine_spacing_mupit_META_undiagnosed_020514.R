# code to combine p values from numbers and spacing of DNMs

# read in p values from clustering analysis, only for genes with >1 mutation

spacing<-read.delim("/Volumes/DDD_meh/Analysis/Exome/For_1133_trio_ms/Clustering_results/de_novo_distance_simulations.increased_set_all.geometric_mean.tsv", header=T)

# read in p values from mupit/daly analyses

probs<-read.delim("/Volumes/DDD_meh/Analysis/Exome/For_1133_trio_ms/Meta_DD_undiagnosed/Mup-it_Daly_020514_META_undiagnosed.output.txt", header=T)

# function to combine p values
fishersMethod = function(x) pchisq(-2 *sum(log(x)),df=2*length(x),lower=FALSE)
# merge files

merged<-merge(probs, spacing, by.x=1, by.y=1, all.x=T )

num.genes<-length(merged[,1])

combined.p.dist<-rep(0,num.genes)

combined.p.dist.cons<-rep(0,num.genes)

for (i in seq(1,num.genes)) {
	
	combined.p.dist[i]<-fishersMethod(c(merged$func_dist_probability[i], merged$daly.p.DNM.func[i] ))
	
#	combined.p.dist.cons[i]<-fishersMethod(c(merged$func_dist_probability[i], merged$func_conservation_probability[i], min(merged$daly.p.DNM.lof[i], merged$daly.p.DNM.func[i] )))

#	combined.p.dist.cons[i]<-fishersMethod(c(merged$func_dist_probability[i], merged$func_conservation_probability[i], merged$daly.p.DNM.func[i]))
	
	}

num.tests<-18500

combined.fdr.dist<-p.adjust(combined.p.dist, method="BH", n=num.tests)

#combined.fdr.dist.cons<-p.adjust(combined.p.dist.cons, method="BH", n=num.tests)

output<-cbind(merged, combined.p.dist, combined.fdr.dist)


write.table(output, file="/Volumes/DDD_meh/Analysis/Exome/For_1133_trio_ms/Meta_DD_undiagnosed/Mup-it_Daly_020514_META_undiagnosed.output.combined.txt", row.names=F, quote=F, sep="\t")

