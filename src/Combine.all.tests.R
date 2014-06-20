# code to combine the different tests statistics together into a single table

CODE_DIR = "/nfs/users/nfs_j/jm33/apps/enrichment_analysis"
DATA_DIR = file.path(CODE_DIR, "data")

# load all the data files in
ddd.all = read.delim("/Volumes/DDD_meh/Analysis/Exome/For_1133_trio_ms/DDD_all/Mup-it_Daly_020514_DDD.output.combined.txt", header=T)
ddd.undiagnosed = read.delim("/Volumes/DDD_meh/Analysis/Exome/For_1133_trio_ms/DDD_undiagnosed/Mup-it_Daly_020514_DDD_undiagnosed.output.combined.txt", header=T)
meta.all = read.delim("/Volumes/DDD_meh/Analysis/Exome/For_1133_trio_ms/Meta_DD_all/Mup-it_Daly_020514_META.output.combined.txt", header=T)
meta.undiagnosed = read.delim(file.path(DATA_DIR, "Mup-it_Daly_020514_META_undiagnosed.output.combined.txt"), header=T)

# need names that are more informative as same across files, add prefix
names(ddd.all) = paste("ddd.all",names(ddd.all), sep=".")
names(ddd.undiagnosed) = paste("ddd.undiagnosed",names(ddd.undiagnosed), sep=".")
names(meta.all) = paste("meta.all",names(meta.all), sep=".")
names(meta.undiagnosed) = paste("meta.undiagnosed",names(meta.undiagnosed), sep=".")

# merge together files, focusing on genes with DNMs in DDD
merged = merge(ddd.all, ddd.undiagnosed, by.x=1, by.y=1, all.x=T)
merged = merge(merged, meta.all, by.x=1, by.y=1, all.x=T)
merged = merge(merged, meta.undiagnosed, by.x=1, by.y=1, all.x=T)

# calculate minimum p value across LoF and func+clustering tests for each dataset
min.ddd.all.p = apply(cbind(merged$ddd.all.daly.p.DNM.lof, merged$ddd.all.combined.p.dist), 1, min, na.rm=T)
min.ddd.undiagnosed.p = apply(cbind(merged$ddd.undiagnosed.daly.p.DNM.lof,  merged$ddd.undiagnosed.combined.p.dist), 1, min, na.rm=T)
min.meta.all.p = apply(cbind(merged$meta.all.daly.p.DNM.lof, merged$meta.all.combined.p.dist), 1, min, na.rm=T)
min.meta.undiagnosed.p = apply(cbind(merged$meta.undiagnosed.daly.p.DNM.lof,  merged$meta.undiagnosed.combined.p.dist), 1, min, na.rm=T)
min.all.p = apply(cbind(min.ddd.all.p, min.meta.all.p), 1, min, na.rm=T)
min.undiagnosed.p = apply(cbind(min.ddd.undiagnosed.p, min.meta.undiagnosed.p), 1, min, na.rm=T)

# write out full table
merged = cbind(merged, min.ddd.all.p, min.ddd.undiagnosed.p, min.meta.all.p, min.meta.undiagnosed.p, min.all.p, min.undiagnosed.p)
write.table(merged, file="/Volumes/DDD_meh/Analysis/Exome/For_1133_trio_ms/All.tests.combined_v4.txt", row.names=F, quote=F, sep="\t")

# write out trimmed table with only useful fields
merged.trim = merged[,c(1,7,8,9,2,3,4,5,6,22,28,41,43,44,45,46,47,63,69,82,84,85,86,87,88,104,110,123,125,126,127,128,129,145,151,164, 166, 167, 168, 169, 170, 171)]
write.table(merged.trim, file="/Volumes/DDD_meh/Analysis/Exome/For_1133_trio_ms/All.tests.combined_v4-trim.txt", row.names=F, quote=F, sep="\t")
