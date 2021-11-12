
library(seqinr)
library(ape)
library(pegas)
#base = "/global/homes/u/ukaraoz/cscratch/mtb/stomp-tb"
base = "/Users/ukaraoz/Work/MTB/uganda/stomp-tb"

dist = read.table(file.path(base, "final.merged.min0.genetic_distances_distmatrix.txt"), check.names = F)
dist2=upper.triangle(as.matrix(dist))
pdf(width = 12, height = 8, file = file.path(base, "final.merged.min0.genetic_distances_distmatrix-SNPdistancedistr.pdf"))
#hist(dist2, breaks=c(0,12,20,50,100,500, 1000, 1500, 1780))
hist(dist2, xlab = "SNP distance", ylab = "count", main = "")
dev.off()

dist.bin = dist
threshold = 12
for(r in 1:nrow(dist)) {
  temp = dist[r,]
  one = which(temp >= threshold)
  zero = which(temp < threshold)
  dist.bin[r, zero] = 0
  dist.bin[r, one] = 1
}
write.table(dist.bin, file=file.path(base, "final.merged.min0.genetic_distances_distmatrix-thresh12.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

h = hclust(as.dist(dist.bin))
hcut = cutree(h, h=0)
sample2cluster = data.frame(sample = names(hcut), cluster = paste0("cluster_", hcut))
write.table(sample2cluster, file=file.path(base, "final.merged.min0.genetic_distances_distmatrix-thresh12.sample2cluster.xls"), quote = F, row.names = F, col.names = T, sep = "\t")

temp = table(sample2cluster[, "cluster"])
cluster2size = data.frame(cluster = names(temp), size = as.numeric(temp))
pdf(width = 8, height = 8, file = file.path(base, "final.merged.min0.genetic_distances_distmatrix-thresh12.clustsizedist.pdf"))
barplot(table(cluster2size[,2]), xlab = "cluster size", ylab = "number of clusters", yaxt='n')
axis(2,at = c(132, 25, 11, 1, 2, 0), labels =  c(132, 25, 11, 1, 2, 0), las = 2)
dev.off()

msa = read.fasta(file.path(base, "final.merged.min0.fasta"))
msa.list = list()
for(i in 1:length(msa)) {
  msa.list[[names(msa)[i]]] = as.character(msa[[i]])
}
write.nexus.data(msa, file = file.path(base, "final.merged.min0.nex"), format="dna")

msa.dnabin <- as.matrix(as.DNAbin(msa.list[1:10]))
(ntz <- mjn(msa.dnabin, 0))


z <- list(c("g", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a"),
          c("a", "g", "g", "a", "a", "a", "a", "a", "a", "a", "a", "a"),
          c("a", "a", "a", "g", "a", "a", "a", "a", "a", "a", "g", "g"),
          c("a", "a", "a", "a", "g", "g", "a", "a", "a", "a", "g", "g"),
          c("a", "a", "a", "a", "a", "a", "a", "a", "g", "g", "c", "c"),
          c("a", "a", "a", "a", "a", "a", "g", "g", "g", "g", "a", "a"))
names(z) <- c("A1", "A2", "B1", "B2", "C", "D")
z <- as.matrix(as.DNAbin(z))
(ntz <- mjn(z, 0))
plotNetMDS(ntz, dist.dna(attr(ntz, "data"), "N"), 3)
