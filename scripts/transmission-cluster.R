library(seqinr)
library(ape)
library(pegas)

plot_clustsize_distrib <- function(clusters, pdf) {
  temp = table(clusters)
  cluster2size = data.frame(cluster = names(temp), size = as.numeric(temp))
  table = table(cluster2size[,2])
  
  pdf(width = 8, height = 8, file = pdf)
  barplot(table, xlab = "cluster size", ylab = "number of clusters", yaxt='n')
  axis(2, at = as.numeric(as.character(table)), labels =  as.character(table), las = 2, cex.axis = 0.7)
  dev.off()
}

plot_distance_distrib <- function(distances, pdf) {
  pdf(width = 12, height = 8, file = pdf)
  hist(distances, xlab = "SNP distance", ylab = "count", main = "")
  dev.off()
}

base = "/Users/ukaraoz/Work/MTB/stomp-tb"
fasta.file = "final.merged.min0.fasta"
nsnps_threshold = 12    # number of nucleotide differences used to define clusters


dist = read.table(file.path(base, "results", sub(".fasta", ".geneticdistances_matrix.txt", fasta.file)), check.names = F)
dist.bin = dist
for(r in 1:nrow(dist)) {
  dist.bin[r, which(dist[r,] < nsnps_threshold)] = 0
  dist.bin[r, which(dist[r,] >= nsnps_threshold)] = 1
}
#write.table(dist.bin, file=file.path(base, "final.merged.min0.genetic_distances_distmatrix-thresh12.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

# cluster
hclust = hclust(as.dist(dist.bin))
hclust.cut = cutree(hclust, h=0)
sample2cluster = data.frame(sample = names(hclust.cut), cluster = paste0("cluster_", hclust.cut))
write.table(sample2cluster, file = file.path(base, "results", paste0(sub(".fasta", "", fasta.file), ".sample2cluster", "-lt", nsnps_threshold, "snps.txt")), 
  quote = F, row.names = F, col.names = T, sep = "\t")

plot_clustsize_distrib(clusters = sample2cluster[, "cluster"], 
                    pdf = file.path(base, "results", paste0(sub(".fasta", "", fasta.file), "clustsize-distrib", "-lt", nsnps_threshold, ".pdf")))

plot_distance_distrib(distances = t(dist)[upper.tri(t(dist))],
                       pdf = file.path(base, "results", paste0(sub(".fasta", "", fasta.file), "distance-distrib", "-lt", nsnps_threshold, ".pdf")))

