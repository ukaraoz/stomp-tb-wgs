# calculate distances between strains as number of SNPs

# on cori, conda environment R4
# install.packages(c("adegenet", "ape", "matrixcalc"), .libPaths()[1])
# install.packages(c("seqinr", "doParallel", "foreach"), .libPaths()[1])

library(adegenet)
library(ape)
library(seqinr)
library(matrixcalc)
library(doParallel)
library(foreach)
library(seqinr)

#fasta.file = "EPI.snp.nodensityfilter.final.merged.min0.fasta"
base = "/Users/ukaraoz/Work/MTB/stomp-tb"
fasta.file = "final.merged.min0.fasta"

#fasta = seqinr::read.fasta(fasta.file, as.string = TRUE, forceDNAtolower = FALSE)
#seqs = unlist(getSequence(fasta, as.string = TRUE))

fasta = adegenet::fasta2DNAbin(file.path(base, "results", fasta.file))
dist = ape::dist.dna(fasta, model="N",as.matrix=TRUE, pairwise.deletion = TRUE)
dist2 = matrixcalc::upper.triangle(as.matrix(dist))

registerDoParallel(cores = floor(detectCores()*0.7))
total_distances=NULL
total_distances<-foreach(i=1:dim(dist2)[1], .errorhandling = 'remove') %dopar%
{
  distsample=NULL
  for (j in i:dim(dist2)[2]){
    if (rownames(dist2)[i]!=colnames(dist2)[j]) distsample=rbind(distsample,c(rownames(dist2)[i],dist2[i,j],colnames(dist2)[j]))
  }
  return(distsample)
}
distances=NULL
for(i in 1:length(total_distances))
{
  distances=rbind(distances, unlist(total_distances[[i]]))
}

#for (i in 1:dim(dist2)[1]){
#  for (j in i:dim(dist2)[2]){
#    if (rownames(dist2)[i]!=colnames(dist2)[j]) total_distances=rbind(total_distances,c(rownames(dist2)[i],dist2[i,j],colnames(dist2)[j]))
#  }
#}
total_ordered=distances[order(as.numeric(distances[,2])),]
colnames(total_ordered) = c("sample1", "dist", "sample2")
write.table(total_ordered, file = file.path(base, "results", sub(".fasta", ".geneticdistances_flat.txt", fasta.file)), quote = F, row.names = FALSE)
write.table(dist, file = file.path(base, "results", sub(".fasta", ".geneticdistances_matrix.txt", fasta.file)), quote = F, sep = "\t", row.names = TRUE, col.names = TRUE)
