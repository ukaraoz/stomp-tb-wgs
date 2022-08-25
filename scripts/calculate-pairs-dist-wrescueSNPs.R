library(dplyr)
library(tidyr)
library(parallel)
library(gtools)
#library(adegenet)
#library(ape)
#library(seqinr)
#library(matrixcalc)
#library(doParallel)
#library(foreach)

pair.dist <- function(id1, id2, data.dir, ext.filtered = ".EPI.snp.nodensityfilter", ext.unfiltered = ".snp.pileup2snp"){
  select_cols = c("Chrom", "Position", "Ref", "Cons",
                "Reads1", "Reads2", "VarFreq", "Pvalue",
                "Qual1", "Qual2", "Strands1", "Strands2")
  data1.unfiltered = read.table(file.path(data.dir, paste0(id1, ext.unfiltered)), header = T) %>% as_tibble() %>% dplyr::select(select_cols) %>% tibble::add_column(type = "unfiltered")
  data1.filtered = read.table(file.path(data.dir, paste0(id1, ext.filtered)), header = T) %>% as_tibble() %>% dplyr::select(select_cols) %>% tibble::add_column(type = "filtered")
  data2.unfiltered = read.table(file.path(data.dir, paste0(id2, ext.unfiltered)), header = T) %>% as_tibble() %>% dplyr::select(select_cols) %>% tibble::add_column(type = "unfiltered")
  data2.filtered = read.table(file.path(data.dir, paste0(id2, ext.filtered)), header = T) %>% as_tibble() %>% dplyr::select(select_cols) %>% tibble::add_column(type = "filtered")

  data1notin2.filtered =  data1.filtered %>% 
    dplyr::anti_join(data2.filtered, by = c("Position" = "Position"))
  # check unfiltered list from 2 to recover SNPs, and add them to data2 list
  data2.filtered.wrescued = data1notin2.filtered %>% 
    dplyr::inner_join(data2.unfiltered, by = c("Position" = "Position"), suffix = c(".data1", ".data2"), keep = T) %>%
    dplyr::select(ends_with(".data2")) %>%
    dplyr::rename_with(.fn = ~gsub(".data2", "", .), .cols = everything())
  if(nrow(data2.filtered) != 0) {
    data2.filtered.wrescued = data2.filtered.wrescued %>%
      dplyr::bind_rows(data2.filtered)
  }

  data2notin1.filtered =  data2.filtered %>% 
    dplyr::anti_join(data1.filtered, by = c("Position" = "Position"))
  # check unfiltered list from 2 to recover SNPs, and add them to data2 list
  data1.filtered.wrescued = data2notin1.filtered %>% 
    dplyr::inner_join(data1.unfiltered, by = c("Position" = "Position"), suffix = c(".data2", ".data1"), keep = T) %>%
    dplyr::select(ends_with(".data1")) %>%
    dplyr::rename_with(.fn = ~gsub(".data1", "", .), .cols = everything())
  if(nrow(data1.filtered) != 0) {
    data1.filtered.wrescued = data1.filtered.wrescued %>%
      dplyr::bind_rows(data1.filtered)
  }


  data1.pos.rescued = data1.filtered.wrescued %>% pull(Position)
  data2.pos.rescued = data2.filtered.wrescued %>% pull(Position)
  data1.nrescued = length(data1.pos.rescued)
  data2.nrescued = length(data2.pos.rescued)
  intersect.nrescued  = length(intersect(data1.pos.rescued, data2.pos.rescued))
  diff.nrescued = length(union(setdiff(data1.pos.rescued, data2.pos.rescued), setdiff(data2.pos.rescued, data1.pos.rescued)))


  data1.pos = data1.filtered %>% pull(Position)
  data2.pos = data2.filtered %>% pull(Position)
  data1.n = length(data1.pos)
  data2.n = length(data2.pos)
  intersect.n = length(intersect(data1.pos, data2.pos))
  diff.n = length(union(setdiff(data1.pos, data2.pos), setdiff(data2.pos, data1.pos)))

  result = data.frame(id1 = id1, 
                      id2 = id2, 
                      id1.n = data1.n,
                      id2.n = data2.n,
                      intersect.n = intersect.n,
                      diff.n = diff.n,
                      id1.nrescued = data1.nrescued,
                      id2.nrescued = data2.nrescued,
                      intersect.nrescued = intersect.nrescued,
                      diff.nrescued = diff.nrescued)
  result
}

check_symetry <- function(matrix) {
  for(i in 1:nrow(matrix)) {
    for(j in (i+1):ncol(matrix)) {
      a = matrix[i, j] - matrix[j, i]
      cat(i, j, matrix[i, j], matrix[j, i],"\n")
      #if(a != 0) {cat(i, j, "\n")}
    }
  }
}

isolate_pairs = cbind(c("5872_S85", "5241_S96", "STP62336_S655_L006", "STP61302_S637_L006", "57875A_S478_L005", "56124_S488_L005", "56124_S488_L005", "STP_56124A_S217_L004", "STP62440_S673_L006"),
  c("STP_35872_S261_L004", "STP_35241_S213_L004", "STP62359_S649_L006", "STP61313_S627_L006", "57875B_S485_L005", "STP_56124A_S217_L004", "STP_56124B_S236_L004", "STP_56124B_S236_L004", "STP62498_S642_L006"))
for(i in 1:nrow(isolate_pairs)) {
  result = pair.dist(isolate_pairs[i, 1], isolate_pairs[i, 2], data.dir = file.path(base, "results"))
  cat(i,  result[1,"diff.nrescued"], "\n")
}

base = "/Users/ukaraoz/Work/MTB/stomp-tb/"
ids = read.table(file.path(base, "stomptb_ids_295.txt")) %>% pull(V1)
ids.wSNPs = sub(".EPI.snp.nodensityfilter", "", list.files(file.path(base, "results"), pattern = "EPI.snp.nodensityfilter"))

base = "/global/homes/u/ukaraoz/cscratch/mtb/stomptb"
ids.wSNPs = sub(".EPI.snp.nodensityfilter.vcf", "", list.files(file.path(base, "MTB_ancestor/out"), pattern = "EPI.snp.nodensityfilter.vcf"))

ids.pairs = expand.grid(ids.wSNPs, ids.wSNPs, stringsAsFactors = F)
#ids.pairs = gtools::combinations(n = 238, r = 2, v = ids, repeats.allowed = F)
#ids.pairs = gtools::permutations(n = 238, r = 2, v = ids, repeats.allowed = T)
distances = parallel::mclapply(1:nrow(ids.pairs),
                    function(r) {
                      result = pair.dist(ids.pairs[r, 1], ids.pairs[r, 2], data.dir = file.path(base, "MTB_ancestor/out"))
                      #saveRDS(result, file = paste0(file.path(base, "out/"), r, ".rds"))
                      result
                    },
                    mc.cores = 40)
distances_rbind = do.call("rbind", distances)
distances_rbind.diff.nrescued = distances_rbind[, c("id1", "id2", "diff.nrescued")] %>% as_tibble() %>% 
  dplyr::arrange(diff.nrescued) %>%
  dplyr::select(id1, diff.nrescued, id2) %>%
  dplyr::rename(sample1 = id1, dist = diff.nrescued, sample2 = id2)

distances_rbind.diff.nrescued.matrix = distances_rbind.diff.nrescued %>%
  tidyr::pivot_wider(names_from = sample1, values_from = dist) %>%
  tibble::column_to_rownames(var = "sample2")

write.table(distances_rbind.diff.nrescued, file = file.path(base, "MTB_ancestor/results", "pairs-dist-wrescueSNPs.txt"), quote = F, row.names = FALSE)
write.table(distances_rbind.diff.nrescued.matrix, file = file.path(base, "MTB_ancestor/results", "distmatrix-wrescueSNPs.txt"), quote = F, row.names = FALSE)




upperTriangle(x, diag=FALSE, byrow=FALSE)



total_ordered=distances1[order(as.numeric(distances1[, "diff.nrescued"])),]
colnames(total_ordered) = c("sample1", "dist", "sample2")
write.table(total_ordered, file = file.path(base, "results", sub(".fasta", ".geneticdistances_flat.txt", fasta.file)), quote = F, row.names = FALSE)
write.table(dist, file = file.path(base, "results", sub(".fasta", ".geneticdistances_matrix.txt", fasta.file)), quote = F, sep = "\t", row.names = TRUE, col.names = TRUE)














