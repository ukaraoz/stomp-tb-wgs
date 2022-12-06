library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)
library(gtools)
#library(adegenet)
#library(ape)
#library(seqinr)
#library(matrixcalc)
#library(doParallel)
#library(foreach)

density.positions <- function(data.file, window = 10, max_density = 2) {
  data = read.table(data.file, header = T, sep = "\t")
  # %>% as_tibble() %>%  dplyr::select(select_cols) %>% tibble::add_column(type = "unfiltered")
  density = 1
  candidates = integer()
  density_pos = integer()
  for(i in 1:(nrow(data)-1)) {
    posA = data[i, "Position"]
    posB = data[i+1, "Position"]
    distAB = posB-posA
    
    if(distAB <= window) {
      density = density + 1
      candidates = c(candidates, posA)
      #cat(i, "\t", "Accumulate:", posA, "\t", distAB, "\t", density, "\n")
    } else {
      if(density >= max_density) {
        #cat(i, "\t", "Density maxed:", posA, "\t", distAB, "\t", density, "\t", paste0(candidates, collapse= ":"), "\n")
        for(c in 1:length(candidates)) {
          density_pos = c(density_pos, candidates[c])
        }
        density_pos = c(density_pos, posA)
        #cat(i, "\t", "Density maxed:", posA, "\t", distAB, "\t", density, "\t", paste0(density_pos, collapse= ":"), "\n")
      }
      candidates = integer()
      density = 1
    }
  }
  return(density_pos)
}

base = "/Users/ukaraoz/Work/MTB/stomp-tb"
data.dir = file.path(base, "results/raw")
ids = read.table(file.path(base, "stomptb_ids_295.txt")) %>% pull(V1)
ids.wSNPs = sub(".EPI.snp.nodensityfilter", "", list.files(file.path(base, "results/raw"), pattern = "EPI.snp.nodensityfilter"))
ext.filtered = ".EPI.snp.nodensityfilter"
ext.density_filtered = ".EPI.snp.withdensityfilter"

data.file = file.path(data.dir, paste0(ids.wSNPs[1], ext.filtered))
density_pos = density.positions(data.file)
data_filtered = read.table(data.file, header = T) %>% as_tibble() %>%
  dplyr::filter(!(Position %in% density_pos))

result = pair.dist(ids.pairs[r, 1], ids.pairs[r, 2], data.dir = file.path(base, "results/raw"))
data.dir = file.path(base, "results/raw")


a = parallel::mclapply(1:length(ids.wSNPs),
                    function(i) {
                      data.file = file.path(data.dir, paste0(ids.wSNPs[i], ext.filtered))
                      density_pos = density.positions(data.file)
                      data_filtered = read.table(data.file, header = T) %>% as_tibble() %>%
                        dplyr::filter(!(Position %in% density_pos))
                      write.table(data_filtered, 
                                  file = file.path(data.dir, paste0(ids.wSNPs[i], ext.density_filtered)),
                                  quote = F, sep = "\t", row.names = FALSE)
                      length(density_pos)
                    },
                    mc.cores = 4)

# compare pre and post filtering
results = parallel::mclapply(1:length(ids.wSNPs),
                    function(i) {
                      data1.file = file.path(data.dir, paste0(ids.wSNPs[i], ext.filtered))
                      data2.file = file.path(data.dir, paste0(ids.wSNPs[i], ext.density_filtered))
                      data1 = read.table(data1.file, header = T)
                      data2 = read.table(data2.file, header = T)
                      result = c(ids.wSNPs[i], nrow(data1), nrow(data2))
                    },
                    mc.cores = 4)
results_rbind = do.call("rbind", results)
colnames(results_rbind) = c("id", "qcfiltered", "qc_and_densityfiltered")
write.table(results_rbind, file = file.path(base, "results", "densityfiltering-window10maxdens2_stats.xls"), quote = F, sep = "\t", row.names = FALSE)

