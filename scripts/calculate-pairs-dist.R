conda activate R4

library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)
library(gtools)


pair.dist <- function(id1, id2, data.dir, ext.filtered, ext.unfiltered){
  select_cols = c("Chrom", "Position", "Ref", "Cons", "VarAllele",
                "Reads1", "Reads2", "VarFreq", "Pvalue",
                "Qual1", "Qual2", "Strands1", "Strands2")

  # rescue SNPs from unfiltered list, by two-way comparisons even if they didn't pass the depth and percent filters
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
  data2.filtered.wrescued = data2.filtered.wrescued %>%
    dplyr::rename_with(.fn = ~paste0(., ".", id2, sep = ""), .cols = -c(Chrom, Position, Ref))
  

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
  data1.filtered.wrescued = data1.filtered.wrescued %>%
    dplyr::rename_with(.fn = ~paste0(., ".", id1, sep = ""), .cols = -c(Chrom, Position, Ref))

  data12 = data1.filtered.wrescued %>% 
    dplyr::full_join(data2.filtered.wrescued, by = c("Position" = "Position")) %>%
    dplyr::select(Chrom.x, Position, Ref.x, 
                  paste0("Cons.", id1),
                  paste0("VarAllele.", id1),
                  paste0("Reads1.", id1),
                  paste0("Reads2.", id1),
                  paste0("VarFreq.", id1),
                  paste0("type.", id1),
                  Chrom.y, Ref.y, 
                  paste0("Cons.", id2),
                  paste0("VarAllele.", id2),
                  paste0("Reads1.", id2),
                  paste0("Reads2.", id2),
                  paste0("VarFreq.", id2),
                  paste0("type.", id2)) %>%
    dplyr::mutate(Chrom.x = case_when(is.na(Chrom.x) ~ "MTB_anc",
                                      TRUE ~ Chrom.x)) %>%
    dplyr::mutate(Ref.x = case_when(is.na(Ref.x) ~ Ref.y,
                                      TRUE ~ Ref.x)) %>%
    dplyr::mutate_at(vars(matches("VarFreq")), ~ as.numeric(sub("%", "", .))) %>%
    dplyr::rename(Chrom = Chrom.x, Ref = Ref.x) %>%
    dplyr::select(-Chrom.y, -Ref.y) %>%
    dplyr::arrange(Position) %>%
    dplyr::mutate(comparison = case_when(!is.na(get(paste0("VarAllele.", id1))) & is.na(get(paste0("VarAllele.", id2))) ~ "unique",
                                         is.na(get(paste0("VarAllele.", id1))) & !is.na(get(paste0("VarAllele.", id2))) ~ "unique",
                                         ((!is.na(get(paste0("VarAllele.", id1)))) & 
                                          (!is.na(get(paste0("VarAllele.", id2)))) &
                                          (get(paste0("VarAllele.", id1)) == get(paste0("VarAllele.", id2)))) ~ "common",
                                         ((!is.na(get(paste0("VarAllele.", id1)))) & 
                                          (!is.na(get(paste0("VarAllele.", id2)))) &
                                          (get(paste0("VarAllele.", id1)) != get(paste0("VarAllele.", id2)))) ~ "unique",
                                      TRUE ~ "other")) %>%
    dplyr::select(Chrom, Position, Ref, comparison, everything())

    data1.pos.rescued = data12 %>% dplyr::filter(get(paste0("type.", id1)) %in% c("filtered", "unfiltered")) %>% dplyr::pull(Position)
    data2.pos.rescued = data12 %>% dplyr::filter(get(paste0("type.", id2)) %in% c("filtered", "unfiltered")) %>% dplyr::pull(Position)
    data1.nrescued = data1.pos.rescued %>% length()
    data2.nrescued = data2.pos.rescued %>% length()
    intersect.nrescued = data12 %>% dplyr::filter(comparison == "common") %>% dplyr::pull(Position) %>% length()
    diff.nrescued = data12 %>% dplyr::filter(comparison == "unique") %>% dplyr::pull(Position) %>% length()

    #result = list()
    #result[[1]] = data.frame(id1 = id1, 
    #                         id2 = id2, 
    #                         id1.nrescued = data1.nrescued,
    #                         id2.nrescued = data2.nrescued,
    #                         intersect.nrescued = intersect.nrescued,
    #                         diff.nrescued = diff.nrescued)
    #result[[2]] = data12
    result = data.frame(id1 = id1, 
                        id2 = id2, 
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


#base = "/Users/ukaraoz/Work/MTB/stomp-tb"
#data.dir = file.path(base, "results/raw")
#ids = read.table(file.path(base, "stomptb_ids_295.txt")) %>% pull(V1)
#ids.wSNPs = sub(".EPI.snp.nodensityfilter", "", list.files(file.path(base, "results/raw"), pattern = "EPI.snp.nodensityfilter"))

base = "/global/homes/u/ukaraoz/cscratch/mtb/stomptb"
data.dir = file.path(base, "MTB_ancestor/out")
ids.wSNPs = sub(".EPI.snp.nodensityfilter.vcf", "", list.files(file.path(base, "MTB_ancestor/out"), pattern = "EPI.snp.nodensityfilter.vcf"))

ids.pairs = expand.grid(ids.wSNPs, ids.wSNPs, stringsAsFactors = F) %>% as_tibble() %>% dplyr::filter(Var1 != Var2)
#ids.pairs = gtools::combinations(n = 238, r = 2, v = ids, repeats.allowed = F)
#ids.pairs = gtools::permutations(n = 238, r = 2, v = ids, repeats.allowed = T)
#ext.filtered = ".EPI.snp.nodensityfilter"
ext.filtered = ".EPI.snp.withdensityfilter" # with density filter
ext.unfiltered = ".snp.pileup2snp"
id1=ids.pairs[r, 1]
id2=ids.pairs[r, 2]
data.dir = data.dir
distances = parallel::mclapply(1:nrow(ids.pairs),
                    function(r) {
                      result = pair.dist(ids.pairs[r, 1], ids.pairs[r, 2], 
                                         data.dir = data.dir,
                                         ext.filtered = ext.filtered,
                                         ext.unfiltered = ext.unfiltered)
                      #saveRDS(result, file = paste0(file.path(base, "out/"), r, ".rds"))
                      result
                    },
                    mc.cores = 40)
distances_rbind = lapply(distances,"[[",1)
distances_rbind = do.call("rbind", distances_rbind)

distances_rbind.diff.nrescued = distances_rbind %>% as_tibble() %>% 
  dplyr::rename(sample1 = Var1, sample2 = Var2, dist = diff.nrescued) %>% 
  dplyr::filter(sample1 != sample2) %>%
  dplyr::arrange(dist) %>%
  dplyr::select(sample1, dist, sample2)

distances_rbind.diff.nrescued.matrix = distances_rbind.diff.nrescued %>%
  tidyr::pivot_wider(names_from = sample1, values_from = dist) %>%
  tibble::column_to_rownames(var = "sample2")

write.table(distances_rbind.diff.nrescued, file = file.path(base, "MTB_ancestor/results", "pairs-dist-densityfiltered_wrescueSNPs.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(distances_rbind.diff.nrescued.matrix, file = file.path(base, "MTB_ancestor/results", "distmatrix-densityfiltered_wrescueSNPs.txt"), quote = F, sep = "\t", row.names = FALSE)

pairs = list()
pairs[[1]] = c("STP63533-S49", "STP63533_S660_L006")
pairs[[2]] = c("STP63747-S53", "STP63409_S658_L006")
pairs[[3]] = c("STP69976-S44", "STP70469-S31")
pairs[[4]] = c("STP68461-S41", "STP68351-S46")
pairs[[5]] = c("57875B_S485_L005", "57875A_S478_L005")
pairs[[6]] = c("57875B_S485_L005", "STP57875-S4")
pairs[[7]] = c("57875A_S478_L005", "STP57875-S4")
pairs[[8]] = c("5872_S85", "STP_35872_S261_L004")
pairs[[9]] = c("5241_S96", "STP_35241_S213_L004")
pairs[[10]] = c("STP62336_S655_L006", "STP62359_S649_L006")
pairs[[11]] = c("STP61302_S637_L006", "STP61313_S627_L006")
pairs[[12]] = c("57875A_S478_L005", "57875B_S485_L005")
pairs[[13]] = c("56124_S488_L005", "STP_56124A_S217_L004")
pairs[[14]] = c("56124_S488_L005", "STP_56124B_S236_L004")
pairs[[15]] = c("STP_56124A_S217_L004", "STP_56124B_S236_L004")
pairs[[16]] = c("STP62440_S673_L006", "STP62498_S642_L006")

for(i in 1:length(pairs)) {
  temp = distances_rbind.diff.nrescued %>%
    dplyr::filter(sample1 == pairs[[i]][1] & sample2 == pairs[[i]][2]) %>%
    as.data.frame()
  cat(temp[1, "sample1"], "\t", temp[1, "dist"], "\t", temp[1, "sample2"], "\n")
 }




