#library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)

#tree.file = "/Users/ukaraoz/Work/MTB/philippines/RAxML_bestTree.result"
#tree = read.tree(tree.file)
#tree = drop.tip(tree, "RO1_18_602_S316_L004")
#write.tree(tree, "/Users/ukaraoz/Work/MTB/philippines/RAxML_bestTree_262.result")

tree = read.tree("/Users/ukaraoz/Work/MTB/philippines/RAxML_bestTree_262.result")
annotation <- read.table("/Users/ukaraoz/Work/MTB/philippines/philippines.lineage.txt", header = T, sep = "\t")
rownames(annotation) = annotation[,1]
annotation$lineage = factor(annotation$lineage)
p <- ggtree(tree,layout='fan', open.angle=5, size=0.2) + 
  geom_tiplab(size = 3, align=TRUE, linetype = "dashed", linesize = 0.0001)
p <- p %<+% annotation

p1 = p + 
  geom_fruit(geom=geom_tile,
         mapping=aes(fill=lineage),
         width=0.01,
         offset=0.6) +
  geom_fruit(geom=geom_tile,
         mapping=aes(fill=lineagesub1),
         width=0.01,
         offset=0.5)
ggsave(p1, filename = "~/Desktop/tree1.png", width = 12, height = 12)
#############


  scale_fill_manual(
      name="Lineage",
      values=yearcolors$year__colour,
      guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=2)
  )

heatmap.colours <- c("red", "darkgreen", "blue")
names(heatmap.colours) <- c("lineage1", "lineage2", "lineage4")
gheatmap(p, annotation, offset = 10, color=NULL, 
         colnames_position="top", 
         colnames_angle=90, colnames_offset_y = 1, 
         hjust=0, font.size=2) +
scale_fill_manual(values=heatmap.colours)

%<+% annotation + xlim(-.1, 4) +
  geom_tiplab(size = 2)
gheatmap(p, annotation, offset=5, width=0.5, font.size=3, colnames_angle=-45, hjust=0) +
    scale_fill_manual(breaks=c("lineage1", "lineage2", "lineage4"), 
        values=c("steelblue", "firebrick", "darkgreen"), name="lineage")

scale_fill_manual(breaks=c("red", "green", "blue"), 
        values=c("lineage1", "lineage2", "lineage4"), name = "lineage")


p2 <- p + geom_tiplab(offset = .6, hjust = .5) +
    geom_tippoint(aes(shape = trophic_habit, color = trophic_habit, size = mass_in_kg)) + 
    theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))
