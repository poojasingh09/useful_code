
## input files

vcf="all.reseq.raw.chr1.snps.filt.vcf"
info="popinfo_edit.txt"


##libraries

library(fastreeR)
options(java.parameters="-Xmx1G")
unloadNamespace("fastreeR")
library(fastreeR)
library(utils)
library(ape)
library(stats)
library(grid)
library(BiocFileCache)
library(WGCNA)

## get stats from vcf

myVcfIstats <- fastreeR::vcf2istats(inputFile = vcf)
plot(myVcfIstats[,7:9])
write.table(myVcfIstats, file = "all.reseq.raw.chr1.snps.filt.stats.txt", sep="\t", row.names = FALSE, quote=FALSE)

pop_code <- read.table(info, header=T)
pop_code <- pop_code[order(match(pop_code[,1],myVcfIstats[,1])),] # this re-orders your pop info to your pca file

myVcfIstats$GENUS <- pop_code$genus
myVcfIstats$LOC <- pop_code$location

## get distance from vcf

myVcfDist <- fastreeR::vcf2dist(inputFile = vcf, threads = 2)


# distance historgram
graphics::hist(myVcfDist, breaks = 100, main=NULL, 
                                xlab = "Distance", xlim = c(0,max(myVcfDist)))


#plot using  fastreeR dist2tree

myVcfTree <- fastreeR::dist2tree(inputDist = myVcfDist)
plot(ape::read.tree(text = myVcfTree), direction = "right", cex = 0.5)
ape::add.scale.bar()
ape::axisPhylo(side = 2)

## plot using hclust
myVcfTreeStats <- stats::hclust(myVcfDist)

pdf("snp_hclust.pdf")
plot(myVcfTreeStats, cex = 0.5)
dev.off()


## hclust plot annotated

anno <-data.frame(genus=factor(pop_code$genus),location=factor(pop_code$location),row.names=pop_code$sample)
anno$all <-paste0(anno$genus,"_",anno$location)

colours = labels2colors(anno)
colnames(colours) <- c("genus", "location")
rownames(colours) = rownames(anno)

pdf("snp_hclust_annotated.pdf")
plotDendroAndColors(myVcfTreeStats, colours, main = "SNP hclust")
dev.off()

myVcfTreeStats$labels <- anno$all

pdf("snp_hclust_annotated_v2.pdf")
plotDendroAndColors(myVcfTreeStats, colours, main = "SNP hclust")
dev.off()



