#pooja.singh09@gmail.com
#feb2022 this script makes a PCA from vcf of SNPs. The data points can be annotated as shapes and colours
# Load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)
library(stringr)
library(ggplot2)

vcf.fn <- "all.reseq.raw.chr1.snps.filt.vcf"

snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")

snpgdsSummary("test.gds")
genofile <- openfn.gds("test.gds")
pca<-snpgdsPCA(genofile)

plot(pca$eigenvect[,1],pca$eigenvect[,2] ,col=as.numeric(substr(pca$sample, 1,3) == 'sample.id:')+3, pch=2)



pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

#plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

# Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

pop_code <- read.table("popinfo") # this is a text file with rows as samples and columns as location, genus, whatever else you want

a <- data.frame(pca$sample.id)
pop_code <- pop_code[order(match(pop_code[,1],a[,1])),] # this re-orders your pop info to your pca file

# assume the order of sample IDs is as the same as population codes
head(cbind(sample.id, pop_code))

colnames(pop_code) <- c("sampleid", "pop_code")

tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(pop_code$pop_code)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

tab$genus <- str_split_fixed(tab$pop, "_", 2)[,1] # since my popinfo file didnt have the location and genus in separate columns, i am doing this here
tab$loc <- str_split_fixed(tab$pop, "_", 2)[,2]

## test run
plot(tab$EV2, tab$EV1, col=as.integer(tab$info.1), col=as.integer(tab$info.2), xlab="PC2 5.6%", ylab="PC1 5.2%")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))

lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)

## real deal

pdf("wgs_test.pdf")

tab$loc <- factor(tab$loc)

ggplot(tab,aes(EV1,EV2,colour=genus, shape=loc), xlab="PC1 5.6%", ylab="PC2 5.2%")+
scale_shape_manual(values=1:nlevels(tab$loc)) +
geom_point()

dev.off()
