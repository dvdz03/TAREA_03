########################################################
############## PREPARACIÓN DE ENTORNO ##################

#instalar paquetes
install.packages("ggplot2")
install.packages("vegan")
install.packages("dplyr")
install.packages("BiocManager")
BiocManager::install("phyloseq",force=TRUE)
BiocManager::install("microbiome",force=TRUE)
install.packages("cluster")
install.packages("foreign")
install.packages("MASS")
install.packages("Matrix")

#cargar librerías
library(phyloseq)
library(ggplot2)
library(vegan)
library(microbiome)
library(dplyr)
library(data.table)
#objeto phyloseq
data("dietswap", package="microbiome")
ps<-dietswap
ps
sample_data(ps)#esto es para ver los metadatos
head(sample_data(ps))

readcount<- data.table(as(sample_data(data.prune), "data.frame"),
                       TotalReads=sample_sums(data.prune),
                       keep.rownames=TRUE)
readcount
setnames(readcount,"rn","SampleID")
ggplot(readcount,aes(TotalReads))+

#curva rarefacción
data.prune=prune_taxa(taxa_sums(ps)>1,ps)
raref<- otu_table(data.prune)
raref<- as.data.frame(t(raref))
sample_names<-rownames(raref)
curaref<- rarecurve(raref, step = 10000)

#diversidad alfa
plot_richness(data.prune, x="group", measures=c("Shannon","Simpson"))

