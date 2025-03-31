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

ggplot(readcount,aes(TotalReads))+geom_histogram()+ggtitle("secuencias")
head(readcount[order(readcount$TotalReads), c("SampleID","TotalReads")])
#curva rarefacción
data.prune=prune_taxa(taxa_sums(ps)>1,ps)
raref<- otu_table(data.prune)
raref<- as.data.frame(t(raref))
sample_names<-rownames(raref)
curaref<- rarecurve(raref, step = 10000)
curaref<- rarecurve(raref[1:20],step=10000, label=FALSE)

#diversidad alfa
head(readcount)
plot_richness(data.prune,x="group")
plot_richness(data.prune,x="timepoint")
plot_richness(data.prune,x="bmi_group")
plot_richness(data.prune,x="TotalReads")
plot_richness(data.prune, x="nationality",measures = c("Shannon","Simpson","Observed"))
plot_richness(data.prune,x="sex",measures=c("Shannon","Simpson","Observed"))
plot_richness(data.prune, x="sex",measures = c("Shannon","Simpson","Observed"))+geom_boxplot()

#gráfica así buena
plot_richness(data.prune,x="sex",color="sex",measure=c("Observed","Shannon","Simpson"))+
  geom_boxplot()+theme_bw()+ggtitle("Gráfica")+
  theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.x=element_text(angle=90,hjust=1))


#gráficas de RANK-ABUNDANCIA
psTopNOTU<-names(sort(taxa_sums(data.prune), TRUE)[1:10])
pstop.prune<-prune_taxa(psTopNOTU, data.prune)
plot_bar(pstop.prune)

plot_bar(data.prune,x="Phylum", y="Abundance")

#Gráficas apiladas
plot_bar(data.prune,x="sex",y="Abundance",fill="Phylum")+geom_bar(aes(color=Phylum, fill=Phylum),stat="identity",position="stack")
plot_bar(data.prune,x="nationality",y="Abundance",fill="Phylum")+geom_bar(aes(color=Phylum, fill=Phylum),stat="identity",position="stack")
plot_bar(data.prune,x="Phylum",y="Abundance",fill="Phylum")+geom_bar(aes(color=Phylum, fill=Phylum),stat="identity",position="stack")#esto es para que la gráfica de abundancias se vea más bonita, así con colores y así

#Exportar los resultados
data("GlobalPatterns")
gp<-GlobalPatterns
gp
otu_table(gp)
tax_table(gp)


#preprocesamiento
filtrado<- prune_taxa(taxa_sums(gp)>=5,gp)
familia<- tax_glom(filtrado,taxrank="Family")
rank_names(gp)
rank_names(filtrado)
familia[ ,c("Soil","Feces","Skin")]
subset(familia,Soil)

