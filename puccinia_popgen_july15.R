#puccinia popgen

library(bigsnpr)

setwd("~/Desktop/puccinia_pool")

snp_readBed("filt_nonmissing_pucc2.bed.bed")

puccsnp<-snp_attach("/Users/Acer/Desktop/puccinia_pool/filt_nonmissing_pucc2.bed.rds")

#trim to good contigs
load("~/Desktop/puccinia_pool/good_novopanici_contigs.rda")
filt_contigs<-puccsnp$map$chromosome %in% good_contigs
snp_good<-snp_subset(puccsnp, ind.col = filt_contigs)
pucc_snp2<-snp_attach(snp_good)

#make pca
svd_pucc<-big_SVD(pucc_snp2$genotypes, k=10)

svd1<-as.data.frame(svd_pucc$u)
svd1$samp<-pucc_snp2$fam$sample.ID
svd1$site<-substr(svd1$samp,1,1)
svd1$col1<-rainbow(8)[match(svd1$site, unique(svd1$site))]

plot(svd_pucc$u[,1], svd_pucc$u[,2], col=svd1$col1)
text(svd_pucc$u[,1], svd_pucc$u[,2], labels=svd1$site)

plot(svd_pucc$u[,1], svd_pucc$u[,3], col=svd1$col1)
text(svd_pucc$u[,1], svd_pucc$u[,3], labels=svd1$site)
library(readxl)
DOE_GWAS_Consolidated_Cultural_Data_1_ <- read_excel("~/Downloads/DOE_GWAS_Consolidated Cultural Data(1).xlsx")
svd1$lat<-DOE_GWAS_Consolidated_Cultural_Data_1_$Latitude[match(svd1$site, substr(DOE_GWAS_Consolidated_Cultural_Data_1_$`Site Code`,1,1))]

sitecode<-data.frame(site=c("B", "L","F","M","C","S","T","P","K"), code=c("BRKG", "LINC","FRMI","KBSM","CLMB","STIL","TMPL","PKLE","KING"))

svd1$SITE<-sitecode$code[match(svd1$site, sitecode$site)]

svd1$lat<-DOE_GWAS_Consolidated_Cultural_Data_1_$Latitude[match(svd1$SITE, DOE_GWAS_Consolidated_Cultural_Data_1_$`Site Code`)]

#plot(svd1$lat, svd1$V2)

library(viridis)

svd1$site<-factor(svd1$site, levels=c("B", "M","F","L","C","S","T","P","K"))

ggplot(svd1)+
  geom_point(aes(x=lat, y=V2, fill=site), shape=21, size=3)+
  scale_fill_manual(values = viridis(9), name="Site")+
  labs(x="Latitude", y="PC2")+
  theme_classic()
ggsave("pucc_PC2_lat_oct23.pdf", height=5, width=5)  

ggplot(svd1)+
  geom_point(aes(x=V1, y=V2, fill=lat), shape=21, size=3)+
  #scale_fill_manual(values = viridis(8))+
  labs(x="PC1", y="PC2")+
  theme_classic()
ggsave("pucc_PCs_lat_oct23.pdf", height=5, width=5)  


