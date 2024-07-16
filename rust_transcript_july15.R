#compare candidate genes RNA

setwd("~/Desktop/GWAS")

library(tidyverse)
#BiocManager::install("DESeq2")
library(DESeq2)

#link lib to genoxtreatment
meta1 <- read.csv("md4geno3site4Acer.csv")
#raw counts
counts1 <- read.csv("cnts4geno3site4Acer.csv")

allhitsout6_NS <- read.csv("~/Desktop/GWAS/allhitsout6_NS_april23.csv")

format_genes<-function(allhitsout6){
  topgenes<-allhitsout6
  
  #sub to focal
  foc<-counts1[unique(topgenes$GeneID),]
  #add ids
  foc$id<-rownames(foc)
  #rm NAs
  foc2<-foc[-which(is.na(foc$BUWGB)),]
  
  rownames(meta1)<-meta1$library
  meta1$eco<-meta1$geno
  meta1$eco[meta1$eco %in% c("AP13", "WBC")]<-"Lowland"
  meta1$eco[meta1$eco %in% c("DAC", "VS16")]<-"Upland"

    #obj construction
  meta1$sizeFactor<-NULL
  foc2$id<-NULL
  ds_test<-DESeqDataSetFromMatrix(foc2, meta1, design = ~eco)
  #ds_test<-DESeqDataSetFromMatrix(foc4, meta1, design = ~eco)
  return(ds_test)
}


allhitsout6_north<-allhitsout6_NS %>% filter(Region=="N")
allhitsout6_south<-allhitsout6_NS %>% filter(Region=="S")

dstest_north<-format_genes(allhitsout6_north)
dstest_south<-format_genes(allhitsout6_south)

smallestGroupSize <- 3
keepN <- rowSums(counts(dstest_north) >= 10) >= smallestGroupSize
keepS <- rowSums(counts(dstest_south) >= 10) >= smallestGroupSize
ddsN <- dstest_north[keepN,]
ddsS <- dstest_south[keepS,]

wd.grp_north <- DESeq(object=ddsN,
  #"Wald" or "LRT
  test="Wald",
   full = ~ eco,
  quiet = F,
  parallel = F)

resnorth<-as.data.frame(results(wd.grp_north))
nrow(resnorth[resnorth$padj<0.05,])

wd.grp_south <- DESeq(object=ddsS,test="Wald",full = ~ eco,quiet = F,parallel = F)

ressouth<-as.data.frame(results(wd.grp_south))
nrow(ressouth[ressouth$padj<0.05,])

#put deseq test onto allhits
allhitsout6_north$DEpadj<-resnorth$padj[match( allhitsout6_north$GeneID, rownames(resnorth))]
allhitsout6_south$DEpadj<-ressouth$padj[match( allhitsout6_south$GeneID, rownames(ressouth))]

write_csv(allhitsout6_north, "allhitsout6_north_wald_may16.csv")
write_csv(allhitsout6_south, "allhitsout6_south_wald_may16.csv")

write_csv(ressouth, "waldtest_south.csv")
write_csv(resnorth, "waldtest_north.csv")

#plotting results####
#plot top 25 N
aho6n<-allhitsout6_north[!duplicated(allhitsout6_north$GeneID),]
length(aho6n$GeneID[aho6n$DEpadj<0.05])
topDE_N<-aho6n %>%  slice_min(DEpadj, n=20) %>% select(DEpadj, Predicted.Function, GeneID, log10pval)
topDE_N$GWAS.log10p<-topDE_N$log10pval*(-1)
topDE_N$DE.log10p<-(-1)*log10(topDE_N$DEpadj)
topDE_N2<-topDE_N %>% select(DE.log10p, Predicted.Function, GeneID, GWAS.log10p)


aho6s<-allhitsout6_south[!duplicated(allhitsout6_south$GeneID),]
length(aho6s$GeneID[aho6s$DEpadj<0.05])
topDE_S<-aho6s %>%  slice_min(DEpadj, n=20) %>% select(DEpadj, Predicted.Function, GeneID, log10pval)
topDE_S$GWAS.log10p<-topDE_S$log10pval*(-1)
topDE_S$DE.log10p<-(-1)*log10(topDE_S$DEpadj)
topDE_S2<-topDE_S %>% select(DE.log10p, Predicted.Function, GeneID, GWAS.log10p)

library(kableExtra)
topDE_N2 %>%
  kbl(caption = "Top DE genes North") %>%
  kable_classic(full_width = F, html_font = "Cambria")
#Use zoom export

topDE_S2 %>%
  kbl(caption = "Top DE genes South") %>%
  kable_classic(full_width = F, html_font = "Cambria")



#plotting####
#use fpkm
resistance<-as.vector(t(counts1["Pavir.1KG382115",]))
names(resistance)<-rownames(t(counts1["Pavir.1KG382115",]))
meta2<-as.data.frame(meta1)
meta2$resist2N<-resistance[match(meta2$library, names(resistance))]

boxplot(meta2$resist2N~meta2$geno)

#get all DE genes
ngenes<-aho6n$GeneID
sgenes<-aho6s$GeneID

#sub counts1 to selected
meta2<-as.data.frame(meta1)

aho6n$GeneID[which(substr(aho6n$Predicted.Function,1,4)=="Terp")]
selgenes<-c("Pavir.1KG382200", "Pavir.1KG382110", "Pavir.1KG382115")

selcounts<-as.data.frame(t(counts1[selgenes,]))
selcounts$library<-rownames(t(counts1["Pavir.1KG382115",]))
meta3<-left_join(meta2, selcounts)
meta4<-pivot_longer(meta3,cols = c(Pavir.1KG382200, Pavir.1KG382110, Pavir.1KG382115))

colnames(meta4)[13]<-"Count"
colnames(meta4)[8]<-"Genotype"
colnames(meta4)[7]<-"Site"

meta4$Site[meta4$Site=="KBSM"]<-"Hickory Corners,\n MI"
meta4$Site[meta4$Site=="CLMB"]<-"Columbia, MO"
meta4$Site[meta4$Site=="PKLE"]<-"Austin, TX"


#meta4$Site<-factor(meta4$Site, levels = c("PKLE", "CLMB", "KBSM"))
fpkm_sel<-as.data.frame(t(fpkm.all[selgenes,]))
fpkm_sel$lib<-rownames(fpkm_sel)
fpkm1<-pivot_longer(fpkm_sel, cols = 1:3)
fpkm1$libgene<-paste0(fpkm1$lib, fpkm1$name)
meta4$libgene<-paste0(meta4$library, meta4$name)

meta4$fpkm<-fpkm1$value[match(meta4$libgene, fpkm1$libgene)]

ggplot(meta4, aes(x=Site, y=(fpkm)))+
  geom_boxplot(aes(fill=Genotype))+
  facet_grid(name~., scales = "free")+
  scale_fill_manual(values = c("firebrick", "blue3", "skyblue", "orange3" ))+
  labs(y="FPKM (Fragments per kilobase per million reads)")+
  theme_classic()
ggsave("terpenoids_transcript_may8.pdf", height=5, width = 5)


#resistance
aho6s$GeneID[which(substr(aho6s$Predicted.Function,1,4)%in%c("Leuc","EF-T"))]
selgenes2<-c("Pavir.9KG462795", "Pavir.3KG551700","Pavir.1KG445800")

selcounts2<-as.data.frame(t(counts1[selgenes2,]))
selcounts2$library<-rownames(t(counts1["Pavir.1KG382115",]))
meta32<-left_join(meta2, selcounts2)
meta42<-pivot_longer(meta32,cols = c(Pavir.9KG462795, Pavir.3KG551700,Pavir.1KG445800))

colnames(meta42)[13]<-"Count"
colnames(meta42)[8]<-"Genotype"
colnames(meta42)[7]<-"Site"

meta42$Site[meta42$Site=="KBSM"]<-"Hickory Corners,\n MI"
meta42$Site[meta42$Site=="CLMB"]<-"Columbia, MO"
meta42$Site[meta42$Site=="PKLE"]<-"Austin, TX"

ggplot(meta42, aes(x=Site, y=Count))+
  #geom_boxplot()+
  geom_boxplot(aes(fill=Genotype))+
  facet_grid(name~., scales = "free")+
  scale_fill_manual(values = c("firebrick", "blue3", "skyblue", "orange3" ))+
  theme_classic()
ggsave("receptors_transcript_april25.pdf", height=5, width = 5)
