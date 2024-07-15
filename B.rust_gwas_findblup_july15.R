#GWAS analysis for rust
#Acer VanWallendael created Oct20 2021

library(tidyverse)
library(switchgrassGWAS)
library(bigsnpr)
library(readxl)

setwd("~/Desktop/GWAS")

#load("audpc_long.rda")
#load("core_long.rda")
#load("core_sum.rda")
#load("csvs_nov17.rda")

snp <- snp_attach("~/Downloads/switch_snps/Pvirgatum_V5_GWAS_630g_33M_tensite_twoyear.rds")

NCORES <- nb_cores()

CHRN<-snp$map$chromosome
chrkey<-data.frame(chr=unique((CHRN)), 
                   num=1:length(unique(CHRN)))
chrnum<-chrkey$num[match(CHRN, chrkey$chr)]
POS<-snp$map$physical.pos
#chrpos<-list(CHRN, chrnum, POS)
#names(chrpos)<-c("CHRN", "chrnum", "POS")
#save(chrpos, file="chrpos.rda")

#run pca
#obj.svd<-snp_autoSVD(snp$genotypes, infos.chr = chrnum, infos.pos = POS)
#save(obj.svd, file="obj.svd_tensite_oct22.rda")
load("obj.svd_tensite_oct22.rda")
snp_ids<-obj.svd$PLANT_ID
#save(snp_ids, file="snp_ids.rda")

DOE <- as.data.frame(read_excel("~/Downloads/DOE_GWAS_Master Plant List (2).xlsx"))


#calculate scaled means####
load("csvs_sum_allsites_april16.rda")

csvs_sum$PLANT_ID<-DOE$PLANT_ID[match(csvs_sum$PLOT_GL, DOE$PLOT_GL)]
csvs_sum$siteyear<-paste0(csvs_sum$SITE, "_", csvs_sum$YEAR)
sy_split_all<-split(csvs_sum, csvs_sum$siteyear)

sy_split2_all<-sy_split_all
for(i in 1:24){
  sy_split2_all[[i]]$sy_scale<-scale(sy_split2_all[[i]]$audpc)[,1]
}
csvs_sum_scale1<-do.call(rbind, sy_split2_all)

csvs_sum_scale2<-csvs_sum_scale1 %>% group_by(PLANT_ID) %>% summarize(mean_audpc=mean(sy_scale, na.rm=T), med_audpc=median(sy_scale, na.rm=T), var_audpc=var(sy_scale, na.rm=T))

#scale NS
sy_N<-sy_split2_all[c(2:10,14:16)]
css_N<-do.call(rbind, sy_N)
css_N2<-css_N%>% group_by(PLANT_ID) %>% summarize(mean_audpc=mean(sy_scale, na.rm=T), med_audpc=median(sy_scale, na.rm=T), var_audpc=var(sy_scale, na.rm=T))
save(css_N, file="css_N_april16.rda")

sy_S<-sy_split2_all[c(11:13,19:24)]
css_S<-do.call(rbind, sy_S)
css_S2<-css_S %>% group_by(PLANT_ID) %>% summarize(mean_audpc=mean(sy_scale, na.rm=T), med_audpc=median(sy_scale, na.rm=T), var_audpc=var(sy_scale, na.rm=T))
save(css_S, file="css_S_april16.rda")

#RUN SOMMER MODELS ON HPCC####
#see shell scripts for LMMs. sommer_findblup.sh B1.sommer_findblup.R

load("blups_S_april16.rda")
load("blups_N_april16.rda")

#calculate correlation between blups and scaled means
blups2<-data.frame(PLANT_ID=snp_ids)
blups2$blupsN<-blupsN[match(blups2$PLANT_ID, names(blupsN))]
blups2$blupsS<-blupsS[match(blups2$PLANT_ID, names(blupsS))]

blups2$scaled_all<-csvs_sum_scale2$mean_audpc[match(blups2$PLANT_ID, csvs_sum_scale2$PLANT_ID)]
blups2$scaled_N<-css_N2$mean_audpc[match(blups2$PLANT_ID, css_N2$PLANT_ID)]
blups2$scaled_S<-css_S2$mean_audpc[match(blups2$PLANT_ID, css_S2$PLANT_ID)]

blups2$scaled_N2<-blups2$scaled_N
blups2$scaled_S2<-blups2$scaled_S

blups2$scaled_N2[is.na(blups2$scaled_N)]<-mean(blups2$scaled_N, na.rm=T)
blups2$scaled_S2[is.na(blups2$scaled_S)]<-mean(blups2$scaled_S, na.rm=T)

#save(blups2, file="blups2_scaledNS_april19.rda")

cor.test(blups2$blupsN, blups2$blupsS)
cor.test(blups2$scaled_N, blups2$scaled_S)
cor.test(blups2$blupsN, blups2$scaled_N)
cor.test(blups2$blupsS, blups2$scaled_S)

#Run GWAS######
#blups with gxe
#load("blups2_scaledNS_april19.rda")

#prep phenos

phenos_blup3<-left_join(data.frame(PLANT_ID=snp$fam$sample.ID), blups2)
phenos_blup4<-phenos_blup3 %>% dplyr::select(PLANT_ID, blupsN)

gwas_blup_N<- pvdiv_gwas(df = phenos_blup4, type = "linear",snp = snp, covar = obj.svd, npcs = 10)
save(gwas_blup_N, file = "gwas_blup_Napril19.rda")

phenos_blup5<-phenos_blup3 %>%dplyr::select(PLANT_ID, finalblupsS)

gwas_blup_S<- pvdiv_gwas(df = phenos_blup5, type = "linear",snp = snp, covar = obj.svd, npcs = 10)
save(gwas_blup_S, file = "gwas_blup_Soct10.rda")

#plotting####
bonf<-(-log10(0.05/(nrow(gwas_blup_N))))
pvals_blup<-stats::predict(gwas_blup_N)
gwas_blup2<-cbind(gwas_blup_N,data.frame(POS=POS,CHR=CHRN, pval=pvals_blup))
gwas_blup2$marker<-paste0(gwas_blup2$CHR,gwas_blup2$POS)
gwas_trim_blup<-    gwas_blup2[which(pvals_blup<(-3)),]
gwas_trim_sigN<-    gwas_blup2[which(pvals_blup<(-bonf)),]

pvals_blupS<-stats::predict(gwas_blup_S)
gwas_blupS2<-cbind(gwas_blup_S,data.frame(POS=POS,CHR=CHRN, pval=pvals_blupS))
gwas_blupS2$marker<-paste0(gwas_blupS2$CHR,gwas_blupS2$POS)
gwas_trim_blupS<-    gwas_blupS2[which(pvals_blupS<(-3)),]
gwas_trim_sigS<-    gwas_blupS2[which(pvals_blupS<(-bonf)),]


#add previous noted regions from VanWallendael et al. 2020, 2022.
qtls<-data.frame(CHR=c("Chr03N","Chr09N"), POSmin=c(45021487,54120780),POSmax=c(44564470, 53325843),ymin=c(0,0),ymax=c(13,13))
phyllos_out<-data.frame(CHR="Chr02N", POS=60928981)

#plot####
library(patchwork)

#pdf(file="gwas_blupN_dec12.pdf", width=15)
A<-ggplot(data=gwas_trim_blup, aes(x=POS, y=(-pval)))+
  geom_rect(data = phyllos_out, aes(xmin=POS,xmax=POS+500, ymin=0,ymax=Inf),col="lightgreen", inherit.aes = F)+
  geom_rect(data = qtls, aes(xmin=POSmin, xmax=POSmax, ymin=0, ymax=Inf), col="skyblue", inherit.aes = F)+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=3, fill=CHR))+
  geom_point(aes(col=CHR), size=.5)+
  geom_hline(aes(yintercept=bonf), col="red3", linewidth=.5)+
  scale_color_manual(values = rep(c("black","grey30"),9), guide="none")+
  scale_fill_manual(values = rep(c("black","grey30"),9), guide="none")+
  facet_grid(.~CHR, scales = "free", space="free")+
  ggtitle(label = "Rust AUDPC North")+
  labs(y="-log p value", x="Genome position")+
  scale_y_continuous(breaks = c(0,5,10,15,20,25)) +
  coord_cartesian(ylim=c(0,25.5))+
  theme_classic()+
  theme(axis.text.x = element_blank())
#dev.off()

#pdf(file="gwas_blupS_dec12.pdf", width=15)
B<-ggplot(data=gwas_trim_blupS, aes(x=POS, y=(-pval)))+
  geom_rect(data = phyllos_out, aes(xmin=POS,xmax=POS+500, ymin=0,ymax=Inf),col="lightgreen", inherit.aes = F)+
  geom_rect(data = qtls, aes(xmin=POSmin, xmax=POSmax, ymin=0, ymax=Inf), col="skyblue", inherit.aes = F)+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=3, fill=CHR))+
  geom_point(aes(col=CHR), size=.5)+
  geom_hline(aes(yintercept=bonf), col="red3", linewidth=.5)+
  scale_color_manual(values = rep(c("black","grey30"),9), guide="none")+
  scale_fill_manual(values = rep(c("black","grey30"),9), guide="none")+
  facet_grid(.~CHR, scales = "free", space="free")+
  ggtitle(label = "Rust AUDPC South")+
  labs(y="-log p value", x="Genome position")+
  scale_y_continuous(breaks = c(0,5,10,15,20,25)) +
  coord_cartesian(ylim=c(0,25.5))+
  theme_classic()+
  theme(axis.text.x = element_blank())
#dev.off()

pdf(file="gwas_blupNS_dec18.pdf", height= 8,width=12)
A/B+plot_annotation(tag_levels = 'A')
dev.off()

#exp for figuremaking
save(gwas_trim_blup, phyllos_out, qtls, gwas_trim_blupS, file = "GWASfigure4_files.rda")

gwas_trim2_blupN<- gwas_trim_blup %>% filter(pval<(-bonf))

northpeaks<-gwas_trim2_blupN %>% filter(pval<(-16.5))
write_csv(northpeaks, "northpeaks_oct12.csv")

gwas_trim2_blupS<- gwas_trim_blupS %>% filter(pval<(-bonf))

southpeaks<-gwas_trim2_blupS %>% filter(pval<(-10.5))
write_csv(southpeaks, "southpeaks_oct12.csv")

#are outliers correlated between NS?####
#take sigN and find corresponding S
gwas_trim_sigN$pvalS<-gwas_blupS2$pval[match(gwas_trim_sigN$marker, gwas_blupS2$marker)]

ggplot(gwas_trim_sigN)+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(x=-pval, y=-pvalS, fill=CHR), shape=21)+
  scale_fill_manual(values=rainbow(18))+
  #facet_grid(.~CHR)+
  theme_classic()

cor.test(gwas_trim_sigN$pval, gwas_trim_sigN$pvalS) 

length(gwas_trim_sigN$pvalS[which(gwas_trim_sigN$pvalS<(-bonf))])
length(gwas_trim_sigN$pval)

gwas_trim_sigS$pvalN<-gwas_blup2$pval[match(gwas_trim_sigS$marker, gwas_blup2$marker)]

ggplot(gwas_trim_sigS)+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(x=-pval, y=-pvalN, fill=CHR), shape=21)+
  scale_fill_manual(values=rainbow(18))+
  #facet_grid(.~CHR)+
  theme_classic()

cor.test(gwas_trim_sigS$pval, gwas_trim_sigS$pvalN) 
length(gwas_trim_sigS$pvalN[which(gwas_trim_sigS$pvalN>(-bonf))])

length(gwas_trim_sigS$pvalN[which(gwas_trim_sigS$pvalN<(-bonf))])
length(gwas_trim_sigS$pval)

#blups by population####
pca1<-as.data.frame(obj.svd$u)
pca1$PLANT_ID<-phenos_blup3$PLANT_ID

blups2$PC1<-pca1$V1[match(blups2$PLANT_ID, pca1$PLANT_ID)]
blups2$PC2<-pca1$V2[match(blups2$PLANT_ID, pca1$PLANT_ID)]
blups2$PC3<-pca1$V3[match(blups2$PLANT_ID, pca1$PLANT_ID)]

blups2$rustN<-NA
blups2$rustS<-NA

blups2$rustN[which(blups2$finalblupsN>0)]<-"Positive"
blups2$rustN[which(blups2$finalblupsN<=0)]<-"Negative"

A1<-blups2 %>% filter(!is.na(rustN)) %>%
  ggplot() + 
  geom_point(aes(x=PC1, y=PC2, fill=rustN), shape=21, size=3)+
  scale_fill_manual(values=c("black", "skyblue"), name="Rust BLUP")+
  annotate(geom="text", x=c(-.04, 0,-0.01,.035), y=c(-.03,0.001,.075,-.05), label=c("Atlantic", "Gulf", "Gulf inland", "Midwest"))+
  theme_classic()
#ggsave("PCA_byrustN_oct13.pdf", height = 5, width = 5)

blups2$rustS[which(blups2$finalblupsS>0)]<-"Positive"
blups2$rustS[which(blups2$finalblupsS<=0)]<-"Negative"

B1<-blups2 %>% filter(!is.na(rustS)) %>%
  ggplot() + 
  geom_point(aes(x=PC1, y=PC2, fill=rustS), shape=21, size=3)+
  scale_fill_manual(values=c("black", "firebrick3"), name="Rust BLUP")+
  annotate(geom="text", x=c(-.04, 0,-0.01,.035), y=c(-.03,0.001,.075,-.05), label=c("Atlantic", "Gulf", "Gulf inland", "Midwest"))+
  theme_classic()
#ggsave("PCA_byrustS_oct13.pdf", height = 5, width = 5)

#stats
#check normality
LAtest1<-aov(value~pop*Region, data=blups3)
TukeyHSD(LAtest1)

shapiro.test(blups3$value)

blups3$popregion<-paste0(blups3$pop, blups3$Region)
dunntest<-dunn.test(blups3$value, blups3$popregion, method="bonferroni", kw=T)

data.frame(comp=dunntest$comparisons, pval=dunntest$P.adjusted)

#pca plot####
PVDIV <- read_excel("~/Downloads/PVDIV_Master Metadata File_5-8-2019.xlsx")

blups2$lat<-as.numeric(PVDIV$LATITUDE[match(blups2$PLANT_ID, PVDIV$PLANT_ID)])
blups2$lon<-as.numeric(PVDIV$LONGITUDE[match(blups2$PLANT_ID, PVDIV$PLANT_ID)])
#install your packages
#install.packages("maps")

blups2$pop<-"Atlantic"
blups2$pop[which(blups2$PC1>.04)]<-"Midwest"
blups2$pop[which(blups2$PC2>.03)]<-"Gulf"

blups3 <- blups2 %>% pivot_longer(cols=c(finalblupsN, finalblupsS))

blups3$Region<-blups3$name
blups3$Region[which(blups3$name=="finalblupsN")]<-"Score_North"
blups3$Region[which(blups3$name=="finalblupsS")]<-"Score_South"

C1<-blups3 %>% 
  #filter(pop != "Atlantic") %>%
  ggplot()+
  geom_boxplot(aes(x=pop, y=value, fill=Region), col="black", alpha=.8)+
  scale_fill_manual(values=c("skyblue","firebrick3" ), name="Region")+
  annotate(geom = "text", x=c("Gulf", "Midwest"), y=140, label="*", size=8)+
  coord_cartesian(ylim=c(-30,145))+
  labs(x="Population", y="Rust Genotypic BLUP")+
  theme_classic()
#ggsave("audpc_bypop_oct13.pdf", height = 5, width=7)

(A1+B1)/C1 + plot_annotation(tag_levels = 'A')
ggsave("GWASfigure3_dec18.pdf", width=9, height = 7)

#exp data for fig
save(blups2, blups3, file="GWASfigure3_files.rda")

#map####
#load libraries
library(maps)
library(sp)

map1<-map_data("usa")
states_map <- map_data("state")

blups2$lon[which(blups2$lon>0)]<-NA
blups2$lat[which(blups2$lat>50)]<-NA

#first just plot points
ggplot(map1) +
  geom_polygon(aes(long, lat, group = group), color = "black", fill="white")+
  geom_point(data=blups2, aes(x=lon, y=lat, fill=rustS), shape=21)+
  scale_fill_manual(values=c("black", "skyblue"), name="Scaled Rust Score")+
  labs(x="Longitude", y="Latitude", title="Susceptible to northern rust")+
  theme_classic()
ggsave("map_south_susceptible.pdf", height=5, width = 9)

ggplot(map1) +
  geom_polygon(aes(long, lat, group = group), color = "black", fill="white")+
  geom_point(data=blups2, aes(x=lon, y=lat, fill=rustN), shape=21)+
  scale_fill_manual(values=c("black", "firebrick"), name="Scaled Rust Score")+
  labs(x="Longitude", y="Latitude", title="Susceptible to southern rust")+
  theme_classic()
ggsave("map_north_susceptible.pdf", height=5, width = 9)



phenos_blup3$Library.ID<-PVDIV$LIBRARY[match(phenos_blup3$PLANT_ID, PVDIV$PLANT_ID)]

subpops <- read.csv("~/Downloads/41586_2020_3127_MOESM18_ESM (2).xlsx - fig_3a.csv")
ecos <- read.csv("~/Downloads/41586_2020_3127_MOESM17_ESM (3).xlsx - fig_2a.csv")

phenos_blup3$pop1<-PVDIV$SUBPOP_SNP[match(phenos_blup3$PLANT_ID, PVDIV$PLANT_ID)]
popkey<-data.frame(pops=unique(phenos_blup3$pop1),
                   newpops=c("Atlantic", "Texas", "Midwest", "Atlantic", "Texas", "Gulf", "Gulf", "Midwest", NA))
phenos_blup3$Subpopulations<-popkey$newpops[match(phenos_blup3$pop1, popkey$pops)]

phenos_blup4<-left_join(phenos_blup3, subpops)

phenos_blup4$Subpopulation<-NA

phenos_blup4$Ecotype<-ecos$ecotype[match(phenos_blup4$PLANT_ID, ecos$plant_ID)]
save(phenos_blup4, file="phenos_blup4_april20.rda")

phenos_blup5<-phenos_blup4 %>%
  filter(!is.na(Ecotype)&!is.na(Subpopulations))

phenos_blup5$Subpopulations[which(phenos_blup5$Subpopulations=="Gulf")]<-"Gulf-coast"

phenos_blup5$Subpopulations[which(phenos_blup5$Subpopulations=="Texas")]<-"Gulf-Texas"

phenos_blup5$Subpopulations<-factor(phenos_blup5$Subpopulations,
                                    levels=c("Gulf-coast", "Gulf-Texas", "Atlantic", "Midwest" ))

ggplot(phenos_blup5, aes(y=finalblupsS))+
  geom_boxplot(aes(x=Subpopulations))+
  #facet_grid(.~Subpopulations, space="free", scales="free")+
  labs(y="Rust BLUP")+
  theme_classic(base_size = 10)

#ggsave("Eco_pop_blup.png", width= 7, height = 5)
#ggsave("Eco_pop_blup.pdf", width= 5, height = 5)

#check overlap regions####
#qtls
#pull qtl region
threen_gwas<-gwas_trim_blup %>% filter(CHR=="Chr03N")
threen_qtl<-threen %>% filter(POS<(qtls$POSmin[1]+5000000) & POS> (qtls$POSmax[1]-5000000))

ggplot(threen_qtl, aes(x=POS, y=-pval))+
  geom_rect(aes(xmin=qtls$POSmax[1], xmax=qtls$POSmin[1], ymin=-Inf, ymax=Inf), fill="skyblue", inherit.aes = F)+
  geom_point(size=.7)+
  theme_classic()

#add in qtl plot
slod_newmap_north <- read.csv("~/Downloads/slod_newmap_north.csv")
#get new pos
markers<-strsplit(slod_newmap_north$marker, split = "_")
markpos<-sapply(markers, FUN=function(x){as.numeric(x[2])})

slod_newmap_north$POS<-markpos
threen_qtls<-slod_newmap_north %>% filter(chr=="3N")
threen_qtls$lev<-1:nrow(threen_qtls)
plot(threen_qtls$POS~threen_qtls$lev)
modpos<-lm(slod_newmap_north$POS~slod_newmap_north$pos)
slod_newmap_north$POS2<-914494*slod_newmap_north$pos+4115371
slod_peak<-slod_newmap_north %>% filter(POS<(qtls$POSmin[1]+5000000) & POS> (qtls$POSmax[1]-5000000))
library(patchwork)



load("~/Desktop/geteffects_oct19 copy.RData")

threen_comb<-threen %>% filter(site %in% c("CLMB", "KBSM", "LINC", "MNHT", "CLMB17", "KBSM17", "MNHT17")) 
#%>% group_by(chrname)%>%
# summarise(POS=mean(pos), mlodsum=sum(mlod))

#use data from newphyt supp

a<-ggplot(threen_comb, aes(x=pos, y=mlod))+
  geom_line(aes(col=site))+
  geom_hline(aes(yintercept=22))+
  theme_classic()

#use data from newphyt supp
qtl_peaks <- read.csv("~/Downloads/Untitled spreadsheet - Sheet1 (5).csv")
peaks_threen<-qtl_peaks %>% filter(Chromosome=="3N")
peakspos<-strsplit(peaks_threen$Marker, split = "_")
peakspos2<-sapply(peakspos, FUN=function(x){as.numeric(x[2])})
peaks_threen$POS<-peakspos2
peaks_threen$siteyear<-paste(peaks_threen$Site, peaks_threen$Year, sep="_")

set.seed(22356)

ggplot(threen_gwas, aes(x=POS, y=-pval))+
  geom_vline(data=peaks_threen,aes(xintercept=POS), col="red2")+
  geom_text(data=peaks_threen, aes(label=siteyear, x=POS, y=sample((150:200)/10, 10, replace=F)), col="firebrick")+
  geom_point(size=.7)+
  theme_classic()
ggsave("threen_qtl_gwas.pdf", height=5, width=7)

#add in genes
peak_3n_genes <- read.csv("~/Downloads/peak_3n - Sheet1 (1).csv")

set.seed(223453)
ggplot(threen_gwas, aes(x=POS, y=-pval))+
  geom_vline(data=peaks_threen,aes(xintercept=POS), col="red2")+
  geom_text(data=peaks_threen, aes(label=siteyear, x=POS, y=sample((150:200)/10, 10, replace=F)), col="firebrick")+
  geom_segment(data=peak_3n_genes, aes(x=POS_start, xend=POS_end, y=16, yend=16), col="blue", 
               size=5,arrow = arrow(length = unit(0.2, "cm")), lineend="butt", linejoin="mitre",inherit.aes = F)+
  geom_point(size=1)+
  coord_cartesian(xlim=c(44400000, 45000000))+
  theme_classic()
ggsave("threen_qtl_genes_zoom.pdf", height=5, width=7)

#QQ####
pdf("snp_qqN_dec1.pdf", height = 10, width = 10)
snp_qq(gwas_blup_N)
dev.off()

pdf("snp_qqS_dec1.pdf", height = 10, width = 10)
snp_qq(gwas_blup_S)
dev.off()

#get North-South correlation
gwas_trimboth<-cbind(gwas_blup2, gwas_blupS2)
gwas_trimboth$combP<-gwas_trimboth[,6] + gwas_trimboth[,13]
gwas_trimboth2<-gwas_trimboth[gwas_trimboth$combP<(-1),]

pdf("gwas_corr_NS.pdf", height=10, width=10)
plot(-gwas_trimboth2[,6], -gwas_trimboth2[,13])
dev.off()

#correlate blups####
plot(phenos_blup3$finalblupsN, phenos_blup3$finalblupsS)
cor.test(phenos_blup3$finalblupsN, phenos_blup3$finalblupsS)
hist(phenos_blup3$finalblupsN)
hist(phenos_blup3$finalblupsS)
