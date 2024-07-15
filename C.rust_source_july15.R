#rust infection by collection location
#started oct 27 2021 by Acer VanWallendael
#data from B
setwd("~/Desktop/GWAS")

library(readxl)
library(tidyverse)

load("obj.svd_tensite_oct22.rda")
csvs_sum_scale1 <- read.csv("~/Desktop/GWAS/csvs_sum_scale1_april19.csv")
DOE <- as.data.frame(read_excel("~/Downloads/DOE_GWAS_Master Plant List (2).xlsx"))

unique(csvs_sum_scale1$SITE)

#max rust infection by ecotype at different sites
#need ecotype
PVDIV<-read_excel("~/Downloads/PVDIV_Master Metadata File_5-8-2019.xlsx")
new_eco <- read_excel("~/Downloads/41586_2020_3127_MOESM17_ESM.xlsx")

csvs_sum_scale1$PLANT_ID[which(substr(csvs_sum_scale1$PLANT_ID,1,4)=="AP13")]<-"AP13"
csvs_sum_scale1$Ecotype<-new_eco$ecotype[match(csvs_sum_scale1$PLANT_ID, new_eco$plant_ID)]
csvs_sum_scale1$Subpopulation<-PVDIV$SUBPOP_SNP[match(csvs_sum_scale1$PLANT_ID, PVDIV$PLANT_ID)]
subpop_key<-data.frame(old=unique(csvs_sum_scale1$Subpopulation), new=c("Gulf Inland", "Gulf", "Gulf Inland", "Midwest","Atlantic","Atlantic", "Octoploid", "Gulf","Midwest"))
csvs_sum_scale1$Subpopulations<-subpop_key$new[match(csvs_sum_scale1$Subpopulation, subpop_key$old)]
csvs_sum_scale2 <- csvs_sum_scale1%>% group_by(PLANT_ID, Subpopulations) %>% summarize(mean_scaled=mean(sy_scale, na.rm=T))

AUDPC_subpop<-ggplot(csvs_sum_scale2, aes(x=Subpopulations, y=mean_scaled))+
  geom_jitter(width=.3, alpha=.2)+
  geom_boxplot(width=.3, alpha=.8,fill="skyblue")+labs(y="Scaled  AUDPC")+
  theme_classic()
save(AUDPC_subpop, file="AUDPC_subpop_fig2b.rda")
ggsave(paste0("rust_source_pop_",substr(date(),5,7),substr(date(),9,10),".pdf"), height=4, width = 3.5)

#add zeroes to brookings for 2020, 2021
csvs_sum_scale1$audpc[which(csvs_sum_scale1$SITE=="BRKG"&csvs_sum_scale1$YEAR!=2019)]
csvs_sum_scale1$SITE<-factor(csvs_sum_scale1$SITE, levels = c("KING", "PKLE", "TMPL", "OVTN", "CLMB", "LINC","FRMI", "KBSM",  "BRKG"))

brkg_dum<-csvs_sum_scale1[csvs_sum_scale1$SITE=="BRKG",]
brkg_dum2<-brkg_dum
brkg_dum$YEAR<-2020
brkg_dum2$YEAR<-2021
brkg_dum$audpc<-0
brkg_dum2$audpc<-0

csvs_sum_scale3<-rbind(csvs_sum_scale1, brkg_dum, brkg_dum2)

ggplot(csvs_sum_scale3, aes(x=SITE, y=audpc))+
  geom_jitter(width=.3, alpha=.5)+
  geom_boxplot(width=.3, fill="skyblue")+labs(y="AUDPC")+
  facet_grid(YEAR~.)+
  theme_classic()

#export data for fig
save(csvs_sum_scale2, csvs_sum_scale3, file="GWAS_figure2ab_files.rda")
