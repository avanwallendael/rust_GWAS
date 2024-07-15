###Basic rust data manipulation for rust GWAS manu
#Acer VanWallendael started Jul14 2021

setwd("~/Desktop/GWAS")

library(readxl)
library(tidyverse)

#READ IN DATA####
setwd("datasheets/csvs_2019")

files_2019<-list.files(pattern = ".csv")
n9<-length(files_2019)
csvs_2019<-vector("list", n9)

for(i in 1:n9){
csvs_2019[[i]]<-read.csv(files_2019[i], na = c("NA", "N/A", "na", "?"))}

names9<-substr(files_2019, 1, 14)
names(csvs_2019)<-names9

setwd("~/Desktop/GWAS/datasheets/csvs_2020")

files_2020<-list.files(pattern = ".csv")
n0<-length(files_2020)
csvs_2020<-vector("list", n0)

for(i in 1:n0){
  csvs_2020[[i]]<-read.csv(files_2020[i], na = c("NA", "N/A", "na","?"))}

names0<-substr(files_2020, 1, 14)
names(csvs_2020)<-names0

#fix CLMB typos
if(csvs_2019$CLMB_GWAS_2019[972,9]==90){
  csvs_2019$CLMB_GWAS_2019[972,9]<-9}

#C6711_CLMB_GWAS_2020
if(csvs_2020$CLMB_GWAS_2020[1236,9]==77){
  csvs_2020$CLMB_GWAS_2020[1236,9]<-7}

setwd("~/Desktop/GWAS/datasheets/csvs_2021")

files_2021<-list.files(pattern = ".csv")
n1<-length(files_2021)
csvs_2021<-vector("list", n1)

for(i in 1:n1){
  csvs_2021[[i]]<-read.csv(files_2021[i], na = c("NA", "N/A", "na","?"))}

names1<-substr(files_2021, 1, 14)
names(csvs_2021)<-names1

#Fix CLMB NAs
#View(csvs_2021$CLMB_GWAS_2021)
CLMB_GWAS_2021_M <- read_excel("~/Downloads/CLMB_GWAS_2021_Master Data_Working.xlsx")
CLMB_GWAS_2021_M$deadones<-0
CLMB_GWAS_2021_M$deadones[which(CLMB_GWAS_2021_M$GR100=="NA")]<-1
csvs_2021$CLMB_GWAS_2021$deadones<-CLMB_GWAS_2021_M$deadones[match(csvs_2021$CLMB_GWAS_2021$PLOT_GL, CLMB_GWAS_2021_M$PLOT_GL)]
clmbrust<-csvs_2021$CLMB_GWAS_2021[,4:9]
for(i in 1:6){
  clmbrust[,i][which(is.na(clmbrust[,i]))]<-0
  clmbrust[,i][which(csvs_2021$CLMB_GWAS_2021$deadones==1)]<-NA
}
csvs_2021$CLMB_GWAS_2021[,4:9]<-clmbrust

setwd("~/Desktop/GWAS/")

csvs<-c((csvs_2019), (csvs_2020), csvs_2021)

#DATA CLEANUP####
#clean out rust_xx columns
badcol<-vector("list", length(csvs))
csvs2<-csvs
for(i in 1:length(csvs)){
     badcol[[i]]<-grep("X", colnames(csvs[[i]]))
     #clip one off the end to stop from losing ones w no grep
     badcol[[i]]<-c(badcol[[i]], ncol(csvs[[i]]))
     csvs2[[i]]<-csvs[[i]][-badcol[[i]]]
     }
lapply(csvs2, colnames)
#save(csvs2, file = "cleaned_csvs_allsite3yr_mar28.rda")


#All Sites Analysis####
#get colnames
csvs2<-lapply(csvs2, as.data.frame)
cols_all<-lapply(csvs2, colnames)

cols_keep1<-lapply(cols_all, FUN=function(x){x[substr(x,1,3) %in% c("SIT","PLO")]})
cols_keep2<-lapply(cols_all, FUN=function(x){x[substr(x,1,3) %in% c("RUS")]})
lapply(cols_all, FUN=function(x){x[substr(x,1,3) %in% c("RUS")]})

csvs3<-csvs2

n<-length(cols_all)

for(i in 1:n){csvs3[[i]]<-csvs2[[i]][,colnames(csvs2[[i]]) %in% cols_keep2[[i]]]}

csvs4 <- lapply(csvs3, function(x) {
  x[] <- lapply(x, as.numeric)
  x
})

for(i in 1:n){
  csvs3[[i]]<-csvs2[[i]][,colnames(csvs2[[i]]) %in% cols_keep1[[i]]]
  csvs2[[i]]<-cbind(csvs3[[i]],csvs4[[i]])
  #rm blank lines
  csvs2[[i]]<-csvs2[[i]][which((csvs2[[i]]$SITE)!=""),]
}

lapply(csvs2, FUN=function(x){sapply(x,class)})

#check
lapply(csvs2,FUN=function(x){colSums(x[-1:-3], na.rm = T)})

#remove frmi zeros
csvs2$FRMI_GWAS_2019<-csvs2$FRMI_GWAS_2019[-4:-8]

#pivot to long form
csvs5<-lapply(csvs2, FUN=function(x){
  pivot_longer(x,  cols = -c(SITE, PLOT_GL, PLOT_LC),  names_to = "rust_DOY", values_to = "score")})

csvs6<-csvs5

namesall<-c(names9, names0, names1)

names(csvs6)<-namesall

for( i in 1:n){
  csvs6[[i]]$DOY<-as.numeric(substr(csvs6[[i]]$rust_DOY,6,8))
  csvs6[[i]]$year<-substr(namesall[i],11,14)
  csvs6[[i]]<-as.data.frame(csvs6[[i]])
}

#fix tmpl zeroes
#View(csvs6$TMPL_GWAS_2021)
TMPL_GWAS_2021_M <- read_excel("gwas_phenos/TMPL_GWAS_2021_Master Data_Working.xlsx")

deadstmpl<-(1:nrow(TMPL_GWAS_2021_M))[which(is.na(TMPL_GWAS_2021_M$FL1))]
deadidstmpl<-TMPL_GWAS_2021_M$PLOT_GL[deadstmpl]

csvs6$TMPL_GWAS_2021$score[is.na(csvs6$TMPL_GWAS_2021$score)]<-0
csvs6$TMPL_GWAS_2021$score[csvs6$TMPL_GWAS_2021$PLOT_GL %in% deadidstmpl]<-NA

#combine all
csvsall_long<-do.call(rbind, csvs6)

#check for weird scores
csvsall_long %>% filter(score>10)
csvsall_long %>% filter(score<0)

csvsall_long %>%
  group_by(SITE, year) %>%
  summarize(meanrust=mean(score, na.rm=T),
            medirust=median(score, na.rm=T),
            times=length(unique(DOY)),
            span=max(DOY)-min(DOY))

#calculate audpc for ALL#####
#need to get ~ same intervals across sites

#vector of dates for each siteyear
dates<-vector('list', n)
for(i in 1:n){dates[[i]]<-csvs6[[i]]$DOY}
dates2<-lapply(dates, unique)
dates3<-lapply(dates2, FUN=function(x){as.numeric(x)})

names(dates3)<-namesall

dates3<-dates3[-c(3,7,12)]

#plot all, choose 4 dates that span exponential growth period
lapply(csvs6, FUN=function(x){
  ggplot(x)+
    geom_line(aes(x=as.character(DOY), y=score, group=PLOT_GL))+
    geom_jitter(aes(x=as.character(DOY), y=score))
})

bestdoysall<-data.frame(BRKG_GWAS_2019=c(219,233,247,261),
CLMB_GWAS_2019=c(213,228,241,259),
FRMI_GWAS_2019=c(217,242,249,255),
KBSM_GWAS_2019=c(205,219,231,248),
KING_GWAS_2019=c(115,129,143,157),
LINC_GWAS_2019=c(218,234,242,254),
PKLE_GWAS_2019=c(104,120,135,149),
TMPL_GWAS_2019=c(107,122,137,150),
CLMB_GWAS_2020=c(203,219,233,246),
FRMI_GWAS_2020=c(233,247,261,275),
KBSM_GWAS_2020=c(217,235,256,270),
KING_GWAS_2020=c(112,126,140,154),
LINC_GWAS_2020=c(233,247,261,275),
PKLE_GWAS_2020=c(128,143,156,175),
TMPL_GWAS_2020=c(153,169,182,196),
CLMB_GWAS_2021=c(207,221,235,249),
FRMI_GWAS_2021=c(237,251,265,279),
KBSM_GWAS_2021=c(197,213,230,250),
KING_GWAS_2021=c(123,137,153,166),
LINC_GWAS_2021=c(193,204,218,235),
OVTN_GWAS_2021=c(139,153,168,183),
PKLE_GWAS_2021=c(201,216,230,244),
TMPL_GWAS_2021=c(182,200,214,229))


#rm ovtn
csvs7<-csvs6[-7]
for(i in 1:length(csvs7)){
  csvs7[[i]]<-csvs7[[i]][csvs7[[i]]$DOY %in% bestdoysall[,i],]
  print(sort(unique(csvs7[[i]]$DOY)))
}

#audpc math
#area=((sall1+sall2)/2)*timediff
trarea<-function(s1,s2,delta){((s1+s2)/2)*delta}

calc_audpc4<-function(ind_df){trarea(ind_df$score[1],ind_df$score[2],ind_df$DOY[2]-ind_df$DOY[1])+
    trarea(ind_df$score[2],ind_df$score[3],ind_df$DOY[3]-ind_df$DOY[2])+
    trarea(ind_df$score[3],ind_df$score[4],ind_df$DOY[4]-ind_df$DOY[3])}

n2<-length(csvs7)
audpc_allsites<-vector('list', n2)

for(i in 1:n2){
  #divide genos into list
  genosplit<-split(csvs7[[i]], as.factor(csvs7[[i]]$PLOT_GL))
  genosplit2<-lapply(genosplit, as.data.frame)
  #sort by DOY
  genosplit3<-lapply(genosplit2, FUN = function(x){x[order(x$DOY,decreasing = F),]})
  audpc_allsites[[i]]<-unlist(lapply(genosplit3, calc_audpc4))
}

names(audpc_allsites)<-names(csvs7)

#find audpc variance across sites & years
audpc_allsites2<-vector('list', n2)

for(i in 1:n2){
  audpc_allsites2[[i]]<-data.frame(PLOT_GL=names(audpc_allsites[[i]]),
                              score=audpc_allsites[[i]],
                              siteyear=names(csvs7)[i])
  audpc_allsites2[[i]]$UID<-paste(audpc_allsites2[[i]]$PLOT_GL,audpc_allsites2[[i]]$siteyear, sep="_")
}

csvsall_long$UID<-paste(csvsall_long$PLOT_GL,csvsall_long$SITE, "GWAS", csvsall_long$year, sep="_")

audpc_allsites_long<-do.call(rbind, audpc_allsites2)

csvs_sum<-csvsall_long %>% 
  group_by(UID) %>%
  summarize(maxrust=max(score, na.rm=T),
            medrust=median(score, na.rm=T),
            sumrust=sum(score, na.rm=T))
#errors because of NAs
csvs_sum$audpc<-audpc_allsites_long$score[match(csvs_sum$UID, audpc_allsites_long$UID)]


#Visualize AUDPC####
csvs_sum<-as.data.frame(csvs_sum)
csvs_sum$PLOT_GL<-substr(csvs_sum$UID,1,5)
csvs_sum$SITE<-substr(csvs_sum$UID,7,10)
csvs_sum$YEAR<-substr(csvs_sum$UID,17,20)

#adjust for number of sites
if(length(unique(csvs_sum$SITE))!=9){print("Revise names")}
csvs_sum$SITE<- factor(csvs_sum$SITE, levels=c("KING", "PKLE","TMPL","OVTN", "CLMB","LINC", "KBSM","FRMI","BRKG"))

ggplot(csvs_sum, aes(x=SITE, y=audpc))+
  geom_jitter(height = 1)+
  geom_boxplot(width=.3)+
  facet_grid(YEAR~.)+
  theme_classic()

csvs_sum$maxrust[which(csvs_sum$maxrust<0)]<-NA
csvs_sum$maxrust[which(csvs_sum$maxrust==77)]<-7

#add zeroes to BRKG for plotting
brkg19<-csvs_sum %>% filter(SITE=="BRKG")
brkg20<-brkg19
brkg20$audpc[-which(is.na(brkg20$audpc))]<-0
brkg21<-brkg20
brkg21$YEAR<-2021
brkg20$YEAR<-2020

csvs_sum2<-rbind(csvs_sum, brkg20, brkg21)

ggplot(csvs_sum2, aes(x=SITE, y=audpc))+
  geom_jitter(width=.3,height = 5, size=1, alpha=.5, col="grey50")+
  #geom_violin()+
  geom_boxplot(width=.3, fill="skyblue", col="black", alpha=.7)+
  facet_grid(YEAR~.)+
  labs(y="AUDPC")+
  theme_classic()
ggsave("audpc_bysiteyear_april16.pdf", height = 5, width=7)

#save(csvs_sum, file="csvs_sum_allsites_april16.rda")

