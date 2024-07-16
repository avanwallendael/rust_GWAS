#puccinia popoolation analysis

setwd("~/Desktop/puccinia_pool")
library(tidyverse)

#loop imports
files<-list.files(pattern=".pi")
n<-length(files)

load("~/Desktop/puccinia_pool/novocor_blast2.rda")
chr_contigs<-novocor_blast2$V1[!is.na(novocor_blast2$Chr)]
novocor_blast2[novocor_blast2$V1=="SMGU01004352.1",]

pis_test<-vector("list", n)
pis_test2<-vector("list", n)

for(i in 1:n){
  pis_test[[i]]<- read.delim(file = files[i], sep =  "\t", 
                                                 stringsAsFactors = F,header = F)
  pis_test[[i]]$pos_dummy<-1:nrow(pis_test[[i]])
  colnames(pis_test[[i]])<-c("locus", "pos", "alleles", "pi", "idk", "pos_dummy")
  pis_test[[i]]$contig_cols<-contig_cols[match(pis_test[[i]]$locus, contigs)]
  # try removing bad contigs
  #pis_test2[[i]]<-pis_test[[i]][(pis_test[[i]]$locus %in% good_contigs),]
  #rm contigs not in Chr
  pis_test2[[i]]<-pis_test[[i]][(pis_test[[i]]$locus %in% chr_contigs),]
  }

names(pis_test)<-substr(files,1,5)
names(pis_test2)<-substr(files,1,5)


#get some summary stats####



#pool across NS
#find shared regions across sites
setwd("~/Desktop/puccinia_pool")
load("~/Desktop/puccinia_pool/puccinia_pool_may26.RData")
library(tidyverse)
for(i in 1:9){
  pis_test[[i]]$locpos<-paste0(pis_test[[i]]$locus,pis_test[[i]]$pos)
  pis_test[[i]]$site<-names(pis_test)[i]
}

pis_long<-do.call(rbind, pis_test)

pis_long$NS<-"N"
pis_long$NS[pis_long$site %in% c("TMPL","PKLE","KING","STIL", "CLMB")]<-"S"

pis_long2<-pis_long %>% group_by(NS, locus, pos, locpos) %>% summarize(meanpi=mean(pi,na.rm=T)) %>% arrange(NS)

poskey<-unique(pis_long$locpos)
poskey1<-data.frame(posall=sort(poskey),
                    posdummyall=1:length(poskey))

pis_long2$posdummy<-poskey1$posdummyall[match(pis_long2$locpos, poskey1$posall)]

#import alignment w triticina
novopanici_to_trit2 <- read.delim("~/Desktop/puccinia_pool/novopanici_to_trit2.blast", header=FALSE)
colnames(novopanici_to_trit2)<-c("novocontig", "tritchr", "score", "start", "end")

#sub to one best hit per contig
novotrit1<-novopanici_to_trit2 %>% group_by(novocontig) %>% arrange(desc(score)) %>% slice_head(n=1)
#how many contigs below 1000?
nrow(novotrit1 %>% filter(score<100))
hist(log(novotrit1$score))
novotrit1$length<-abs(novotrit1$end-novotrit1$start)
plot(novotrit1$score~novotrit1$length)
#novotrit1 %>% filter(novocontig=="SMGU01004352.1")
#outlier
novotrit1 %>% filter(score>8000)

chrkey<-data.frame(cmchr=sort(unique(novotrit1$tritchr)),
                   tritchr=c(paste0(rep("Chr0",18),c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9),rep(c("A","B"),9)), "JARPQD010000021.1", "JARPQD010000030.1","JARPQD010000031.1"))
novotrit1$Chr<-chrkey$tritchr[match(novotrit1$tritchr, chrkey$cmchr)]

pis_long2$Chr<-novotrit1$Chr[match(pis_long2$locus, novotrit1$novocontig)]
pis_long2$mapscore_trit<-novotrit1$score[match(pis_long2$locus, novotrit1$novocontig)]
pis_long2$end_trit<-novotrit1$end[match(pis_long2$locus, novotrit1$novocontig)]

pis_long3a<-as.data.frame(pis_long2 %>% filter(!is.na(Chr)))

#split by chr, sort by end, give posdummy and posreal, rejoin

pi_split<-split(pis_long3a, pis_long3a$NS)
pis_split2<-lapply(pi_split, FUN=function(x){split(x, x$Chr)})

pisn<-pis_split2[[1]]
pis_s<-pis_split2[[2]]

for(i in 1:length(pisn)){
  pisn[[i]]<-pisn[[i]] %>% arrange(end_trit)
  pisn[[i]]$posdummy1<-1:nrow(pisn[[i]])
  pisn[[i]]$posreal<-pisn[[i]]$end_trit-pisn[[i]]$pos
}

for(i in 1:length(pis_s)){
  pis_s[[i]]<-pis_s[[i]] %>% arrange(end_trit)
  pis_s[[i]]$posdummy1<-1:nrow(pis_s[[i]])
  pis_s[[i]]$posreal<-pis_s[[i]]$end_trit-pis_s[[i]]$pos
}

#rejoin
pis_long4a<-rbind(do.call(rbind, pisn), do.call(rbind, pis_s))
pis_long4a$Subgenome<-substr(pis_long4a$Chr,6,6)
library(viridis)

pis_long4a %>%
ggplot( aes(x=posreal, y=meanpi))+
  geom_point(aes(col=Subgenome))+
  #scale_color_manual(values=cc2)+
  facet_grid(NS~Chr, space="free", scales="free_x" )+
  labs(x="Approximate position", y="Mean Pi")+
  #coord_cartesian(ylim = c(0,.1))+
  scale_color_manual(values=c('black', 'grey70', "skyblue"), guide="none")+
  theme_classic()+
  theme(axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"))
ggsave(paste0("pi_ns2_",substr(date(),5,7),substr(date(),9,10),".pdf"), height=5, width=12)

pis_long4a$NS<-c("North", "South")[match(pis_long4a$NS, c("N","S")) ]

#rm outlier
p1<-pis_long4a %>% filter(meanpi<.5)%>%
  ggplot( aes(x=posreal, y=meanpi))+
  geom_point(aes(col=Subgenome))+
  #scale_color_manual(values=cc2)+
  facet_grid(NS~Chr, space="free" ,scales = "free_x" )+
  #geom_hline(yintercept=0.05)+
  labs(x="Approximate position", y="Mean Pi")+
  #coord_cartesian(ylim = c(0,.1))+
  scale_color_manual(values=c('chocolate4', 'burlywood3', "skyblue"), guide="none")+
  theme_classic()+
  theme(axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"))

p1 + plot_annotation(tag_levels = "A")

ggsave(paste0("pi_ns_nooutlier_",substr(date(),5,7),substr(date(),9,10),".pdf"), height=5, width=12)


nrow(pis_long2%>%filter(NS=="S"))
nrow(pis_long2%>%filter(NS=="N"))

pis_long4a %>% filter(meanpi<.5)%>% group_by(NS) %>% summarise(meanpi=mean(meanpi))


pis_long4<-pis_long2 %>% filter(!is.na(Chromosome))
chrs<-unique(pis_long4$Chromosome)
chrs2<-NA
for(i in 1:36){chrs2[i]<-paste0(substr("Chr0",1,6-nchar(chrs[i])),chrs[i])}
chrkey<-data.frame(chrs=chrs, newchr=chrs2)
pis_long4$Chrom<-chrkey$newchr[match(pis_long4$Chromosome, chrkey$chrs)]
pis_long4$NS[which(pis_long4$NS=="N")]<-"North"
pis_long4$NS[which(pis_long4$NS=="S")]<-"South"

pis_long4$Subgenome<-substr(pis_long4$Chrom,6,6)
pis_long4$Chrom2<-substr(pis_long4$Chrom,4,6)

library(patchwork)

ggplot(pis_long4, aes(x=posdummy, y=meanpi))+
  geom_point(aes(col=Subgenome))+
  #scale_color_manual(values=cc2)+
  facet_grid(NS~Chrom2, space="free", scales="free", )+
  labs(x="Contig position", y="Mean Pi")+
  #coord_cartesian(ylim = c(0,.1))+
  scale_color_manual(values=c('black', 'grey70'), guide="none")+
  theme_classic()+
  theme(axis.text.x = element_blank(), panel.spacing.x=unit(0, "lines"))
ggsave(paste0("pi_ns_",substr(date(),5,7),substr(date(),9,10),".pdf"), height=5, width=12)


#stats without big outlier
meanN<-(pis_long2 %>% filter(NS=="N"&meanpi<.5) %>% select(meanpi))
mean(meanN$meanpi, na.rm=T)

meanS<-(pis_long2 %>% filter(NS=="S"&meanpi<.5) %>% select(meanpi))
mean(meanS$meanpi, na.rm=T)

#peaks
peakn<-pis_long2 %>%filter(NS=="N" & meanpi>.5)
peaks<-pis_long2 %>%filter(NS=="S" & meanpi>.5)

#number of outliers
countN<-meanN %>% filter(meanpi>0.025)
length(unique(countN$locus))

countS<-meanS %>% filter(meanpi>0.025)
length(unique(countS$locus))


