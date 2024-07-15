#analyze gwas peaks table for top 100 outliers

setwd("~/Desktop/GWAS")

southpeaks_april19 <- read.csv("~/Desktop/GWAS/southpeaks_april19.csv")
northpeaks_april19 <- read.csv("~/Desktop/GWAS/northpeaks_april19.csv")

southpeaks<-southpeaks_april19 %>% arrange(pval)%>% slice_head(n = 100)
northpeaks<-northpeaks_april19 %>% arrange(pval)%>% slice_head(n = 100)

#Find gene overlaps with top 100 hits####
#auto add fxn from phytozyme
#https://data.jgi.doe.gov/refine-download/phytozome?organism=Pvirgatum&expanded=516&_gl=1*1qor4xr*_ga*MTA4Mjk2NDY4Mi4xNjc4OTA0Mjc5*_ga_YBLMHYR3C2*MTY5ODc2MzMzNi4xMS4wLjE2OTg3NjMzMzYuMC4wLjA.

pv_anno <- read.delim("Pvirgatum_516_v5.1.annotation_info.txt", header=FALSE, comment.char="#")
pv_gff  <- read.delim("Pvirgatum_516_v5.1.gene.gff3", header=FALSE, comment.char="#")

#crop gff to mrna
pv_gff2<-pv_gff %>% filter(V3 == "mRNA")

#match to anno
pv_gff2$gene<-substr(pv_gff2$V9,4,18)
pv_gff2$gene_scaf<-substr(pv_gff2$V9,4,16)

pv_anno$start<-pv_gff2$V4[match(pv_anno$V2, pv_gff2$gene)]

#scaffolds don't match. add here
pv_anno$start1<-pv_gff2$V4[match(pv_anno$V2, pv_gff2$gene_scaf)]
pv_anno$start[is.na(pv_anno$start)]<-pv_anno$start1[is.na(pv_anno$start)]
pv_anno$end<-pv_gff2$V5[match(pv_anno$V2, pv_gff2$gene)]
pv_anno$end1<-pv_gff2$V4[match(pv_anno$V2, pv_gff2$gene_scaf)]
pv_anno$end[is.na(pv_anno$end)]<-pv_anno$end1[is.na(pv_anno$end)]

#iranges to find overlaps
library(IRanges)
pv_anno2<-pv_anno
pv_anno2$CHR<-paste0("Chr0",substr(pv_anno2$V2,7,8))
#pv_anno2<-pv_anno2[which(substr(pv_anno2$CHR,2,2) %in% c("N","K")),]

ranges_gff<-split(IRanges(pv_anno2$start, pv_anno2$end), pv_anno2$CHR)

northpeaks$region<-"N"
southpeaks$region<-"S"

GWAS_peaks<-rbind(northpeaks, southpeaks)
peaks1<-GWAS_peaks[!is.na(GWAS_peaks$POS),]
peaks1$POS<-as.numeric(peaks1$POS)
#peaks1<-pv_anno2[which(substr(pv_anno2$CHR,2,2) %in% c("N","K")),]
ranges_hits<-split(IRanges(peaks1$POS-10000, peaks1$POS+10000), peaks1$CHR)

overlaps1<-findOverlaps(ranges_gff, ranges_hits)
overlaps2<-vector('list', 18)
for(i in 1:18){overlaps2[[i]]<-as.matrix(overlaps1[[i]])}
names(overlaps2)<-names(overlaps1)[1:18]

#get best hit for each locus. query is pv_anno, subject is peaks
pv_annolist<-split(pv_anno2, pv_anno2$CHR)
peaks_list<-split(peaks1, peaks1$CHR)
overlaps3<-overlaps2[names(peaks_list)]
#check for ones with no overlaps
lapply(overlaps3, nrow)
length(overlaps3)

pv_annolist2<-pv_annolist[names(overlaps3)]

allhitsout<-overlaps3

for(i in 1:length(overlaps3)){
  peaks_list[[i]]$tophit<-"no_match"
  peaks_list[[i]]$tophit[unique(overlaps3[[i]][,2])]<-"tophit"
  #split overlaps3 by subjecthits
  overdf<-as.data.frame(overlaps3[[i]])
  overlaps_temp<-split(overdf, overdf[,2])
  #sub first query hit into peakslist
  peaks_list[[i]]$hitmatch<-NA
  #for each hit in overlaps, pick the first query, then sub pv_anno$v2 to this rownum
  for(j in 1:length(overlaps_temp)){
  peaks_list[[i]]$hitmatch[overlaps_temp[[j]][1,2]]<-pv_annolist2[[i]][,2][overlaps_temp[[j]][1,1]]
  }
  #print(table(peaks_list2[[i]]$tophit))
  #need to output all also. better to link peaks_list onto hits?
  overdf$genename<-pv_annolist2[[i]][overdf$queryHits,2]
  overdf$fxn<-pv_annolist2[[i]][overdf$queryHits,13]
  overdf$GO<-pv_annolist2[[i]][overdf$queryHits,10]
  overdf$panther<-pv_annolist2[[i]][overdf$queryHits,6]
  #add in peak info
  overdf<-cbind(overdf, peaks_list[[i]][overdf$subjectHits,])
  allhitsout[[i]]<-overdf
}

allhitsout2<-do.call(rbind, allhitsout)

#clean up
allhitsout3<- allhitsout2 %>% select(POS, CHR, pval, genename, fxn, tophit, region)
allhitsout4<-allhitsout3[!duplicated(allhitsout3),]

#add back in nonhits
pl1<-do.call(rbind, peaks_list)
pl2<-pl1 %>% filter(tophit=="no_match")
pl2$genename<-NA
pl2$fxn<-NA
pl3<-pl2 %>%select(POS, CHR, pval, genename, fxn, tophit, region)

allhitsout5<-rbind(allhitsout4, pl3)
allhitsout6<-allhitsout5 %>% group_by(CHR) %>%arrange(POS,.by_group = T)

paste0("there are ",length(unique(allhitsout6$POS)), " outlier loci linked to ", nrow(allhitsout6[!is.na(allhitsout6$genename),]), " genes")

colnames(allhitsout6)<-c("Position", "Chromosome", "log10pval", "GeneID", "Predicted Function", "Match", "Region")
allhitsout6$Match[allhitsout6$Match=="tophit"]<-""
allhitsout6$Match[allhitsout6$Match=="no_match"]<-"No gene"

#write_csv(allhitsout6, "allhitsout6_NS_april23.csv")
allhitsout6N<-allhitsout6 %>% filter(Region=="N")
allhitsout6S<-allhitsout6 %>% filter(Region=="S")


unique_genesN<-data.frame(genes=unique(allhitsout6N$GeneID))
#write_csv(unique_genesN, "uniquegenes_north_april19.csv")

unique_genesS<-data.frame(genes=unique(allhitsout6S$GeneID))
#write_csv(unique_genesS, "uniquegenes_south_april19.csv")

#find overrepresented annotations####
#format for global vs local comparison
#need to dedup prots. 
allhitsout6_north<- allhitsout6N
protsN<-unique(allhitsout6_north$GeneID)

sample_set<-allhitsout6_north$`Predicted Function`[match(protsN, allhitsout6_north$GeneID)]
#split by space, if last element can be coerced into numeric, drop it. 
ss1<-strsplit(sample_set,split = " ")
ss2<-vector("list", length(ss1))
for(i in 1:length(ss1)){
  n1<-length(ss1[[i]])
  tryCatch({
    skip_to_next<-F  
    
    if (is.na(suppressWarnings(as.numeric(ss1[[i]][n1])))){
      ss2[[i]]<-paste0(ss1[[i]], collapse = " ")}
    else {
      ss2[[i]]<-paste0(ss1[[i]][-n1], collapse = " ")}
    
  },error=function(x){skip_to_next<<-T})
  if(skip_to_next) { next }   
}

#get  global
global_set<-pv_anno2$V13

gs1<-strsplit(global_set,split = " ")
gs2<-vector("list", length(gs1))
for(i in 1:length(gs1)){
  n1<-length(gs1[[i]])
  tryCatch({
    skip_to_next<-F  
    
    if (is.na(suppressWarnings(as.numeric(gs1[[i]][n1])))){
      gs2[[i]]<-paste0(gs1[[i]], collapse = " ")}
    else {
      gs2[[i]]<-paste0(gs1[[i]][-n1], collapse = " ")}
    
  },error=function(x){skip_to_next<<-T})
  if(skip_to_next) { next }   
}

testdiff<-function(gs2, ss2){# Create a vector of all unique phrases from both sets
  all_phrases <- unique(c(gs2, ss2))
  
  # Create frequency tables for global and sample sets
  global_freq <- table(factor(gs2, levels = all_phrases))
  sample_freq <- table(factor(ss2, levels = all_phrases))
  
  # Create a data frame for the comparison
  comparison_data <- data.frame(Global = global_freq, Sample = sample_freq)
  
  # Perform a chi-squared test
  #chisq_result <- chisq.test(x = comparison_data$Global.Freq, p=comparison_data$Sample.Freq/sum(comparison_data$Sample.Freq))
  
  #chisq_result$p.value
  
  #add in freqs
  comparison_data$glob_freq<-comparison_data$Global.Freq/(sum(comparison_data$Global.Freq))
  comparison_data$samp_freq<-comparison_data$Sample.Freq/(sum(comparison_data$Sample.Freq))
  comparison_data$samp_se<-sqrt(comparison_data$samp_freq * (1-comparison_data$samp_freq)/sum(sum(comparison_data$Sample.Freq)))
  comparison_data$zscore <- (comparison_data$glob_freq-comparison_data$samp_freq)/comparison_data$samp_se
  
  #get diff, then sort
  comparison_data$freqdiff<-comparison_data$samp_freq-comparison_data$glob_freq
  return(comparison_data)
}

alpha=.05
#zscore cutoff
cutoff=qnorm(1 - alpha)

comparison_dataN<-testdiff(gs2,ss2)
all_phrases <- unique(c(gs2, ss2))

#south
allhitsout_south<-allhitsout6S
protsS<-unique(allhitsout_south$GeneID)
sample_setS<-allhitsout_south$`Predicted Function`[match(protsS, allhitsout_south$GeneID)]

ss1S<-strsplit(sample_setS,split = " ")
ss2S<-vector("list", length(ss1S))
for(i in 1:length(ss1S)){
  n1<-length(ss1S[[i]])
  tryCatch({
    skip_to_next<-F  
    
    if (is.na(suppressWarnings(as.numeric(ss1S[[i]][n1])))){
      ss2S[[i]]<-paste0(ss1S[[i]], collapse = " ")}
    else {
      ss2S[[i]]<-paste0(ss1S[[i]][-n1], collapse = " ")}
    
  },error=function(x){skip_to_next<<-T})
  if(skip_to_next) { next }   
}

comparison_dataS<-testdiff(gs2, ss2S)


#bootstrap test####
#1. resample from global
#2. compute table
#3. compute chi-sq value sum((o-e)^2/e)
#how many reps had a greater chi-sq?

#first get expected based on n
comparison_dataN$exp<-comparison_dataN$glob_freq*sum(comparison_dataN$Sample.Freq)
gs3<-unlist(gs2)
comparison_dataN2<-comparison_dataN 

nchi<-10000
chisums<-vector("list", nchi)
freqdiffN<-vector("list", nchi)

for( i in 1:nchi){
  #resample from global
  samp<-sample(gs3, nrow(northpeaks), replace=F)
  #create table
  tab_freq <- as.data.frame(table(factor(samp, levels = all_phrases)))
  #add to comp (overwriting)
  comparison_dataN2$sampX<-tab_freq$Freq[match(comparison_dataN2$Global.Var1, tab_freq$Var1)]
  #calc X2
  chisums[[i]]<-((comparison_dataN2$sampX-comparison_dataN2$exp)^2)/comparison_dataN2$exp
  freqdiffN[[i]]<-comparison_dataN2$sampX/nrow(northpeaks)
}

chisums2<-sapply(chisums, FUN=function(x){sum(x, na.rm=T)})

#what is actual chisum
comparison_dataN2$chireal<-((comparison_dataN2$Sample.Freq-comparison_dataN2$exp)^2)/comparison_dataN2$exp
chireal<-sum(comparison_dataN2$chireal[which(!is.infinite(comparison_dataN2$chireal))], na.rm = T)

hist(chisums2)
#in how many of the sims was chisum greater than chireal?
length(chisums2[chisums2>chireal])
#7786/10000

#check how many had higher X2 for each line
comparison_dataN2$testfeature<-NA

for(i in 1:5190){
  nulldist<-NA
  for(k in 1:10000){
    nulldist[k]<-chisums[[k]][i]
  }
  comparison_dataN2$testfeature[i]<-(length(nulldist[nulldist>comparison_dataN2$chireal[i]])/10000)
}
comparison_dataN2$Global.Var1[comparison_dataN2$testfeature>.95]

comparison_dataN2$testfreq<-NA
for(i in 1:5190){
  nulldist<-NA
  for(k in 1:10000){
    nulldist[k]<-freqdiffN[[k]][i]
  }
  comparison_dataN2$testfreq[i]<-(length(nulldist[abs(nulldist)>abs(comparison_dataN2$samp_freq[i])])/10000)
}
comparison_dataN2$Global.Var1[comparison_dataN2$testfreq>.95]

#for south
comparison_dataS$exp<-comparison_dataS$glob_freq*sum(comparison_dataS$Sample.Freq)
comparison_dataS2<-comparison_dataS 

nchi<-10000
chisumsS<-vector("list", nchi)
freqdiffS<-vector("list", nchi)
set.seed(487573)
for( i in 1:nchi){
  #resample from global
  samp<-sample(gs3, nrow(southpeaks), replace=F)
  #create table
  tab_freq <- as.data.frame(table(factor(samp, levels = all_phrases)))
  #add to comp (overwriting)
  comparison_dataS2$sampX<-tab_freq$Freq[match(comparison_dataS2$Global.Var1, tab_freq$Var1)]
  #calc X2
  chisumsS[[i]]<-((comparison_dataS2$sampX-comparison_dataS2$exp)^2)/comparison_dataS2$exp
  freqdiffS[[i]]<-comparison_dataS2$sampX/nrow(southpeaks)
}

chisumsS2<-sapply(chisumsS, FUN=function(x){sum(x, na.rm=T)})

#what is actual chisum
comparison_dataS2$chireal<-((comparison_dataS2$Sample.Freq-comparison_dataS2$exp)^2)/comparison_dataS2$exp
chirealS<-sum(comparison_dataS2$chireal[which(!is.infinite(comparison_dataS2$chireal))], na.rm = T)

hist(chisumsS2)
#in how many of the sims was chisum greater than chireal?
length(chisumsS2[chisumsS2>chirealS])
#11/nchi

#check how many had higher X2 for each line
comparison_dataS2$testfeature<-NA

#for each fxn
for(i in 1:5190){
  nulldist<-NA
  #count across all bootstraps
  for(k in 1:10000){
    nulldist[k]<-chisumsS[[k]][i]
  }
  comparison_dataS2$testfeature[i]<-(length(nulldist[nulldist>comparison_dataS2$chireal[i]])/10000)
}
comparison_dataS2$Global.Var1[comparison_dataS2$testfeature>.95]

comparison_dataS2$testfreq<-NA
for(i in 1:5190){
  nulldist<-NA
  for(k in 1:10000){
    nulldist[k]<-freqdiffS[[k]][i]
  }
  comparison_dataS2$testfreq[i]<-(length(nulldist[abs(nulldist)>abs(comparison_dataS2$samp_freq[i])])/10000)
}

#plotting####

top_overN<-comparison_dataN %>% arrange(desc(freqdiff)) %>% filter(Global.Freq>0)%>% slice_head(n=25)
varsN<-as.character(top_overN$Sample.Var1)
top_overN$Sample.Var1<- varsN
top_overN$Sample.Var1<-as.factor(top_overN$Sample.Var1)

ggplot(top_overN, aes(x=reorder(Sample.Var1,freqdiff,mean), y=freqdiff))+geom_col(fill="forestgreen", col="black")+coord_flip()+theme_classic()+
  ggtitle("Top 25 Functional annotations North")  +
  labs(y="Frequency Difference",x= "Annotation")
ggsave("Overrep_N_april19.pdf", height=7, width=10)

top_underN<-comparison_dataN %>% arrange(desc(freqdiff)) %>% slice_tail(n=25)
#edit long names
top_underN$Global.Var1<-as.character(top_underN$Global.Var1)
top_underN$Global.Var1[nchar(as.character(top_underN$Global.Var1))>70]<-c("Bifunctional inhibitor/lipid-transfer superfamily protein",
                                                                          "S-adenosyl-L-methionine-dependent methyltransferases",
                                                                          "P-loop containing nucleoside triphosphate hydrolases")

ggplot(top_underN, aes(x=reorder(Global.Var1,freqdiff,mean), y=freqdiff))+geom_col(fill="firebrick", col="black")+coord_flip()+theme_classic()+
  ggtitle("Underrepresented Functional annotations North")  +
  labs(y="Frequency Difference",x= "Annotation")
ggsave("Underrep_N_april19.pdf", height=7, width=10)

#south####


#run bootstrap for all N####
count<-1
max = 100000

top25n<-vector('list', 25)
date()
while(count < max){
  #sample
  samp1<-sample(global_full, 113, replace = T)
  #get full table w zeroes
  sampall<-factor(samp1, levels = uniq_glob)
  #make table
  tab1<-as.data.frame(table(sampall))
  
  #only do 1
  for(i in 1:25){
    #if samp is greater than target, append
    target<-top_overN$Sample.Freq[i]
    samp<-tab1$Freq[as.character(tab1$sampall)==as.character(top_overN$Sample.Var1)[i]]
    if(samp>= target){
      top25n[[i]]<-append(top25n[[i]],1)
    }
  }
  count=count+1
}
date()
save.image("reformat_peaks_table_dec18.RData")

top_overN2<-top_overN[1:20,]
top_overN2$Global.Var1<-as.character(top_overN2$Global.Var1)
top_overN2$Global.Var1[nchar(as.character(top_overN2$Global.Var1))>63]<-c("Terpenoid cyclases/Protein prenyltransferases","2-oxoglutarate (2OG) and Fe(II)-dependent oxygenase"," COBRA-like extracellular glycosyl-phosphatidyl inositol-anchored")


ggplot(top_overN2, aes(x=reorder(Global.Var1,freqdiff,mean), y=freqdiff))+geom_col(fill="forestgreen", col="black")+coord_flip()+theme_classic()+
  ggtitle("Top 20 Functional annotations North")  +
  annotate(geom="text", x=c(14.8,17.8,13.8,19.8), y=.001, label="*", size=11)+
  labs(y="Frequency Difference",x= "Annotation")
ggsave("Overrep_N_dec18.pdf", height=7, width=10)

bonf_n<-.05/113


top_overS<-comparison_dataS %>% arrange(desc(freqdiff)) %>% slice_head(n=25)
varsS<-as.character(top_overS$Sample.Var1)
top_overS$Sample.Var1<- varsS
#top_overS$Sample.Var1[2]<-"Unknown"
top_overS$Sample.Var1<-as.factor(top_overS$Sample.Var1)

top_underS<-comparison_dataS %>% arrange(desc(freqdiff)) %>% slice_tail(n=25)

ggplot(top_overS, aes(x=reorder(Sample.Var1,freqdiff,mean), y=freqdiff))+geom_col(fill="forestgreen", col="black")+coord_flip()+theme_classic()+
  ggtitle("Overrepresented Functional annotations South")  +
  labs(y="Frequency Difference",x= "Annotation")
ggsave("Overrep_S_nov9.pdf", height=7, width=10)

ggplot(top_underS, aes(x=reorder(Sample.Var1,freqdiff,mean), y=freqdiff))+geom_col(fill="firebrick", col="black")+coord_flip()+theme_classic()+
  ggtitle("Underrepresented Functional annotations South")  +
  labs(y="Frequency Difference",x= "Annotation")
ggsave("Underrep_S_nov9.pdf", height=7, width=10)

