setwd("~/Desktop/P1 Pri-miRNA scaffold/Analysis cleavage Bioinformatic/")
miRBase <- read.table(file="high_confidence_miRBase21-master.tsv", header = TRUE, fill = TRUE, sep =  " ")
arm <- unique(miRBase[,c(1,16)])

setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/")
require(data.table)
iso_summary <- fread(file = "TUT_KO_Qiagen_rerun_2_5.isomir.sequence_info.tsv", header = T)
iso_summary$SAMPLE <- gsub(".fastq_ready","",iso_summary$SAMPLE)
unique(iso_summary$SAMPLE)
iso_summary <- unique(iso_summary)
iso_summary$MIRNA <- gsub(x = iso_summary$MIRNA,"hsa-miR-375","hsa-miR-375-3p")
arm$MIRNA <- gsub(x = arm$MIRNA, "hsa-miR-378c","hsa-miR-378c-5p")
arm$MIRNA <- gsub(x = arm$MIRNA, "hsa-miR-548s","hsa-miR-548s-3p")
iso_summary$MIRNA <- gsub(x = iso_summary$MIRNA,"mir","miR")


agg_sample_miRNA <- aggregate(iso_summary$READS, by=list(SAMPLE=iso_summary$SAMPLE,MIRNA=iso_summary$MIRNA),sum)
selected_samples <- c("1_WT_293T_1","2_WT_293T_2")

all_reads <- aggregate(agg_sample_miRNA$x, by=list(MIRNA=agg_sample_miRNA$MIRNA),sum)
all_reads <- merge(all_reads, arm,by="MIRNA")
all_reads <- all_reads[order(-all_reads$x),]
all_5P <- all_reads[which(all_reads$STRAND=="5P"),]
all_3P <- all_reads[which(all_reads$STRAND=="3P"),]
all_5P <- all_5P[c(1:100),]
all_3P <- all_3P[c(1:100),]
all_5P3P <- rbind(all_5P,all_3P)

anytail <- iso_summary[which(iso_summary$LEN_TAIL>0),]
anytail <- anytail[anytail$SAMPLE %in% selected_samples,]
anytail <- anytail[anytail$MIRNA %in% all_5P3P$MIRNA,]

tail_test <- unique(anytail[,c(2,3)])
tail_test$templated <- NA

i <- 1
for (i in 1:nrow(tail_test)) {
  seq_test <- tail_test$SEQUENCE[i]
  miRNA_test <- tail_test$MIRNA[i]
  miRBase_test <- miRBase
  miRBase_test <- miRBase_test[which(miRBase_test$MIRNA==miRNA_test),]
  miRBase_test$found <- regexpr(seq_test, miRBase_test$EXTENDED.SEQUENCE)
  miRBase_test <- miRBase_test[which(miRBase_test$found!=-1),]
  if(nrow(miRBase_test)>0){
    tail_test$templated[i] <- T
  }
  if(nrow(miRBase_test)==0){
    tail_test$templated[i] <- F
  }
  rm(miRBase_test,seq_test,miRNA_test)
  print(i/nrow(tail_test)*100)
}

tail_test <- tail_test[which(tail_test$templated==F),]

anytail2 <- anytail[anytail$SEQUENCE %in% tail_test$SEQUENCE,]
anytail2$pureA <- FALSE
anytail2$pureU <- FALSE

i <- 1
for (i in 1:nrow(anytail2)) {
  seq_test <- anytail2$SEQ_TAIL[i]
  U_test <- gsub("T","", seq_test)
  A_test <- gsub("A","", seq_test)

  if(nchar(U_test)==0){
    anytail2$pureU[i] <- T
  }
  if(nchar(A_test)==0){
    anytail2$pureA[i] <- T
  }
  rm(seq_test,U_test,A_test)
  print(i/nrow(anytail2)*100)
}
#saved <- anytail2
anytail2 <- saved
#anytail2 <- anytail2[anytail2$MIRNA %in% all_5P$MIRNA,]
anytail2 <- anytail2[anytail2$MIRNA %in% all_3P$MIRNA,]


Atail <- anytail2[which(anytail2$pureA==T),]
Atail <- aggregate(Atail$READS, by=list(SAMPLE=Atail$SAMPLE,MIRNA=Atail$MIRNA,LEN_TAIL=Atail$LEN_TAIL),sum)
Utail <- anytail2[which(anytail2$pureU==T),]
Utail <- aggregate(Utail$READS, by=list(SAMPLE=Utail$SAMPLE,MIRNA=Utail$MIRNA,LEN_TAIL=Utail$LEN_TAIL),sum)
Mtail <- anytail2[which(anytail2$pureU==F & anytail2$pureA==F),]
Mtail <- aggregate(Mtail$READS, by=list(SAMPLE=Mtail$SAMPLE,MIRNA=Mtail$MIRNA,LEN_TAIL=Mtail$LEN_TAIL),sum)

Atail <- merge(Atail,agg_sample_miRNA,by=c("SAMPLE", "MIRNA"))
Utail <- merge(Utail,agg_sample_miRNA,by=c("SAMPLE", "MIRNA"))
Mtail <- merge(Mtail,agg_sample_miRNA,by=c("SAMPLE", "MIRNA"))
Atail$percent_len_tail <- Atail$x.x/Atail$x.y*100
Utail$percent_len_tail <- Utail$x.x/Utail$x.y*100
Mtail$percent_len_tail <- Mtail$x.x/Mtail$x.y*100

Atail_2 <- aggregate(Atail$percent_len_tail, by=list(MIRNA=Atail$MIRNA,LEN_TAIL=Atail$LEN_TAIL),mean)
Atail_mean <- aggregate(Atail_2$x, by=list(LEN_TAIL=Atail_2$LEN_TAIL),mean)
colnames(Atail_mean) <- c("LEN_TAIL","Atail_mean")
Atail_sd <- aggregate(Atail_2$x, by=list(LEN_TAIL=Atail_2$LEN_TAIL),sd)
colnames(Atail_sd) <- c("LEN_TAIL","Atail_sd")
Atail_2$count <- 1
Atail_N <- aggregate(Atail_2$count, by=list(LEN_TAIL=Atail_2$LEN_TAIL),sum)
colnames(Atail_N) <- c("LEN_TAIL","Atail_N")


Utail_2 <- aggregate(Utail$percent_len_tail, by=list(MIRNA=Utail$MIRNA,LEN_TAIL=Utail$LEN_TAIL),mean)
Utail_mean <- aggregate(Utail_2$x, by=list(LEN_TAIL=Utail_2$LEN_TAIL),mean)
colnames(Utail_mean) <- c("LEN_TAIL","Utail_mean")
Utail_sd <- aggregate(Utail_2$x, by=list(LEN_TAIL=Utail_2$LEN_TAIL),sd)
colnames(Utail_sd) <- c("LEN_TAIL","Utail_sd")
Utail_2$count <- 1
Utail_N <- aggregate(Utail_2$count, by=list(LEN_TAIL=Utail_2$LEN_TAIL),sum)
colnames(Utail_N) <- c("LEN_TAIL","Utail_N")

Mtail_2 <- aggregate(Mtail$percent_len_tail, by=list(MIRNA=Mtail$MIRNA,LEN_TAIL=Mtail$LEN_TAIL),mean)
Mtail_mean <- aggregate(Mtail_2$x, by=list(LEN_TAIL=Mtail_2$LEN_TAIL),mean)
colnames(Mtail_mean) <- c("LEN_TAIL","Mtail_mean")
Mtail_sd <- aggregate(Mtail_2$x, by=list(LEN_TAIL=Mtail_2$LEN_TAIL),sd)
colnames(Mtail_sd) <- c("LEN_TAIL","Mtail_sd")
Mtail_2$count <- 1
Mtail_N <- aggregate(Mtail_2$count, by=list(LEN_TAIL=Mtail_2$LEN_TAIL),sum)
colnames(Mtail_N) <- c("LEN_TAIL","Mtail_N")

anytail5 <- merge(Atail_mean, Atail_sd, by="LEN_TAIL", all=T)
anytail5 <- merge(anytail5, Atail_N, by="LEN_TAIL", all=T)
anytail5 <- merge(anytail5, Utail_mean, by="LEN_TAIL", all=T)
anytail5 <- merge(anytail5, Utail_sd, by="LEN_TAIL", all=T)
anytail5 <- merge(anytail5, Utail_N, by="LEN_TAIL", all=T)
anytail5 <- merge(anytail5, Mtail_mean, by="LEN_TAIL", all=T)
anytail5 <- merge(anytail5, Mtail_sd, by="LEN_TAIL", all=T)
anytail5 <- merge(anytail5, Mtail_N, by="LEN_TAIL", all=T)
anytail5[is.na(anytail5)] <- 0


sum(anytail5$Atail_mean)+sum(anytail5$Utail_mean)+sum(anytail5$Mtail_mean)

write.table(anytail5, "tail_by_length_composition_3p.txt", sep="\t", append = FALSE, row.names = F)
