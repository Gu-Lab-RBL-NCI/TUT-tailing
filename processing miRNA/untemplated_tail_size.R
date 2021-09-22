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
anytail2 <- anytail[anytail$MIRNA %in% all_5P$MIRNA,]
anytail2 <- anytail[anytail$MIRNA %in% all_3P$MIRNA,]



anytail2 <- aggregate(anytail2$READS, by=list(SAMPLE=anytail2$SAMPLE,MIRNA=anytail2$MIRNA,LEN_TAIL=anytail2$LEN_TAIL),sum)
anytail3 <- merge(anytail2,agg_sample_miRNA,by=c("SAMPLE", "MIRNA"))
anytail3$percent_len_tail <- anytail3$x.x/anytail3$x.y*100
anytail4 <- aggregate(anytail3$percent_len_tail, by=list(MIRNA=anytail3$MIRNA,LEN_TAIL=anytail3$LEN_TAIL),mean)
anytail5 <- aggregate(anytail4$x, by=list(LEN_TAIL=anytail4$LEN_TAIL),mean)
colnames(anytail5) <- c("LEN_TAIL","mean")
anytail5_sd <- aggregate(anytail4$x, by=list(LEN_TAIL=anytail4$LEN_TAIL),sd)
colnames(anytail5_sd) <- c("LEN_TAIL","sd")
anytail4$count <- 1
anytail5_N <- aggregate(anytail4$count, by=list(LEN_TAIL=anytail4$LEN_TAIL),sum)
colnames(anytail5_N) <- c("LEN_TAIL","N")
anytail5 <- merge(anytail5, anytail5_sd, by="LEN_TAIL")
anytail5 <- merge(anytail5, anytail5_N, by="LEN_TAIL")
sum(anytail5$mean)
write.table(anytail5, "tail_by_length_3p.txt", sep="\t", append = FALSE, row.names = F)
