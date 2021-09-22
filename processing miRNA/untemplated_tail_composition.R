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

canonical <- iso_summary[which(iso_summary$LEN_TRIM==0 & iso_summary$LEN_TAIL==0),]
canonical <- aggregate(canonical$READS, by=list(SAMPLE=canonical$SAMPLE,MIRNA=canonical$MIRNA),sum)
canonical <- merge(canonical,agg_sample_miRNA,by=c("MIRNA","SAMPLE"))
canonical$percent_can <- canonical$x.x/canonical$x.y*100
canonical <- canonical[,c(1,2,4,5)]

trimmed <- iso_summary[which(iso_summary$LEN_TRIM>0 & iso_summary$LEN_TAIL==0),]
trimmed <- aggregate(trimmed$READS, by=list(SAMPLE=trimmed$SAMPLE,MIRNA=trimmed$MIRNA),sum)
trimmed <- merge(trimmed,agg_sample_miRNA,by=c("MIRNA","SAMPLE"))
trimmed$percent_trim <- trimmed$x.x/trimmed$x.y*100
trimmed <- trimmed[,c(1,2,5)]

trimtail <- iso_summary[which(iso_summary$LEN_TRIM>0 & iso_summary$LEN_TAIL>0),]
trimtail <- aggregate(trimtail$READS, by=list(SAMPLE=trimtail$SAMPLE,MIRNA=trimtail$MIRNA),sum)
trimtail <- merge(trimtail,agg_sample_miRNA,by=c("MIRNA","SAMPLE"))
trimtail$percent_trimtail <- trimtail$x.x/trimtail$x.y*100
trimtail <- trimtail[,c(1,2,5)]

tail <- iso_summary[which(iso_summary$LEN_TRIM==0 & iso_summary$LEN_TAIL>0),]
#tail <- aggregate(tail$READS, by=list(SAMPLE=tail$SAMPLE,MIRNA=tail$MIRNA),sum)
#tail <- merge(tail,agg_sample_miRNA,by=c("MIRNA","SAMPLE"))
#tail$percent_tail <- tail$x.x/tail$x.y*100
#tail <- tail[,c(1,2,5)]
tail_test <- unique(tail[,c(2,3)])
tail_test$templated <- NA
tail_test <- unique(tail_test[,c(2,3)])
i <- 1
for (i in 1:nrow(tail_test)) {
  seq_test <- tail_test$SEQUENCE[i]
  miRBase_test <- miRBase
  miRBase_test$found <- regexpr(seq_test, miRBase_test$EXTENDED.SEQUENCE)
  miRBase_test <- miRBase_test[which(miRBase_test$found!=-1),]
  if(nrow(miRBase_test)>0){
    tail_test$templated[i] <- T
  }
  if(nrow(miRBase_test)==0){
    tail_test$templated[i] <- F
  }
  rm(miRBase_test,seq_test)
  print(i/nrow(tail_test)*100)
}

tail <- merge(tail, tail_test, by="SEQUENCE")
templated <- tail[which(tail$templated==T),]
templated <- aggregate(templated$READS, by=list(SAMPLE=templated$SAMPLE,MIRNA=templated$MIRNA),sum)
templated <- merge(templated,agg_sample_miRNA,by=c("MIRNA","SAMPLE"))
templated$percent_templated <- templated$x.x/templated$x.y*100
templated <- templated[,c(1,2,5)]

untemplated <- tail[which(tail$templated==F),]
untemplated$tail_events <-untemplated$LEN_TAIL*untemplated$READS
untemplated$fract_A <- (nchar(untemplated$SEQ_TAIL)-nchar(gsub("A","",untemplated$SEQ_TAIL)))*untemplated$READS
untemplated$fract_T <- (nchar(untemplated$SEQ_TAIL)-nchar(gsub("T","",untemplated$SEQ_TAIL)))*untemplated$READS
untemplated$fract_C <- (nchar(untemplated$SEQ_TAIL)-nchar(gsub("C","",untemplated$SEQ_TAIL)))*untemplated$READS
untemplated$fract_G <- (nchar(untemplated$SEQ_TAIL)-nchar(gsub("G","",untemplated$SEQ_TAIL)))*untemplated$READS

all_tailing <- aggregate(untemplated$tail_events, by=list(SAMPLE=untemplated$SAMPLE,MIRNA=untemplated$MIRNA),sum) 
A_tailing <- aggregate(untemplated$fract_A, by=list(SAMPLE=untemplated$SAMPLE,MIRNA=untemplated$MIRNA),sum) 
A_tailing <- merge(A_tailing,all_tailing,by=c("MIRNA","SAMPLE"))
A_tailing$A_frac <- A_tailing$x.x/A_tailing$x.y
A_tailing <- A_tailing[,c(1,2,5)]
T_tailing <- aggregate(untemplated$fract_T, by=list(SAMPLE=untemplated$SAMPLE,MIRNA=untemplated$MIRNA),sum) 
T_tailing <- merge(T_tailing,all_tailing,by=c("MIRNA","SAMPLE"))
T_tailing$T_frac <- T_tailing$x.x/T_tailing$x.y
T_tailing <- T_tailing[,c(1,2,5)]
C_tailing <- aggregate(untemplated$fract_C, by=list(SAMPLE=untemplated$SAMPLE,MIRNA=untemplated$MIRNA),sum) 
C_tailing <- merge(C_tailing,all_tailing,by=c("MIRNA","SAMPLE"))
C_tailing$C_frac <- C_tailing$x.x/C_tailing$x.y
C_tailing <- C_tailing[,c(1,2,5)]
G_tailing <- aggregate(untemplated$fract_G, by=list(SAMPLE=untemplated$SAMPLE,MIRNA=untemplated$MIRNA),sum) 
G_tailing <- merge(G_tailing,all_tailing,by=c("MIRNA","SAMPLE"))
G_tailing$G_frac <- G_tailing$x.x/G_tailing$x.y
G_tailing <- G_tailing[,c(1,2,5)]
sum_tailing <- merge(A_tailing,T_tailing,by=c("MIRNA","SAMPLE"))
sum_tailing <- merge(sum_tailing,G_tailing,by=c("MIRNA","SAMPLE"))
sum_tailing <- merge(sum_tailing,C_tailing,by=c("MIRNA","SAMPLE"))
sum_tailing$sum <- sum_tailing$A_frac+sum_tailing$T_frac+sum_tailing$G_frac+sum_tailing$C_frac

selected_samples <- c("1_WT_293T_1","2_WT_293T_2","3_TUT2KO_S4","4_TUT2KO_S6","5_TUT4KO_S1","6_TUT4KO_S2","7_TUT7KO_combined","8_TUT7KO_S2",
                      "10_DKO-24-2-IP-miseq","10_DKO_24_2","11_TKO_7A","11_TKO-7A-IP_miseq")
agg_sample_miRNA <- agg_sample_miRNA[agg_sample_miRNA$SAMPLE %in% selected_samples,]

all2 <- sum_tailing[sum_tailing$SAMPLE %in% selected_samples,]
all2$SAMPLE <- gsub("1_WT_293T_1","A_WT_293T_1",all2$SAMPLE)
all2$SAMPLE <- gsub("2_WT_293T_2","B_WT_293T_2",all2$SAMPLE)
all2$SAMPLE <- gsub("3_TUT2KO_S4","C_TUT2KO_S4",all2$SAMPLE)
all2$SAMPLE <- gsub("4_TUT2KO_S6","D_TUT2KO_S6",all2$SAMPLE)
all2$SAMPLE <- gsub("5_TUT4KO_S1","E_TUT4KO_S1",all2$SAMPLE)
all2$SAMPLE <- gsub("6_TUT4KO_S2","F_TUT4KO_S2",all2$SAMPLE)
all2$SAMPLE <- gsub("7_TUT7KO_combined","G_TUT7KO_combined",all2$SAMPLE)
all2$SAMPLE <- gsub("8_TUT7KO_S2","H_TUT7KO_S2",all2$SAMPLE)
all2$SAMPLE <- gsub("10_DKO","I_DKO",all2$SAMPLE)
all2$SAMPLE <- gsub("10_DKO_24_2","J_DKO_24_2",all2$SAMPLE)
all2$SAMPLE <- gsub("11_TKO_7A","K_TKO_7A",all2$SAMPLE)
all2$SAMPLE <- gsub("11_TKO-7A-IP_miseq","L_TKO-7A-IP_miseq",all2$SAMPLE)

all_reads <- aggregate(agg_sample_miRNA$x, by=list(MIRNA=agg_sample_miRNA$MIRNA),sum)
all_reads <- merge(all_reads, arm,by="MIRNA")
all_reads <- all_reads[order(-all_reads$x),]
all_5P <- all_reads[which(all_reads$STRAND=="5P"),]
all_3P <- all_reads[which(all_reads$STRAND=="3P"),]
all_5P <- all_5P[c(1:100),]
all_3P <- all_3P[c(1:100),]
all_5P3P <- rbind(all_5P,all_3P)

all2 <- all2[all2$MIRNA %in% all_5P3P$MIRNA,]
all2 <- all2[all2$MIRNA %in% all_5P$MIRNA,]
all2 <- all2[all2$MIRNA %in% all_3P$MIRNA,]

A_tail <- aggregate(all2$A_frac, by=list(SAMPLE=all2$SAMPLE),mean)
U_tail <- aggregate(all2$T_frac, by=list(SAMPLE=all2$SAMPLE),mean)
C_tail <- aggregate(all2$C_frac, by=list(SAMPLE=all2$SAMPLE),mean)
G_tail <- aggregate(all2$G_frac, by=list(SAMPLE=all2$SAMPLE),mean)

