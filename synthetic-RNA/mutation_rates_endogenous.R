setwd("~/Desktop/P1 Pri-miRNA scaffold/Analysis cleavage Bioinformatic/")
miRBase <- read.table(file="high_confidence_miRBase21-master.tsv", header = TRUE, fill = TRUE, sep =  " ")

setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/")
require(data.table)
iso_summary <- fread(file = "TUT_KO_Qiagen_rerun_2_5.isomir.sequence_info.tsv", header = T)
iso_summary$SAMPLE <- gsub(".fastq_ready","",iso_summary$SAMPLE)

unique(iso_summary$SAMPLE)
iso_summary <- iso_summary[which(iso_summary$SAMPLE=="1_WT_293T_1"|iso_summary$SAMPLE=="2_WT_293T_2"),]
iso_summary <- unique(iso_summary)

unmodified <- iso_summary[which(iso_summary$LEN_TRIM==0 & iso_summary$LEN_TAIL==0 & iso_summary$DISTANCE==0),]
max_mirna <- aggregate(unmodified$READS, by=list(MIRNA=unmodified$MIRNA,SEQUENCE=unmodified$SEQUENCE),sum)
max_mirna2 <- aggregate(max_mirna$x, by=list(MIRNA=max_mirna$MIRNA),max)
max_mirna <- merge(max_mirna, max_mirna2, by=c("MIRNA","x"))
max_mirna$LEN_READ <- nchar(max_mirna$SEQUENCE)

max_mirna$noA <- NA
max_mirna$noG <- NA
i <- 13
for(i in 1:nrow(max_mirna)){
  mir_test <- max_mirna$MIRNA[i]
  seq_test <- max_mirna$SEQUENCE[i]

  miRBase_test <- miRBase[which(miRBase$MIRNA==mir_test),]
  miRBase_test$findA <- regexpr(paste0(seq_test,"A"), miRBase_test$EXTENDED.SEQUENCE)> -1
  miRBase_test <- miRBase_test[which(miRBase_test$findA==T),]
  if(nrow(miRBase_test)==0){
    max_mirna$noA[i] <- T
  }
  if(nrow(miRBase_test)>0){
    max_mirna$noA[i] <- F
  }
  
  miRBase_test <- miRBase[which(miRBase$MIRNA==mir_test),]
  miRBase_test$findG <- regexpr(paste0(seq_test,"G"), miRBase_test$EXTENDED.SEQUENCE)> -1
  miRBase_test <- miRBase_test[which(miRBase_test$findG==T),]
  if(nrow(miRBase_test)==0){
    max_mirna$noG[i] <- T
  }
  if(nrow(miRBase_test)>0){
    max_mirna$noG[i] <- F
  }
  print(i/nrow(max_mirna)*100)
  rm(mir_test,seq_test,miRBase_test)
}

max_mirna <- max_mirna[which(max_mirna$noA==T & max_mirna$noG==T),]

i <- 1
seq_test <- max_mirna$SEQUENCE[1]
A_tail <- iso_summary[which(iso_summary$SEQUENCE==paste0(seq_test,"A")),c(1,2,5)]
colnames(A_tail) <- c("SAMPLE","MIRNA","READS_A")
G_tail <- iso_summary[which(iso_summary$SEQUENCE==paste0(seq_test,"G")),c(1,2,5)]
colnames(G_tail) <- c("SAMPLE","MIRNA","READS_G")
AtoG_final <- merge(A_tail,G_tail, by=c("SAMPLE","MIRNA"))

for(i in 1:nrow(max_mirna)){
  seq_test <- max_mirna$SEQUENCE[i]
  A_tail <- iso_summary[which(iso_summary$SEQUENCE==paste0(seq_test,"A")),c(1,2,5)]
  colnames(A_tail) <- c("SAMPLE","MIRNA","READS_A")
  G_tail <- iso_summary[which(iso_summary$SEQUENCE==paste0(seq_test,"G")),c(1,2,5)]
  colnames(G_tail) <- c("SAMPLE","MIRNA","READS_G")
  AtoG <- merge(A_tail,G_tail, by=c("SAMPLE","MIRNA"))
  AtoG_final <- rbind(AtoG_final,AtoG)
  print(i/nrow(max_mirna)*100)
}

AtoG_final <- unique(AtoG_final)
AtoG_final$rateAtoG <- AtoG_final$READS_G/(AtoG_final$READS_A+AtoG_final$READS_G)*100
AtoG_final$total <- AtoG_final$READS_A+AtoG_final$READS_G

AtoG_final2 <- aggregate(AtoG_final$rateAtoG, by=list(MIRNA=AtoG_final$MIRNA),mean)
AtoG_final3 <- aggregate(AtoG_final$total, by=list(MIRNA=AtoG_final$MIRNA),mean)

AtoG_final4 <- merge(AtoG_final2, AtoG_final3, by="MIRNA")
colnames(AtoG_final4) <- c("MIRNA","rateAtoG","total")
AtoG_final4 <- AtoG_final4[order(-AtoG_final4$total),]
AtoG_final4 <- AtoG_final4[c(1:100),]

setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/synthetic-RNA/")
write.table(AtoG_final4, "AtoG_non-templated-A-to-G-miRNA.tsv", sep="\t", append = FALSE, row.names = F)
