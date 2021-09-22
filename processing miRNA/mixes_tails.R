setwd("~/Desktop/P1 Pri-miRNA scaffold/Analysis cleavage Bioinformatic/")
miRBase <- read.table(file="high_confidence_miRBase21-master.tsv", header = TRUE, fill = TRUE, sep =  " ")
arm <- unique(miRBase[,c(1,16)])

setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/")
require(data.table)
iso_summary <- fread(file = "TUT_KO_Qiagen_rerun_2_5.isomir.sequence_info.tsv", header = T)
miRNA_summary <- fread(file = "TUT_KO_Qiagen_rerun_2_5.isomir.tsv", header = T)
miRNA_summary$CPM <- miRNA_summary$TOTAL_READS/miRNA_summary$TOTAL_READS_IN_SAMPLE*1000000
agg_mirna <- aggregate(iso_summary$READS, by=list(SAMPLE=iso_summary$SAMPLE, MIRNA=iso_summary$MIRNA),sum)
iso_summary <- merge(iso_summary, agg_mirna, by=c("MIRNA","SAMPLE"))
agg_sample <- aggregate(miRNA_summary$TOTAL_READS, by=list(SAMPLE=miRNA_summary$SAMPLE),sum)
miRNA_summary <- merge(miRNA_summary, agg_sample, by=c("SAMPLE"))
miRNA_summary$CPM_miRNA <- miRNA_summary$TOTAL_READS/miRNA_summary$x*1000000

iso_summary$ratio2 <- iso_summary$READS/iso_summary$x*100
iso_summary$SAMPLE <- gsub(".fastq_ready","",iso_summary$SAMPLE)
unique(iso_summary$SAMPLE)


true_false <- TRUE
agg_mirseq <- aggregate(iso_summary$READS, by=list(MIRNA=iso_summary$MIRNA,SEQUENCE=iso_summary$SEQUENCE),sum)
description <- unique(iso_summary[,c(1,3,4,7,8,9)])
all <- merge(agg_mirseq, description, by=c("MIRNA","SEQUENCE"), all=true_false)
colnames(all) <- c("MIRNA","SEQUENCE","COUNTS","LEN_READ","LEN_TRIM","LEN_TAIL","SEQ_TAIL")

WT <- unique(iso_summary[which(iso_summary$SAMPLE=="1_WT_293T_1"),c(1,3,16)])
colnames(WT) <- c("MIRNA","SEQUENCE","WT")
all <- merge(all, WT, by=c("MIRNA","SEQUENCE"), all=true_false)
WT2 <- unique(iso_summary[which(iso_summary$SAMPLE=="2_WT_293T_2"),c(1,3,16)])
colnames(WT2) <- c("MIRNA","SEQUENCE","WT2")
all <- merge(all, WT2, by=c("MIRNA","SEQUENCE"), all=true_false)
all[is.na(all)] <- 0

TUT2 <- unique(iso_summary[which(iso_summary$SAMPLE=="3_TUT2KO_S4"),c(1,3,16)])
colnames(TUT2) <- c("MIRNA","SEQUENCE","TUT2")
all <- merge(all, TUT2, by=c("MIRNA","SEQUENCE"), all=true_false)
TUT2_2 <- unique(iso_summary[which(iso_summary$SAMPLE=="4_TUT2KO_S6"),c(1,3,16)])
colnames(TUT2_2) <- c("MIRNA","SEQUENCE","TUT2_2")
all <- merge(all, TUT2_2, by=c("MIRNA","SEQUENCE"), all=true_false)

TUT4 <- unique(iso_summary[which(iso_summary$SAMPLE=="5_TUT4KO_S1"),c(1,3,16)])
colnames(TUT4) <- c("MIRNA","SEQUENCE","TUT4")
all <- merge(all, TUT4, by=c("MIRNA","SEQUENCE"), all=true_false)
TUT4_2 <- unique(iso_summary[which(iso_summary$SAMPLE=="6_TUT4KO_S2"),c(1,3,16)])
colnames(TUT4_2) <- c("MIRNA","SEQUENCE","TUT4_2")
all <- merge(all, TUT4_2, by=c("MIRNA","SEQUENCE"), all=true_false)

TUT7 <- unique(iso_summary[which(iso_summary$SAMPLE=="7_TUT7KO_combined"),c(1,3,16)])
colnames(TUT7) <- c("MIRNA","SEQUENCE","TUT7")
all <- merge(all, TUT7, by=c("MIRNA","SEQUENCE"), all=true_false)
TUT7_2 <- unique(iso_summary[which(iso_summary$SAMPLE=="8_TUT7KO_S2"),c(1,3,16)])
colnames(TUT7_2) <- c("MIRNA","SEQUENCE","TUT7_2")
all <- merge(all, TUT7_2, by=c("MIRNA","SEQUENCE"), all=true_false)

DKO <- unique(iso_summary[which(iso_summary$SAMPLE=="10_DKO-24-2-IP-miseq"),c(1,3,16)])
colnames(DKO) <- c("MIRNA","SEQUENCE","DKO")
all <- merge(all, DKO, by=c("MIRNA","SEQUENCE"), all=true_false)
DKO2 <- unique(iso_summary[which(iso_summary$SAMPLE=="10_DKO_24_2"),c(1,3,16)])
colnames(DKO2) <- c("MIRNA","SEQUENCE","DKO2")
all <- merge(all, DKO2, by=c("MIRNA","SEQUENCE"), all=true_false)

TKO <- unique(iso_summary[which(iso_summary$SAMPLE=="11_TKO_7A"),c(1,3,16)])
colnames(TKO) <- c("MIRNA","SEQUENCE","TKO")
all <- merge(all, TKO, by=c("MIRNA","SEQUENCE"), all=true_false)
TKO2 <- unique(iso_summary[which(iso_summary$SAMPLE=="11_TKO-7A-IP_miseq"),c(1,3,16)])
colnames(TKO2) <- c("MIRNA","SEQUENCE","TKO2")
all <- merge(all, TKO2, by=c("MIRNA","SEQUENCE"), all=true_false)

rm(WT, WT2,TUT2,TUT2_2,TUT4,TUT4_2,TUT7,TUT7_2,DKO,DKO2, TKO,TKO2, agg_mirna,agg_sample, agg_mirseq, description)
all[is.na(all)] <- 0
all$avg <- (all$WT+all$WT2)/2

all$name <- paste0(all$MIRNA, " " ,"[",all$SEQUENCE,"]")
all$rWT <- all$WT/all$avg*100
all$rWT2 <- all$WT2/all$avg*100
all$rTUT2 <- all$TUT2/all$avg*100
all$rTUT22 <- all$TUT2_2/all$avg*100
all$rTUT4 <- all$TUT4/all$avg*100
all$rTUT42 <- all$TUT4_2/all$avg*100
all$rTUT7 <- all$TUT7/all$avg*100
all$rTUT72 <- all$TUT7_2/all$avg*100
all$rDKO <- all$DKO/all$avg*100
all$rDKO2 <- all$DKO2/all$avg*100
all$rTKO <- all$TKO/all$avg*100
all$rTKO2 <- all$TKO2/all$avg*100
all <- all[order(-all$COUNTS),]

all$templated <- NA
#templated non-templated
i <- 1
for(i in 1:nrow(all)){
  mir_test <- all$MIRNA[i]
  seq_test <- all$SEQUENCE[i]
  miRBase_test <- miRBase[which(miRBase$MIRNA==mir_test),]
  miRBase_test$find <- regexpr(seq_test, miRBase_test$EXTENDED.SEQUENCE)> -1
  miRBase_test <- miRBase_test[which(miRBase_test$find==T),]
  if(nrow(miRBase_test)>0){
    all$templated[i] <- T
  }
  if(nrow(miRBase_test)==0){
    all$templated[i] <- F
  }
  print(i/nrow(all)*100)
  rm(mir_test,seq_test,miRBase_test)
}

mixed_tails <- all
mixed_tails <- mixed_tails[which(mixed_tails$templated==F),]
mixed_tails <- mixed_tails[which(mixed_tails$LEN_TAIL>=2 & mixed_tails$LEN_TRIM<=1),]
mixed_tails <- mixed_tails[,c(1:7,21:33)]
mixed_tails <- mixed_tails[which(mixed_tails$COUNTS>=10),]
mixed_tails$t1 <- NA
#templated non-templated
i <- 1
for(i in 1:nrow(mixed_tails)){
  mir_test <- mixed_tails$MIRNA[i]
  seq_test <- mixed_tails$SEQUENCE[i]
  seq_test <- substr(seq_test, 1, nchar(seq_test)-mixed_tails$LEN_TAIL[i])
  seq_test <- paste0(seq_test,substr(mixed_tails$SEQ_TAIL[i], 1, 1))
  miRBase_test <- miRBase[which(miRBase$MIRNA==mir_test),]
  miRBase_test$find <- regexpr(seq_test, miRBase_test$EXTENDED.SEQUENCE)> -1
  miRBase_test <- miRBase_test[which(miRBase_test$find==T),]
  if(nrow(miRBase_test)>0){
    mixed_tails$t1[i] <- T
  }
  if(nrow(miRBase_test)==0){
    mixed_tails$t1[i] <- F
  }
  print(i/nrow(mixed_tails)*100)
  rm(mir_test,seq_test,miRBase_test)
}
mixed_tails <- mixed_tails[which(mixed_tails$t1==F),]
tail_AA <- mixed_tails[which(mixed_tails$SEQ_TAIL=="AA"),]
tail_AA <- tail_AA[c(1:100),]
tail_UU <- mixed_tails[which(mixed_tails$SEQ_TAIL=="TT"),]
tail_UU <- tail_UU[c(1:100),]
tail_AU <- mixed_tails[which(mixed_tails$SEQ_TAIL=="AT"|mixed_tails$SEQ_TAIL=="TA"),]
tail_AU <- tail_AU[c(1:100),]
write.table(tail_AA, "tail_AA.txt", sep="\t", append = FALSE, row.names = F)
write.table(tail_UU, "tail_UU.txt", sep="\t", append = FALSE, row.names = F)
write.table(tail_AU, "tail_AU.txt", sep="\t", append = FALSE, row.names = F)
