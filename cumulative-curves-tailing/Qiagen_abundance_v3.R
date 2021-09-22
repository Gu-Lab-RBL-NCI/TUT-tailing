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

WT <- unique(iso_summary[which(iso_summary$SAMPLE=="1_WT_293T_1"),c(1,3,5)])
colnames(WT) <- c("MIRNA","SEQUENCE","WT")
all <- merge(all, WT, by=c("MIRNA","SEQUENCE"), all=true_false)
WT2 <- unique(iso_summary[which(iso_summary$SAMPLE=="2_WT_293T_2"),c(1,3,5)])
colnames(WT2) <- c("MIRNA","SEQUENCE","WT2")
all <- merge(all, WT2, by=c("MIRNA","SEQUENCE"), all=true_false)

TUT2 <- unique(iso_summary[which(iso_summary$SAMPLE=="3_TUT2KO_S4"),c(1,3,5)])
colnames(TUT2) <- c("MIRNA","SEQUENCE","TUT2")
all <- merge(all, TUT2, by=c("MIRNA","SEQUENCE"), all=true_false)
TUT2_2 <- unique(iso_summary[which(iso_summary$SAMPLE=="4_TUT2KO_S6"),c(1,3,5)])
colnames(TUT2_2) <- c("MIRNA","SEQUENCE","TUT2_2")
all <- merge(all, TUT2_2, by=c("MIRNA","SEQUENCE"), all=true_false)

TUT4 <- unique(iso_summary[which(iso_summary$SAMPLE=="5_TUT4KO_S1"),c(1,3,5)])
colnames(TUT4) <- c("MIRNA","SEQUENCE","TUT4")
all <- merge(all, TUT4, by=c("MIRNA","SEQUENCE"), all=true_false)
TUT4_2 <- unique(iso_summary[which(iso_summary$SAMPLE=="6_TUT4KO_S2"),c(1,3,5)])
colnames(TUT4_2) <- c("MIRNA","SEQUENCE","TUT4_2")
all <- merge(all, TUT4_2, by=c("MIRNA","SEQUENCE"), all=true_false)

TUT7 <- unique(iso_summary[which(iso_summary$SAMPLE=="7_TUT7KO_combined"),c(1,3,5)])
colnames(TUT7) <- c("MIRNA","SEQUENCE","TUT7")
all <- merge(all, TUT7, by=c("MIRNA","SEQUENCE"), all=true_false)
TUT7_2 <- unique(iso_summary[which(iso_summary$SAMPLE=="8_TUT7KO_S2"),c(1,3,5)])
colnames(TUT7_2) <- c("MIRNA","SEQUENCE","TUT7_2")
all <- merge(all, TUT7_2, by=c("MIRNA","SEQUENCE"), all=true_false)

DKO <- unique(iso_summary[which(iso_summary$SAMPLE=="10_DKO-24-2-IP-miseq"),c(1,3,5)])
colnames(DKO) <- c("MIRNA","SEQUENCE","DKO")
all <- merge(all, DKO, by=c("MIRNA","SEQUENCE"), all=true_false)
DKO2 <- unique(iso_summary[which(iso_summary$SAMPLE=="10_DKO_24_2"),c(1,3,5)])
colnames(DKO2) <- c("MIRNA","SEQUENCE","DKO2")
all <- merge(all, DKO2, by=c("MIRNA","SEQUENCE"), all=true_false)

TKO <- unique(iso_summary[which(iso_summary$SAMPLE=="11_TKO_7A"),c(1,3,5)])
colnames(TKO) <- c("MIRNA","SEQUENCE","TKO")
all <- merge(all, TKO, by=c("MIRNA","SEQUENCE"), all=true_false)
TKO2 <- unique(iso_summary[which(iso_summary$SAMPLE=="11_TKO-7A-IP_miseq"),c(1,3,5)])
colnames(TKO2) <- c("MIRNA","SEQUENCE","TKO2")
all <- merge(all, TKO2, by=c("MIRNA","SEQUENCE"), all=true_false)

rm(WT, WT2,TUT2,TUT2_2,TUT4,TUT4_2,TUT7,TUT7_2,DKO,DKO2, TKO,TKO2, agg_mirna,agg_sample, agg_mirseq, description)

all[is.na(all)] <- 0.5

all$WT <- all$WT/sum(all$WT)*1000000
all$WT2 <- all$WT2/sum(all$WT2)*1000000
all$TUT2 <- all$TUT2/sum(all$TUT2)*1000000
all$TUT2_2 <- all$TUT2_2/sum(all$TUT2_2)*1000000
all$TUT4 <- all$TUT4/sum(all$TUT4)*1000000
all$TUT4_2 <- all$TUT4_2/sum(all$TUT4_2)*1000000
all$TUT7 <- all$TUT7/sum(all$TUT7)*1000000
all$TUT7_2 <- all$TUT7_2/sum(all$TUT7_2)*1000000
all$DKO <- all$DKO/sum(all$DKO)*1000000
all$DKO2 <- all$DKO2/sum(all$DKO2)*1000000
all$TKO <- all$TKO/sum(all$TKO)*1000000
all$TKO2 <- all$TKO2/sum(all$TKO2)*1000000
all$avgWT <- (all$WT+all$WT2)/2

all$FC_TUT2 <- all$TUT2/all$avgWT
all$FC_TUT2_2 <- all$TUT2_2/all$avgWT
all$FC_TUT4 <- all$TUT4/all$avgWT
all$FC_TUT4_2 <- all$TUT4_2/all$avgWT
all$FC_TUT7 <- all$TUT7/all$avgWT
all$FC_TUT7_2 <- all$TUT7_2/all$avgWT
all$FC_DKO <- all$DKO/all$avgWT
all$FC_DKO2 <- all$DKO2/all$avgWT
all$FC_TKO <- all$TKO/all$avgWT
all$FC_TKO2 <- all$TKO2/all$avgWT

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

#templated non-templated mono
all$templated_mono <- NA
i <- 1
for(i in 1:nrow(all)){
  mir_test <- all$MIRNA[i]
  seq_test <- all$SEQUENCE[i]
  seq_test <- substr(seq_test,1,nchar(seq_test)-1)
  miRBase_test <- miRBase[which(miRBase$MIRNA==mir_test),]
  miRBase_test$find <- regexpr(seq_test, miRBase_test$EXTENDED.SEQUENCE)> -1
  miRBase_test <- miRBase_test[which(miRBase_test$find==T),]
  if(nrow(miRBase_test)>0){
    all$templated_mono[i] <- T
  }
  if(nrow(miRBase_test)==0){
    all$templated_mono[i] <- F
  }
  print(i/nrow(all)*100)
  rm(mir_test,seq_test,miRBase_test)
}


all <- all[order(-all$COUNTS),]
all$MIRNA <- gsub(x = all$MIRNA,"hsa-miR-375","hsa-miR-375-3p")
arm$MIRNA <- gsub(x = arm$MIRNA, "hsa-miR-378c","hsa-miR-378c-5p")
arm$MIRNA <- gsub(x = arm$MIRNA, "hsa-miR-548s","hsa-miR-548s-3p")
all$MIRNA <- gsub(x = all$MIRNA,"mir","miR")
all <- merge(all,arm,by="MIRNA")

all$log2FC_TUT2 <- log2((all$FC_TUT2+all$FC_TUT2_2)/2)
all$log2FC_TUT4 <- log2((all$FC_TUT4+all$FC_TUT4_2)/2)
all$log2FC_TUT7 <- log2((all$FC_TUT7+all$FC_TUT7_2)/2)
all$log2FC_DKO <- log2((all$FC_DKO+all$FC_DKO2)/2)
all$log2FC_TKO <- log2((all$FC_TKO+all$FC_TKO2)/2)

all_nontemplated <- all[which(all$templated==F & all$templated_mono==T & all$COUNTS>=10),]
all_templated <- all[which(all$templated==T & all$COUNTS>=10),]

breaks <-  seq(-10, 10, by=0.01)

canonical<- all_templated[which(all_templated$LEN_TRIM==0 & all_templated$SEQ_TAIL=="-"),]
#canonical<- all_templated[which(all_templated$SEQ_TAIL=="-"),]
canonical <- canonical[order(-canonical$COUNTS), ]

#Adenylation
##TUT2KO
monoA_nontemplated <- all_nontemplated[which(all_nontemplated$LEN_TRIM==0 & all_nontemplated$SEQ_TAIL=="A"),]
monoA_nontemplated <- monoA_nontemplated[order(-monoA_nontemplated$COUNTS), ]

monoA_templated <- all_templated[which(all_templated$LEN_TRIM==0 & all_templated$SEQ_TAIL=="A"),]
monoA_templated <- monoA_templated[order(-monoA_templated$COUNTS), ]

par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_TUT2)
monoA_nontemplated_TUT2KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
plot(breaks, monoA_nontemplated_TUT2KO, type="n", xlab = "mono-adenylated isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoA_nontemplated_TUT2KO, col="red")

datLog2FC <-  as.numeric(monoA_templated$log2FC_TUT2)
monoA_templated_TUT2KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_templated))) 
lines(breaks, monoA_templated_TUT2KO, col="black")

datLog2FC <-  as.numeric(canonical$log2FC_TUT2)
canonical_TUT2KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_TUT2KO, col="black")

wilcox.test(monoA_nontemplated$log2FC_TUT2,canonical$log2FC_TUT2)[3]
wilcox.test(monoA_templated$log2FC_TUT2,canonical$log2FC_TUT2)[3]
wilcox.test(monoA_templated$log2FC_TUT2,monoA_nontemplated$log2FC_TUT2)[3]
setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/cumulative-curves-tailing/")
write.table(rbind(breaks, monoA_templated_TUT2KO, monoA_nontemplated_TUT2KO,canonical_TUT2KO), "cumulative-Adenylation-TU2KO.txt", sep="\t", append = FALSE)

##TKO
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_TKO)
monoA_nontemplated_TKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
plot(breaks, monoA_nontemplated_TKO, type="n", xlab = "mono-adenylated isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoA_nontemplated_TKO, col="red")

datLog2FC <-  as.numeric(monoA_templated$log2FC_TKO)
monoA_templated_TKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_templated))) 
lines(breaks, monoA_templated_TKO, col="black")

datLog2FC <-  as.numeric(canonical$log2FC_TKO)
canonical_TKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_TKO, col="black")

wilcox.test(monoA_nontemplated$log2FC_TKO,canonical$log2FC_TKO)[3]
wilcox.test(monoA_templated$log2FC_TKO,canonical$log2FC_TKO)[3]
wilcox.test(monoA_nontemplated$log2FC_TKO,monoA_templated$log2FC_TKO)[3]
write.table(rbind(breaks, monoA_nontemplated_TKO, monoA_templated_TKO, canonical_TKO), "cumulative-Adenylation-TKO.txt", sep="\t", append = FALSE)

##TUT4-Andenylation
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_TUT4)
monoA_nontemplated_TUT4KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
plot(breaks, monoA_nontemplated_TUT4KO, type="n", xlab = "mono-adenylated isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoA_nontemplated_TUT4KO, col="red")

datLog2FC <-  as.numeric(monoA_templated$log2FC_TUT4)
monoA_templated_TUT4KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_templated))) 
lines(breaks, monoA_templated_TUT4KO, col="black")

datLog2FC <-  as.numeric(canonical$log2FC_TUT4)
canonical_TUT4KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_TUT4KO, col="black")

wilcox.test(monoA_nontemplated$log2FC_TUT4,canonical$log2FC_TUT4)[3]
wilcox.test(monoA_templated$log2FC_TUT4,canonical$log2FC_TUT4)[3]
wilcox.test(monoA_nontemplated$log2FC_TUT4,monoA_templated$log2FC_TUT4)[3]
write.table(rbind(breaks, monoA_nontemplated_TUT4KO, monoA_templated_TUT4KO, canonical_TUT4KO), "cumulative-Adenylation-TUT4KO.txt", sep="\t", append = FALSE)

##TUT7-Andenylation
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_TUT7)
monoA_nontemplated_TUT7KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
plot(breaks, monoA_nontemplated_TUT7KO, type="n", xlab = "mono-adenylated isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoA_nontemplated_TUT7KO, col="red")

datLog2FC <-  as.numeric(monoA_templated$log2FC_TUT7)
monoA_templated_TUT7KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_templated))) 
lines(breaks, monoA_templated_TUT7KO, col="black")

datLog2FC <-  as.numeric(canonical$log2FC_TUT7)
canonical_TUT7KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_TUT7KO, col="black")

wilcox.test(monoA_nontemplated$log2FC_TUT7,canonical$log2FC_TUT7)[3]
wilcox.test(monoA_templated$log2FC_TUT7,canonical$log2FC_TUT7)[3]
wilcox.test(monoA_nontemplated$log2FC_TUT7,monoA_templated$log2FC_TUT7)[3]
write.table(rbind(breaks, monoA_nontemplated_TUT7KO, monoA_templated_TUT7KO, canonical_TUT7KO), "cumulative-Adenylation-TUT7KO.txt", sep="\t", append = FALSE)

##DKO-Andenylation
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_DKO)
monoA_nontemplated_DKOKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
plot(breaks, monoA_nontemplated_DKOKO, type="n", xlab = "mono-adenylated isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoA_nontemplated_DKOKO, col="red")

datLog2FC <-  as.numeric(monoA_templated$log2FC_DKO)
monoA_templated_DKOKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_templated))) 
lines(breaks, monoA_templated_DKOKO, col="black")

datLog2FC <-  as.numeric(canonical$log2FC_DKO)
canonical_DKOKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_DKOKO, col="black")

wilcox.test(monoA_nontemplated$log2FC_DKO,canonical$log2FC_DKO)[3]
wilcox.test(monoA_templated$log2FC_DKO,canonical$log2FC_DKO)[3]
wilcox.test(monoA_nontemplated$log2FC_DKO,monoA_templated$log2FC_DKO)[3]
write.table(rbind(breaks, monoA_nontemplated_DKOKO, monoA_templated_DKOKO, canonical_DKOKO), "cumulative-Adenylation-DKOKO.txt", sep="\t", append = FALSE)







#Uridylation
monoU_nontemplated <- all_nontemplated[which(all_nontemplated$LEN_TRIM==0 & all_nontemplated$SEQ_TAIL=="T"),]
monoU_nontemplated <- monoU_nontemplated[order(-monoU_nontemplated$COUNTS), ]

monoU_templated <- all_templated[which(all_templated$LEN_TRIM==0 & all_templated$SEQ_TAIL=="T"),]
monoU_templated <- monoU_templated[order(-monoU_templated$COUNTS), ]

#monoU inTUT2KO
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoU_nontemplated$log2FC_TUT2)
monoU_nontemplated_TUT2KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_nontemplated))) 
plot(breaks, monoU_nontemplated_TUT2KO, type="n", xlab = "mono-adenylated isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoU_nontemplated_TUT2KO, col="red")

datLog2FC <-  as.numeric(monoU_templated$log2FC_TUT2)
monoU_templated_TUT2KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_templated))) 
lines(breaks, monoU_templated_TUT2KO, col="black")

datLog2FC <-  as.numeric(canonical$log2FC_TUT2)
canonical_TUT2KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_TUT2KO, col="black")

wilcox.test(monoU_nontemplated$log2FC_TUT2,canonical$log2FC_TUT2)[3]
wilcox.test(monoU_templated$log2FC_TUT2,canonical$log2FC_TUT2)[3]
wilcox.test(monoU_templated$log2FC_TUT2,monoU_nontemplated$log2FC_TUT2)[3]
setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/cumulative-curves-tailing/")
write.table(rbind(breaks, monoU_templated_TUT2KO, monoU_nontemplated_TUT2KO,canonical_TUT2KO), "cumulative-Uridylation-TU2KO.txt", sep="\t", append = FALSE)


##TUT4KO
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoU_nontemplated$log2FC_TUT4) #canonical FC
monoU_nontemplated_TUT4KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_nontemplated))) 
plot(breaks, monoU_nontemplated_TUT4KO, type="n", xlab = "mono-U isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoU_nontemplated_TUT4KO, col="red")

datLog2FC <-  as.numeric(monoU_templated$log2FC_TUT4) #canonical FC
monoU_templated_TUT4KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_templated))) 
lines(breaks, monoU_templated_TUT4KO, col="black")
wilcox.test(monoU_nontemplated$log2FC_TUT4,monoU_templated$log2FC_TUT4)[3]

datLog2FC <-  as.numeric(canonical$log2FC_TUT4)
canonical_TUT4KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_TUT4KO, col="black")

wilcox.test(monoU_nontemplated$log2FC_TUT4,canonical$log2FC_TUT4)[3]
wilcox.test(monoU_templated$log2FC_TUT4,canonical$log2FC_TUT4)[3]
wilcox.test(monoU_templated$log2FC_TUT4,monoU_nontemplated$log2FC_TUT4)[3]
setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/cumulative-curves-tailing/")
write.table(rbind(breaks, monoU_templated_TUT4KO, monoU_nontemplated_TUT4KO,canonical_TUT4KO), "cumulative-Uridylation-TUT4KO.txt", sep="\t", append = FALSE)

##TUT4KO monoA
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_TUT4) #canonical FC
monoU_nontemplated_TUT4KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
plot(breaks, monoU_nontemplated_TUT4KO, type="n", xlab = "mono-U isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoU_nontemplated_TUT4KO, col="red")

datLog2FC <-  as.numeric(monoA_templated$log2FC_TUT4) #canonical FC
monoU_templated_TUT4KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
lines(breaks, monoU_templated_TUT4KO, col="black")
wilcox.test(monoU_nontemplated$log2FC_TUT4,monoU_templated$log2FC_TUT4)[3]

datLog2FC <-  as.numeric(canonical$log2FC_TUT4)
canonical_TUT4KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_TUT4KO, col="black")

wilcox.test(monoU_nontemplated$log2FC_TUT4,canonical$log2FC_TUT4)[3]
wilcox.test(monoU_templated$log2FC_TUT4,canonical$log2FC_TUT4)[3]
wilcox.test(monoU_templated$log2FC_TUT4,monoU_nontemplated$log2FC_TUT4)[3]
setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/cumulative-curves-tailing/")
write.table(rbind(breaks, monoU_templated_TUT4KO, monoU_nontemplated_TUT4KO,canonical_TUT4KO), "cumulative-Uridylation-TUT4KO.txt", sep="\t", append = FALSE)


##TUT7KO
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoU_nontemplated$log2FC_TUT7) #canonical FC
monoU_nontemplated_TUT7KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_nontemplated))) 
plot(breaks, monoU_nontemplated_TUT7KO, type="n", xlab = "mono-U isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoU_nontemplated_TUT7KO, col="red")

datLog2FC <-  as.numeric(monoU_templated$log2FC_TUT7) #canonical FC
monoU_templated_TUT7KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_templated))) 
lines(breaks, monoU_templated_TUT7KO, col="black")
wilcox.test(monoU_nontemplated$log2FC_TUT7,monoU_templated$log2FC_TUT7)[3]

datLog2FC <-  as.numeric(canonical$log2FC_TUT7)
canonical_TUT7KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_TUT7KO, col="black")

wilcox.test(monoU_nontemplated$log2FC_TUT7,canonical$log2FC_TUT7)[3]
wilcox.test(monoU_templated$log2FC_TUT7,canonical$log2FC_TUT7)[3]
wilcox.test(monoU_templated$log2FC_TUT7,monoU_nontemplated$log2FC_TUT7)[3]
setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/cumulative-curves-tailing/")
write.table(rbind(breaks, monoU_templated_TUT7KO, monoU_nontemplated_TUT7KO,canonical_TUT7KO), "cumulative-Uridylation-TUT7KO.txt", sep="\t", append = FALSE)

##DKO
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoU_nontemplated$log2FC_DKO) #canonical FC
monoU_nontemplated_DKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_nontemplated))) 
plot(breaks, monoU_nontemplated_DKO, type="n", xlab = "mono-U isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoU_nontemplated_DKO, col="red")

datLog2FC <-  as.numeric(monoU_templated$log2FC_DKO) #canonical FC
monoU_templated_DKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_templated))) 
lines(breaks, monoU_templated_DKO, col="black")
wilcox.test(monoU_nontemplated$log2FC_D,monoU_templated$log2FC_D)[3]

datLog2FC <-  as.numeric(canonical$log2FC_DKO)
canonical_DKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_DKO, col="black")

wilcox.test(monoU_nontemplated$log2FC_DKO,canonical$log2FC_DKO)[3]
wilcox.test(monoU_templated$log2FC_DKO,canonical$log2FC_DKO)[3]
wilcox.test(monoU_templated$log2FC_DKO,monoU_nontemplated$log2FC_DKO)[3]
setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/cumulative-curves-tailing/")
write.table(rbind(breaks, monoU_templated_DKO, monoU_nontemplated_DKO,canonical_DKO), "cumulative-Uridylation-DKO.txt", sep="\t", append = FALSE)

##TKO
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoU_nontemplated$log2FC_TKO) #canonical FC
monoU_nontemplated_TKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_nontemplated))) 
plot(breaks, monoU_nontemplated_TKO, type="n", xlab = "mono-U isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoU_nontemplated_TKO, col="red")

datLog2FC <-  as.numeric(monoU_templated$log2FC_TKO) #canonical FC
monoU_templated_TKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_templated))) 
lines(breaks, monoU_templated_TKO, col="black")
wilcox.test(monoU_nontemplated$log2FC_D,monoU_templated$log2FC_D)[3]

datLog2FC <-  as.numeric(canonical$log2FC_TKO)
canonical_TKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(canonical))) 
lines(breaks, canonical_TKO, col="black")

wilcox.test(monoU_nontemplated$log2FC_TKO,canonical$log2FC_TKO)[3]
wilcox.test(monoU_templated$log2FC_TKO,canonical$log2FC_TKO)[3]
wilcox.test(monoU_templated$log2FC_TKO,monoU_nontemplated$log2FC_TKO)[3]
setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/cumulative-curves-tailing/")
write.table(rbind(breaks, monoU_templated_TKO, monoU_nontemplated_TKO,canonical_TKO), "cumulative-Uridylation-TKO.txt", sep="\t", append = FALSE)







##TUT7KO
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_TUT7) #canonical FC
monoA_nontemplated_TUT7KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
plot(breaks, monoA_nontemplated_TUT7KO, type="n", xlab = "mono-adenylated isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoA_nontemplated_TUT7KO, col="red")

datLog2FC <-  as.numeric(monoA_templated$log2FC_TUT7) #canonical FC
monoA_templated_TUT7KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_templated))) 
lines(breaks, monoA_templated_TUT7KO, col="black")
wilcox.test(monoA_nontemplated$log2FC_TUT7,monoA_templated$log2FC_TUT7)[3]

##DKO
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_DKO) #canonical FC
monoA_nontemplated_DKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
plot(breaks, monoA_nontemplated_DKO, type="n", xlab = "mono-adenylated isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoA_nontemplated_DKO, col="red")

datLog2FC <-  as.numeric(monoA_templated$log2FC_DKO) #canonical FC
monoA_templated_DKO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_templated))) 
lines(breaks, monoA_templated_DKO, col="black")
wilcox.test(monoA_nontemplated$log2FC_DKO,monoA_templated$log2FC_DKO)[3]


#Canonical
##TUT2KO
monoA_nontemplated <- all_nontemplated[which(all_nontemplated$LEN_TRIM==0 & all_nontemplated$SEQ_TAIL=="A"),]
monoA_nontemplated <- monoA_nontemplated[order(-monoA_nontemplated$COUNTS), ]




par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_TUT2)
monoA_nontemplated_TUT2KO = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
plot(breaks, monoA_nontemplated_TUT2KO, type="n", xlab = "mono-adenylated isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, monoA_nontemplated_TUT2KO, col="red")













datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_TUT4) #canonical FC
notargets2 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
lines(breaks, notargets2, col="blue")

datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_TUT7) #canonical FC
notargets2 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
lines(breaks, notargets2, col="green")

datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_DKO) #canonical FC
notargets2 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
lines(breaks, notargets2, col="brown")

datLog2FC <-  as.numeric(monoA_nontemplated$log2FC_TKO) #canonical FC
notargets2 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoA_nontemplated))) 
lines(breaks, notargets2, col="black")

#Uridylation
monoU_nontemplated <- all_nontemplated[which(all_nontemplated$LEN_TRIM==0 & all_nontemplated$SEQ_TAIL=="T"),]
monoU_nontemplated <- monoU_nontemplated[order(-monoU_nontemplated$COUNTS), ]
monoU_nontemplated <- monoU_nontemplated[c(1:500),]

par(mfrow=c(1,1))
datLog2FC <-  as.numeric(monoU_nontemplated$log2FC_TUT2) #canonical FC
notargets1 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_nontemplated))) 
plot(breaks, notargets1, type="n", xlab = "mono-Uridylation isomiRs fold-change Log2 (KO/WT)", ylab = "Cumulative fraction",ylim=c(0,1), xlim=c(-5,5))
lines(breaks, notargets1, col="red")

datLog2FC <-  as.numeric(monoU_nontemplated$log2FC_TUT4) #canonical FC
notargets2 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_nontemplated))) 
lines(breaks, notargets2, col="blue")

datLog2FC <-  as.numeric(monoU_nontemplated$log2FC_TUT7) #canonical FC
notargets2 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_nontemplated))) 
lines(breaks, notargets2, col="green")

datLog2FC <-  as.numeric(monoU_nontemplated$log2FC_DKO) #canonical FC
notargets2 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_nontemplated))) 
lines(breaks, notargets2, col="brown")

datLog2FC <-  as.numeric(monoU_nontemplated$log2FC_TKO) #canonical FC
notargets2 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(monoU_nontemplated))) 
lines(breaks, notargets2, col="black")
