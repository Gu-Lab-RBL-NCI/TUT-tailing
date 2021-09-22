setwd("~/Desktop/P1 Pri-miRNA scaffold/Analysis cleavage Bioinformatic/")
miRBase <- read.table(file="high_confidence_miRBase21-master.tsv", header = TRUE, fill = TRUE, sep =  " ")
arm <- unique(miRBase[,c(1,16)])
arm[which(arm$MIRNA=="hsa-miR-548bc-3p"),] <- "3P"
rm(miRBase)

setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/")
require(data.table)
miRNA_summary <- fread(file = "TUT_KO_Qiagen_rerun_2_5.isomir.tsv", header = T)
miRNA_summary <- unique(miRNA_summary)
miRNA_summary$CPM <- miRNA_summary$TOTAL_READS/miRNA_summary$TOTAL_READS_IN_SAMPLE*1000000
miRNA_summary$CPM_miRNA <- miRNA_summary$TOTAL_READS/sum(miRNA_summary$TOTAL_READS)*1000000

miRNA_summary$SAMPLE <- gsub(".fastq_ready","",miRNA_summary$SAMPLE)

miRNA_summary$MIRNA <- gsub(x = miRNA_summary$MIRNA,"hsa-miR-375","hsa-miR-375-3p")
arm$MIRNA <- gsub(x = arm$MIRNA, "hsa-miR-378c","hsa-miR-378c-5p")
arm$MIRNA <- gsub(x = arm$MIRNA, "hsa-miR-548s","hsa-miR-548s-3p")
miRNA_summary$MIRNA <- gsub(x = miRNA_summary$MIRNA,"mir","miR")


all <- unique(miRNA_summary[,c(2,4)])
all <- merge(all, arm, by="MIRNA")

true_false <- TRUE
WT1 <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="1_WT_293T_1"),c(2,19)])
colnames(WT1) <- c("MIRNA","WT1")
all <- merge(all, WT1, by="MIRNA", all=true_false)
WT2 <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="2_WT_293T_2"),c(2,19)])
colnames(WT2) <- c("MIRNA","WT2")
all <- merge(all, WT2, by="MIRNA", all=true_false)

TUT2 <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="3_TUT2KO_S4"),c(2,19)])
colnames(TUT2) <- c("MIRNA","TUT2")
all <- merge(all, TUT2, by="MIRNA", all=true_false)
TUT2_2 <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="4_TUT2KO_S6"),c(2,19)])
colnames(TUT2_2) <- c("MIRNA","TUT2_2")
all <- merge(all, TUT2_2, by="MIRNA", all=true_false)

TUT4 <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="5_TUT4KO_S1"),c(2,19)])
colnames(TUT4) <- c("MIRNA","TUT4")
all <- merge(all, TUT4, by="MIRNA", all=true_false)
TUT4_2 <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="6_TUT4KO_S2"),c(2,19)])
colnames(TUT4_2) <- c("MIRNA","TUT4_2")
all <- merge(all, TUT4_2, by="MIRNA", all=true_false)

TUT7 <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="7_TUT7KO_combined"),c(2,19)])
colnames(TUT7) <- c("MIRNA","TUT7")
all <- merge(all, TUT7, by="MIRNA", all=true_false)
TUT7_2 <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="8_TUT7KO_S2"),c(2,19)])
colnames(TUT7_2) <- c("MIRNA","TUT7_2")
all <- merge(all, TUT7_2, by="MIRNA", all=true_false)

DKO <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="10_DKO-24-2-IP-miseq"),c(2,19)])
colnames(DKO) <- c("MIRNA","DKO")
all <- merge(all, DKO, by="MIRNA", all=true_false)
DKO2 <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="10_DKO_24_2"),c(2,19)])
colnames(DKO2) <- c("MIRNA","DKO2")
all <- merge(all, DKO2, by="MIRNA", all=true_false)

TKO <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="11_TKO_7A"),c(2,19)])
colnames(TKO) <- c("MIRNA","TKO")
all <- merge(all, TKO, by="MIRNA", all=true_false)
TKO2 <- unique(miRNA_summary[which(miRNA_summary$SAMPLE=="11_TKO-7A-IP_miseq"),c(2,19)])
colnames(TKO2) <- c("MIRNA","TKO2")
all <- merge(all, TKO2, by="MIRNA", all=true_false)

rm(WT1, WT2,TUT2,TUT2_2,TUT4,TUT4_2,TUT7,TUT7_2,DKO,DKO2, TKO,TKO2, arm)

all[is.na(all)] <- 0
#all[is.na(all)] <- 0.001

all$mWT <- (all$WT1+all$WT2)/2
all$mTUT2 <- (all$TUT2+all$TUT2_2)/2
all$mTUT4 <- (all$TUT4+all$TUT4_2)/2
all$mTUT7 <- (all$TUT7+all$TUT7_2)/2
all$mDKO <- (all$DKO+all$DKO2)/2
all$mTKO <- (all$TKO+all$TKO2)/2


all$FcTUT2 <- all$mTUT2/all$mWT
all$FcTUT4 <- all$mTUT4/all$mWT
all$FcTUT7 <- all$mTUT7/all$mWT
all$FcDKO <- all$mDKO/all$mWT
all$FcTKO <- all$mTKO/all$mWT

all <- all[order(-all$mWT),]

n <- 1800
n <- nrow(all)

final <- all[c(1:n),]
plot(log10(final$mWT),log2(final$FcDKO))
plot(log10(final$mWT),log2(final$FcTKO))


write.table(final, "expression.tsv", append = F, row.names = F, sep = "\t")
