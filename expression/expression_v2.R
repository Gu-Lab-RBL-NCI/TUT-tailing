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


all$FcTUT2 <- (all$mTUT2/all$mWT)
all$FcTUT4 <- (all$mTUT4/all$mWT)
all$FcTUT7 <- (all$mTUT7/all$mWT)
all$FcDKO <- (all$mDKO/all$mWT)
all$FcTKO <- (all$mTKO/all$mWT)

all <- all[order(-all$mWT),]

n <- 1800
n <- nrow(all)

final <- all[c(1:n),]
final <- final[which(final$FcTUT2!=0),]
final <- final[which(final$FcTUT4!=0),]
final <- final[which(final$FcTUT7!=0),]
final <- final[which(final$FcDKO!=0),]
final <- final[which(final$FcTKO!=0),]

final$log2FcTUT2 <- log2(final$mTUT2/final$mWT)
final$log2FcTUT4 <- log2(final$mTUT4/final$mWT)
final$log2FcTUT7 <- log2(final$mTUT7/final$mWT)
final$log2FcDKO <- log2(final$mDKO/final$mWT)
final$log2FcTKO <- log2(final$mTKO/final$mWT)

plot(log10(final$mWT),log2(final$FcDKO))
plot(log10(final$mWT),(final$log2FcDKO))

plot(log10(final$mWT),log2(final$FcTKO))

corrTUT2 <- 0-median(final$log2FcTUT2[c(1:320)]) #upper quartile
final$log2FcTUT2 <- final$log2FcTUT2+corrTUT2

plot(log10(final$mWT),log2(final$FcTUT2))
abline(h=0)
plot(log10(final$mWT),(final$log2FcTUT2))
abline(h=0)

corrTUT4 <- 0-median(final$log2FcTUT4[c(1:320)])
final$log2FcTUT4 <- final$log2FcTUT4+corrTUT4

corrTUT7 <- 0-median(final$log2FcTUT7[c(1:320)])
final$log2FcTUT7 <- final$log2FcTUT7+corrTUT7

corrDKO <- 0-median(final$log2FcDKO[c(1:320)])
final$log2FcDKO <- final$log2FcDKO+corrDKO

corrTKO <- 0-median(final$log2FcTKO[c(1:320)])
final$log2FcTKO <- final$log2FcTKO+corrTKO

TUT2model <- final[,c(16,27)]
TUT2model <- TUT2model[which(TUT2model$mWT>0),]
TUT2model <- unique(TUT2model)
TUT2model$upper <- NA
i <- 1
for (i in 1:nrow(TUT2model)) {
  a <- i+200
  sd_test <- sd(TUT2model$log2FcTUT2[c(i:a)])
  sd_test <- sd_test*2
  TUT2model$upper[i] <- sd_test
  rm(sd_test)
  print(i)
}
TUT2model$lower <- 0-TUT2model$upper


set.seed(123)
donwsample <- TUT2model[complete.cases(TUT2model),]
index <- sample(x = 1:nrow(donwsample), size = 500)
index <- c(1,nrow(donwsample),index)
index <- unique(index)
donwsample <- donwsample[index, ]
donwsample <- donwsample[order(-donwsample$mWT),]
plot(log10(TUT2model$mWT),TUT2model$log2FcTUT2, cex=0.7, ylim = c(-10,10))
abline(h=0)
lines(log10(donwsample$mWT),donwsample$upper, col="red")
lines(log10(donwsample$mWT),donwsample$lower, col="red")

donwsample <- donwsample[,c(1,3,4)]

plot(log10(final$mWT),final$log2FcDKO, cex=0.7, ylim = c(-10,10))
abline(h=0)
lines(log10(donwsample$mWT),donwsample$upper, col="red")
lines(log10(donwsample$mWT),donwsample$lower, col="red")

TUT2model <- TUT2model[complete.cases(TUT2model),]
TUT2model <- aggregate(TUT2model$upper, by=list(mWT=TUT2model$mWT),median)
TUT2model$lower <- 0-TUT2model$x
TUT2model$upper <- 0+TUT2model$x
TUT2model <- TUT2model[,c(1,4,3)]
TUT2model <- TUT2model[which(TUT2model$mWT>=1),]
TUT2model <- merge(TUT2model, final, by="mWT")
TUT2_change <- TUT2model[which(TUT2model$log2FcTUT2>=TUT2model$upper | TUT2model$log2FcTUT2<=TUT2model$lower),]
TUT4_change <- TUT2model[which(TUT2model$log2FcTUT4>=TUT2model$upper | TUT2model$log2FcTUT4<=TUT2model$lower),]
TUT7_change <- TUT2model[which(TUT2model$log2FcTUT7>=TUT2model$upper | TUT2model$log2FcTUT7<=TUT2model$lower),]
DKO_change <- TUT2model[which(TUT2model$log2FcDKO>=TUT2model$upper | TUT2model$log2FcDKO<=TUT2model$lower),]
TKO_change <- TUT2model[which(TUT2model$log2FcTKO>=TUT2model$upper | TUT2model$log2FcTKO<=TUT2model$lower),]
TUT2_change$change <- "TUT2KO" 
TUT4_change$change <- "TUT4KO" 
TUT7_change$change <- "TUT7KO" 
DKO_change$change <- "DKO" 
TKO_change$change <- "TKO" 
all_changes <- rbind(TUT2_change,TUT4_change,TUT7_change,DKO_change,TKO_change)
all_changes <- all_changes[,c(4,34)]
write.table(all_changes, "venn.tsv", append = F, row.names = F, sep = "\t")
write.table(final, "expression2.tsv", append = F, row.names = F, sep = "\t")




write.table(TUT2model, "TUT2model.tsv", append = F, row.names = F, sep = "\t")



