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
all$mWT <- (all$WT1+all$WT2)/2

all$log2FcWT_1 <- log2(all$WT1/all$mWT)
all$log2FcWT_2 <- log2(all$WT2/all$mWT)
all$log2FcTUT2_1 <- log2(all$TUT2/all$mWT)
all$log2FcTUT2_2 <- log2(all$TUT2_2/all$mWT)
all$log2FcTUT4_1 <- log2(all$TUT4/all$mWT)
all$log2FcTUT4_2 <- log2(all$TUT4_2/all$mWT)
all$log2FcTUT7_1 <- log2(all$TUT7/all$mWT)
all$log2FcTUT7_2 <- log2(all$TUT7_2/all$mWT)
all$log2FcDKO_1 <- log2(all$DKO/all$mWT)
all$log2FcDKO_2 <- log2(all$DKO2/all$mWT)
all$log2FcTKO_1 <- log2(all$TKO/all$mWT)
all$log2FcTKO_2 <- log2(all$TKO2/all$mWT)

all <- all[order(-all$mWT),]
plot(log10(all$mWT),(all$log2FcTUT2_1))
abline(h=0)
corr <- 0-median(all$log2FcWT_1[c(1:320)]) #upper quartile
all$log2FcWT_1 <- all$log2FcWT_1+corr
corr <- 0-median(all$log2FcWT_2[c(1:320)]) #upper quartile
all$log2FcWT_2 <- all$log2FcWT_2+corr

corr <- 0-median(all$log2FcTUT2_1[c(1:320)]) #upper quartile
all$log2FcTUT2_1 <- all$log2FcTUT2_1+corr
corr <- 0-median(all$log2FcTUT2_2[c(1:320)]) #upper quartile
all$log2FcTUT2_2 <- all$log2FcTUT2_2+corr

corr <- 0-median(all$log2FcTUT4_1[c(1:320)]) #upper quartile
all$log2FcTUT4_1 <- all$log2FcTUT4_1+corr
corr <- 0-median(all$log2FcTUT4_2[c(1:320)]) #upper quartile
all$log2FcTUT4_2 <- all$log2FcTUT4_2+corr

corr <- 0-median(all$log2FcTUT7_1[c(1:320)]) #upper quartile
all$log2FcTUT7_1 <- all$log2FcTUT7_1+corr
corr <- 0-median(all$log2FcTUT7_2[c(1:320)]) #upper quartile
all$log2FcTUT7_2 <- all$log2FcTUT7_2+corr

corr <- 0-median(all$log2FcDKO_1[c(1:320)]) #upper quartile
all$log2FcDKO_1 <- all$log2FcDKO_1+corr
corr <- 0-median(all$log2FcDKO_2[c(1:320)]) #upper quartile
all$log2FcDKO_2 <- all$log2FcDKO_2+corr

corr <- 0-median(all$log2FcTKO_1[c(1:320)]) #upper quartile
all$log2FcTKO_1 <- all$log2FcTKO_1+corr
corr <- 0-median(all$log2FcTKO_2[c(1:320)]) #upper quartile
all$log2FcTKO_2 <- all$log2FcTKO_2+corr

n <- 1800
n <- nrow(all)

final <- all[c(1:n),]
final <- final[complete.cases(final),]
final <- final[!is.infinite(rowSums(final[,c(17:28)])),]

final$log2FcWT <- (final$log2FcWT_1+final$log2FcWT_2)/2
final$log2FcTUT2 <- (final$log2FcTUT2_1+final$log2FcTUT2_2)/2
final$log2FcTUT4 <- (final$log2FcTUT4_1+final$log2FcTUT4_2)/2
final$log2FcTUT7 <- (final$log2FcTUT7_1+final$log2FcTUT7_2)/2
final$log2FcDKO <- (final$log2FcDKO_1+final$log2FcDKO_2)/2
final$log2FcTKO <- (final$log2FcTKO_1+final$log2FcTKO_2)/2

TUT2model <- final[,c(16,19,20)]
TUT2model <- TUT2model[which(TUT2model$mWT>0),]
TUT2model <- unique(TUT2model)
TUT2model$avg <- (TUT2model$log2FcTUT2_1+TUT2model$log2FcTUT2_2)/2
plot(log10(TUT2model$mWT),(TUT2model$avg), cex=0.7, ylim = c(-10,10))
abline(h=0)
TUT2model$upper <- NA
i <- 1
for (i in 1:nrow(TUT2model)) {
  a <- i+50
  sd_test <- sd(TUT2model$avg[c(i:a)])
  sd_test <- sd_test*2
  TUT2model$upper[i] <- sd_test
  rm(sd_test)
  print(i)
}
TUT2model$lower <- 0-TUT2model$upper
lines(log10(TUT2model$mWT),TUT2model$upper, col="red")
lines(log10(TUT2model$mWT),TUT2model$lower, col="red")

plot(log10(final$mWT),final$log2FcDKO, cex=0.7, ylim = c(-10,10))
abline(h=0)
lines(log10(TUT2model$mWT),TUT2model$upper, col="red")
lines(log10(TUT2model$mWT),TUT2model$lower, col="red")

TUT2model <- TUT2model[complete.cases(TUT2model),]
TUT2model <- aggregate(TUT2model$upper, by=list(mWT=TUT2model$mWT),median)
TUT2model$lower <- 0-TUT2model$x
TUT2model$upper <- 0+TUT2model$x
TUT2model <- TUT2model[,c(1,4,3)]
TUT2model <- TUT2model[which(TUT2model$mWT>=1),]

model <- merge(TUT2model, final, by="mWT")

TUT2_change <- model[which(model$log2FcTUT2_1>=model$upper),]
TUT2_change2 <- model[which(model$log2FcTUT2_2>=model$upper),]
TUT2_change <- TUT2_change[TUT2_change$MIRNA %in% TUT2_change2$MIRNA,]
TUT2_change3 <- model[which(model$log2FcTUT2_1<=model$lower),]
TUT2_change4 <- model[which(model$log2FcTUT2_2<=model$lower),]
TUT2_change3 <- TUT2_change3[TUT2_change3$MIRNA %in% TUT2_change4$MIRNA,]
TUT2_change <- rbind(TUT2_change,TUT2_change3)

TUT4_change <- model[which(model$log2FcTUT4_1>=model$upper),]
TUT4_change2 <- model[which(model$log2FcTUT4_2>=model$upper),]
TUT4_change <- TUT4_change[TUT4_change$MIRNA %in% TUT4_change2$MIRNA,]
TUT4_change3 <- model[which(model$log2FcTUT4_1<=model$lower),]
TUT4_change4 <- model[which(model$log2FcTUT4_2<=model$lower),]
TUT4_change3 <- TUT4_change3[TUT4_change3$MIRNA %in% TUT4_change4$MIRNA,]
TUT4_change <- rbind(TUT4_change,TUT4_change3)

TUT7_change <- model[which(model$log2FcTUT7_1>=model$upper),]
TUT7_change2 <- model[which(model$log2FcTUT7_2>=model$upper),]
TUT7_change <- TUT7_change[TUT7_change$MIRNA %in% TUT7_change2$MIRNA,]
TUT7_change3 <- model[which(model$log2FcTUT7_1<=model$lower),]
TUT7_change4 <- model[which(model$log2FcTUT7_2<=model$lower),]
TUT7_change3 <- TUT7_change3[TUT7_change3$MIRNA %in% TUT7_change4$MIRNA,]
TUT7_change <- rbind(TUT7_change,TUT7_change3)

DKO_change <- model[which(model$log2FcDKO_1>=model$upper),]
DKO_change2 <- model[which(model$log2FcDKO_2>=model$upper),]
DKO_change <- DKO_change[DKO_change$MIRNA %in% DKO_change2$MIRNA,]
DKO_change3 <- model[which(model$log2FcDKO_1<=model$lower),]
DKO_change4 <- model[which(model$log2FcDKO_2<=model$lower),]
DKO_change3 <- DKO_change3[DKO_change3$MIRNA %in% DKO_change4$MIRNA,]
DKO_change <- rbind(DKO_change,DKO_change3)

TKO_change <- model[which(model$log2FcTKO_1>=model$upper),]
TKO_change2 <- model[which(model$log2FcTKO_2>=model$upper),]
TKO_change <- TKO_change[TKO_change$MIRNA %in% TKO_change2$MIRNA,]
TKO_change3 <- model[which(model$log2FcTKO_1<=model$lower),]
TKO_change4 <- model[which(model$log2FcTKO_2<=model$lower),]
TKO_change3 <- TKO_change3[TKO_change3$MIRNA %in% TKO_change4$MIRNA,]
TKO_change <- rbind(TKO_change,TKO_change3)

TUT2_change$change <- "TUT2KO" 
TUT4_change$change <- "TUT4KO" 
TUT7_change$change <- "TUT7KO" 
DKO_change$change <- "DKO" 
TKO_change$change <- "TKO" 
all_changes <- rbind(TUT2_change,TUT4_change,TUT7_change,DKO_change,TKO_change)
all_changes <- all_changes[,c(4,37)]
write.table(all_changes, "venn.tsv", append = F, row.names = F, sep = "\t")
unique(all_changes$MIRNA)
write.table(all, "expression3.tsv", append = F, row.names = F, sep = "\t")

write.table(model, "model.tsv", append = F, row.names = F, sep = "\t")



