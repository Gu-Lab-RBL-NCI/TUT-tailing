setwd("~/Desktop/P0 TUT specificty/Qiagen_Dec2020/synthetic-RNA/")
require(data.table)
iso_summary <- fread(file = "11miRNA_4methods.isomir.sequence_info.tsv", header = T)
iso_summary$SAMPLE <- gsub(".fastq_ready","",iso_summary$SAMPLE)

unique(iso_summary$SAMPLE)
iso_summary <- iso_summary[which(iso_summary$SAMPLE=="12miRNA-Qia_S2_L001_R1_001"),]
iso_summary <- unique(iso_summary)

unmodified <- iso_summary[which(iso_summary$LEN_TRIM==0 & iso_summary$LEN_TAIL==0),]
max_mirna <- aggregate(unmodified$READS, by=list(MIRNA=unmodified$MIRNA,SEQUENCE=unmodified$SEQUENCE),sum)
max_mirna2 <- aggregate(max_mirna$x, by=list(MIRNA=max_mirna$MIRNA),max)
max_mirna2 <- max_mirna2[which(max_mirna2$x>=1000),]
max_mirna <- merge(max_mirna, max_mirna2, by=c("MIRNA","x"))
max_mirna$LEN_READ <- nchar(max_mirna$SEQUENCE)

all_unmod <- aggregate(unmodified$READS, by=list(MIRNA=unmodified$MIRNA,SEQUENCE=unmodified$SEQUENCE),sum)
all_unmod <- all_unmod[all_unmod$MIRNA %in% max_mirna$MIRNA, ]
all_unmod$LEN_READ <- nchar(all_unmod$SEQUENCE)

all_ratained <- all_unmod
all_ratained$dist <- NA
all_ratained <- all_ratained[c(0),]

#string distance calculation
require('stringdist')

i <- 1
for(i in 1:nrow(max_mirna)){
  miRNA_test <- max_mirna$MIRNA[i]
  len_test <- max_mirna$LEN_READ[i]
  seq_test <- max_mirna$SEQUENCE[i]
  reads_test <- all_unmod[which(all_unmod$MIRNA==miRNA_test & all_unmod$LEN_READ==len_test),]
  reads_test$dist <- NA
  a <- 1
  for (a in 1:nrow(reads_test)) {
    dist_test <- stringdist(reads_test$SEQUENCE[a], seq_test, method="dl")
    reads_test$dist[a] <- dist_test
    rm(dist_test)
    print(a)
  }
  reads_test <- reads_test[which(reads_test$dist<=1),]
  all_ratained <- rbind(all_ratained,reads_test)
  rm(miRNA_test,len_test,seq_test,reads_test)
  print(i/nrow(max_mirna)*100)
}

max_mirna$second <- substr(max_mirna$SEQUENCE, 2, 2)
all_ratained$second <- substr(all_ratained$SEQUENCE, 2, 2)

#max_mirna$second <- substr(max_mirna$SEQUENCE, 3, 3) #3rd pos
#all_ratained$second <- substr(all_ratained$SEQUENCE, 3, 3)

A_ref <- max_mirna[which(max_mirna$second=="A"),]
A_test <- all_ratained[which(all_ratained$second!="A"),]
A_test <- A_test[A_test$MIRNA %in% A_ref$MIRNA,]
mut_A <- sum(A_test$x)/(sum(A_test$x)+sum(A_ref$x))*100
mut_AtoU <- sum(A_test[which(A_test$second=="T"),3])/(sum(A_test$x)+sum(A_ref$x))*100
mut_AtoC <- sum(A_test[which(A_test$second=="C"),3])/(sum(A_test$x)+sum(A_ref$x))*100
mut_AtoG <- sum(A_test[which(A_test$second=="G"),3])/(sum(A_test$x)+sum(A_ref$x))*100
(sum(A_test$x)+sum(A_ref$x))

U_ref <- max_mirna[which(max_mirna$second=="T"),]
U_test <- all_ratained[which(all_ratained$second!="T"),]
U_test <- U_test[U_test$MIRNA %in% U_ref$MIRNA,]
mut_U <- sum(U_test$x)/(sum(U_test$x)+sum(U_ref$x))*100
mut_UtoA <- sum(U_test[which(U_test$second=="A"),3])/(sum(U_test$x)+sum(U_ref$x))*100
mut_UtoC <- sum(U_test[which(U_test$second=="C"),3])/(sum(U_test$x)+sum(U_ref$x))*100
mut_UtoG <- sum(U_test[which(U_test$second=="G"),3])/(sum(U_test$x)+sum(U_ref$x))*100
(sum(U_test$x)+sum(U_ref$x))

C_ref <- max_mirna[which(max_mirna$second=="C"),]
C_test <- all_ratained[which(all_ratained$second!="C"),]
C_test <- C_test[C_test$MIRNA %in% C_ref$MIRNA,]
mut_C <- sum(C_test$x)/(sum(C_test$x)+sum(C_ref$x))*100
mut_CtoA <- sum(C_test[which(C_test$second=="A"),3])/(sum(C_test$x)+sum(C_ref$x))*100
mut_CtoU <- sum(C_test[which(C_test$second=="T"),3])/(sum(C_test$x)+sum(C_ref$x))*100
mut_CtoG <- sum(C_test[which(C_test$second=="G"),3])/(sum(C_test$x)+sum(C_ref$x))*100
(sum(C_test$x)+sum(C_ref$x))

G_ref <- max_mirna[which(max_mirna$second=="G"),]
G_test <- all_ratained[which(all_ratained$second!="G"),]
G_test <- G_test[G_test$MIRNA %in% G_ref$MIRNA,]
mut_G <- sum(G_test$x)/(sum(G_test$x)+sum(G_ref$x))*100
mut_GtoA <- sum(G_test[which(G_test$second=="A"),3])/(sum(G_test$x)+sum(G_ref$x))*100
mut_GtoU <- sum(G_test[which(G_test$second=="T"),3])/(sum(G_test$x)+sum(G_ref$x))*100
mut_GtoC <- sum(G_test[which(G_test$second=="C"),3])/(sum(G_test$x)+sum(G_ref$x))*100
(sum(G_test$x)+sum(G_ref$x))
sum(all_ratained$x)

