# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome")
# library(BSgenome)
# GRCh38.dna <- readDNAStringSet("/data/references/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa", "fasta") ## load the file.
# save(GRCh38.dna,file="/home/sguelfi/projects/R/hipp/data/general/GRCh38.dna.rda")

library(tidyverse)
## load the fasta sequences
load("/home/sguelfi/projects/R/hipp/data/general/GRCh38.dna.rda")
## load the split reads
load("/home/sguelfi/projects/R/splicing_tolerance/ensembl_annotated_split_reads.rda")

## remove weird chromosome scaffolds and change the name of the mitochondrial one
annotated_split_reads <- annotated_split_reads %>% filter(as.character(chr) %in% c(as.character(1:22),"X","Y","M"))
annotated_split_reads$chr[annotated_split_reads$chr=="M"] <- "MT"


library(Biostrings)
## function that get the sequences from the coordinates to use as input for maxent score
getss3ss5 <- function(GRCh38.dna,coordinate){
  
  selChr <- GRCh38.dna[grep(paste0("GRCh38:",as.character(coordinate$chr),":"),names(GRCh38.dna))]
  if (as.character(coordinate$strand) =="-")
  {
    ss3 <- as.data.frame(reverseComplement(subseq(selChr, start=as.numeric(coordinate$start)-3, end=as.numeric(coordinate$start)+19)))[,1]
    ss5 <- as.data.frame(reverseComplement(subseq(selChr, start=as.numeric(coordinate$stop)-5, end=as.numeric(coordinate$stop)+3)))[,1]
  }else{
    ss5 <- as.data.frame(subseq(selChr, start=as.numeric(coordinate$start)-3, end=as.numeric(coordinate$start)+5))[,1]
    ss3 <- as.data.frame(subseq(selChr, start=as.numeric(coordinate$stop)-19, end=as.numeric(coordinate$stop)+3))[,1]
  }
  return(c(junID=coordinate$junID, ss5=ss5,ss3=ss3))
}

## for testing 
#annotated_split_reads <-  annotated_split_reads[sample(1:nrow(annotated_split_reads), 2000, replace=F),]

#parallelise the operation
library(doParallel)
library(foreach)

cl <- makeCluster(10)
clusterExport(cl, c("getss3ss5","reverseComplement","subseq","as.data.frame"))
registerDoParallel(cl)

ss5ss3 <- NULL
ss5ss3 <- foreach(i=1:nrow(annotated_split_reads),.combine=rbind,.verbose=F)%dopar%getss3ss5(GRCh38.dna = GRCh38.dna,annotated_split_reads[i,])
stopCluster(cl)
rm(cl)

# Save the output
save(ss5ss3,file="/home/sguelfi/projects/R/tmp/ss5ss3_GTEx.rda")
load("/home/sguelfi/projects/R/tmp/ss5ss3_GTEx.rda")

tmp.file <- tempfile()
## get the maxentscan for the 5' splice site
write.table(gsub("N","A",ss5ss3[,2]),file=tmp.file,row.names=F,col.names=F,quote=F)
setwd("/tools/maxEntScan/fordownload/")
ss5score <- read.delim(pipe(paste0("perl /tools/maxEntScan/fordownload/score5.pl ", tmp.file)),header = F)

## get the maxentscan for the 3' splice site
write.table(gsub("N","A",ss5ss3[,3]),file=tmp.file,row.names=F,col.names=F,quote=F)
ss3score <- read.delim(pipe(paste0("perl /tools/maxEntScan/fordownload/score3.pl ", tmp.file)),header = F)


ss5score[,1] <- ss5ss3[,2]
ss3score[,1] <- ss5ss3[,3]

ss5score[,3] <- ss5ss3[,1]
ss3score[,3] <- ss5ss3[,1]


## Assign NA score to sequences that contain Ns
ss5score[grep("N",as.character(ss5score[,1])),2] <- "NA"
colnames(ss5score) <- c("ss5sequence","ss5score","junID" )

ss3score[grep("N",as.character(ss3score[,1])),2] <- "NA"
colnames(ss3score) <- c("ss3sequence","ss3score","junID" )

head(ss5score)
#   ss5sequence ss5score junID
# 1   CCAGTAAGT     9.09   538
# 2   AAGGTGAGC      9.6  2641
# 3   CAGGTGGGA     6.71  4373
# 4   AAGGTAGGT    10.29  6011
# 5   CAGGTGGGA     6.71  6457
# 6   CCTGTATGG     1.54  7128


head(ss3score)
#               ss3sequence ss3score junID
# 1 AGGCTCCTGTCTCCCCCCAGGTG     11.9   538
# 2 CTGTGGCTTTCCCGTTGCAGTGA     9.17  2641
# 3 TCCCCACCTCCCGGCTCCAGTCC     6.51  4373
# 4 CCCCACACTTTGTGTTTCAGACC     8.69  6011
# 5 GTCCCTGCCTAATCTTGCAGGTC    10.64  6457
# 6 CCCCCATGTCGCCTCTGTAGGTA       10  7128

tmp.table <- cbind(tmp.table,ss5score,ss3score)
head(tmp.table)

save(tmp.table,file="/home/sguelfi/projects/R/SplicingProject/data/.rda")
load(file="~/projects/R/hipp/data/expression/splitReads/jun.table.ann.GTEx.maxEntScan.rda")
tmp.table.maxENTSCAN <- tmp.table

load(file="~/projects/R/hipp/data/expression/splitReads/jun.table.ann.GTEx.NABEC.rda")

identical(tmp.table.maxENTSCAN$junId,tmp.table$junId)

tmp.table.maxENTSCAN$NABEC <- tmp.table$NABEC

tmp.table <- tmp.table.maxENTSCAN[which(tmp.table.maxENTSCAN$`Brain-Hippocampus`>0 | tmp.table.maxENTSCAN$NABEC>0),]



#### Perform analysis on the sequence and the scores


tmp <- ss5score %>% inner_join(ss3score,by = "junID")

tmp_2 <- tmp %>% inner_join(annotated_split_reads,by = "junID")


tmp_2_annotated <- tmp_2 %>% filter(junction != "")
head(tmp_2_annotated)

mean(as.numeric(tmp_2_annotated$ss5score))

tmp_2_unannotated <- tmp_2 %>% filter(junction == "", donor == "",  acceptor == "")


tmp_2_unannotated <- tmp_2 %>% filter(junction == "", donor != "",  acceptor == "")




strsplit(x, "")

  library(data.table)

plotSeqLogo <- function(sequences){
  sequences <- do.call(rbind, str_split(sequences,  ""))
  df <- t(as.data.frame(c(A=c(0),T=c(0),C=c(0),G=c(0))))
  for(i in 1:ncol(sequences)){
      df <- rbind.fill(as.data.frame(t(as.matrix(table(sequences[,i])))),as.data.frame(df))
  }

  df <- df[1:(nrow(df)-1),]

  proportion <- function(x){
  rs <- sum(x,na.rm = T);
  return(x / rs);
  }

  #create position weight matrix
  pwm <- apply(df, 1, proportion)
  pwm[is.na(pwm)] <- 0

  library(seqLogo)
  pwm <- makePWM(pwm)
  seqLogo(pwm)
}


tmp_2_unannotated <- 


tmp_2_unannotated <- tmp_2 %>% filter(junction == "", donor != "",  acceptor == "")

## both ends are annotated
plotSeqLogo(tmp_2 %>% filter(junction != "") %>% select (ss5sequence) %>% unlist(c()))
## donor is annotatated
plotSeqLogo(tmp_2 %>% filter(junction == "", donor != "",  acceptor == "") %>% select (ss5sequence) %>% unlist(c()))

## both ends are annotated
plotSeqLogo(tmp_2 %>% filter(junction == "", donor == "",  acceptor == "") %>% select (ss5sequence) %>% unlist(c()))





            
            

library(devtools)

install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/DECIPHER_2.10.2.tar.gz", repos=NULL)

https://omarwagih.github.io/ggseqlogo/

library(DECIPHER)

data("TrainingSet_16S")
# import test sequences
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)

# remove any gaps in the sequences
dna <- RemoveGaps(dna)

# classify the test sequences
ids <- IdTaxa(dna, TrainingSet_16S, strand="top")
ids

# view the results
plot(ids, TrainingSet_16S)



