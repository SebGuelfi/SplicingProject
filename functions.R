### Functions to manipulate the counts of split


### get only samples for GTEx


getGTEx <- function()
{
    GTEx.info <- read.delim("/data/recount/GTEx_SRP012682/sample_ids.tsv",header = F)
    recount.info <- read.delim("/data/recount/GTEx_SRP012682/SRP012682.tsv",header = T)
    
    library(dplyr)
    
    GTEx.info <- GTEx.info %>% 
        filter(V2=="SRP012682") %>%  
        mutate(V3_chr = as.character(V3)) %>%  
        left_join(recount.info %>% 
                      mutate(run_chr = as.character(run)) %>% 
                      select(run_chr,smts,smtsd),
                  by=c("V3_chr"="run_chr"))
    return(GTEx.info)

}


## Get the frequency of junc for a tissue
## inputs: 
## junc, the junc to calculate the frequency
## tissue, the tissue to calculate the frequency
## countsDetection, threshold of split reads detection, by default 1
getFreqTissueJunc <- function(junc, tissue,countsDetection=1)
{
    ## load the data
    #GTEx.info <- getGTEx()
    #tissues <- names(table(GTEx.info$smtsd))
    # tb <- fread(paste0("/data/recount/GTEx_SRP012682/gtex_split_read_table_annotated/",tolower(tissue),"_split_read_table_annotated.csv"))
    # get the line of the junction
    nLine <- read.delim(pipe(paste0("cut -d ',' -f 1 /data/recount/GTEx_SRP012682/gtex_full_split_read_count_table/",tissue,"/",
                           tissue,".csv | grep -nw ",junc," | cut -f1 -d:")),header = F) 
    # load the line from nLine
    tmp <- read.csv(pipe(paste0("sed '",nLine,"q;d' /data/recount/GTEx_SRP012682/gtex_full_split_read_count_table/",tissue,"/",
                                  tissue,".csv")),header = F)
    
    return((rowSums(tmp[,2:ncol(tmp)]>countsDetection,na.rm = T)/(ncol(tmp)-1))*100)
}

## example
## getFreqTissueJunc(2859,"Bladder")

## Get frequency for tissue 
getFreqTissue <- function(tissue,countsDetection=1)
{
    library(data.table)
    tb <- fread(paste0("/data/recount/GTEx_SRP012682/gtex_full_split_read_count_table/",tissue,"/",tissue,".csv"))
    tb <- tb[which(rowSums(tb[,2:ncol(tb)]>=countsDetection,na.rm = T)>0),]
    return((rowSums(tb[,2:ncol(tb)]>=countsDetection,na.rm = T)/(ncol(tb)-1))*100)
    
}
tissue <- "Brain-FrontalCortex_BA9"

## Get mean per junction
## returns a vector containing the mean per junction
getMeanJuncTissue <- function(tissue,countsDetection=1)
{
  library(data.table)
  tb <- fread(paste0("/data/recount/GTEx_SRP012682/gtex_full_split_read_count_table/",tissue,"/",tissue,".csv"))
  tb <- tb[which(rowSums(tb[,2:ncol(tb)]>=countsDetection,na.rm = T)>0),]
  #tb[is.na(tb)==T] <- 0
  message("Calculating the mean...")
  return(apply(tb[,2:ncol(tb)],1,function(x){mean(x,na.rm = T)}))
}

## get example
# head(getFreqTissue("Bladder"))
# 
# recount.info <- read.delim("/data/recount/GTEx_SRP012682/SRP012682.tsv",header = T)
# 
# table(table(unlist(lapply(strsplit(as.character(recount.info$sampid),split = "-"),function(x){return(paste0(x[1],"-",x[2]))}))))
# 
# recount.info[grep("K-562",recount.info$sampid),]



# recount.info <- read.delim("/data/recount/GTEx_SRP012682/SRP012682.tsv",header = T)
# subjID <- read.delim("/data/recount/GTEx_SRP012682/GTEx_v7_Annotations_SubjectPhenotypesDS.txt",header = T)
# head(subjID)
# 
# recount.info$subj_id <- unlist(lapply(strsplit(as.character(recount.info$sampid),"-",fixed = T),function(x){paste(x[1:2],collapse = "-")}))
# 
# 
# table(recount.info$subj_id)
# head(recount.info)


#recount.info$age <- subjID$AGE[match(recount.info$subj_id,subjID$SUBJID)]

#write.table(recount.info,"/home/sguelfi/projects/data/splicing_tolerance/GTEX_info.txt")




## Get tissue per indivual
## This function returns the tissues for the individual given by the input
get_individuals_id_per_sample <- function(samples_id,database_path="/data/splicing_tolerance/splicing_tolerance.sqlite")
{
    library(RSQLite)
    db <-  dbConnect(SQLite(), dbname=database_path)
    res <- dbSendQuery(db, paste0("SELECT subj_id,sample_recount_id FROM GTEX_info WHERE sample_recount_id IN (",
                                  paste(samples_id,collapse = ","),")"))
    return(dbFetch(res))
}

## Get tissue per indivual
## This function returns the tissues for the individual given by the input
get_tissues_individual <- function(individual_id,database_path="/data/splicing_tolerance/splicing_tolerance.sqlite")
{
    library(RSQLite)
    db <-  dbConnect(SQLite(), dbname=database_path)
    res <- dbSendQuery(db, paste0("SELECT smtsd FROM GTEX_info WHERE subj_id = '",individual_id,"'"))
    return(dbFetch(res)$smtsd)
}


## Get tissue and sample_id per indivual
## This function returns the tissues for the individual given by the input
get_recount_id_individual <- function(individual_id,database_path="/data/splicing_tolerance/splicing_tolerance.sqlite")
{
    library(RSQLite)
    db <-  dbConnect(SQLite(), dbname=database_path)
    res <- dbSendQuery(db, paste0("SELECT sample_recount_id,smtsd FROM GTEX_info WHERE subj_id = '",individual_id,"'"))
    return(dbFetch(res))
}

## example get_tissues_individual("GTEX-QMR6")
get_junc_freq_tissue <- function(junc_id,tissue,detection_threshold,database_path="/data/splicing_tolerance/splicing_tolerance.sqlite")
{
    library(RSQLite)
    db <-  dbConnect(SQLite(), dbname=database_path)
    res <- dbSendQuery(db, paste0("SELECT * FROM ",tissue,"  WHERE juncID IN (",paste(junc_id,collapse = ","),")"))
    res <- as.data.frame(dbFetch(res))
    ## we get from index 2, because the first value is the juncID
    rownames(res) <- res$juncID
    res$juncID <- NULL
    apply(res,1,function(x){
     (table(as.numeric(x)>=detection_threshold)[TRUE]/(ncol(res)))*100     
    })
    
}


## get the samples that detect the junc for specific tissue
get_samples_per_junc_detected <- function(junc_id,individual_id,tissue,detection_threshold,database_path="/data/splicing_tolerance/splicing_tolerance.sqlite")
{
    library(RSQLite)
    db <-  dbConnect(SQLite(), dbname=database_path)
    res <- dbSendQuery(db, paste0("SELECT * FROM ",tissue,"  WHERE juncID IN (",paste(junc_id,collapse = ","),")"))
    res <- as.data.frame(dbFetch(res))
    ## we get from index 2, because the first value is the juncID
    rownames(res) <- res$juncID
    res$juncID <- NULL
    res <- t(res)
    res <- unique(rownames(res)[as.numeric(res)>=detection_threshold])
    return(res[!is.na(res)])
}


## Get the samples junc_id
get_samples_per_junc_detected <- function(junc_id,individual_id, detection_threshold, database_path="/data/splicing_tolerance/splicing_tolerance.sqlite")
{
    sample_ids <- get_recount_id_individual(individual_id = individual_id )
    
    library(RSQLite)
    db <-  dbConnect(SQLite(), dbname=database_path)
    res <- apply(sample_ids,1,function(x){
        # print(paste0("SELECT * FROM '",gsub("\\)","",
        #                                    gsub("\\(","",gsub(" ","",x[2]))),"' WHERE juncID = '",junc_id,"'"))
        
        return(dbFetch(dbSendQuery(db,     paste0("SELECT * FROM '",gsub("\\)","",
        gsub("\\(","_",gsub(" ","",x[2]))),"' WHERE juncID = '",junc_id,"'"))))
    })
    
    names(res) <- sample_ids[,2]
    
    res_freq <- apply(sample_ids,1,function(x){
        #res[[x[2]]][,x[1]]
        #print(!is.na(res[[x[2]]][,x[1]]))
        if(nrow(res[[x[2]]])>0){
           if (!is.na(res[[x[2]]][,x[1]])){
            return((table(as.numeric(res[[x[2]]][,x[1]])>=detection_threshold))[TRUE])/(ncol(res[[x[2]]]))*100
           }
        }
        return(0)
    })

    names(res_freq) <- sample_ids
    as.numeric(res[[1]][],na.rm)
    res <- dbSendQuery(db,     paste0("SELECT * FROM ",gsub("\\)","",
        gsub("\\(","",gsub(" ","",sample_ids[1,2])))," WHERE juncID = '",junc_id,"'"))
    
    
    get_tissues_individual()
    res <- dbSendQuery(db, paste0("SELECT * FROM ",tissue,"  WHERE juncID IN (",paste(junc_id,collapse = ","),")"))
    res <- as.data.frame(dbFetch(res))
    ## we get from index 2, because the first value is the juncID
    rownames(res) <- res$juncID
    res$juncID <- NULL
    res <- t(res)
    res <- unique(rownames(res)[as.numeric(res)>=detection_threshold])
    return(res[!is.na(res)])
}




## this function filters the split read
#' @param GTFPath The path to the GTF file
#' @param splitReadTable a genomic ranges table with cordinates of split reads to be annotated
#' @return minCounts the minimum number of counts that are going to be used as
#' threshold (1 by default)
annotateSplitReads <- function(GTFPath,splitReadTable,geneList=NULL)
{
    library(data.table)
    library(tidyverse)
    library(refGenome)
    # change to directory where you downloaded GTF file
    #setwd("/home/sguelfi/neuroscience/WT_BRAINEAC/hipp/salmon/salmonReferences/GTF/")
    setwd(dirname(GTFPath))
    # create ensemblGenome object for storing Ensembl genomic annotation data
    ens <- ensemblGenome()
    # read GTF file into ensemblGenome object
    basedir(ens) <- dirname(GTFPath)
    message(paste(Sys.time(),"Loading the GTF..."))
    read.gtf(ens, "Homo_sapiens.GRCh38.87.gtf")
    message(paste(Sys.time(),"GTF Loaded"))
    
    jens <- getSpliceTable(ens)
    jens <- as.data.frame(jens@ev$gtf)
    
    if(!is.null(geneList))
    {
        jens <- jens %>% filter(as.character(gene_id)%in%as.character(geneList))
        message(paste(Sys.time(),length(unique(jens$gene_id)),"genes retained"))
        rm(geneList)
    }
    
    jens$strand <- as.character(jens$strand)
    
    
    jens$seqid <- gsub("chr","",as.character(jens$seqid))
    splitReadTable$chr <- gsub("chr","",as.character(splitReadTable$chr))
    
    # head(jens)
    # id seqid    lstart      lend    rstart      rend         gene_id gene_name strand   transcript_id   lexid   rexid transcript_biotype
    # 1     7 127588345 127588565 127589083 127589163 ENSG00000004059      ARF5      + ENST00000000233 1007522 1007525     protein_coding
    # 2     7 127589083 127589163 127589485 127589594 ENSG00000004059      ARF5      + ENST00000000233 1007525 1007527     protein_coding
    # 3     7 127589485 127589594 127590066 127590137 ENSG00000004059      ARF5      + ENST00000000233 1007527 1007529     protein_coding
    # 4     7 127590066 127590137 127590963 127591088 ENSG00000004059      ARF5      + ENST00000000233 1007529 1007531     protein_coding
    # 5     7 127590963 127591088 127591213 127591705 ENSG00000004059      ARF5      + ENST00000000233 1007531 1007533     protein_coding
    # 6    12   8940365   8941940   8942416   8942542 ENSG00000003056      M6PR      - ENST00000000412 1566634 1566632     protein_coding
    
    
    ## -1 the start and +1 the stop of the junctions to match the splice site positions
    ## get the coordinates
    splitReadTable$start <- splitReadTable$start -1
    splitReadTable$stop <- splitReadTable$stop +1
    
    ## convert the strand from numeric to "+,-,*" because this is what the ensembl reference uses.
    if(any(names(table(splitReadTable$strand))==1) | any(names(table(splitReadTable$strand))==2) | any(names(table(splitReadTable$strand))==0))
    {
        splitReadTable$strand[splitReadTable$strand==1] <- "+"
        splitReadTable$strand[splitReadTable$strand==2] <- "-"
        splitReadTable$strand[splitReadTable$strand==0] <- "*"
    }
    
    
    message(paste(Sys.time(), 'Getting acceptor sites annotation'))
    ## get the acceptor for the forward strand
    splitReadTable.forward <- splitReadTable
    splitReadTable.forward.tmp <- left_join(splitReadTable.forward,
                                            jens[as.character(jens$strand)=="+",c("seqid","lstart","strand","transcript_id")],
                                            by=c("chr"="seqid","stop"="lstart","strand"="strand"))
    
    splitReadTable.forward <- left_join(splitReadTable.forward,
                                        jens[as.character(jens$strand)=="+",c("seqid","rstart","strand","transcript_id")],
                                        by=c("chr"="seqid","stop"="rstart","strand"="strand"))
    
    ## get the acceptor for the reverse strand
    splitReadTable.reverse <- splitReadTable
    splitReadTable.reverse.tmp <- left_join(splitReadTable.reverse,
                                            jens[as.character(jens$strand)=="-",c("seqid","lend","strand","transcript_id")],
                                            by=c("chr"="seqid","start"="lend","strand"="strand"))
    
    splitReadTable.reverse <- left_join(splitReadTable.reverse,
                                        jens[as.character(jens$strand)=="-",c("seqid","rend","strand","transcript_id")],
                                        by=c("chr"="seqid","start"="rend","strand"="strand"))
    
    
    ## merge the acceptor for forward and reverse strand
    splitReadTable.acceptor <- rbind(splitReadTable.forward,splitReadTable.forward.tmp,splitReadTable.reverse.tmp,splitReadTable.reverse)
    rm(splitReadTable.forward,splitReadTable.forward.tmp,splitReadTable.reverse.tmp,splitReadTable.reverse)
    ## filter duplicates
    splitReadTable.acceptor <- splitReadTable.acceptor[!duplicated(splitReadTable.acceptor[,c("junID","transcript_id")]),]
    setDT(splitReadTable.acceptor)
    ## collapse the transcripts that have the same location of the split site
    splitReadTable.acceptor <- splitReadTable.acceptor %>% group_by(junID) %>% summarise(transcript_id = paste0(na.omit(transcript_id), collapse=","))
    # head(splitReadTable.acceptor)
    # junId   transcript_id
    # 1              NA
    # 2              NA
    # 3 ENST00000488147
    # 4 ENST00000488147
    # 5 ENST00000488147
    # 6 ENST00000488147
    
    ## set back to data.frame
    setDF(splitReadTable.acceptor)
    
    ## assign the transcript to the acceptor column
    splitReadTable[match(splitReadTable.acceptor$junID,splitReadTable$junID),"acceptor"] <- splitReadTable.acceptor$transcript_id
    rm(splitReadTable.acceptor)
    
    ## get the donor
    message(paste(Sys.time(), 'Getting donor sites annotation'))
    splitReadTable.forward <- splitReadTable
    
    ## get the donor for the forward strand
    splitReadTable.forward.tmp <- left_join(splitReadTable.forward,
                                            jens[as.character(jens$strand)=="+",c("seqid","lend","strand","transcript_id")],
                                            by=c("chr"="seqid","start"="lend","strand"="strand"))
    
    splitReadTable.forward <- left_join(splitReadTable.forward,
                                        jens[as.character(jens$strand)=="+",c("seqid","rend","strand","transcript_id")],
                                        by=c("chr"="seqid","start"="rend","strand"="strand"))
    
    ## get the donor for the reverse strand
    splitReadTable.reverse <- splitReadTable
    splitReadTable.reverse.tmp <- left_join(splitReadTable.reverse,
                                            jens[as.character(jens$strand)=="-",c("seqid","lstart","strand","transcript_id")],
                                            by=c("chr"="seqid","stop"="lstart","strand"="strand"))
    
    splitReadTable.reverse <- left_join(splitReadTable.reverse,
                                        jens[as.character(jens$strand)=="-",c("seqid","rstart","strand","transcript_id")],
                                        by=c("chr"="seqid","stop"="rstart","strand"="strand"))
    
    splitReadTable.donor <- rbind(splitReadTable.forward,splitReadTable.forward.tmp,splitReadTable.reverse.tmp,splitReadTable.reverse)
    rm(splitReadTable.forward,splitReadTable.forward.tmp,splitReadTable.reverse.tmp,splitReadTable.reverse)
    ## filter duplicates
    splitReadTable.donor <- splitReadTable.donor[!duplicated(splitReadTable.donor[,c("junID","transcript_id")]),]
    setDT(splitReadTable.donor)
    splitReadTable.donor <- splitReadTable.donor %>% group_by(junID) %>% summarise(transcript_id = paste0(na.omit(transcript_id), collapse=","))
    setDF(splitReadTable.donor)
    
    splitReadTable[match(splitReadTable.donor$junID,splitReadTable$junID),"donor"] <- splitReadTable.donor$transcript_id
    rm(splitReadTable.donor)
    # head(splitReadTable)
    #   junId chr start  stop strand intronMotif inAnnotation countsSamples        acceptor           donor
    # 1     1   1 14829 14930      -           2            1            12
    # 2     2   1 14829 14970      -           2            1           100
    # 3     3   1 15038 15796      -           2            1            59 ENST00000488147 ENST00000488147
    # 4     4   1 15947 16607      -           2            1            15 ENST00000488147 ENST00000488147
    # 5     5   1 16765 16858      -           2            1            27 ENST00000488147 ENST00000488147
    # 6     6   1 17055 17233      -           2            1            46 ENST00000488147 ENST00000488147
    
    message(paste(Sys.time(), 'Getting junction sites annotation'))
    
    splitReadTable.junction <-left_join(splitReadTable,
                                        jens[,c("seqid","lend","rstart","strand","transcript_id")],
                                        by=c("chr"="seqid","start"="lend", "stop"="rstart","strand"="strand"))
    
    
    ## filter duplicates
    splitReadTable.junction <- splitReadTable.junction[!duplicated(splitReadTable.junction[,c("junID","transcript_id")]),]
    setDT(splitReadTable.junction)
    splitReadTable.junction <- splitReadTable.junction %>% group_by(junID) %>% summarise(transcript_id = paste0(na.omit(transcript_id), collapse=","))
    setDF(splitReadTable.junction)
    splitReadTable[match(splitReadTable.junction$junID,splitReadTable$junID),"junction"] <- splitReadTable.junction$transcript_id
    rm(splitReadTable.junction)
    #head(splitReadTable)
    
    # junId chr start  stop strand intronMotif inAnnotation countsSamples        acceptor           donor        junction
    # 1     1   1 14829 14930      -           2            1            12
    # 2     2   1 14829 14970      -           2            1           100
    # 3     3   1 15038 15796      -           2            1            59 ENST00000488147 ENST00000488147 ENST00000488147
    # 4     4   1 15947 16607      -           2            1            15 ENST00000488147 ENST00000488147 ENST00000488147
    # 5     5   1 16765 16858      -           2            1            27 ENST00000488147 ENST00000488147 ENST00000488147
    # 6     6   1 17055 17233      -           2            1            46 ENST00000488147 ENST00000488147 ENST00000488147
    
    #save(splitReadTable,file="~/projects/R/hipp/data/expression/splitReads/jun.table.ann.tmp.rda")
    library(GenomicRanges)
    message(paste(Sys.time(),"Annotating split reads with no exact boundaries"))
    splitReadTable$precBoundDonor <- FALSE
    splitReadTable$precBoundAcceptor <- FALSE
    ## add information about correct boundaries
    splitReadTable[which(splitReadTable$junction !=""),c("precBoundDonor","precBoundAcceptor")] <- TRUE
    splitReadTable[which(splitReadTable$donor!=""),c("precBoundDonor")] <- TRUE
    splitReadTable[which(splitReadTable$acceptor!=""),c("precBoundAcceptor")] <- TRUE
    
    ## In this section I check whether the split reads lie in coding regions of the transcript, but not on splice site.
    ## forward donor
    
    splitReadTable.notPreBoun.forward <- splitReadTable[which((!splitReadTable$precBoundDonor) & splitReadTable$strand=="+"),]
    #head(splitReadTable.notPreBoun)
    
    jens.forward <- jens[jens$strand=="+",]
    
    juncti.GR <- GRanges(c(as.character(splitReadTable.notPreBoun.forward$chr)),
                         IRanges(start=c(splitReadTable.notPreBoun.forward$start),
                                 end = c(splitReadTable.notPreBoun.forward$start)),
                         strand = c(splitReadTable.notPreBoun.forward$strand),
                         juncId=c(splitReadTable.notPreBoun.forward$junID))
    
    jens.GR <- GRanges(c(as.character(jens.forward$seqid),as.character(jens.forward$seqid)),
                       IRanges(start=c(jens.forward$lstart,jens.forward$rstart),
                               end = c(jens.forward$lend,jens.forward$rend)),
                       strand = c(as.character(jens.forward$strand),as.character(jens.forward$strand)),
                       juncId=c(jens.forward$id,jens.forward$id))
    
    tabOver <- as.data.frame(findOverlaps(juncti.GR,jens.GR,ignore.strand=FALSE))
    tabOver <- cbind(tabOver,transcripts=jens[match(mcols(jens.GR[tabOver$subjectHits])[,"juncId"],jens$id),"transcript_id"])
    ## remove duplicates
    tabOver <- tabOver[!duplicated(tabOver[,c("queryHits","transcripts")]),]
    if(nrow(tabOver)>0)
    {
        setDT(tabOver)
        tabOver <- tabOver %>% group_by(queryHits) %>% mutate(transcripts = str_c(na.omit(transcripts), collapse=","))
        setDF(tabOver)
        splitReadTable[match(mcols(juncti.GR[tabOver$queryHits])[,"juncId"],splitReadTable$junID),"donor"] <- tabOver$transcripts
    }
    
    rm(splitReadTable.notPreBoun.forward,tabOver,jens.GR,juncti.GR)
    
    ## donor reverse
    splitReadTable.notPreBoun.reverse <- splitReadTable[which((!splitReadTable$precBoundDonor) &splitReadTable$strand=="-"),]
    jens.reverse <- jens[jens$strand=="-",]
    
    juncti.GR <- GRanges(c(as.character(splitReadTable.notPreBoun.reverse$chr)),
                         IRanges(start=c(splitReadTable.notPreBoun.reverse$stop),
                                 end = c(splitReadTable.notPreBoun.reverse$stop)),
                         strand = c(splitReadTable.notPreBoun.reverse$strand),
                         juncId=c(splitReadTable.notPreBoun.reverse$junID))
    
    jens.GR <- GRanges(c(as.character(jens.reverse$seqid),as.character(jens.reverse$seqid)),
                       IRanges(start=c(jens.reverse$lstart,jens.reverse$rstart),
                               end = c(jens.reverse$lend,jens.reverse$rend)),
                       strand = c(as.character(jens.reverse$strand),as.character(jens.reverse$strand)),
                       juncId=c(jens.reverse$id,jens.reverse$id))
    
    tabOver <- as.data.frame(findOverlaps(juncti.GR,jens.GR,ignore.strand=FALSE))
    tabOver <- cbind(tabOver,transcripts=jens[match(mcols(jens.GR[tabOver$subjectHits])[,"juncId"],jens$id),"transcript_id"])
    ## remove duplicates
    tabOver <- tabOver[!duplicated(tabOver[,c("queryHits","transcripts")]),]
    if(nrow(tabOver)>0)
    {
        setDT(tabOver)
        tabOver <- tabOver %>% group_by(queryHits) %>% summarise(transcripts = str_c(na.omit(transcripts), collapse=","))
        setDF(tabOver)
        splitReadTable[match(mcols(juncti.GR[tabOver$queryHits])[,"juncId"],splitReadTable$junID),"donor"] <- tabOver$transcripts
    }
    rm(splitReadTable.notPreBoun.reverse,tabOver,jens.GR,juncti.GR)
    
    ## Acceptor forward
    
    splitReadTable.notPreBoun.forward <- splitReadTable[which((!splitReadTable$precBoundAcceptor) & splitReadTable$strand=="+"),]
    
    juncti.GR <- GRanges(c(as.character(splitReadTable.notPreBoun.forward$chr)),
                         IRanges(start=c(splitReadTable.notPreBoun.forward$stop),
                                 end = c(splitReadTable.notPreBoun.forward$stop)),
                         strand = c(splitReadTable.notPreBoun.forward$strand),
                         juncId=c(splitReadTable.notPreBoun.forward$junID))
    
    jens.GR <- GRanges(c(as.character(jens.forward$seqid),as.character(jens.forward$seqid)),
                       IRanges(start=c(jens.forward$lstart,jens.forward$rstart),
                               end = c(jens.forward$lend,jens.forward$rend)),
                       strand = c(as.character(jens.forward$strand),as.character(jens.forward$strand)),
                       juncId=c(jens.forward$id,jens.forward$id))
    
    tabOver <- as.data.frame(findOverlaps(juncti.GR,jens.GR,ignore.strand=FALSE))
    tabOver <- cbind(tabOver,transcripts=jens[match(mcols(jens.GR[tabOver$subjectHits])[,"juncId"],jens$id),"transcript_id"])
    ## remove duplicates
    
    tabOver <- tabOver[!duplicated(tabOver[,c("queryHits","transcripts")]),]
    if(nrow(tabOver)>0)
    {
        setDT(tabOver)
        tabOver <- tabOver %>% group_by(queryHits) %>% summarise(transcripts = str_c(na.omit(transcripts), collapse=","))
        setDF(tabOver)
        splitReadTable[match(mcols(juncti.GR[tabOver$queryHits])[,"juncId"],splitReadTable$junID),"acceptor"] <- tabOver$transcripts
    }
    
    rm(splitReadTable.notPreBoun.forward,tabOver,jens.GR,juncti.GR,jens.forward)
    
    ##Acceptor reverse
    splitReadTable.notPreBoun.reverse <- splitReadTable[which((!splitReadTable$precBoundDonor) &splitReadTable$strand=="-"),]
    
    juncti.GR <- GRanges(c(as.character(splitReadTable.notPreBoun.reverse$chr)),
                         IRanges(start=c(splitReadTable.notPreBoun.reverse$start),
                                 end = c(splitReadTable.notPreBoun.reverse$start)),
                         strand = c(splitReadTable.notPreBoun.reverse$strand),
                         juncId=c(splitReadTable.notPreBoun.reverse$junID))
    
    jens.GR <- GRanges(c(as.character(jens.reverse$seqid),as.character(jens.reverse$seqid)),
                       IRanges(start=c(jens.reverse$lstart,jens.reverse$rstart),
                               end = c(jens.reverse$lend,jens.reverse$rend)),
                       strand = c(as.character(jens.reverse$strand),as.character(jens.reverse$strand)),
                       juncId=c(jens.reverse$id,jens.reverse$id))
    
    tabOver <- as.data.frame(findOverlaps(juncti.GR,jens.GR,ignore.strand=FALSE))
    tabOver <- cbind(tabOver,transcripts=jens[match(mcols(jens.GR[tabOver$subjectHits])[,"juncId"],jens$id),"transcript_id"])
    
    ## remove duplicates
    tabOver <- tabOver[!duplicated(tabOver[,c("queryHits","transcripts")]),]
    if(nrow(tabOver)>0)
    {
        setDT(tabOver)
        tabOver <- tabOver %>% group_by(queryHits) %>% summarise(transcripts = str_c(na.omit(transcripts), collapse=","))
        setDF(tabOver)
        splitReadTable[match(mcols(juncti.GR[tabOver$queryHits])[,"juncId"],splitReadTable$junID),"acceptor"] <- tabOver$transcripts
    }
    
    rm(tabOver,jens.GR,juncti.GR,jens.reverse)
    
    ## get original coordinates
    splitReadTable$start <- splitReadTable$start +1
    splitReadTable$stop <- splitReadTable$stop -1
    
    return(splitReadTable)
}
