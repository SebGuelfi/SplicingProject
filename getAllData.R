library("stringr")
library(data.table)
library("RSQLite")
con = dbConnect(SQLite(), dbname="/data/splicing_tolerance/splicing_tolerance.sqlite")


#######################################
######## GET ALL TABLES ###############
#######################################
query <- dbSendQuery(con, "select * from sqlite_master;")
all_tables <- dbFetch(query)
dbClearResult(query)
table_names <- all_tables$name
gtex_table <- table_names[56]
table_names <- table_names[-1] #Delete "Junc_coordinates" name-table
table_names <- table_names[-55] #Delete "GTEX" name-table
########################################
########## GET GTEX DATA ###############
########################################
# Get all 'USE ME' people from GTEX project 
query <- dbSendQuery(con, paste0("SELECT DISTINCT subj_id FROM '",gtex_table,"' WHERE smafrze = 'USE ME';"))
people <- dbFetch(query)
dbClearResult(query)
# Get all samplesID from 'USE ME' people from GTEX project
samplesID <- list()
for(person in people$subj_id){
  query<-dbSendQuery(con, paste0("SELECT sample_recount_id FROM '",gtex_table,"' WHERE subj_id = '",person,"';"))
  samplesID[[person]] <- dbFetch(query)$sample_recount_id;
  dbClearResult(query)
}
###############################################
######## GET COLNAMES FROM EVERY TISSUE #######
###############################################
col_names <- list()
for(table_name in table_names){
  #Only select non-sex-specific tissues
  if(!str_detect(tolower(table_name), "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast")){
    query<-dbSendQuery(con, paste0("PRAGMA table_info('",table_name,"');"))
    columns <- dbFetch(query)$name
    dbClearResult(query)
    #Only select those tissues with more than 75 samples
    if(length(columns) >= 75)
      col_names[[table_name]] <- columns[-1] #Delete juncID column
  }
}
###################################################
######## GET ALL JUNCID PER PERSON PER TISSUE #####
###################################################
juncIDs_person_tissue <- list()
juncIDs_tissue <- list()

for(person in people$subj_id){
  print(paste0("###### ", person, " ######"))
  
  for(table_name in table_names){
    i<-1
    for(sample in samplesID[[person]]){
      
      if(length(which(col_names[[table_name]] == sample)) == 1){ #every sampleID is unique, so should be equal to 1
        
        print(paste0("                 ###### ", table_name, " ######"))
        query<-dbSendQuery(con, paste0("SELECT juncID FROM '",table_name,"' WHERE `",sample,"` >= 1 ;"))
        juncIDs <- dbFetch(query)[[1]];
        dbClearResult(query)
        if(table_name %in% names(juncIDs_tissue))
        {
          name <- paste0(table_name, "_",i)
          juncIDs_tissue[[name]]<- juncIDs
          
        }
        else
          juncIDs_tissue[[table_name]]<- juncIDs
        
      }
    }
  }
  juncIDs_person_tissue[[person]] <- juncIDs_tissue
  juncIDs_tissue <- NULL
}
all_data <- list(samplesID_person=samplesID,col_names=col_names,juncIDs_person_tissue=juncIDs_person_tissue)
saveRDS(all_data, "~/R/splicing_project/SplicingProject/Results/nonSexSpecificJunctions.rda")
print("File saved!")

