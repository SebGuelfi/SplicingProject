library(data.table)
library("RSQLite")
con = dbConnect(SQLite(), dbname="/data/splicing_tolerance/splicing_tolerance.sqlite")

test <-F

#Get table names
all_tables <- dbFetch(dbSendQuery(con, "select * from sqlite_master;"))
table_names <- all_tables$name

gtex_table <- table_names[56]
table_names <- table_names[-1] #Delete "Junc_coordinates" name-table
table_names <- table_names[-55] #Delete "GTEX" name-table

query <- dbSendQuery(con, paste0("SELECT DISTINCT subj_id FROM '",gtex_table,"';"))
people <- dbFetch(query)
dbClearResult(query)

person_samples <- list()
person_tissues <- list()
tissues_col_oper <- list()


for(person in people$subj_id){
  
  
  print(paste0("#### ",person, " ####"))
  query<-dbSendQuery(con, paste0("SELECT sample_recount_id FROM '",gtex_table,"' WHERE subj_id = '",person,"';"))
  samplesID <- dbFetch(query)$sample_recount_id;
  dbClearResult(query)
  
  
  for(table_name in table_names){
    
    if(grepl("FrontalCortex", table_name)){
      
      query<-dbSendQuery(con, paste0("PRAGMA table_info('",table_name,"')"))
      col_names_tissue <- dbFetch(query)$name;
      dbClearResult(query)
      x<-match(col_names_tissue,samplesID, nomatch = FALSE)
      
      if(length(which(x>0)) >= 2){
        col_names_tissue <- col_names_tissue[which(x>0)]
        
        for(col_name in col_names_tissue){
          
          query<-dbSendQuery(con, paste0("SELECT sum(`",col_name,"`) FROM '",table_name,"';"))
          total_junt_sample <- dbFetch(query)[[1]]
          dbClearResult(query)
          
          query<-dbSendQuery(con, paste0("SELECT `", col_name, "` FROM '",table_name,"';"))
          junc_sample <- dbFetch(query)[[1]]
          dbClearResult(query)
          
          junc_sample[is.na(junc_sample)] <- -1
          junc_sample <- junc_sample/total_junt_sample
          tissues_col_oper[[col_name]] <- junc_sample
          
          print(paste0(table_name," - ", col_name, " - ", total_junt_sample))
        }
        person_tissues[[table_name]]<-tissues_col_oper
      }
    }
  }
  if(length(person_tissues) > 0){
    person_samples[[person]] <- list(name=person,samples_id=samplesID,junc_tissues=person_tissues,total=total_junt_sample)
    person_tissues <- list()
    tissues_col_oper <- list()
    if(test)
      break;
  }
}
dbDisconnect(con)