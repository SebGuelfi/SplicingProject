library("stringr")
library(data.table)
library("RSQLite")

getAllJunctionsAndCounts <- function(numberOfCounts=1){
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
  tissue_names <- NULL
  for(tissue in table_names){
    #Only select non-sex-specific tissues
    if(!str_detect(tolower(tissue), "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast")){
      query<-dbSendQuery(con, paste0("PRAGMA table_info('",tissue,"');"))
      columns <- dbFetch(query)$name
      dbClearResult(query)
      #Only select those tissues with more than 70 samples (71 because later I remove the juncID column)
      if(length(columns) >= 71){
        col_names[[tissue]] <- columns[-1] #Delete juncID column
        tissue_names <- c(tissue_names,tissue)
      }
    }
  }
  ###################################################
  ######## GET ALL JUNCID PER PERSON PER TISSUE #####
  ###################################################
  juncIDs_person_tissue <- list()
  juncIDs_tissue <- list()
  
  for(person in people$subj_id){
    print(paste0("###### ", person, " ######"))
    
    for(tissue in tissue_names){
      i<-1
      for(sample in samplesID[[person]]){
        
        if(length(which(col_names[[tissue]] == sample)) == 1 ){ #every sampleID is unique, so should be equal to 1
          
          print(paste0("                 ###### ", tissue, " ######"))
          
          query<-dbSendQuery(con, paste0("SELECT juncID,`",sample,"` FROM '",tissue,"' WHERE `",sample,"` >= ",numberOfCounts," ;"))
          juncIDs <- dbFetch(query);
          dbClearResult(query)
          if(tissue %in% names(juncIDs_tissue))
          {
            name <- paste0(tissue, "_",i)
            juncIDs_tissue[[name]]<- list(juncID = juncIDs[[1]], number = juncIDs[[2]])
          }
          else{
            juncIDs_tissue[[tissue]]<- list(juncID = juncIDs[[1]], number = juncIDs[[2]])
          }
        }
      }
    }
    juncIDs_person_tissue[[person]] <- juncIDs_tissue
    juncIDs_tissue <- NULL
  }
  fileName <- paste0("~/R/splicing_project/SplicingProject/Results/allJunctions_NonSexSpecific_Counts",numberOfCounts,".rda")
  saveRDS(list(juncIDs=juncIDs_person_tissue,tissues=tissue_names), fileName)
  print("File saved!")
}

getAllJunctionsAndCountsAdjusted <- function(sourceData="~/R/splicing_project/SplicingProject/Results/allJunctions_NonSexSpecific_Counts1.rda",
                                             saveResults = T,
                                             fileName = "allJunctions_NonSexSpecific_Counts1_Normalised.rda") {
  
  ## READ JUNCTIONS FILE 
  all.JunctionsCount1 <-readRDS(file = sourceData)
  
  ## declare variables
  allJunctionsCount1 <- all.JunctionsCount1$juncIDs
  
  i<-1
  threshold <- 0
  minimumCount<-NULL
  maximumCount<-0

  
  for(person in names(allJunctionsCount1)){
    #print(paste0("###### ", person, " ######"))
    
    for(tissue in names(allJunctionsCount1[[person]])){
      
      #divide each cell by its column total
      totalCounts <- sum(allJunctionsCount1[[person]][[tissue]][[2]])
      allJunctionsCount1[[person]][[tissue]][[2]] <- allJunctionsCount1[[person]][[tissue]][[2]]/totalCounts
      
      ## GET THE MINIMUM/MAXIMUM NUMBER OF TOTAL COUNTS 
      if(i == 1 || totalCounts < minimumCount){
        minimumCount <- totalCounts
        print(paste0("New minimum: ",minimumCount))
      }
      if(i == 1 || maximumCount < totalCounts){
        maximumCount <- totalCounts
        print(paste0("New maximum: ",maximumCount))
      }
      i<-i+1
    }
  }
  if(saveResults){
      saveRDS(object = list(juncID=allJunctionsCount1,minimum=minimumCount,maximum=maximumCount),
              file = paste0("~/R/splicing_project/SplicingProject/Results/",fileName))
      print("File saved!")
  }
}

getJunctionsAndCountsAdjustedByThreshold <- function(sourceData="~/R/splicing_project/SplicingProject/Results/allJunctions_NonSexSpecific_Counts1_Normalised.rda",
                                                     saveResults = T,
                                                     isMinimum = T) {
  
  
  ###########################################################
  ######## FILTER ALL JUNCTIONS BY USING A THRESHOLD ########
  ###########################################################
  
  ## READ JUNCTIONS FILE 
  all.junctions_normalised <- readRDS(file = sourceData)
  
  ## DECLARE VARIABLES
  threshold = NULL
  fileName = NULL
  
  if(isMinimum){
    threshold <- 3/all.junctions_normalised$minimum
    fileName = "allJunctions_NonSexSpecific_Counts3_Filtered_Minimum.rda"
  }else{
    threshold <- 3/all.junctions_normalised$maximum
    fileName = "allJunctions_NonSexSpecific_Counts3_Filtered_Maximum.rda"
  }
  
  all_junctions_adjusted <- all.junctions_normalised$juncID
 
  
  ## START LOOPING
  adjustedJuctionsToSave <- list()
  tissuesToSave <- list()
  
  for(person in names(all_junctions_adjusted)){
    
    for(tissue in names(all_junctions_adjusted[[person]])){
      
      #only first tissue samples
      if(str_sub(tissue,-2,-2) != "_"){
        
        #get all junctions which have more than three adjusted counts
        indexAdjustedJunctionsCounts_GreaterThan3 <- which(all_junctions_adjusted[[person]][[tissue]]$number >= threshold)
       
        if(length(indexAdjustedJunctionsCounts_GreaterThan3) > 0){
          print(paste0("#################### ", person, " ########### ", tissue, " #################### ", length(indexAdjustedJunctionsCounts_GreaterThan3)))
          #save results
          tissuesToSave[[tissue]] <- list(juncID=all_junctions_adjusted[[person]][[tissue]]$juncID[indexAdjustedJunctionsCounts_GreaterThan3],
                                          counts=all_junctions_adjusted[[person]][[tissue]]$number[indexAdjustedJunctionsCounts_GreaterThan3])
        }
      }
    }
    if(!is.null(tissuesToSave)){
      adjustedJuctionsToSave[[person]] <- tissuesToSave
      tissuesToSave <- list()
    }
  }
  
  if(saveResults){
    saveRDS(adjustedJuctionsToSave, paste0("~/R/splicing_project/SplicingProject/Results/", fileName))
    print("file saved!")
  }
}

getUniqueAdjustedJunctions <- function(sourceData="allJunctions_NonSexSpecific_Counts3_Filtered_Minimum.rda",
                                       fileName="uniqueJunctions_NonSexSpecific_Counts3_Filtered_Minimum.rda",
                                       saveResults = T){
  
  ## READ JUNCTIONS FILE 
  all_junctions_adjusted_GT3 <-readRDS(file = paste0("~/R/splicing_project/SplicingProject/Results/",sourceData))
  
  
  adjustedUniqueJunctions <- list()
  
  for(person in names(all_junctions_adjusted_GT3)){
    
    adjustedUniqueJunctions[[person]] <- list()
    
    
    for(tissue in names(all_junctions_adjusted_GT3[[person]])){
      
      if(str_sub(tissue,-2,-2) != "_"){
        
        print(paste0("#################### ", person, " ########### ", tissue, " ####################"))
        
        #set the initial data
        uniqueJuncIDs <- all_junctions_adjusted_GT3[[person]][[tissue]]$juncID
        
        adjustedUniqueJunctions[[person]][[tissue]] <- list()
        adjustedUniqueJunctions[[person]][[tissue]][[person]] <- list(juncID = uniqueJuncIDs, percentage = 100)
        
        ###########################################################
        ######## CROSS PERSON1 DATA WITH THE OTHER PEOPLE
        ######## COMPARE RESULTS --> INTERSECTION BETWEEN BOTH
        ###########################################################
        
        for(personN in names(all_junctions_adjusted_GT3)){
          
          if(person != personN && length(all_junctions_adjusted_GT3[[personN]][[tissue]]) > 0){
            
            outerLeft <- setdiff(uniqueJuncIDs,all_junctions_adjusted_GT3[[personN]][[tissue]]$juncID)
            
            if(length(outerLeft) > 0){
              uniqueJuncIDs <- outerLeft
              
              percentage <- length(uniqueJuncIDs)*100/length(all_junctions_adjusted_GT3[[person]][[tissue]]$juncID)
              adjustedUniqueJunctions[[person]][[tissue]][[personN]] <- list(juncID = uniqueJuncIDs, percentage = percentage)
              
            }
          }
        }
      }
    }
  }
  if(saveResults){
    saveRDS(adjustedUniqueJunctions, paste0("~/R/splicing_project/SplicingProject/Results/", fileName))
    print("file saved!")
  }
}


# getRandomizedUniqueJunctionAdjusted <- function(sourceData="allJunctions_NonSexSpecific_Counts3_AdjustedMaximum.rda",
#                                                 fileName="randomized_UniqueJunctions_NonSexSpecific_Counts3_AdjustedMaximum.rda",
#                                                 saveResults = T){
#   ## READ JUNCTIONS FILE 
#   all_junctions_adjusted_GT3 <-readRDS(file = paste0("~/R/splicing_project/SplicingProject/Results/",sourceData))
#    
# 
#   dataRandomizedPeople <- list()
#   for(item in 1:30){
#      
#     #Randomize people
#     people <- sample(names(all_junctions_adjusted_GT3),replace = F)
#     adjustedUniqueJunctions <- list()
#       
#     for(person in people){
#         
#         adjustedUniqueJunctions[[person]] <- list()
#         
#         
#         for(tissue in names(all_junctions_adjusted_GT3[[person]])){
#           
#           if(str_sub(tissue,-2,-2) != "_"){
#             
#             print(paste0("#################### ", person, " ########### ", tissue, " ####################"))
#             
#             #set the initial data
#             uniqueJuncIDs <- all_junctions_adjusted_GT3[[person]][[tissue]]$juncID
#             
#             adjustedUniqueJunctions[[person]][[tissue]] <- list()
#             adjustedUniqueJunctions[[person]][[tissue]][[person]] <- list(juncID = uniqueJuncIDs, percentage = 100)
#             
#             ###########################################################
#             ######## CROSS PERSON1 DATA WITH THE OTHER PEOPLE
#             ######## COMPARE RESULTS --> INTERSECTION BETWEEN BOTH
#             ###########################################################
#             
#             for(personN in people){
#               
#               if(person != personN && length(all_junctions_adjusted_GT3[[personN]][[tissue]]) > 0){
#                 
#                 outerLeft <- setdiff(uniqueJuncIDs,all_junctions_adjusted_GT3[[personN]][[tissue]]$juncID)
#                 
#                 if(length(outerLeft) > 0){
#                   uniqueJuncIDs <- outerLeft
#                   
#                   percentage <- length(uniqueJuncIDs)*100/length(all_junctions_adjusted_GT3[[person]][[tissue]]$juncID)
#                   adjustedUniqueJunctions[[person]][[tissue]][[personN]] <- list(juncID = uniqueJuncIDs, percentage = percentage)
#                   
#                 }
#               }
#             }
#           }
#         }
#       }
#      dataRandomizedPeople[[item]] <- adjustedUniqueJunctions
#    }
#    if(saveResults){
#      saveRDS(dataRandomizedPeople, file = paste0("~/R/splicing_project/SplicingProject/Results/",fileName))
#      print("File saved!")
#    }
#  }
# 

getUniqueAdjustedJunctions(sourceData="allJunctions_NonSexSpecific_Counts3_Filtered_Maximum.rda",
                           fileName="uniqueJunctions_NonSexSpecific_Counts3_Filtered_Maximum.rda")