library("stringr")
library(data.table)
library("RSQLite")
library(plotly)
library("processx")
library("webshot")



con = dbConnect(SQLite(), dbname="/data/splicing_tolerance/splicing_tolerance.sqlite")


#color <- colors()[sample(1:650,40,replace=F)]
color <- c("red","blue","green","purple","orange","pink","brown", "grey", "black", "grey",'chartreuse3', "deeppink", 'cornflowerblue', 
           'darkgoldenrod1', 'peachpuff3','mediumorchid2', 'turquoise3', 'wheat4', 'slategray2',"salmon4","tomato","deepskyblue3","brown4",
           "greenyellow")


########################################
######## READ JUNCTIONS FILE ###########
########################################
all_data <-readRDS(file = "~/R/splicing_project/SplicingProject/Results/nonSexSpecificJunctions.rda")
samplesID <- all_data[[1]]
col_names <- all_data[[2]]
juncIDs_person_tissue <- all_data[[3]]



#############################################
######## GET ONLY BLOOD & BRAIN #############
#############################################
personsData <- list()

for(person in names(juncIDs_person_tissue)){
 
  i<-1
  
  uniqueJuncIDs <- NULL
  allUniqueJuncIDs <- NULL
  tissuesToPlot <- NULL
  coloursTissuesToPlot <- NULL
  
  initialNumberUniqueJuncID <- NULL
  finalNumberUniqueJuncID <- NULL
  
  #only people with both brain & blood samples
  if(length(which(str_detect(tolower(names(juncIDs_person_tissue[[person]])),"blood")) == T) > 0 &&
     length(which(str_detect(tolower(names(juncIDs_person_tissue[[person]])),"brain")) == T) > 0)
  {
  
    for(tableName in names(juncIDs_person_tissue[[person]])){
      
      #only blood and brain tissues
      if(str_detect(tolower(tableName), "blood") || str_detect(tolower(tableName), "brain")){
        
        print(paste0("#################### ", person, " ########### ", tableName, " ####################"))
        
        #set the initial data
        uniqueJuncIDs <- juncIDs_person_tissue[[person]][[tableName]]
        allUniqueJuncIDs <- length(uniqueJuncIDs)
        initialNumberUniqueJuncID <- c(initialNumberUniqueJuncID,allUniqueJuncIDs)
        tissuesToPlot <- c(tissuesToPlot,tableName)
        coloursTissuesToPlot <- c(coloursTissuesToPlot, color[i])
        
        
        ###########################################################
        ######## CROSS PERSON1 DATA WITH THE OTHER PEOPLE
        ######## COMPARE RESULTS --> INTERSECTION BETWEEN BOTH
        ###########################################################
        
        for(personN in names(juncIDs_person_tissue)){
          
          if(person != personN && length(juncIDs_person_tissue[[personN]][[tableName]]) > 0){
 
            outerLeft <- setdiff(uniqueJuncIDs,juncIDs_person_tissue[[personN]][[tableName]])

            if(length(outerLeft) > 0){
              uniqueJuncIDs <- outerLeft
              #print(paste0(" --> ",length(uniqueJuncIDs)))
              allUniqueJuncIDs <- c(allUniqueJuncIDs, length(uniqueJuncIDs))
            }
          }
        }
        finalNumberUniqueJuncID <- c(finalNumberUniqueJuncID,length(uniqueJuncIDs))
        
        ######################################################
        ######## PLOT RESULTS (COLOURED LINES GRAPH) #########
        ######################################################
        if(i==1){
          pdf(paste0("~/R/splicing_project/SplicingProject/Results/Graphics/BloodBrain/",person,"-lines.pdf"))
          op <- par(cex = 1)
          plot(allUniqueJuncIDs, col = coloursTissuesToPlot[i], ylim=c(0,275000), xlim = c(0,250), type = "l", xlab="individuals", ylab="unique juncID", main = paste0("Individual ",person," - Brain&Blood tissues"))
        }
        else{
          lines(allUniqueJuncIDs, type = "l", col = coloursTissuesToPlot[i])
        }
        i <- i + 1
        
      }
    }
    
    op <- par(cex = 0.6)
    legend(50,220000, legend=tissuesToPlot, lty = 1, col=coloursTissuesToPlot)
    dev.off()
    
    ###################################
    ######## PERCENTAGE GRAPH #########
    ###################################
    print(paste0(initialNumberUniqueJuncID, " - ", finalNumberUniqueJuncID))
    
    plotData <- data.frame(tissuesToPlot, initialNumberUniqueJuncID, finalNumberUniqueJuncID)
    p <- plot_ly(plotData, x = ~tissuesToPlot, y = ~initialNumberUniqueJuncID, type = 'bar', name = 'Total JuncID') %>%
      add_trace(y = ~finalNumberUniqueJuncID, name = 'Unique JuncID') %>%
      layout(yaxis = list(title = 'Num junctions'),xaxis = list(title = 'Tissues'),barmode = 'overlay') %>% 
      add_annotations(y = ~finalNumberUniqueJuncID, text = paste0(round((plotData$finalNumberUniqueJuncID*100)/plotData$initialNumberUniqueJuncID),"%"),showarrow = T,ax=0)
    export(p, file = paste0("Results/Graphics/BloodBrain/",person,"-percentage.pdf"))
    
    
  }
}

########################################
######## GET UNIQUE JUNCID #############
########################################
# personsData <- list()
# 
# for(person in people$subj_id){
#   
#   print("")
#   print(paste0("####################################################"))
#   print(paste0("#################### ", person, " ####################"))
#   print(paste0("####################################################"))
#   print("")
#   i<-1
#   
#   for(table_name in table_names){
#     
#     if(!str_detect(tolower(table_name), "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast")){
#       initial_juncIDs <- NULL
#       unique <- NULL
#       brain_tissues <- NULL
#       last_intersection <- NULL
#   
#       
#       ######################################################
#       ######## If that person has junctions on that tissue
#       ######################################################
#       if(length(juncIDs_person_tissue[[person]][[table_name]]) > 0){
# 
#         print(paste0("#################### ", table_name, " ####################"))
#         
#         #set the initial data
#         last_intersection <- juncIDs_person_tissue[[person]][[table_name]]
#         all_intersections <- length(last_intersection)
#         initial_juncIDs <- c(initial_juncIDs, length(last_intersection))
#         brain_tissues<-c(brain_tissues,table_name)
#         brain_colours<-c(brain_colours, color[i])
#   
#         ###########################################################
#         ######## Cross person1's data with the rest of individuals
#         ###########################################################
#   
#         for(personN in people$subj_id){
#   
#           if(person != personN && length(juncIDs_person_tissue[[personN]][[table_name]]) > 0){
#   
#             #######################################################################
#             ######## COMPARE RESULTS --> INTERSECTION BETWEEN BOTH
#             #######################################################################
#             if(length(intersect(last_intersection,juncIDs_person_tissue[[personN]][[table_name]])) > 0){
#               last_intersection <- intersect(last_intersection,juncIDs_person_tissue[[personN]][[table_name]])
#               all_intersections <- c(all_intersections, length(last_intersection))
#             }
#           }
#         }
#   
#         if(i==1){
#           pdf(paste0("~/R/splicing_project/SplicingProject/Results/Graphics/BloodBrain/",person,".pdf"))
#           i<-i+1
#         }
#         
#         #plot the data
#         op <- par(cex = 1)
#         plot(all_intersections, col="red", ylim=c(0,300000), xlim = c(0,500), type = "l", xlab="individuals", ylab="unique juncID", main = paste0("Individual ",person," - ", table_name))
#         op <- par(cex = 0.6)
#         legend(300, 300000, legend=brain_tissues, lty = 1, col="red")
#   
#         DT[[table_name]] <- all_intersections
#         unique <- c(unique, length(last_intersection))
#   
#       }
#     }
#   }
#   dev.off()
#   personsData[[person]] <- list(barplot = data.frame(brain_tissues, initial_juncIDs, unique), boxplot = DT)
# }
# 
# saveRDS(personsData, "~/R/splicing_project/SplicingProject/Results/allUniquesJuncID.rda")
