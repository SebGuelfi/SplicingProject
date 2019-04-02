library("stringr")
library(data.table)
library("RSQLite")
library(plotly)


plotMeanAcrossTissuesBrainBlood <- function(data="~/R/splicing_project/SplicingProject/Results/uniqueJunctions_NonSexSpecific_Counts3_AdjustedMinimum.rda",
                                            plotsFolderName = "~/R/splicing_project/SplicingProject/Results/Graphics/AllTissues/GreaterThan3/Adjusted",
                                            isMinimum = T,
                                            savePlot = T,
                                            linesPot = T,
                                            boxPlot = T,
                                            barPlot = T,
                                            zoomed = F,
                                            meanMinMax = T) {
  
  ## ########################################
  ## ###### PLOT RESULTS ALL TISSUES ########
  ## ########################################
  
  tissueNames <- readRDS(file ="~/R/splicing_project/SplicingProject/Results/tissuesUsed.rda")
  print(tissueNames)
  
  ## READ JUNCTIONS FILE 
  uniqueAdjustedJunctions_GT3 <- readRDS(file = data)
  
  ## DECLARE VARIABLES
  if(isMinimum){
    nameFile <- "Minimum"
    plotsFolderName <- paste0(plotsFolderName,"/Minimum")
  }else{
    nameFile <- "Maximum"
    plotsFolderName <- paste0(plotsFolderName,"/Maximum")
  }
  if(zoomed)
    nameFile <- paste0(nameFile,"-ZOOMED")
  nameFile
  plotsFolderName
  
  ## ###### GET MEAN ACROSS INDIVIDUALS ############
  tissuePoints <- list()
  points <- list()
  
  for(tissue in tissueNames){

      points <- list()
      for(person in names(uniqueAdjustedJunctions_GT3)){
        print(person)
        
        if(length(uniqueAdjustedJunctions_GT3[[person]][[tissue]]) > 0){
          
          for(item in 1:length(uniqueAdjustedJunctions_GT3[[person]][[tissue]])){
            
            if(item<=length(points)){
              points[[item]] <- c(points[[item]],uniqueAdjustedJunctions_GT3[[person]][[tissue]][[item]]$percentage)
            }else
              points[[item]] <- uniqueAdjustedJunctions_GT3[[person]][[tissue]][[item]]$percentage
          }
        }
      }
      if(length(points)>0)
        tissuePoints[[tissue]] <- points
  }
  max<-list()
  min<-list()
  mean<-list()
  for(tissue in names(tissuePoints)){
    for(item in 1:length(tissuePoints[[tissue]])){
      if(tissue %in% names(mean)){
        mean[[tissue]] <- c(mean[[tissue]],mean(tissuePoints[[tissue]][[item]]))
        max[[tissue]] <- c(max[[tissue]],max(tissuePoints[[tissue]][[item]]))
        min[[tissue]] <- c(min[[tissue]],min(tissuePoints[[tissue]][[item]]))
      }else{
        mean[[tissue]] <- mean(tissuePoints[[tissue]][[item]])
        max[[tissue]] <- max(tissuePoints[[tissue]][[item]])
        min[[tissue]] <- min(tissuePoints[[tissue]][[item]])
      }
    }
  }
  
  
  ## ############ PLOT RESULTS  ##############
  if(linesPot){
    i<-1
    allColours <- NULL
    color <- c("red","blue","green","purple","orange","pink","brown", "grey", "black", "grey",'chartreuse3', "deeppink", 'cornflowerblue', 
               'darkgoldenrod1', 'peachpuff3','mediumorchid2', 'turquoise3', 'wheat4', 'slategray2',"salmon4","tomato","deepskyblue3","brown4",
               "greenyellow")

    
    if(savePlot){
      dir.create(file.path(plotsFolderName), showWarnings = T)
      pdf(paste0(plotsFolderName,"/MeanAcrossIndividuals-AllTissues-Lines-Adjusted",nameFile,".pdf"))
    }
    
    for(tissue in names(mean)){
      print(tissue)
      if(str_detect(tolower(tissue), "brain")){
        color <- "red"
      }
      else if(str_detect(tolower(tissue), "blood")){
        color <- "blue"
      }
      else
        color <- "black"
      allColours <- c(allColours,color)

      if(i==1){
        op <- par(cex = 1) #ylim=c(0,275000), xlim = c(0,250),
        if(zoomed){
          plot(mean[tissue][[1]], col = color, ylim=c(0,25),xlim = c(0,25), type = "l", xlab="nº individuals",ylab="mean of unique juncID (in %)", main = paste0("Mean of unique junctions across individuals - All Tissues"))
        }
        else{
          plot(mean[tissue][[1]], col = color, type = "l",xlim = c(0,450),xlab="nº individuals",ylab="mean of unique juncID (in %)", main = paste0("Mean of unique junctions across individuals - All Tissues"))
        }
        abline(h=0, col="red", lwd=2, lty=2)
        abline(v=0, col="red", lwd=2, lty=2)
      }else{
        lines(mean[tissue][[1]], type = "l", col = color)
      }
      i<-i+1
    }
    op <- par(cex = 0.5)
    if(zoomed){
      legend(15,25, legend=names(mean), lty = 1, col=allColours)
    }else{
      legend(250,100, legend=names(mean), lty = 1, col=allColours)
    }
    if(savePlot)
      dev.off()
  }
  if(boxPlot){
    dir.create(file.path(paste0(plotsFolderName,"/boxPlot")), showWarnings = T)
    for(tissue in names(tissuePoints)){
      if(savePlot)
        pdf(paste0(plotsFolderName,"/boxPlot/MeanAcrossIndividuals-", tissue,"-", nameFile,"-BoxPlot.pdf"))
      boxplot(tissuePoints[[tissue]][2:10], ylim=c(0,100),col = "red", main = tissue, xlab="nº individuals",ylab="unique juncID")
      if(savePlot)
        dev.off()
    }
  }
  if(barPlot){
    percentage <- NULL
    minLength<- 1000
    for(tissue in names(mean)){
      print(minLength)
      print(tissue)
      if(length(mean[[tissue]]) < minLength)
        minLength <- length(mean[[tissue]])
    }
    for(tissue in names(mean)){
      percentage<-c(percentage,mean[[tissue]][[minLength]])
    }

    
    plotData <- data.frame(names(mean), percentage)
    
    plotData<-plotData[order(plotData$percentage, decreasing = T),] 
    
    p <- plot_ly(plotData, x = ~names(mean), y = ~percentage, type = 'bar', name = '% Unique JuncID') %>%
     layout(yaxis = list(title = 'mean % Unique JuncID'),xaxis = list(title = 'Tissues')) %>%
     add_annotations(y = ~percentage, text = paste0(round(x = percentage,digits = 2),"%"),showarrow = T,ax=0)
    p
    graphName <- paste0("Mean-UniquenessPercentage-Counts3-",nameFile,".png")
    graphName
    #if(savePlot)
    #  export(p, file = paste0(plotsFolderName,"/Mean-UniquenessPercentage-Counts3-",nameFile,".pdf"))
  }
  
  if(meanMinMax){
    for(tissue in names(mean)){
      
      if(savePlot){
        dir.create(file.path(paste0(plotsFolderName,"/MeanMinMax")), showWarnings = T)
        pdf(paste0(plotsFolderName,"/MeanMinMax","/MeanMinMaxAcrossIndividuals-Lines-Adjusted",tissue,"-",nameFile,".pdf"))
      }
      total <- length(mean[tissue][[1]])
      plot(mean[tissue][[1]], col = "red", type = "l", xlab="nº individuals",ylab="unique % of juncID", main = paste0("Mean/Min/Max % of unique junctions across individuals\n", tissue))
      text(x = (total-15), y = 45, paste0(round(x = mean[tissue][[1]][total],digits = 2),"%"), col = "red")
      
      lines(max[tissue][[1]], type = "l", col = "black")
      text(x = (total-15), y = 50, paste0(round(x = max[tissue][[1]][total],digits = 2),"%"), col = "black")
      
      lines(min[tissue][[1]], type = "l", col = "blue")
      text(x = (total-15), y = 40, paste0(round(x = min[tissue][[1]][total],digits = 2),"%"), col = "blue")
      
      legend((total-60),70, legend=c("max","mean","min"), lty = 1, col=c("black","red","blue"))
      if(savePlot)
        dev.off()
    }
  }
}



plotHistograms <- function(allJunctionsData="~/R/splicing_project/SplicingProject/Results/allJunctions_NonSexSpecific_Counts3_AdjustedMinimum.rda",
                           uniqueJunctionsData="~/R/splicing_project/SplicingProject/Results/uniqueJunctions_NonSexSpecific_Counts3_AdjustedMinimum.rda",
                           folderName = "~/R/splicing_project/SplicingProject/Results/Graphics/AllTissues/GreaterThan3/Adjusted/",
                           generatePlot = T,
                           isMinimum = T,
                           savePlot = T,
                           saveResults = T){

  #############################################
  ########### PLOT HISTOGRAM LINE #############
  #############################################
  
  ######## Read file
  allJuncIDs <-readRDS(file = allJunctionsData)
  allUniqueJuncIDs <-readRDS(file = uniqueJunctionsData)
  
 
  ## DECLARE VARIABLES
  if(isMinimum){
    folderName <- paste0(folderName,"Minimum/Histogram")
  }else
    folderName <- paste0(folderName,"Maximum/Histogram")
  folderName
  ######## Plot histogram
  all.nonunique <- list()
  all.unique <- list()
   
  for(person in names(allUniqueJuncIDs)){
    
   uniqueJuncIDs <- NULL
   all.JuncIDs <- NULL
   
  
   for(tissue in names(allUniqueJuncIDs[[person]])){
     
     #only first samples 
     if(str_sub(tissue,-2,-2) != "_"){
       
       print(paste0("#################### ", person, " ########### ", tissue, " ####################"))
       
       all.JuncIDs <- allUniqueJuncIDs[[person]][[tissue]][[1]]$juncID
       total <- length(allUniqueJuncIDs[[person]][[tissue]])
       uniqueJuncIDs <- allUniqueJuncIDs[[person]][[tissue]][[total]]$juncID
       
       indexes <- match(uniqueJuncIDs, all.JuncIDs)
       
       #Get the counts from the non-unique juncID
       non.unique <- allJuncIDs[[person]][[tissue]]$counts[-indexes]
       if(tissue %in% names(all.nonunique)){
        all.nonunique[[tissue]] <- c(all.nonunique[[tissue]], non.unique)
       }
       else
        all.nonunique[[tissue]] <- non.unique
       
       
       #Get the counts from the unique juncID
       unique <- allJuncIDs[[person]][[tissue]]$counts[indexes]
       if(tissue %in% names(all.unique)){
         all.unique[[tissue]] <- c(all.unique[[tissue]], unique)
       }
       else
         all.unique[[tissue]] <- unique
       
       #Print some results
       print(length(all.JuncIDs))
       print(length(uniqueJuncIDs)) 
       print(length(non.unique))
       print(length(unique))

     }
   }
  }
  if(generatePlot){
    for(tissue in names(all.nonunique)){
      if(savePlot){
        dir.create(file.path(folderName), showWarnings = FALSE)
        if(isMinimum)
          pdf(paste0(folderName,"/Histogram-ConcatenatedData-Minimum-",tissue,".pdf"))
        else
          pdf(paste0(folderName,"/Histogram-ConcatenatedData-Maximum-",tissue,".pdf"))
      }
      
      plot(density(log(all.nonunique[[tissue]])), col = "black", type = "l",  main = paste0("Concatenated Counts Across Individuals - ",tissue), xlab = "log(counts) concatenated across all individuals")
      lines(density(log(all.unique[[tissue]])), type = "l", col = "red")
      
      t<-ks.test(all.nonunique[[tissue]],all.unique[[tissue]])
      text(x = -8, y = 0.6, paste0("p-value = ",t$p.value))
      legend(x = -8, y = 0.5, legend = c("non.unique junctions", "unique junctions"), lty = 1, col=c("black", "red"))
      
      if(savePlot)
        dev.off()
    }
  }
}


getRandomizedDataPerPersonAndTissue <- function(data="~/R/splicing_project/SplicingProject/Results/uniqueJunctions_Randomized_GT3.rda",
                                                    saveResults = T){
  
 
  
  ##################################################################################################################################
  ##################################################################################################################################
  
  allRandomizedData <-readRDS(file = data)
  randomizedIndividualList <- allRandomizedData$randomizedIndividualList
  
  firstRandomList <- randomizedIndividualList[[1]]
  randomizedIndividualList <- randomizedIndividualList[-1]
  
  points <- list()

  for(person in names(firstRandomList)){
    myperson <- firstRandomList[[person]]
    points[[person]] <- list()
    print(person)
    for(tissue in names(myperson)){
      i<-1
      points[[person]][[tissue]] <- list()
    
        
      dataToMean <- myperson[[tissue]]
      for(item in 1:length(dataToMean)){
        if(item<=length(points[[person]][[tissue]]))
          points[[person]][[tissue]][[item]] <- c(points[[person]][[tissue]][[item]],dataToMean[[item]]) 
        else
          points[[person]][[tissue]][[item]] <- dataToMean[[item]]
      }
      
      for(randomList in randomizedIndividualList){
        
        personN <- randomList[[person]]
        dataToMean <- personN[[tissue]]
       
        for(item in 1:length(dataToMean)){
          if(item<=length(points[[person]][[tissue]]))
            points[[person]][[tissue]][[item]] <- c(points[[person]][[tissue]][[item]],dataToMean[[item]]) 
          else
            points[[person]][[tissue]][[item]] <- dataToMean[[item]]
        }
      }
    }
  }
  #data <- list(randomizedIndividualList=dataRandomizedPeople)
  if(saveResults)
    saveRDS(points, "~/R/splicing_project/SplicingProject/Results/randomizedJunctionsJoinedPerPersonAndTissue.rda")
}

plotMeanRandomizedUniqueJunctionAllTissues <- function(data="~/R/splicing_project/SplicingProject/Results/randomizedJunctionsJoinedPerPersonAndTissue.rda",
                                                       saveResults = T){
  
  
  randomizedData <-readRDS(file = data)
  
  #######################################################
  ########### GET THE MEAN PER TISSUE ###################
  #######################################################
  finalMeanData<-list()
  
  for(person in names(randomizedData)){
    finalMeanData[[person]]<-list()
    print(person)

    for(tissue in names(randomizedData[[person]])){
      
      meanDataPerTissue <-NULL
      for(item in 1:length(randomizedData[[person]][[tissue]])){
        meanDataPerTissue <- c(meanDataPerTissue, mean(randomizedData[[person]][[tissue]][[item]]))
      }
      finalMeanData[[person]][[tissue]] <- meanDataPerTissue
      
    }

  }
  # uniqueJunction <- uniqueJunctions$uniqueJuncID
  if(saveResults)
    saveRDS(finalMeanData, "~/R/splicing_project/SplicingProject/Results/randomizedJunctionsMeanPerTissue.rda")
  
  
  
  
  ###############################################################
  ########### GET THE MEAN ACROSS INDIVIDUALS ###################
  ###############################################################
  

  
  finalMeanData <-readRDS(file = "~/R/splicing_project/SplicingProject/Results/randomizedJunctionsMeanPerTissue.rda")
  
  ############## GET TISSUES' NAME ######################
  query <- dbSendQuery(con, "select * from sqlite_master;")
  all_tables <- dbFetch(query)
  dbClearResult(query)
  table_names <- all_tables$name
  gtex_table <- table_names[56]
  table_names <- table_names[-1] #Delete "Junc_coordinates" name-table
  table_names <- table_names[-55] #Delete "GTEX" name-table

  ############## GET MEAN ACROSS INDIVIDUALS ######################
  tissuePoints <- list()
  points <- list()

  for(tissue in table_names){
    if(!str_detect(tolower(tissue), "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast")){
      points <- list()
      for(person in names(finalMeanData)){
        print(person)

        if(length(finalMeanData[[person]][[tissue]]) > 0){

          for(item in 1:length(finalMeanData[[person]][[tissue]])){

            if(item<=length(points)){
              points[[item]] <- c(points[[item]],finalMeanData[[person]][[tissue]][[item]])
            }
            else
              points[[item]] <- finalMeanData[[person]][[tissue]][[item]]

          }
        }
      }
      if(length(points)>0)
      tissuePoints[[tissue]] <- points
    }
  }
  mean<-list()
  for(tissue in names(tissuePoints)){

  
    for(item in 1:length(tissuePoints[[tissue]])){
      if(tissue %in% names(mean))
        mean[[tissue]] <- c(mean[[tissue]],mean(tissuePoints[[tissue]][[item]]))
      else
        mean[[tissue]] <- mean(tissuePoints[[tissue]][[item]])
    }
  }
  
  
  
  color <- c("red","blue","green","purple","orange","pink","brown", "grey", "black", "grey",'chartreuse3', "deeppink", 'cornflowerblue', 
             'darkgoldenrod1', 'peachpuff3','mediumorchid2', 'turquoise3', 'wheat4', 'slategray2',"salmon4","tomato","deepskyblue3","brown4",
             "greenyellow")
  tissues<-NULL
  tissuesColour<-NULL
  i<-1
  pdf(paste0("~/R/splicing_project/SplicingProject/Results/Graphics/AllTissues/GreaterThan3/Randomized/ALLTISSUES-AllIndividuals-RandomizedMean-ZOOMED.pdf"))
  for(tissue in names(mean))
  {
    if(str_detect(tolower(tissue), "brain"))
      tissueColor <- "red"
    else if(str_detect(tolower(tissue), "blood"))
      tissueColor <- "blue"
    else
      tissueColor <- "black"
    
    
    #if(str_detect(tolower(tissue), "brain") || str_detect(tolower(tissue), "blood")){
      tissuesColour <- c(tissuesColour,tissueColor)
      tissues<-c(tissues,tissue)
      if(i==1){
        
        #ylim=c(0,5000), xlim = c(0,150),
        plot(unlist(mean[[tissue]], use.names = F), type = "l",  ylim=c(0,6000), xlim = c(0,150),col=tissueColor,xlab="individuals", ylab="Mean unique juncID across randomized indiv", main = paste0("AllTissues - All Individuals - Randomized Mean"))
        
      }
      else{
        print(mean[[tissue]])
        lines(unlist(mean[[tissue]], use.names = F), type = "l", col=tissueColor)
      }
      i<-i+1
    #}
    
  }
  op <- par(cex = 0.5)
  legend(x = 100, y = 6000, legend = tissues, lty = 1, col = tissuesColour)
  dev.off()
}


