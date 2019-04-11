
##########################################
## Get the unique junctions per sample ###
##########################################

all_data_3 <- readRDS("/home/sruiz/R/splicing_project/SplicingProject/Results/nonSexSpecificJunctionsGreaterThan3.rda")
#all_data <- readRDS("/home/sruiz/R/splicing_project/SplicingProject/Results/nonSexSpecificJunctions_NumEvents_GT1.rda")
all_data_1 <- readRDS("/home/sruiz/R/splicing_project/SplicingProject/Results/nonSexSpecificJunctionsGreaterThan1.rda")

all_data_3$tissueIndividuals <- all_data_1$col_names

all_data_3$tissueIndividuals <- sapply(names(all_data_3$tissueIndividuals),function(tissue){
  print(tissue)
  tmp <- c()
  for(i in names(all_data_3$juncIDs_person_tissue))
  {
    
    if(any((tissue == names(all_data_3$juncIDs_person_tissue[[i]]))))
    {
      tmp <- c(tmp,i)
    }
  }
  return(tmp)
})

tmp_unique <- list()
for(tissue in names(all_data_3$tissueIndividuals)){
  print(tissue)
  tmp <- list()
  for(i in 1:length(all_data_3$tissueIndividuals[[tissue]]))
  {
    print(i)
    
    all_other <- lapply(all_data_3$juncIDs_person_tissue[all_data_3$tissueIndividuals[[tissue]][-i]],function(x){
      return(x[[tissue]])
    })
    
    tmp[[all_data_3$tissueIndividuals[[tissue]][i]]] <- Reduce(setdiff,
                                                             c(all_data_3$juncIDs_person_tissue[[all_data_3$tissueIndividuals[[tissue]][i]]][tissue],
                                                               all_other))
  }
  tmp_unique[[tissue]] <- tmp
  rm(tmp)
  
}

all_data_3$tissueIndividuals <- tmp_unique
rm(tmp_unique)

tmp <- list()
for (tissue in names(all_data_3$tissueIndividuals))
{
  for(j in names(all_data_3$tissueIndividuals[[tissue]]))
  {
    tmp[[tissue]][[j]] <- length(all_data_3$tissueIndividuals[[tissue]][[j]])/length(all_data_3$juncIDs_person_tissue[[j]][[tissue]]) 
  }
}


all_data_3$proportionUnique <- tmp

# save(all_data_3, file="/home/sguelfi/projects/R/SplicingProject/data/JunctionsGreaterThan3UniqueSplitReadPerSample.rda")
load(file="/home/sguelfi/projects/R/SplicingProject/data/JunctionsGreaterThan3UniqueSplitReadPerSample.rda")

color_code <- names(all_data_3$proportionUnique)
names(color_code) <- names(all_data_3$proportionUnique)
color_code[1:length(color_code)] <- "grey"
color_code[grep("Brain",names(color_code))] <- "yellow"

par(mar=c(12,4,2,2))

ordered_unique <- sort(unlist(lapply(all_data_3$proportionUnique,mean)),decreasing = T)

boxplot(all_data_3$proportionUnique[names(ordered_unique)],las=2,col = color_code[names(ordered_unique)],main="Proportion of uniqueness",
        ylab="Proportion of unique junctions per sample")



library(RSQLite)
database_path="/data/splicing_tolerance/splicing_tolerance.sqlite"
db <-  dbConnect(SQLite(), dbname=database_path)
res <- dbSendQuery(db, paste0("SELECT age,mapped_read_count,smtsd,subj_id,sample_recount_id  FROM GTEX_info WHERE smafrze=='USE ME'"))

GTEx <- dbFetch(res)
head(GTEx)

mean_mapped <- GTEx %>% group_by(smtsd) %>% summarize(mean_size = mean(mapped_read_count, na.rm = TRUE))
GTEx$age <- as.numeric(as.factor(GTEx$age))
mean_age <- GTEx %>% group_by(smtsd) %>% summarize(mean_age = mean(age, na.rm = TRUE))

mean_mapped$smtsd <- gsub('\\(','_',mean_mapped$smtsd)
mean_mapped$smtsd <- gsub('\\)','',mean_mapped$smtsd)
mean_mapped$smtsd <- gsub(" ","",mean_mapped$smtsd)

mean_proportion_unique <- unlist(lapply(all_data_3$proportionUnique,mean))

mean_mapped <- mean_mapped %>% filter(smtsd %in% names(mean_proportion_unique))

par(mar=c(4,4,2,2))
color_code[color_code == "grey"] <- "black"
color_code[color_code == "yellow"] <- "#FFD700"
plot(mean_mapped$mean_size,mean_proportion_unique[gsub(" ","",mean_mapped$smtsd)],pch=16,
     xlab="Mean mapped reads", ylab="Proportion of splicing uniqueness",col=color_code[names(mean_proportion_unique)])
cor.test(mean_mapped$mean_size,mean_proportion_unique[gsub(" ","",mean_mapped$smtsd)])


plot(mean_age$mean_age,mean_proportion_unique[gsub(" ","",mean_age$smtsd)],
     xlab="Mean age", ylab="Proportion of splicing uniqueness",xaxt="n",
     col=color_code[names(mean_proportion_unique)],pch=16)
axis(1, at=1:6, labels=c("20-29", "30-39", "40-49","50-59", "60-69" ,"70-79"))
cor.test(mean_age$mean_age,mean_proportion_unique[gsub(" ","",mean_age$smtsd)])

##############################################################################################################
### Obtain the proportion of unique splice junctions that an individual expresses in all the brain tissues ###
##############################################################################################################

Proportion_personalised_splicing <- list()
numb_personal_spicing <- list()
prop_personal_spicing_at_least_one <- list()

for(x in names(all_data_3$juncIDs_person_tissue))
{
  tmp <- all_data_3$juncIDs_person_tissue[[x]]
  if(length(tmp)>5)
  {
    tmp_per_ind <- lapply(all_data_3$tissueIndividuals[names(tmp)],function(y)
    {
      return(y[[x]])
    })
    Proportion_personalised_splicing[x] <- length(Reduce(intersect, tmp_per_ind ))/length(Reduce(union, tmp_per_ind ))
    numb_personal_spicing[x] <- length(Reduce(intersect, tmp_per_ind ))
    prop_personal_spicing_at_least_one[x] <- sum(table(unlist(tmp_per_ind))>1)/length(unique(unlist(tmp_per_ind)))
  }
  
}



par(mar=c(4,4,2,2))
hist(unlist(prop_personal_spicing_at_least_one)*100,breaks=100,main="Percentage of personal splicing",
     xlab = "% of personal splicing",ylab="# individuals")

unlist(prop_personal_spicing_at_least_one)[unlist(prop_personal_spicing_at_least_one)>0.05]

numTissuePerSample <- lapply(all_data_3$juncIDs_person_tissue,function(x){
  return(length(x))
})

plot((unlist(prop_personal_spicing_at_least_one)*100),
     unlist(numTissuePerSample[names(unlist(prop_personal_spicing_at_least_one))]),
     xlab="% of personal splicing",ylab=" Number of samples per individual")
cor((unlist(prop_personal_spicing_at_least_one)*100),
    unlist(numTissuePerSample[names(unlist(prop_personal_spicing_at_least_one))]))



library(RSQLite)
database_path="/data/splicing_tolerance/splicing_tolerance.sqlite"
db <-  dbConnect(SQLite(), dbname=database_path)
res <- dbSendQuery(db, paste0("SELECT age,mapped_read_count,smtsd,subj_id,sample_recount_id  FROM GTEX_info WHERE smafrze=='USE ME'"))

GTEx <- dbFetch(res)
head(GTEx)

tmp <- as.data.frame(unlist(prop_personal_spicing_at_least_one)) %>% rownames_to_column("subj_id")

tmp_2 <- GTEx %>% filter(!duplicated(subj_id)) %>% select(age,subj_id) %>% 
  inner_join(tmp)

plot(factor(tmp_2[,1]),as.numeric(tmp_2[,3]*100),
     ylab="% of personal splicing", xlab="Age ranges")

cor.test(as.numeric(factor(tmp_2[,1])),as.numeric(tmp_2[,3]*100))
GTEx$age_num <- as.numeric(factor(GTEx$age))

abline(h=-log10(0.05),col="red")




cor_per_tissue <- lapply(all_data_3$proportionUnique,function(x)
{
  tmp <- as.data.frame(unlist(x)) %>% rownames_to_column("subj_id")
  tmp_2 <- GTEx %>% filter(!duplicated(subj_id)) %>% select(age_num,subj_id) %>% 
    inner_join(tmp)
  return(cor.test(as.numeric(factor(tmp_2[,1])),as.numeric(tmp_2[,3]*100)))
})


## correlation
tmp <- lapply(cor_per_tissue, function(x){
  x$estimate
})

tmp_2 <- sort(unlist(tmp),decreasing = T)
color_code <- c()
color_code[1:length(tmp_2)] <- "grey"
names(color_code) <- names(tmp_2)
color_code[grep("Brain",names(color_code))] <- "yellow"

par(mar=c(12,4,2,2))

names(tmp_2) <- gsub(".cor", "",names(tmp_2))
barplot(tmp_2,las=2,ylim = c(-0.30,0.30),col=color_code,ylab="Correlation")




tmp <- lapply(cor_per_tissue, function(x){
  x$p.value
})

tmp_2 <- sort(-log10(unlist(tmp)),decreasing = T)
color_code <- c()
color_code[1:length(tmp_2)] <- "grey"
names(color_code) <- names(tmp_2)
color_code[grep("Brain",names(color_code))] <- "yellow"

par(mar=c(12,4,2,2))

names(tmp_2) <- gsub(".cor", "",names(tmp_2))
barplot(tmp_2,las=2,col=color_code,ylab="-log10(p.value)")
abline(h=-log10(0.05),col="red")



### Here we compare 

## first we compare the personal splicing 




# personal_splincing_perc <- list()
# for(i in names(all_data_3$juncIDs_person_tissue))
# {
#   print(i)
#   per_tissue_ind <- list()
#   ## we now get all tissue per sample
#   for(t in names(all_data_3$juncIDs_person_tissue[[i]])){
#     per_tissue <- list()
#     ## we compare against all tissues
#     for(t_to_compare in names(all_data_3$juncIDs_person_tissue[[i]])){
#       
#       per_tissue[t_to_compare] <- (length(intersect(all_data_3$tissueIndividuals[[t]][[i]],all_data_3$juncIDs_person_tissue[[i]][[t_to_compare]]))/
#                           length(all_data_3$juncIDs_person_tissue[[i]][[t_to_compare]]))*100
#       
#     }
#     per_tissue_ind[[t]] <- per_tissue
#     rm(per_tissue)
#   }
#   personal_splincing_perc[[i]] <- per_tissue_ind
#   rm(per_tissue_ind)
# }


## Get the comparison for individuals and tissues 
inter_personal_splincing_perc <- list()
for(i in names(all_data_3$juncIDs_person_tissue))
{
  print(paste(i,Sys.time()))
  per_tissue_ind <- list()
  ## we now get all tissue per sample
  for(t in names(all_data_3$juncIDs_person_tissue[[i]])){
    ## we compare against all tissues
    per_ind <- list()
    for(ind_to_compare in names(all_data_3$tissueIndividuals[[t]]))
    {
      per_tissue <- list()
      for(t_to_compare in names(all_data_3$juncIDs_person_tissue[[ind_to_compare]])){
        
        per_tissue[t_to_compare] <- (length(intersect(all_data_3$tissueIndividuals[[t]][[i]],all_data_3$juncIDs_person_tissue[[ind_to_compare]][[t_to_compare]]))/
                                       length(all_data_3$juncIDs_person_tissue[[ind_to_compare]][[t_to_compare]]))*100
        
      }
      per_ind[[ind_to_compare]] <- per_tissue
      rm(per_tissue)
    }
    
    per_tissue_ind[[t]] <- per_ind
    rm(per_ind)
  }
  inter_personal_splincing_perc[[i]] <- per_tissue_ind
  rm(per_tissue_ind)
}

save(inter_personal_splincing_perc, file="/home/sguelfi/projects/R/SplicingProject/data/interPersonalSimilarityScore.rda")






per_ind <- list()

for(i in names(inter_personal_splincing_perc))
{
  ## loop through the tissues
  per_tissue_ind <- list()
  for (t in names(inter_personal_splincing_perc[[i]]))
  {
    tmp_t <- names(inter_personal_splincing_perc[[i]])
    tmp_t <- tmp_t[which(tmp_t !=  t)]
    per_tissue <- list()
    for (t_to_compare in tmp_t)
    {
      tmp_i <- names(inter_personal_splincing_perc[[i]][[t]])
      val_per <- inter_personal_splincing_perc[[i]][[t]][[i]][[t_to_compare]]
      tmp_percentage_distr <- lapply(inter_personal_splincing_perc[[i]][[t]][which(tmp_i != i)],function(x)
      {
        return(x[[t_to_compare]])
      })
      
      per_tissue[t_to_compare] <- (1+sum(unlist(tmp_percentage_distr) >= val_per))/(length(tmp_percentage_distr)+1)
      rm(tmp_i,val_per,tmp_percentage_distr)
    }
    per_tissue_ind[[t]] <- per_tissue
    rm(per_tissue,tmp_t)
  }
  per_ind[[i]] <- per_tissue_ind
  rm(per_tissue_ind)
}



lapply(per_ind$`GTEX-QMR6`$WholeBlood,function(x){
  x[[]]
})

par(mar=c(12,4,2,2))
barplot(sort(-log10(unlist(per_ind$`GTEX-QMR6`$`Brain-FrontalCortex_BA9`)),decreasing = T),las=2,ylab="-log10(p-value)",
        main="Ind GTEX-QMR6 frontalCortex")
abline(h=-log10(0.05),col="red")

par(mar=c(12,4,2,2))
barplot(sort(-log10(unlist(per_ind$`GTEX-QMR6`$WholeBlood)),decreasing = T),las=2,ylab="-log10(p-value)",
        main="Ind GTEX-QMR6 whole blood")
abline(h=-log10(0.05),col="red")


par(mar=c(12,4,2,2))
barplot(sort(-log10(unlist(per_ind$`GTEX-POMQ`$WholeBlood)),decreasing = T),las=2)
abline(h=-log10(0.05),col="red")



-log10(unlist(per_ind$`GTEX-QMR6`$`Brain-FrontalCortex_BA9`$`Brain-Caudate_basalganglia`))
par(mar=c(12,4,2,2))
barplot(sort(-log10(unlist(per_ind$`GTEX-QMR6`$`Brain-FrontalCortex_BA9`)),decreasing = T),las=2)
abline(h=-log10(0.05),col="red")




hist(unlist(lapply(inter_personal_splincing_perc[[i]][[t]][which(tmp != i)],function(x)
{
  print(x[["Brain-Cortex"]])
})),breaks=40)

abline(v=inter_personal_splincing_perc[[i]][[t]][[i]][["Brain-Cortex"]],col="red")


inter_personal_splincing_perc$`GTEX-QMR6`$`Brain-Amygdala`$`GTEX-OHPN`$`Brain-Cortex`

for(t in names(all_data_3$tissueIndividuals)){

  for(s in names(all_data_3$tissueIndividuals[[t]])){
      personal <- lapply(all_data_3$tissueIndividuals[[t]],function(j){
             return((length(intersect(j[[s]],all_data_3$juncIDs_person_tissue[[s]][[t]]))/
             length(all_data_3$juncIDs_person_tissue[[s]][[t]]))*100)
  })
  }
  
}





for(t in names(all_data_3$tissueIndividuals)){

  t <- names(all_data_3$tissueIndividuals)[1]
    
  for(s in names(all_data_3$tissueIndividuals[[t]])){
    
  ## take all other tissues to compare but t  
  personal <- lapply(all_data_3$tissueIndividuals[[-t]],function(j){
    return(length(intersect(j[[s]],all_data_3$juncIDs_person_tissue[s][[))/
                  length(all_data_3$juncIDs_person_tissue$`GTEX-QMR6`$`Brain-Cortex`))*100
    
    
  })
    
                          ) 
val_per <- (length(intersect(all_data_3$tissueIndividuals$`Brain-FrontalCortex_BA9`$`GTEX-QMR6`,
                             all_data_3$juncIDs_person_tissue$`GTEX-QMR6`$`Brain-Cortex`))/
              length(all_data_3$juncIDs_person_tissue$`GTEX-QMR6`$`Brain-Cortex`))*100
  
      
    lapply(all_data_3$tissueIndividuals[[t]])
  {
      
  })
  
  val_per <- (length(intersect(all_data_3$tissueIndividuals[[t]]
    all_data_3$tissueIndividuals$`Brain-FrontalCortex_BA9`$`GTEX-QMR6`,
                                 all_data_3$juncIDs_person_tissue$`GTEX-QMR6`$`Brain-Cortex`))/
                  length(all_data_3$juncIDs_person_tissue$`GTEX-QMR6`$`Brain-Cortex`))*100
    
    
    
  }
    
  
  
  
  distributionOfInter <- lapply(all_data_3$juncIDs_person_tissue[names(all_data_3$tissueIndividuals$`Brain-Cortex`)[-1]],function(x)
  {
    (length(intersect(all_data_3$tissueIndividuals$`Brain-FrontalCortex_BA9`$`GTEX-QMR6`,
                      x[["Brain-Cortex"]]))/length(x[["Brain-Cortex"]]))*100
  })
  
  
  round((1+sum(unlist(distributionOfInter) >= val_per))/(length(distributionOfInter)+1),3)
  
  
    
}



mean(unlist(distributionOfInter))


hist(unlist(distributionOfInter),breaks = 40,xlim=c(0,0.15))
abline(v=length(intersect(all_data_3$tissueIndividuals$`Brain-FrontalCortex_BA9`$`GTEX-QMR6`,
                          all_data_3$juncIDs_person_tissue$`GTEX-QMR6`$`Brain-Cortex`))/
         length(all_data_3$juncIDs_person_tissue$`GTEX-QMR6`$`Brain-Cortex`)*100,
       col="red")

       
       


all_data_3$tissueIndividuals$`Brain-FrontalCortex_BA9`$`GTEX-QMR6`













