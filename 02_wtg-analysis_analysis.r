# Set working directory
setwd ("Q:/1PROJECTS/W. T. Grant-Exploring Treatment Effectiveness/8 Datasets and Analyses/1 Analysis/Alicia Update 2016.03.21/Archive")
setwd ("Q:/1PROJECTS/W. T. Grant-Exploring Treatment Effectiveness/8 Datasets and Analyses/1 Analysis/Alicia Update 2016.03.21/WTGAnalysis")
setwd ("Q:/1PROJECTS/W. T. Grant-Exploring Treatment Effectiveness/8 Datasets and Analyses")


# Activate packages
library(foreign)
library(metafor)
library(dplyr)

# Disable scientific notation in output
options(scipen = 999)

# For reproducability (in random ES selection)
set.seed(62784)

data <- 

##### Prep the dataframes. Some are required and others are optional. Read the comments.

# REQUIRED: FULL DATASET DATAFRAME
data <- dataFinal


# REQUIRED: DATA AND MODEL INFORMATION DATAFRAME
# This block of code creates a data frame that contains information about the dataset.
# This information will be used when when running models and generating model output .csv files
# DO NOT CHANGE unless to add/remove analyses or change dv2_new/prev_strat coding.
key <- data.frame(model=rep("",38), dv=rep("", 38), ps=rep("", 38), es=rep("", 38), stringsAsFactors=FALSE) # Change (increase) 821 if add more meas. chars to examine
n <- 0
# Model Output Name (semi-descriptive), dv2_new, prev_strat (1=universal, 2=indicated), effect size type (SMD, OR)
n <- n+1; key[n,]<-c("adjustmentcoping_u_SMD",1,1,"SMD")
n <- n+1; key[n,]<-c("adjustmentcoping_i_SMD",1,2,"SMD")
n <- n+1; key[n,]<-c("externcrimeSMD_u_SMD",27,1,"SMD")
n <- n+1; key[n,]<-c("externcrimeSMD_i_SMD",27,2,"SMD")
n <- n+1; key[n,]<-c("externotherSMD_u_SMD",29,1,"SMD")
n <- n+1; key[n,]<-c("externotherSMD_i_SMD",29,2,"SMD")
n <- n+1; key[n,]<-c("schoolprogSMD_i_SMD",31,2,"SMD")
n <- n+1; key[n,]<-c("attknowcvd_u_SMD",4,1,"SMD")
n <- n+1; key[n,]<-c("socprobsolve_u_SMD",5,1,"SMD")
n <- n+1; key[n,]<-c("socprobsolve_i_SMD",5,2,"SMD")
n <- n+1; key[n,]<-c("empathy_u_SMD",7,1,"SMD")
n <- n+1; key[n,]<-c("empathy_i_SMD",7,2,"SMD")
n <- n+1; key[n,]<-c("internalizing_u_SMD",13,1,"SMD")
n <- n+1; key[n,]<-c("internalizing_i_SMD",13,2,"SMD")
n <- n+1; key[n,]<-c("loc_u_SMD",15,1,"SMD")
n <- n+1; key[n,]<-c("loc_i_SMD",15,2,"SMD")
n <- n+1; key[n,]<-c("parentfamily_u_SMD",16,1,"SMD")
n <- n+1; key[n,]<-c("parentfamily_i_SMD",16,2,"SMD")
n <- n+1; key[n,]<-c("schoolother_u_SMD",19,1,"SMD")
n <- n+1; key[n,]<-c("schoolother_i_SMD",19,2,"SMD")
n <- n+1; key[n,]<-c("schoolperform_u_SMD",21,1,"SMD")
n <- n+1; key[n,]<-c("schoolperform_i_SMD",21,2,"SMD")
n <- n+1; key[n,]<-c("selfcontrol_u_SMD",22,1,"SMD")
n <- n+1; key[n,]<-c("selfcontrol_i_SMD",22,2,"SMD")
n <- n+1; key[n,]<-c("selfesteem_u_SMD",23,1,"SMD")
n <- n+1; key[n,]<-c("selfesteem_i_SMD",23,2,"SMD")
n <- n+1; key[n,]<-c("peeraccept_u_SMD",24,1,"SMD")
n <- n+1; key[n,]<-c("peeraccept_i_SMD",24,2,"SMD")
n <- n+1; key[n,]<-c("socialskills_u_SMD",25,1,"SMD")
n <- n+1; key[n,]<-c("socialskills_i_SMD",25,2,"SMD")
n <- n+1; key[n,]<-c("schoolcomplete_u_OR",20,1,"OR")
n <- n+1; key[n,]<-c("schoolcomplete_i_OR",20,2,"OR")
n <- n+1; key[n,]<-c("externcrimeOR_u_OR",28,1,"OR")
n <- n+1; key[n,]<-c("externcrimeOR_i_OR",28,2,"OR")
n <- n+1; key[n,]<-c("externotherOR_i_OR",30,2,"OR")
n <- n+1; key[n,]<-c("schoolprogOR_i_OR",32,2,"OR")
n <- n+1; key[n,]<-c("schoolattendance_i_SMD",33,1,"SMD")
n <- n+1; key[n,]<-c("schoolattendance_i_SMD",33,2,"SMD")

##### OPTIONAL: RANDOMLY SELECTED EFFECT SIZE DATA FRAME
# Run this block if you are running single-level analyses in which the effect size is randomly selected
# (However, if you accidentally run this block, it will cause no problems)
data$rand_rnd <- sample(nrow(data), size=nrow(data), replace=FALSE) # Add a column of random numbers to the dataset. replace=FALSE, so no duplicate values.
agg_rnd <- aggregate(data$rand_rnd, by=list(data$studygroup, data$dv2_new, data$prev_strat), min) # Find the smallest randomly generated number within each studygroup+dv+ps combo.
colnames(agg_rnd)[c(1:4)] <- c("studygroup","dv2_new", "prev_strat","min_rand_rnd")
data <- merge(data, agg_rnd, by=c("studygroup","dv2_new","prev_strat"))
data$select_rnd <- ifelse(data$rand_rnd == data$min_rand_rnd, 1, 0) # We'll keep the 1 in our analysis and drop the 0
data <- subset(data, select = -c(rand_rnd, min_rand_rnd))
rm(agg_rnd)

##### OPTIONAL: WEIGHTED MEAN EFFECT SIZE AND WEIGHTED MEAN VARIANCE DATA FRAME
# Run this block if you are running single level analyses in which the effect size and variance are weighted means of the study group
# (However, if you accidentally run this block, it will cause no problems)
agg_wtmeans <-data %>%
  group_by(studygroup, dv2_new, prev_strat) %>%
  summarise(wmean = weighted.mean(tukey1, es3), #Change back to tukey 2 when done with winsorizing experiment! (6/23)
            wvar = weighted.mean(esvar, es3))
agg_wtmeans <- as.data.frame(agg_wtmeans)
data <- merge(data, agg_wtmeans, by=c("studygroup","dv2_new","prev_strat"))
# Below is really an overly complated way to select one of the rows - it doesn't need to be random because they should all be the same. But it works.
data$rand_wt <- sample(nrow(data), size=nrow(data), replace=FALSE) # Add a column of random numbers to the dataset. replace=FALSE, so no duplicate values.
agg_wtmeans <- aggregate(data$rand_wt, by=list(data$studygroup, data$dv2_new, data$prev_strat), min)# Smallest randomly generated number within each studygroup+dv+ps combo.
colnames(agg_wtmeans)[c(1:4)] <- c("studygroup","dv2_new", "prev_strat","min_rand_wt")
data <- merge(data, agg_wtmeans, by=c("studygroup","dv2_new","prev_strat")) # Merge the smallest randomly generated number (by studygroup+dv+ps) back into the dataset
data$select_wt <- ifelse(data$rand_wt == data$min_rand_wt, 1, 0) # We'll keep the 1 in our analysis and drop the 0
data <- subset(data, select = -c(rand_wt, min_rand_wt))
rm(agg_wtmeans)

##### OPTIONAL: NO STUDY GROUP SINGLETONS DATA FRAME
# Run this block if you are running multi-level analyses and wish to drop study groups with only one effect size (within each analysis)
# (However, if you accidentally run this block, it will cause no problems)
data$rand_wt <- sample(nrow(data), size=nrow(data), replace=FALSE)
agg_singletons <- aggregate(data$rand_wt, by=list(data$studygroup, data$dv2_new, data$prev_strat), length)
colnames(agg_singletons)[c(1:4)] <- c("studygroup","dv2_new", "prev_strat","EScount")
agg_singletons$single <- ifelse(agg_singletons$EScount==1, 1, 0) # single = 1 if it's a study group singleton, 0 if not. We will drop the 1, keep the 0.
data <- merge(data, agg_singletons, by=c("studygroup","dv2_new","prev_strat"))
data <- subset(data, select = -c(rand_wt, EScount))
rm(agg_singletons)

##### ANALYSIS ################################################################################################
# The following blocks of code perform different types of analyses.
# Do not run all blocks - only run the block for the analysis you want to perform.
# If multiple analyses are to be run - run one analysis, then export it (see next section) then run the next analysis, then export it, etc.
# Each analysis will be preceeded by a brief description. Be sure to read the rma.mv/uni.mv syntax closely, as you may wish to make small changes (e.g. method, moderators)

##### NO MODERATORS

##### ANALYSIS: Multilevel analysis - full dataset, no moderators
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=es_preadj, V=esvar, intercept=TRUE, 
           data=data, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
    ))
  print(summary(get(key$model[i]))) # These print lines can be commented out if you don't want to clutter the console.
}

### ANALYSIS: using wmean and wvar from espicker syntax
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
           data=data, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
  ))
  print(summary(get(key$model[i]))) # These print lines can be commented out if you don't want to clutter the console.
}

##### ANALYSIS: Multilevel analysis - weights=1
data$wone <- 1
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=es_preadj, 
           V=esvar, 
           W=wone, 
           #digits=20,
           intercept=TRUE, 
           data=data, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
  ))
  print(summary(get(key$model[i]))) # These print lines can be commented out if you don't want to clutter the console.
}

##### ANALYSIS: Multilevel analysis - no singletons (single==0), no moderators 
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=tukey2, V=esvar, intercept=TRUE, 
           data=data, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i] & single==0, 
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
  ))
  print(summary(get(key$model[i])))
}

##### ANALYSIS: Single level analysis - Weighted means and variances (wmean, wvar, select_wt==1), no moderators 
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.uni(yi=wmean, vi=wvar, intercept=TRUE, 
            data=data,
            subset=dv2_new==key$dv[i] & prev_strat==key$ps[i] & select_wt==1,
            method="REML")
  ))
  print(summary(get(key$model[i])))
}

##### ANALYSIS: Single level analysis - Randomly selected means and variances (select_rnd==1), no moderators 
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.uni(yi=es_preadj, vi=esvar, intercept=TRUE, 
            data=data,
            subset=dv2_new==key$dv[i] & prev_strat==key$ps[i] & select_rnd==1,
            method="REML")
  ))
  print(summary(get(key$model[i])))
}


##### ANALYSIS: Single level analysis - Full dataset data as univariate (this is for diagnostic purposes only - not a valid analysis) 
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.uni(yi=tukey2, vi=esvar, intercept=TRUE, 
            data=data,
            subset=dv2_new==key$dv[i] & prev_strat==key$ps[i],
            method="REML")
  ))
  print(summary(get(key$model[i])))
}


##### WITH MODERATORS

# Moderator Set Up
data$mod_crt <- ifelse(data$design_rct==1,1,0) # 1 if crt
data$mod_design_miss <- ifelse(data$design_rct==1 | data$design_crt==1 | data$design_qed==1, 0, 1) # 1 if missing
data$mod_selfreport <- ifelse(data$dv5==1, 1, 0) # 1 if self-report
data$mod_selfreport[is.na(data$mod_selfreport)] <-0 # deal with the NAs (if statement doesn't work on them)
data$mod_informant_miss <- ifelse(is.na(data$dv5) | data$dv5==99, 1, 0)
data$mod_avgage <- data$ageavg # Average age (continuous)
data$mod_avgage_miss <-ifelse(is.na(data$ageavg), 1, 0) #

# List of moderators
mod_list <- c("mod_crt", "mod_selfreport","mod_avgage")


### MODERATOR SECTION IN PROGRESS

##### EXPORT ################################################################################################
# After running each analysis, export the output.
output <-data.frame()
# Add a brief description of the output for the filename
fn <- "Multi_Level_WtEq"

##### Run one of the two code blocks. One works on rma.mv (1 - Run if Multi Level). The other works on rma.uni (2 - Run if Single Level)
##### Run if Multi-Level
for(i in 1:nrow(key)){
  temp_ci <- confint(eval(as.name(key$model[i])))
  output[i,1] <- key$model[i] # Model Name
  output[i,2] <- round(eval(as.name(key$model[i]))[[8]][1], digits=4) # Study Group Level Tau^2
  output[i,3] <- round(sqrt(eval(as.name(key$model[i]))[[8]][1]), digits=4) # Study Group Level Sqrt(Tau^2)
  output[i,4] <- round(temp_ci[[1]][[1]][3], digits=4) # Study group tau^2, CI lower bound
  output[i,5] <- round(temp_ci[[1]][[1]][5], digits=4) # Study group tau^2, CI upper bound
  output[i,6] <- eval(as.name(key$model[i]))$s.nlevels[1] # Study Group Levels
  output[i,7] <- round(eval(as.name(key$model[i]))[[8]][2], digits=4) # Effect Size Level Tau^2
  output[i,8] <- round(sqrt(eval(as.name(key$model[i]))[[8]][2]), digits=4) # Effect Size Level Sqrt(Tau^2)
  output[i,9] <- round(temp_ci[[2]][[1]][3], digits=4) # Effect Size tau^2, CI lower bound
  output[i,10] <- round(temp_ci[[2]][[1]][5], digits=4) # Effect Size tau^2, CI upper bound
  output[i,11] <- eval(as.name(key$model[i]))$s.nlevels[2] # Effect Size Levels
  output[i,12] <- round(eval(as.name(key$model[i]))[[8]][1] + eval(as.name(key$model[i]))[[8]][2], digits=4) # Total variance
  output[i,13] <- round(eval(as.name(key$model[i]))[[8]][1]/(eval(as.name(key$model[i]))[[8]][1] + eval(as.name(key$model[i]))[[8]][2])*100, digits=4) # % Variance at study group level
  output[i,14] <- round(eval(as.name(key$model[i]))[[8]][2]/(eval(as.name(key$model[i]))[[8]][1] + eval(as.name(key$model[i]))[[8]][2])*100, digits=4) # % Variance at effect size level
  output[i,15] <- round(eval(as.name(key$model[i]))[[20]], digits=4) # Q
  output[i,16] <- eval(as.name(key$model[i]))[[14]]-1 # Q df (in output is always = Effect Size levels - 1)
  output[i,17] <- ifelse( key$es[i]=="SMD", 
                              round(eval(as.name(key$model[i]))[[1]][1], digits=4), #IF SMD (Model effect size mean)
                              exp(round(eval(as.name(key$model[i]))[[1]][1], digits=4))) #IF OR (Model effect size mean)
  output[i,18] <- round(eval(as.name(key$model[i]))[[2]][1], digits=4) 
  output[i,19] <- toString(eval(as.name(key$model[1]))$call)
  #mod_length <- length(eval(as.name(key$model[i]))$b)
  #if( mod_length > 1 ){ # If there are moderators (1=intercept)
   # for(j in 2:mod_length) ){
    #  output[i,19 + mod_length - 1] <- paste(row.names(eval(as.name(key$model[i]))$b)[j], round(eval(as.name(key$model[i]))$b[j], digits=4))
   # }
  #}
}
colnames(output) <- c("Model", 
                          "SG Tau^2", "SG Sqrt(Tau^2)","SG CI lower", "SG CI upper", "SG k", 
                          "ES Tau^2", "ES Sqrt(Tau^2)", "ES CI lower", "ES CI upper","ES n", 
                          "Total Variance", "% at SG level", "% at ES level",
                          "Q", "Q df", "Mean ES", "SE", "Model")
write.csv(output, file= paste0("RMAoutput_",fn,"_",format(Sys.time(), "%Y%m%d_%H%M%S_"),".csv") )

##### Run if Single Level
for(i in 1:nrow(key)){
  temp_ci <- confint(eval(as.name(key$model[i])))
  output[i,1] <- key$model[i] # Model Name
  output[i,2] <- round(eval(as.name(key$model[i]))[[8]], digits=4) # Tau^2
  output[i,3] <- round(sqrt(eval(as.name(key$model[i]))[[8]]), digits=4) # Sqrt(Tau^2) 
  output[i,4] <- round(temp_ci[[1]][5], digits=4) # tau^2, CI lower bound
  output[i,5] <- round(temp_ci[[1]][9], digits=4) # tau^2, CI upper bound
  output[i,6] <- eval(as.name(key$model[i]))[[11]]  # Study Group Levels (k)
  output[i,7] <- round(eval(as.name(key$model[i]))[[18]], digits=4) # Q
  output[i,8] <- eval(as.name(key$model[i]))[[11]]-1 # Q df (in output is always = Effect Size levels - 1)
  output[i,9] <- ifelse( i<35 , 
                             round(eval(as.name(key$model[i]))[[1]][1], digits=4),
                             exp(round(eval(as.name(key$model[i]))[[1]][1], digits=4)))# Model effect size mean
  output[i,10] <- round(eval(as.name(key$model[i]))[[2]], digits=4) # Model effect size standard error 
  output[i,11] <- toString(eval(as.name(key$model[1]))$call)
}
colnames(output) <- c("Model", 
                          "Tau^2", "Sqrt(Tau^2)","CI lower", "CI upper", "SG k", 
                          "Q", "Q df", "Mean ES", "SE","Model")
write.csv(output, file= paste0("RMAoutput_",fn,"_",format(Sys.time(), "%Y%m%d_%H%M%S_"),".csv") )












#########################################################################################################################
#########################################################################################################################
# TESTING AND IN PROGRESS
confint(rma.uni(yi=tukey2, vi=esvar, intercept=TRUE, 
        data=data,
        subset=dv2_new==1 & prev_strat==1,
        method="REML"))



# TO DO: HANDLE EXTRACTION OF MODERATOR COEFFICIENTS AND STANDARD ERRORS



mod_list <- c("design_rct", "mod_selfreport")
length(multi2$b) # 1. if >1, we know there are moderators in the model
length(eval(as.name(key$model[i]))$b)
paste(row.names(multi2$b)[2], round(multi2$b[2], digits=4))
paste(row.names(multi2$b)[3], round(multi2$b[3], digits=4))
paste(row.names(multi2$b)[4], round(multi2$b[4], digits=4))

# CREATE MODERATORS
data$mod_crt <- ifelse(data$design_rct==1,1,0) # 1 if crt
data$mod_design_miss <- ifelse(data$design_rct==1 | data$design_crt==1 | data$design_qed==1, 0, 1) # 1 if missing
data$mod_selfreport <- ifelse(data$dv5==1, 1, 0) # 1 if self-report
data$mod_selfreport[is.na(data$mod_selfreport)] <-0 # deal with the NAs (if statement doesn't work on them)
data$mod_informant_miss <- ifelse(is.na(data$dv5) | data$dv5==99, 1, 0)
data$mod_avgage <- data$ageavg # Average age (continuous)
data$mod_avgage_miss <-ifelse(is.na(data$ageavg), 1, 0) #

multi0 <- rma.mv(yi=tukey2, V=esvar, intercept=TRUE, 
                 data=data, 
                 subset = dv2_new==1 & prev_strat==1, 
                 sparse = TRUE, 
                 random = list(~1|studygroup,~1|esid), 
                 method = "REML")
multi1 <- rma.mv(yi=tukey2, V=esvar, intercept=TRUE, 
                 data=data, 
                 subset = dv2_new==1 & prev_strat==1,
                 mods=cbind(design_rct),
                 sparse = TRUE, 
                 random = list(~1|studygroup,~1|esid), 
                 method = "REML")
multi2 <- rma.mv(yi=tukey2, V=esvar, intercept=TRUE, 
                 data=data, 
                 subset = dv2_new==1 & prev_strat==1, 
                 mods=cbind(design_rct, mod_selfreport),
                 sparse = TRUE, 
                 random = list(~1|studygroup,~1|esid), 
                 method = "REML")
uni0 <- rma.uni(yi=wmean, vi=wvar, intercept=TRUE, 
                data=data,
                subset=dv2_new==1 & prev_strat==1 & select_rnd==1,
                method="REML")
uni1 <- rma.uni(yi=wmean, vi=wvar, intercept=TRUE, 
                data=data,
                mods=cbind(design_rct),
                subset=dv2_new==1 & prev_strat==1 & select_rnd==1,
                method="REML")
uni2 <- rma.uni(yi=wmean, vi=wvar, intercept=TRUE, 
                data=data,
                mods=cbind(design_rct, mod_selfreport),
                subset=dv2_new==1 & prev_strat==1 & select_rnd==1,
                method="REML")
listtest <-c("multi0", "multi1", "multi2", "uni0", "uni1","uni2") 




mod_list <- c("design_rct", "mod_selfreport")
length(multi0$b) # 1. if >1, we knwo there are moderators in the model
length(multi1$b) # 2. 
length(multi2$b) # 3. 
row.names(multi2$b)[1] == "intrcpt" #TRUE
row.names(multi2$b)[2] == mod_list[2-1] #TRUE
row.names(multi2$b)[3] == mod_list[3-1] #TRUE. stop here, length= mod_list+1. Have to go off of the list, as we don't know if one model is missing a moderator
# What if a moderator is missing from one of a set of models (due to not varying)

i=1
for(i in 1:length(uni2$b)){
  print(  row.names(uni2$b)[i]  )
  print( uni2$b[i] )
  if ( identical( row.names(uni2$b)[i+1] , mod_list[i]) ){
    print("Match!")
  } else { 
    print ("Missing!") 
  }
}

i<-4
uni2$b
length(uni2$b)
row.names(uni2$b)
row.names(uni2$b)[i] 
is.na(row.names(uni2$b)[i]) # RETURNS TRUE

multi2$b
row.names(multi2$b)
row.names(multi2$b)[i]
is.na(row.names(multi2$b)[i]) # RETURNS TRUE

toString(multi2$b[,1][1])

paste(row.names(multi2$b)[1], round(multi2$b[1], digits=4))


# EXPORT ESID FOR FILE WITH COMPLETE EFFECT SIZES TO MERGE INTO SPSS.
# THIS PERMITS FILTERING FOR DESCRIPTIVES IN SPSS.
# Set working directory
setwd ("Q:/1PROJECTS/W. T. Grant-Exploring Treatment Effectiveness/8 Datasets and Analyses")
dataexp<-subset(dataFinal, select=c(studygroup, esid))
write.csv(dataexp, file="datafinal.csv")

# meta-regression models

universalSMD<-subset(dataFinal, dv2_new==)

key1 <- data.frame(model=rep("",4), ps=rep("", 4), es=rep("", 4), stringsAsFactors=FALSE) # Change (increase) 821 if add more meas. chars to examine
n <- 0
# Model Output Name (semi-descriptive), prev_strat (1=universal, 2=indicated), effect size type (SMD, OR)
n <- n+1; key1[n,]<-c("universal smd",1,"SMD")
n <- n+1; key1[n,]<-c("universal OR",1,"OR")
n <- n+1; key1[n,]<-c("selected smd",2,"SMD")
n <- n+1; key1[n,]<-c("selected OR",2,"OR")


saveRDS(dataFinal, file="R_ANALYTIC_DATASETv2.rds")
dataFinal<-readRDS(file="R_ANALYTIC_DATASETv2.rds")
ls(dataFinal)


##### ANALYSIS: Multilevel analysis - full dataset, no moderators
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
           data=dataFinal, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
  ))
  print(summary(get(key$model[i]))) # These print lines can be commented out if you don't want to clutter the console.
}

##### ANALYSIS: Multilevel analysis - full dataset, no moderators, 2-level
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
           data=dataFinal, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           sparse = TRUE, 
           random = list(~1|studygroup), 
           method = "ML")
  ))
  print(summary(get(key$model[i]))) # These print lines can be commented out if you don't want to clutter the console.
}

summary()

##### ANALYSIS: Multilevel analysis - full dataset, no moderators, univariate
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.uni(yi=wmean, vi=wvar, intercept=TRUE, 
           data=dataFinal, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           sparse = TRUE, 
           method = "REML")
  ))
  print(summary(get(key$model[i])),transf=exp) # These print lines can be commented out if you don't want to clutter the console.
}

summary(externcrimeOR_i_OR)
summary(externcrimeOR_u_OR)
summary(externotherOR_i_OR)


predict(externcrimeOR_i_OR, transf=exp)
predict(externcrimeOR_u_OR, transf=exp)
predict(externotherOR_i_OR, transf=exp)
predict(schoolprogOR_i_OR, transf=exp)
predict(schoolcomplete_i_OR, transf=exp)
predict(schoolcomplete_u_OR, transf=exp)

# CREATE MODERATORS
dataFinal$decade<-dataFinal$h1/1000
dataFinal$ssize<-log(dataFinal$es3)
dataFinal$dv5_sr <- ifelse(dataFinal$dv5==1,1,0) # 1 if self-report
dataFinal$dv5_teach <- ifelse(dataFinal$dv5==4,1,0) # 1 if teacher report
dataFinal$dv5_rec <- ifelse(dataFinal$dv5==7|dataFinal$dv5==8|dataFinal$dv5==10 ,1,0) # 1 if records
dataFinal$dv5_obs <- ifelse(dataFinal$dv5==13,1,0) # 1 if observation
dataFinal$dv5_par <- ifelse(dataFinal$dv5==2,1,0) # 1 if parent report

table(dataFinal$h2alt_d1)
table(dataFinal$h8h9_d3)
table(dataFinal$h38r_d1)
table(dataFinal$h38r_d2)
table(dataFinal$h38r_d3)
table(dataFinal$h38r_d4)


#create subsets of final dataset for prev_strat and type of es
adjustmentcoping_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==1)
adjustmentcoping_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==1)
externcrime_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==27)
externcrime_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==27)
externother_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==29)
externother_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==29)
schoolprog_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==31)
attknowcvd_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==4)
socprobsolve_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==5)
socprobsolve_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==5)
empathy_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==7)
empathy_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==7)
internalizing_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==13)
internalizing_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==13)
loc_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==15)
loc_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==15)
parentfamily_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==16)
parentfamily_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==16)
schooleng_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==19)
schooleng_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==19)
schoolperform_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==21)
schoolperform_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==21)
selfcontrol_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==22)
selfcontrol_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==22)
selfesteem_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==23)
selfesteem_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==23)
peeraccept_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==24)
peeraccept_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==24)
socialskills_u_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==25)
socialskills_i_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==25)
schoolcomplete_u_OR<-subset(dataFinal, prev_strat==1 & dv2_new==20)
schoolcomplete_i_OR<-subset(dataFinal, prev_strat==2 & dv2_new==20)
externcrime_u_OR<-subset(dataFinal, prev_strat==1 & dv2_new==28)
externcrime_i_OR<-subset(dataFinal, prev_strat==2 & dv2_new==28)
externother_i_OR<-subset(dataFinal, prev_strat==2 & dv2_new==30)
schoolprog_i_OR<-subset(dataFinal, prev_strat==2 & dv2_new==32)
schoolattendance_i_SMD<-subset(dataFinal, prev_strat==1 & dv2_new==33)
schoolattendance_u_SMD<-subset(dataFinal, prev_strat==2 & dv2_new==33)


table(adjustmentcoping_u_SMD$dv5)
table(externother_u_SMD$dv5)
table(externother_i_SMD$dv5)
table(schoolperform_u_SMD$dv5)
table(schoolperform_i_SMD$dv5)




table(dataFinal$h35_d1)
table(externother_u_SMD$routine_eval_dev)

#first pass at moderator analysis - method and context
extusmd1<-rma.uni(yi=wmean, vi=wvar, intercept=TRUE, data=externother_u_SMD, sparse = TRUE, method = "DL",
                  mods =~ decade1 + h2alt_d1 + h2alt_d2+dv5_sr + dv5_teach)
summary(extusmd1)





##univariate
extusmd6a<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                 data=externother_u_SMD, 
                 mods= ~ decade ,
                 sparse = TRUE, 
                 random = list(~1|studygroup,~1|esid), 
                 method = "REML")
summary(extusmd6a)

extusmd6b<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~ h2alt_d1 + h2alt_d2 ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6b)

extusmd6c<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~  h8h9_d1 + h8h9_d2 ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6c)

extusmd6d<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~ h14_d1 + h14_d2 ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6d)

extusmd6e<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~  h33 ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6e)
extusmd6f<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~ h35_d1 ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6f)
extusmd6g<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~ dv5_sr + dv5_teach ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6g)
extusmd6h<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~ es3,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6h)
extusmd6i<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~ factor(routine_eval_dev) -1 ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6i)
extusmd6j<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~  h25alt_d1 ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6j)
extusmd6k<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~  h22_d2 ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6k)

extusmd6l<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~  h26h85_d1 ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6l)

extusmd6m<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~  h37r_d1 ,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6m)

extusmd6n<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~  h68alt_d1,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd6n)

##multivariate
extusmd2<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                  data=externother_u_SMD, 
                  mods= ~ decade + h2alt_d1 + h2alt_d2 + h8h9_d1 + h8h9_d2 + h14_d1 + h14_d2 + h33 + h35_d1 + dv5_sr + dv5_teach + es3,
                  sparse = TRUE, 
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(extusmd2)

extusmd3<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                 data=externother_u_SMD, 
                 mods= ~ decade + h2alt_d1 + h2alt_d2 + h8h9_d1 + h8h9_d2 + h14_d1 + h14_d2 + h33 + h35_d1 + dv5_sr + dv5_teach + es3
                        + factor(routine_eval_dev) -1,
                 sparse = TRUE, 
                 random = list(~1|studygroup,~1|esid), 
                 method = "REML")
summary(extusmd3)

extusmd4<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                 data=externother_u_SMD, 
                 mods= ~ decade + h2alt_d1 + h2alt_d2 + h8h9_d1 + h8h9_d2 + h14_d1 + h14_d2 + h33 + h35_d1 + dv5_sr + dv5_teach + es3
                 + factor(routine_eval_dev) -1 + h25alt_d1,
                 sparse = TRUE, 
                 random = list(~1|studygroup,~1|esid), 
                 method = "REML")
summary(extusmd4)

table(externother_u_SMD$h68alt_d1)

extusmd5<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                 data=externother_u_SMD, 
                 mods= ~ decade + h2alt_d1 + h2alt_d2 + h8h9_d1 + h8h9_d2 + h14_d1 + h14_d2 + h33 + h35_d1 + dv5_sr + dv5_teach + es3
                 + factor(routine_eval_dev) -1 + h25alt_d1 + h22_d2 + h26h85_d1 + h37r_d1,
                 sparse = TRUE, 
                 random = list(~1|studygroup,~1|esid), 
                 method = "REML")
summary(extusmd5)

table(externother_u_SMD$h68alt_d1)

extusmd6<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                 data=externother_u_SMD, 
                 mods= ~ decade + h2alt_d1 + h2alt_d2 + h8h9_d1 + h8h9_d2 + h14_d1 + h14_d2 + h33 + h35_d1 + dv5_sr + dv5_teach + es3
                 + factor(routine_eval_dev) -1 + h25alt_d1 + h22_d2 + h26h85_d1 + h37r_d1 + h68alt_d1,
                 sparse = TRUE, 
                 random = list(~1|studygroup,~1|esid), 
                 method = "REML")
summary(extusmd6)




adjsmd1<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                data=adjustmentcoping_u_SMD, 
                mods= ~ decade + h2alt_d1 + h8h9_d3 + h33 + h35_d1 + dv5_sr + es3,
                sparse = TRUE, 
                random = list(~1|studygroup,~1|esid), 
                method = "REML")
summary(adjsmd1)

adjsmd1<-rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
                data=adjustmentcoping_u_SMD, 
                mods= ~ decade + h2alt_d1 + h8h9_d3 + h33 + h35_d1 + dv5_sr + es3,
                sparse = TRUE, 
                random = list(~1|studygroup,~1|esid), 
                method = "REML")
summary(adjsmd1)



##############################################################################################################################################
###code graveyard and saving account
##### ANALYSIS: Multilevel analysis - full dataset, no moderators
for (i in 1:nrow(key1)){
  assign(eval(key1$model[i]),eval(
    rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
           data=dataFinal, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           mods=cbind(design_rct, mod_selfreport),
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
  ))
  print(summary(get(key$model[i]))) # These print lines can be commented out if you don't want to clutter the console.
}

##### ANALYSIS: Multilevel analysis - full dataset, Study and Method Mods
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=wmean, V=wvar, intercept=TRUE, 
           data=dataFinal, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           mods=cbind(decade, design_rct, mod_selfreport),
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
  ))
  print(summary(get(key$model[i]))) # These print lines can be commented out if you don't want to clutter the console.
}



table(universalSMD$h35)
table(universalSMD$ct_diffcontext)
table(universalSMD$ct_samecontext)
table(universalSMD$evaluator_role)


indicated<-subset(data, prev_strat==2)
universalSMD<-subset(universal, dv2_new!=28 & dv2_new!=20)
universalOR<-subset(universal, dv2_new==28 | dv2_new==20)
indicatedSMD<-subset(indicated, dv2_new!=28 & dv2_new!=30 & dv2_new!=32 & dv2_new!=20)
indicatedOR<-subset(indicated,dv2_new==28|dv2_new==30|dv2_new==32|dv2_new==20)
unidv2_29<-subset(universalSMD, dv2_new==29)
unidv2_25<-subset(universalSMD, dv2_new==25)
unidv2_21<-subset(universalSMD, dv2_new==21)
unidv2_5<-subset(universalSMD, dv2_new==5)

uni25smd2<-rma.mv(yi=es_preadj, V=esvar, intercept=TRUE, 
                  data=unidv2_25, 
                  sparse = TRUE, 
                  mods= ~ h2_d1 + h2_d2 + design_rct + design_crt + h33 + ct_samecontext + dv5_sr + dv5_teach,
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(uni25smd2)

uni21smd2<-rma.mv(yi=es_preadj, V=esvar, intercept=TRUE, 
                  data=unidv2_21, 
                  sparse = TRUE, 
                  mods= ~ h2_d1 + h2_d2 + design_rct + design_crt + h33 + ct_samecontext + dv5_sr + dv5_teach,
                  random = list(~1|studygroup,~1|esid), 
                  method = "REML")
summary(uni21smd2)

uni5smd2<-rma.mv(yi=es_preadj, V=esvar, intercept=TRUE, 
                 data=unidv2_5, 
                 sparse = TRUE, 
                 mods= ~ h2_d1 + h2_d2 + design_rct + design_crt + h33 + ct_samecontext + dv5_sr + dv5_teach,
                 random = list(~1|studygroup,~1|esid), 
                 method = "REML")
summary(uni5smd2)

ls(dataFinal)

dataexp<-subset(dataFinal, select=c(studygroup, esid, dv5, dv5_sr, dv5_teach, dv5_rec, h8h9_d1, h8h9_d2, h33, h35_d1, ssize,
                                    adjes, ageavg, agerange, analyticestype, assign_crandom, assign_imatched, assign_irandom, assign_qed, 
                                    attrition1, avg_units_ct_cluster, avg_units_tx_cluster, cluster_err_ct, 
))
write.csv(dataexp, file="datafinal.csv")



