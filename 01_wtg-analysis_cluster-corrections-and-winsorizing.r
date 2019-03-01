# Set working directory
setwd ("Q:/1PROJECTS/W. T. Grant-Exploring Treatment Effectiveness/8 Datasets and Analyses")

# Activate packages
library(foreign)
library(metafor)
library(dplyr)
library(robustbase)

# Disable scientific notation in output
options(scipen = 999)

# For reproducability (in random ES selection)
set.seed(62784)

# IMPORT DATAFRAMES
pre <- read.spss("06 PreTest.sav",use.value.labels=F,to.data.frame=T)
post <- read.spss("06 PostTestAnalytic.sav",use.value.labels=F,to.data.frame=T) #unimputed version
post <- read.spss("06 PostTestAnalyticImputed.sav",use.value.labels = F,to.data.frame = T) #imputed

# Inner and outer tukey fence functions
# Note that these versions will bring the value to the fence, not the highest/lowest non-outlier
# (should be able to condense all of these versions into one function, but work on that for fun later)
inner_tukey <- function(x, data){
  arguments <- as.list(match.call()) # See: http://www.r-bloggers.com/passing-columns-of-a-dataframe-to-a-function-without-quotes/
  x = eval(arguments$x, data) # Same as above.
  print(round(quantile(x)[1:5], digits=4))
  low_tukey <- quantile(x, names=FALSE)[2] - (1.5 * (quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2])) # Get low fence (Q1 - 1.5*IQR)
  high_tukey <- quantile(x, names=FALSE)[4] + (1.5 * (quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2])) # Get high fence (Q3 + 1.5*IQR)
  high_outlier <- ifelse(x >= high_tukey, "yes", "no") # Is this effect size a high outlier? Yes/no
  low_outlier <- ifelse(x <= low_tukey, "yes", "no") # Is this effect size a low outlier? Yes/no
  if(all(x==x[1])){ # Special case for when all values are the same
    es_tukey <- x
  } else {
    es_tukey <- ifelse(high_outlier=="no" & low_outlier=="no", x, NA) # Copy values for non-outliers or if all the same, NAs for outliers
    es_tukey <- ifelse(high_outlier=="yes", high_tukey, es_tukey) # 
    es_tukey <- ifelse(low_outlier=="yes", low_tukey, es_tukey) # 
  }
  return(es_tukey) # Return the effect sizes
}
inner_tukey_skew <- function(x, data){
  arguments <- as.list(match.call()) # See: http://www.r-bloggers.com/passing-columns-of-a-dataframe-to-a-function-without-quotes/
  x = eval(arguments$x, data) # Same as above.
  low_tukey <- adjboxStats(x, coef=1.5)$fence[1]  # Get low fence 
  high_tukey <-  adjboxStats(x, coef=1.5)$fence[2] # Get high fence (Q3 + 1.5*IQR)
  high_outlier <- ifelse(x >= high_tukey, "yes", "no") # Is this effect size a high outlier? Yes/no
  low_outlier <- ifelse(x <= low_tukey, "yes", "no") # Is this effect size a low outlier? Yes/no
  if(all(x==x[1])){ # Special case for when all values are the same
    es_tukey <- x
  } else {
    es_tukey <- ifelse(high_outlier=="no" & low_outlier=="no", x, NA) # Copy values for non-outliers or if all the same, NAs for outliers
    es_tukey <- ifelse(high_outlier=="yes", high_tukey, es_tukey) # 
    es_tukey <- ifelse(low_outlier=="yes", low_tukey, es_tukey) # 
  }
  return(es_tukey) # Return the effect sizes
}
outer_tukey <- function(x, data){
  arguments <- as.list(match.call()) # See: http://www.r-bloggers.com/passing-columns-of-a-dataframe-to-a-function-without-quotes/
  x = eval(arguments$x, data) # Same as above.
  print(round(quantile(x)[1:5], digits=4))
  low_tukey <- quantile(x, names=FALSE)[2] - (3 * (quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2])) # Get low fence (Q1 - 1.5*IQR)
  high_tukey <- quantile(x, names=FALSE)[4] + (3 * (quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2])) # Get high fence (Q3 + 1.5*IQR)
  high_outlier <- ifelse(x >= high_tukey, "yes", "no") # Is this effect size a high outlier? Yes/no
  low_outlier <- ifelse(x <= low_tukey, "yes", "no") # Is this effect size a low outlier? Yes/no
  if(all(x==x[1])){ # Special case for when all values are the same
    es_tukey <- x
  } else {
    es_tukey <- ifelse(high_outlier=="no" & low_outlier=="no", x, NA) # Copy values for non-outliers or if all the same, NAs for outliers
    es_tukey <- ifelse(high_outlier=="yes", high_tukey, es_tukey) # If high outlier, substitute high fence
    es_tukey <- ifelse(low_outlier=="yes", low_tukey, es_tukey) # If low outlier, substitute low fence
  }
  return(es_tukey) # Return the effect sizes
}

outer_tukey_val <- function(x, data){
  arguments <- as.list(match.call()) # See: http://www.r-bloggers.com/passing-columns-of-a-dataframe-to-a-function-without-quotes/
  x = eval(arguments$x, data) # Same as above.
  print(round(quantile(x)[1:5], digits=4))
  low_tukey <- quantile(x, names=FALSE)[2] - (3 * (quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2])) # Get low fence (Q1 - 1.5*IQR)
  high_tukey <- quantile(x, names=FALSE)[4] + (3 * (quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2])) # Get high fence (Q3 + 1.5*IQR)
  high_outlier <- ifelse(x >= high_tukey, "yes", "no") # Is this effect size a high outlier? Yes/no
  low_outlier <- ifelse(x <= low_tukey, "yes", "no") # Is this effect size a low outlier? Yes/no
  if(all(x==x[1])){ # Special case for when all values are the same
    es_tukey <- x
  } else {
    es_tukey <- ifelse(high_outlier=="no" & low_outlier=="no", x, NA) # Copy values for non-outliers or if all the same, NAs for outliers
    es_tukey <- ifelse(high_outlier=="yes", max(es_tukey, na.rm=TRUE), es_tukey) # If high outlier, substitute high fence
    es_tukey <- ifelse(low_outlier=="yes", min(es_tukey, na.rm=TRUE), es_tukey) # If low outlier, substitute low fence
  }
  return(es_tukey) # Return the effect sizes
}

outer_tukey_skew <- function(x, data){
  arguments <- as.list(match.call()) # See: http://www.r-bloggers.com/passing-columns-of-a-dataframe-to-a-function-without-quotes/
  x = eval(arguments$x, data) # Same as above.
  low_tukey <- adjboxStats(x, coef=3)$fence[1]  # Get low fence 
  high_tukey <-  adjboxStats(x, coef=3)$fence[2] # Get high fence (Q3 + 1.5*IQR)
  high_outlier <- ifelse(x > high_tukey, "yes", "no") # Is this effect size a high outlier? Yes/no
  low_outlier <- ifelse(x < low_tukey, "yes", "no") # Is this effect size a low outlier? Yes/no
  if(all(x==x[1])){ # Special case for when all values are the same
    es_tukey <- x
  } else {
    es_tukey <- ifelse(high_outlier=="no" & low_outlier=="no", x, NA) # Copy values for non-outliers or if all the same, NAs for outliers
    es_tukey <- ifelse(high_outlier=="yes", high_tukey, es_tukey) # If high outlier, substitute max non-outlier value
    es_tukey <- ifelse(low_outlier=="yes", low_tukey, es_tukey) # If low outlier, substitute min non-outlier value
  }
  return(es_tukey) # Return the effect sizes
}

# REQUIRED: DATA AND MODEL INFORMATION DATAFRAME
# This block of code creates a data frame that contains information about the dataset.
# This information will be used when when we need to do things WITHIN analytic datasets
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

################################################################################################
# Do this once with the pre-test and once with the post-test

dataset <- post

dataset <- pre

##################################################################################################
## Calculate cluster size(s) when missing and the study used cluster-level assignment
# Create the variables we need
dataset$cluster_err_tx <- 0
dataset$cluster_err_ct <- 0
dataset$h10_new <- NA
dataset$h12_new <- NA
dataset$is_cluster <- 0

dataset<- within(dataset,{ # Dummy code the cases that have missing cluster size (and a cluster size is needed)
  cluster_err_tx[( h10 %in% c(999,9999,888,8888) | is.na(h10) ) & h8 %in% c(2,3)] <- 1 #h10 = cluster size, h8 = assignment (2 & 3 are cluster)
  cluster_err_ct[( h12 %in% c(999,9999,888,8888) | is.na(h12) ) & h8 %in% c(2,3)] <- 1 #h12 = cluster size, h8 = assignment (2 & 3 are cluster)
})

# Find the median treatment cluster size within each analytic dataset
agg_cluster_sz_tx <-dataset %>%
  filter( h8 %in% c(2,3) & cluster_err_tx==0 ) %>% # Only interested in records with valid cluster size (no valid exceed 888)
  group_by(dv2_new, prev_strat) %>%
  summarise(med_num_cluster_tx = round(median(h10), digits=0))
agg_cluster_sz_tx <- as.data.frame(agg_cluster_sz_tx)

# Find the median control cluster size within each analytic dataset
agg_cluster_sz_ct <-dataset %>%
  filter( h8 %in% c(2,3) & cluster_err_ct==0 ) %>% # Only interested in records with valid cluster size (no valid exceed 888)
  group_by(dv2_new, prev_strat) %>%
  summarise(med_num_cluster_ct = round(median(h12), digits=0))
agg_cluster_sz_ct <- as.data.frame(agg_cluster_sz_ct)

# Merge the median cluster sizes into the original dataset and remove the unneeded objects
dataset <- merge(dataset, agg_cluster_sz_tx, by=c("dv2_new","prev_strat"), all.x = TRUE) # All x = keep rows (in orig dataset) with no match
dataset <- merge(dataset, agg_cluster_sz_ct, by=c("dv2_new","prev_strat"), all.x = TRUE) # All x = keep rows (in orig dataset) with no match
rm(agg_cluster_sz_tx)
rm(agg_cluster_sz_ct)

# If the cluster size was missing, use the median for that analytic dataset, else use the old cluster size
dataset <- within(dataset,{
  h10_new <- ifelse(cluster_err_tx == 1, med_num_cluster_tx, h10)
  h12_new <- ifelse(cluster_err_ct == 1, med_num_cluster_ct, h12)
})

##################################################################################################
## Cluster correct the sample sizes
dataset$rho <- NA
dataset$eff_es3 <- NA
dataset <- within(dataset,{
  rho <- ifelse(dv2_new == 21, 0.2, 0.1) # ICC = .2 for achievement, .1 for behavioral/attitudinal
  avg_units_tx_cluster <- es1/h10_new
  avg_units_ct_cluster <- es2/h12_new
  tx_cluster_inv <- 1/avg_units_tx_cluster
  ct_cluster_inv <- 1/avg_units_ct_cluster
  cluster_inv_wtavg <- ( tx_cluster_inv * h10_new + ct_cluster_inv * h12_new ) / ( h10_new + h12_new )
  cluster_size_synth <- 1/cluster_inv_wtavg
  total_sample_size <- cluster_size_synth * ( h10_new + h12_new )
  eff_es3 <- ifelse(h8==1, es3,  # If individual assignment, old sample size is effective sample size (h8 has no missing values/9's)
                    total_sample_size / ( 1 + (cluster_size_synth-1) * rho)) # Compute the new effective sample size
})

# test <- dataset[dataset$h8 %in% c(2,3),]
# test <- dataset[dataset$eff_es3>100000,]
# test <- dataset
# test2 <- cbind(test$h8, test$h10, test$h12, test$h10_new, test$h12_new, test$es1, test$es2, test$es3, test$rho, test$avg_units_tx_cluster, test$avg_units_ct_cluster, test$tx_cluster_inv, test$ct_cluster_inv, test$cluster_inv_wtavg, test$cluster_size_synth, test$total_sample_size, test$eff_es3)
# test2 <- as.data.frame(test2)
# colnames(test2) <- c("Assignment", "tx clusters", "ct clusters", "tx clusters NEW", "ct clusters NEW", "TXN", "CTN", "N", "Rho", "Avg Units Tx Cluster", "Avg Units Ct Cluster", "Tx cluster inv", "ct cluster inv", "cluster inv wtavg", "cluster size synth", "total sample size", "eff es3")
# write.csv(test2, file="CheckAgain.csv")
#rm(test)
#rm(test2)
##################################################################################################
## Winsorize the Sample Size Outliers (Inner Fence, Bring to Fence, use the version that accounts for skewness)
dataset$eff_es3_tukey <- NA
for(i in 1:nrow(key)){
  dataset[dataset$dv2_new==key$dv[i] & dataset$prev_strat==key$ps[i], ]$eff_es3_tukey <- inner_tukey_skew(eff_es3, dataset[dataset$dv2_new==key$dv[i] & dataset$prev_strat==key$ps[i], ])
}

# Take a peek
#winsor_check <- dataset %>% group_by(dv2_new, prev_strat) %>% summarise(MINeff_es3_tukey = min(eff_es3_tukey), MAXeff_es3_tukey = max(eff_es3_tukey), MEANeff_es3_tukey = mean(eff_es3_tukey))
#winsor_check <- as.data.frame(winsor_check)
#rm(winsor_check)
##################################################################################################
## Recalculate the intervention and comparison group sample sizes based on our adjustments
dataset$eff_es1_tukey <- NA
dataset$eff_es2_tukey <- NA
dataset <- within(dataset,{
  eff_es1_tukey <- (es1/es3)*eff_es3_tukey
  eff_es2_tukey <- (es2/es3)*eff_es3_tukey 
})

##################################################################################################
## For Odds Ratios with 2x2 frequency tables, re-compute the cell frequencies using the updated sample sizes
# Step 1: Fill in incomplete 2x2 tables
dataset <- within(dataset,{
  xa <- ifelse(is.na(xa) & !is.na(xb), es1-xb, xa) # If xa is missing and xb is not, xa = tx group n - xb. We're using the old sample sizes here - NOT an error.
  xb <- ifelse(is.na(xb) & !is.na(xa), es1-xa, xb)
  xc <- ifelse(is.na(xc) & !is.na(xd), es2-xd, xc)
  xd <- ifelse(is.na(xd) & !is.na(xc), es2-xc, xd)
  xa_pr <- ifelse(is.na(xa_pr) & !is.na(xb_pr), 1-xb_pr, xa_pr) # If xa proportion is missing and xb proportion is not, xa proportion = 1 - xb proportion
  xb_pr <- ifelse(is.na(xb_pr) & !is.na(xa_pr), 1-xa_pr, xb_pr)
  xc_pr <- ifelse(is.na(xc_pr) & !is.na(xd_pr), 1-xd_pr, xc_pr)
  xd_pr <- ifelse(is.na(xd_pr) & !is.na(xc_pr), 1-xc_pr, xd_pr)
})

# Check whether 2x2s wtih missing data are all percents or proportions (if mix, we will have to deal with it via code)
#test <- cbind(dataset$xa_pr, dataset$xb_pr, dataset$xc_pr, dataset$xd_pr)
#test <- as.data.frame(test)
#test$filter <- ifelse(is.na(test$V1) & is.na(test$V2) & is.na(test$V3) & is.na(test$V4), 1, 0) 
#test <- test[test$filter==0,]
#test$filter <- ifelse(!is.na(test$V1) & !is.na(test$V2) & !is.na(test$V3) & !is.na(test$V4), 1, 0) 
#test <- test[test$filter==0,]
#test$sum <- rowSums(cbind(test$V1, test$V2, test$V3, test$V4), na.rm=TRUE)
#test$na <- apply(cbind(test$V1, test$V2, test$V3, test$V4), 1, function (z) sum(is.na(z)))
# All proportions!


# Checking that either of the 2x2 tables are complete DO NOT NEED TO RERUN
#testtest <- cbind(dataset[dataset$es4==2,]$es23,
#                  dataset[dataset$es4==2,]$xa,dataset[dataset$es4==2,]$xb,dataset[dataset$es4==2,]$xc,dataset[dataset$es4==2,]$xd,
#                  dataset[dataset$es4==2,]$xa_pr,dataset[dataset$es4==2,]$xb_pr,dataset[dataset$es4==2,]$xc_pr,dataset[dataset$es4==2,]$xd_pr)
#testtest<-as.data.frame(testtest)
#colnames(testtest)<- c("es23", "xa","xb","xc","xd","xa_pr","xb_pr","xc_pr","xd_pr")
#testtest$complete <-NA
#testtest <- within(testtest, {
#  complete <- ifelse(es23==2 & !is.na(xa) & !is.na(xb) & !is.na(xc) & !is.na(xd), 1,
#                     ifelse(es23==3 & !is.na(xa_pr) & !is.na(xb_pr) & !is.na(xc_pr) & !is.na(xd_pr), 2, 
#                            ifelse(es23==17,3,0))
#  )
#})
#rm(testtest)

# Step 2: Calculate percentage cells for all cases - will be needed to redistribute the updated sample size among cells. Convert any proportions to percentages.
dataset[dataset$es23==2,] <- within(dataset[dataset$es23==2,],{ # Only need to deal with the frequency ORs
  xa_pr <- xa/(xa+xb)*100
  xb_pr <- xb/(xa+xb)*100
  xc_pr <- xc/(xc+xd)*100
  xd_pr <- xd/(xc+xd)*100
})

dataset$convertORPR <- NA
dataset[dataset$es23==2 | dataset$es23==3,] <- within(dataset[dataset$es23==2 | dataset$es23==3,],{
  convertORPR <- ifelse(xa_pr + xb_pr + xc_pr + xd_pr < 2.5, 1, 0)
  xa_pr <- ifelse(convertORPR == 1, xa_pr*100, xa_pr)
  xb_pr <- ifelse(convertORPR == 1, xb_pr*100, xb_pr)
  xc_pr <- ifelse(convertORPR == 1, xc_pr*100, xc_pr)
  xd_pr <- ifelse(convertORPR == 1, xd_pr*100, xd_pr)
})

# Step 3: Use the percentages and updated sample size to calculate the cell counts for all odds ratios (that are not author reported ORs)
dataset$xa_eff <- NA
dataset$xb_eff <- NA
dataset$xc_eff <- NA
dataset$xd_eff <- NA
dataset[dataset$es23==2 | dataset$es23==3,] <- within(dataset[dataset$es23==2 | dataset$es23==3,],{ 
  xa_eff <- (xa_pr/100) * eff_es1_tukey
  xb_eff <- (xb_pr/100) * eff_es1_tukey
  xc_eff <- (xc_pr/100) * eff_es2_tukey
  xd_eff <- (xd_pr/100) * eff_es2_tukey
})

# Note about the OR calculation:
# The OR calculation was checked against the FileMaker output and appears correct. HOWEVER, FileMaker incorrectly computes Odds Ratios for N/count ORs with 
# missing value(s) in the 2x2 table. If I changed the 2x2 cell NAs to 0s, the odds ratios matched, however, at that point they are both incorrect!

# Step 4: Recompute those odds ratios!
dataset$es81_eff <- NA
dataset <- within(dataset,{
  es81_eff <- # OR uncorrected for favored group AND not including author calculated ORs
    ifelse(es23 %in% c(2,3) & xa_eff>0 & xb_eff>0 & xc_eff>0 & xd_eff>0, (xa_eff*xd_eff)/(xb_eff*xc_eff),
    ifelse(es23 %in% c(2,3) & xa_eff==0 | xb_eff==0 | xc_eff==0 | xd_eff==0, ((xa_eff+.5)*(xd_eff+.5))/((xb_eff+.5)*(xc_eff+.5)),
    NA))
})

# Step 5: Change the OR based on favored group
dataset$es81_new <- NA
dataset <- within(dataset, {
  es81_new <- # Syntax taken from Filemaker. "log" in R is natural log.
    ifelse (es23==17 & es4==2 & es17 %in% c(1,2), es81, # If author calculated and OR, just use the originally computed FM effect size (ES60 isn't exported by default)
    ifelse (es23 %in% c(2, 3, 11, 17) & es17==1 & es81_eff > 1, es81_eff,
    ifelse (es23 %in% c(2, 3, 11, 17) & es17==1 & es81_eff < 1, exp(abs(log(es81_eff))),
    ifelse (es23 %in% c(2, 3, 11, 17) & es17==2 & es81_eff < 1, es81_eff,
    ifelse (es23 %in% c(2, 3, 11, 17) & es17==2 & es81_eff > 1, exp(-1*(log( es81_eff ))),
    ifelse (es23 %in% c(2, 3, 11, 17) & es17==3, 1, 
    ifelse (es23 %in% c(2, 3, 11, 17) & es81_eff==1, 1, 
    NA)))))))
})

# Step 6: Recompute SMDs (and change them based on favored group @ the same time)
dataset$es21_new <- NA
# Code creates same results as FM (checked to 4 decimal places) when using the non-corrected sample size variables. (Extremely small differences w/ no rounding). Good to go!
dataset <- within(dataset, {
  es21_new <- # Create new FINAL SMD variable. Will paste this variable back into SPSS dataset.
    # for author-reported SMDs*/ es23=17
    ifelse(es4==1 & es23==17, es21,
    # for SMDs from means and standard deviations
    ifelse(es4==1 & eff_es3_tukey>0 & !is.na(es9) & !is.na(es10) & !is.na(es12) & !is.na(es13) & es17==1, abs(es11/sqrt(((es12*es12)*(eff_es1_tukey-1)+(es13*es13)*(eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    ifelse(eff_es3_tukey>0 & !is.na(es9) & !is.na(es10) & !is.na(es12) & !is.na(es13) & es17==2, 0-abs(es11/sqrt(((es12*es12)*(eff_es1_tukey-1)+(es13*es13)*(eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    # for SMDs from means and variances
    ifelse(es4==1 & eff_es3_tukey>0 & !is.na(es9) & !is.na(es10) & !is.na(es40) & !is.na(es41) & es17==1, abs(es11/sqrt(((es40)*(eff_es1_tukey-1)+(es41)*(eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    ifelse(eff_es3_tukey>0 & !is.na(es9) & !is.na(es10) & !is.na(es40) & !is.na(es41) & es17==2, 0-abs(es11/sqrt(((es40)*(eff_es1_tukey-1)+(es41)*(eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    # for SMDs from means and standard errors
    ifelse(es4==1 & eff_es3_tukey>0 & !is.na(es9) & !is.na(es10) & !is.na(es42) & !is.na(es43) & es17==1, abs(es11/sqrt((((es42*(sqrt(eff_es1_tukey)))*(es42*(sqrt(eff_es1_tukey))))*(eff_es1_tukey-1)+((es43*(sqrt(eff_es2_tukey)))*(es43*(sqrt(eff_es2_tukey))))*(eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    ifelse(eff_es3_tukey>0 & !is.na(es9) & !is.na(es10) & !is.na(es42) & !is.na(es43) & es17==2, 0-abs(es11/sqrt((((es42*(sqrt(eff_es1_tukey)))*(es42*(sqrt(eff_es1_tukey))))*(eff_es1_tukey-1)+((es43*(sqrt(eff_es2_tukey)))*(es43*(sqrt(eff_es2_tukey))))*(eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    # for SMDs from t-tests (es23==5)
    ifelse(es4==1 & eff_es1_tukey>1 & eff_es2_tukey>1 & es23==5 & !is.na(es5) & es17==1, abs(es5*sqrt((eff_es1_tukey+eff_es2_tukey)/(eff_es1_tukey*eff_es2_tukey))),
    ifelse(es4==1 & eff_es1_tukey>1 & eff_es2_tukey>1 & es23==5 & !is.na(es5) & es17==2, 0-abs(es5*sqrt((eff_es1_tukey+eff_es2_tukey)/(eff_es1_tukey*eff_es2_tukey))),
    # for SMDs from F-tests (es23==6)
    ifelse(es4==1 & eff_es1_tukey>1 & eff_es2_tukey>1 & es23==6 & !is.na(es7) & es17==1, abs((sqrt(es7))*sqrt((eff_es1_tukey+eff_es2_tukey)/(eff_es1_tukey*eff_es2_tukey))),
    ifelse(es4==1 & eff_es1_tukey>1 & eff_es2_tukey>1 & es23==6 & !is.na(es7) & es17==2, 0-abs((sqrt(es7))*sqrt((eff_es1_tukey+eff_es2_tukey)/(eff_es1_tukey*eff_es2_tukey))),
    # for SMDs from chi-squares and no marginals or cell info (es23==9)
    ifelse(es4==1 & !is.na(es6) & es23==9 & es17==1, abs(2*sqrt(es6/(eff_es3_tukey-es6))),
    ifelse(es4==1 & !is.na(es6) & es23==9 & es17==2, 0-abs(2*sqrt(es6/(eff_es3_tukey-es6))),
    # for hand computed effect sizes
    ifelse(!is.na(es70) & es17==1, abs(es70), 
    ifelse(!is.na(es70) & es17==2, 0-abs(es70), 
    # for zero effect sizes
    ifelse(es4==1 & es17==3, 0,NA))))))))))))))))
})


##################################################################################################
# Cluster correct the SMDs. August 15 - Decided to cluster correct the sample size only
#dataset[dataset$es4==1 & dataset$h8 %in% c(2,3),] <- within(dataset[dataset$es4==1 & dataset$h8 %in% c(2,3),], {
#  es21_new <- es21_new * (sqrt(1-((2* (cluster_size_synth-1) * rho )/(total_sample_size-2))))
#})



# Recompute the effect sizes and variances using the updated sample sizes
dataset$adjes <- NA
dataset$V <- NA
dataset$W <- NA
dataset[dataset$es4==1,] <- within(dataset[dataset$es4==1,], {
  adjes <- (1-(3/((4*(eff_es3_tukey))-9)))*es21_new
  V <- ((eff_es1_tukey+eff_es2_tukey)/(eff_es1_tukey*eff_es2_tukey))+((adjes**2)/(2*(eff_es1_tukey+eff_es2_tukey)))
  W <- 1/V
})


# Adjust es81_new=0 to a very small number, so we dont' have to deal with -Inf logs. This will spit out an error message if there are no cases. Do not worry about it.
dataset[dataset$es81_new==0 & dataset$es4==2,] <- within(dataset[dataset$es81_new==0 & dataset$es4==2,], {
  es81_new <- .00001
})

dataset$LOR <- NA
dataset$seLOR <- NA
dataset$varLOR <- NA
# Compute LOR (Need to do the +.5 if 0 here as well, else will not compute for some cases)
dataset[dataset$es4==2,] <- within(dataset[dataset$es4==2,], {
  LOR <- log(es81_new)
  seLOR <- ifelse(xa_eff==0 | xb_eff==0 | xc_eff==0 | xd_eff==0, 
                  sqrt( (1/(xa_eff+.5))+(1/(xb_eff+.5))+(1/(xc_eff+.5))+(1/(xd_eff+.5)) ),
                  sqrt( (1/xa_eff)+(1/xb_eff)+(1/xc_eff)+(1/xd_eff) ) 
                  )
  varLOR <- seLOR * seLOR
})

# Check out the problematic OR SE/Vars
# test3 <- cbind(dataset[dataset$es4==2 & is.na(dataset$seLOR),]$universalid, 
#                dataset[dataset$es4==2 & is.na(dataset$seLOR),]$es23,
#                dataset[dataset$es4==2 & is.na(dataset$seLOR),]$seLOR,
#                dataset[dataset$es4==2 & is.na(dataset$seLOR),]$xa,
#                dataset[dataset$es4==2 & is.na(dataset$seLOR),]$xb,
#                dataset[dataset$es4==2 & is.na(dataset$seLOR),]$xc,
#                dataset[dataset$es4==2 & is.na(dataset$seLOR),]$xd,
#                dataset[dataset$es4==2 & is.na(dataset$seLOR),]$xa_eff,
#                dataset[dataset$es4==2 & is.na(dataset$seLOR),]$xb_eff,
#                dataset[dataset$es4==2 & is.na(dataset$seLOR),]$xc_eff,
#                dataset[dataset$es4==2 & is.na(dataset$seLOR),]$xd_eff)
# test3 <- as.data.frame(test3)


# Set ADJES to NA if OR with non-missing LOR. SET OR to NA if SMD with not missing ADJES.
dataset[dataset$es4==2 & !is.na(dataset$LOR),]$adjes <- NA #If OR and there is an OR effect size, set adjes to NA
dataset[dataset$es4==1 & !is.na(dataset$adjes),]$LOR <- NA #If SMD and there is an adjes, set OR to NA
# Remove cases where the effect size can't be computed (we'll handle variacnes later, as we might be able to get around it thorugh conversion). This effects 0 SMDs and 2 ORs
dataset<- dataset[(dataset$es4==2 & !is.na(dataset$LOR)) | (dataset$es4==1 & !is.na(dataset$adjes)),]

# Convert ES types for all cases.
dataset$LORconvert <- NA
dataset$LORconvertV <- NA
dataset$SMDconvert <- NA
dataset$SMDconvertV <- NA
dataset <- within(dataset,{
  LORconvert <- es21_new*(pi/sqrt(3))
  LORconvertV <- V*((pi*pi)/3)
  SMDconvert <- LOR*(sqrt(3)/pi)
  SMDconvertV <- varLOR*(3/(pi*pi))
  # Copy over the converted values (see above where we set to NA)
  varLOR <- ifelse(is.na(LOR), LORconvertV, varLOR)
  LOR <- ifelse(is.na(LOR), LORconvert, LOR)
  V <- ifelse(is.na(adjes), SMDconvertV, V)
  adjes <- ifelse(is.na(adjes), (1-(3/((4*(eff_es3_tukey))-9)))*SMDconvert, adjes)
})

# Standardize EStype use within analytic datasets
dataset$isSMD <- ifelse(dataset$es4==1, 1, 0)
dataset$isLOR <- ifelse(dataset$es4==2, 1, 0)


########### POST TEST ONLY. FIGURE OUT ANALYTIC ES TYPE. PRETEST WILL USE THE POST-TEST DECISION. 
escount <- dataset %>% group_by(dv2_new, prev_strat) %>% summarize(NumSMD=sum(isSMD), NumOR=sum(isLOR))
escount <- as.data.frame(escount)
escount$analyticestype <- ifelse(escount$NumSMD>=escount$NumOR, 1, 2)
escount <- escount[-c(3,4)] #drop columns we don't want to merge (counts)
# Posttest only over #############################################################################

dataset <- merge(dataset, escount, by=c("dv2_new","prev_strat"), all.x = TRUE) # All x = keep rows (in orig dataset) with no match

##################################################################################################
## Convert to OR/SMD as necessary.
dataset$es <- NA
dataset$esvar <- NA
dataset$ORtoSMD <- NA
dataset$SMDtoOR <- NA
dataset$ESsame <- NA
dataset <- within(dataset,{
  es <- ifelse(analyticestype==1, adjes, LOR)
  esvar <- ifelse(analyticestype==1, V, varLOR)
  ORtoSMD <- ifelse(isLOR==1 & analyticestype==1, 1, 0)
  SMDtoOR <- ifelse(isSMD==1 & analyticestype==2, 1, 0)
  ESsame <- ifelse( (isSMD==1 & analyticestype==1) | (isLOR==1 & analyticestype==2), 1, 0)
})

# Drop cases with missing variances (These are all odds ratios - the 10 expected posttest cases)
dataset <- dataset[!is.na(dataset$esvar),]

##################################################################################################
## Winsorize the Effect Size Outliers within analytic datasets(Outer Fence, don't account for skewness)
# REQUIRED: WINSORIZE THE EFFECT SIZES
dataset$estukey1 <- NA
for(i in 1:nrow(key)){
  dataset[dataset$dv2_new==key$dv[i] & dataset$prev_strat==key$ps[i], ]$estukey1 <- outer_tukey(es, dataset[dataset$dv2_new==key$dv[i] & dataset$prev_strat==key$ps[i], ])
}


##################################################################################################
## Re-save dataset as pre/post test
## BE CAREFUL TO DO THIS CORRECTLY!
postFINAL <- dataset
rm(post)

preFINAL <- dataset
rm(pre)



##################################################################################################
##################################################################################################
##################################################################################################
## REPEAT ALL OF THE  ABOVE STEPS WITH THE PRETEST DATASET PRIOR TO MERGING 
## (Except for creating the escount dataframe!)
##################################################################################################
##################################################################################################
##################################################################################################



##################################################################################################
## Merge pre-test and post-test datasets

# Keep only the  variables needed for the merge
preFinalMerge <- subset(preFINAL, select = c(studyid, substudy, source, breakid, subgrp, typevar, varno, estukey1))
colnames(preFinalMerge)[8] <- "estukey1_pre"
dataFinal <- merge(postFINAL, preFinalMerge, by=c("studyid", "substudy", "source", "breakid", "subgrp", "typevar", "varno"), all.x = TRUE) # All x = keep rows (in orig dataset) with no match

##################################################################################################
## Adjust post-test scores using pre-test scores (when appropriate)
dataFinal$es_preadj <- NA
dataFinal$es_preadj <- ifelse(is.na(dataFinal$estukey1_pre) | dataFinal$es50 %in% c(2, 4), dataFinal$estukey1, # If there is no pretest OR already adjusted on pretest, use the posttest
                              dataFinal$estukey1 - dataFinal$estukey1_pre) # Else, subtract the pretest from the posttest

##################################################################################################
## Export file for use in later analyses

save(dataFinal,file="R_ANALYTIC_DATASET2.Rda")
# This file can be used later. Or just move on to the analysis now!


########## EFFECT SIZE PICKER ########## 
## not fully tested yet, but makes sense. Need to run the weighted mean effect size section in the Analysis.R file
## (as some study groups will still have > 1 ES)
## STILL TO DO: 
##### 1) SPOT CHECK OUTPUT IN REAL DATASET, 
##### 2) ADD VARIABLE FOR: 1) ALL OR TO SMD / SMD TO OR, 2) ALL SAME, 3) MIX OF SAME AND CONVERTED

rm(escount)

# Function to find the mode (why?????)
Mode <- function(x) {
  ux <- unique(x) # Get list of unique values
  ux[which.max #Get location of maximum
     (tabulate # Get counts of number of appearances
      (match(x, ux)) # Return vector of matches (list, comparison table)
     )]
}

# Create a construct variable that combines the macro and micro construct (micro construct IDs are not unique across construct families ("macro" constructs))
dataFinal$macro_micro <- dataFinal$dv1 + dataFinal$dv2/100 #  Max is 99, so can divide by 100.

###### Creating a helpful dataframe with useful information
# Get list of the studygroups in each analytic dataset
escount <- dataFinal %>% group_by(dv2_new, prev_strat, studygroup) %>% 
  summarize(ds_num=length(esid))
escount <- as.data.frame(escount)
# The loop will take a few minutes. Could make first half more efficient but then we'd have to merge twice which seems like a lot of effort to me right now
escount$d_inf <- NA
escount$d_con <- NA
escount$ds_inf <- NA
escount$ds_con <- NA
for(i in 1:nrow(escount)){
  # Get the most common informant and construct across all effect sizes within each analytic dataset
  escount[i,5] <- Mode(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat==escount$prev_strat[i],]$dv5) # Informant
  escount[i,6] <- Mode(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat==escount$prev_strat[i],]$macro_micro) # Construct 
  # Get the most common informant and construct across all effect sizes within each studygroup within each analytic dataset
  escount[i,7] <- Mode(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat==escount$prev_strat[i] & dataFinal$studygroup==escount$studygroup[i],]$dv5) # Informant
  escount[i,8] <- Mode(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat==escount$prev_strat[i] & dataFinal$studygroup==escount$studygroup[i],]$macro_micro) # Construct 
}

#### Merge the helpful dataframe into the main dataset (we will delete these variables once we're done with them)
dataFinal <- merge(dataFinal, escount, by=c("dv2_new", "prev_strat", "studygroup"))

####### Set up our decision flags!
####### This is an inefficient way to do this, but all of my other ideas were bad and this works
# This first set of variables isn't too necessary, but they will greatly simply the code I write next
dataFinal$is_single <- ifelse(dataFinal$ds_num==1, 1, 0)
dataFinal$is_adj <- ifelse(!is.na(dataFinal$estukey1_pre) | dataFinal$es50 %in% c(2, 4), # Flag adjusted effect sizes (this might be better placed above / in the main cleaning code)
                           1, 0)# If there is a pretest OR already adjusted on pretest, 1, else = 0
dataFinal$is_d_inf <- ifelse(dataFinal$dv5 == dataFinal$d_inf, 1, 0) # Does the informant match the most popular for the analytic dataset?
dataFinal$is_d_con <- ifelse(dataFinal$macro_micro == dataFinal$d_con, 1, 0) # Does the construct match the most popular for the analytic dataset?
dataFinal$is_ds_inf <- ifelse(dataFinal$dv5 == dataFinal$ds_inf, 1, 0) # Does the informant match the most popular for the analytic dataset and studygroup?
dataFinal$is_ds_con <- ifelse(dataFinal$macro_micro == dataFinal$ds_con, 1, 0) # Does the construct match the most popular for the analytic dataset and studygroup?

dataFinal$keepFinal <- 0

# This loop will take a VERY VERY LONG TIME to run. (Loops are pretty slow in R in general, but if you're running a loop, you're probably already Doing It Wrong (inefficiently)). Do as I say not as I do.
for(i in 1:nrow(escount)){
  if (nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 1,]) > 0) {
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 1,]$keepFinal <- 1 # First priority, keep if singleton
  } 
  else if (nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 1 & dataFinal$is_d_inf == 1 & dataFinal$is_d_con == 1,]) > 0) {
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 1 & dataFinal$is_d_inf == 1 & dataFinal$is_d_con == 1,]$keepFinal <- 1 # Else, Keep if adjusted and match common analytic dataset informant and construct
  } 
  else if(nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 1 & dataFinal$is_d_inf == 1 & dataFinal$is_d_con == 0,]) > 0){
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 1 & dataFinal$is_d_inf == 1 & dataFinal$is_d_con == 0,]$keepFinal <- 1 # Else, Keep if adjusted and match common analytic dataset informant 
  } 
  else if (nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 1 & dataFinal$is_ds_inf == 1 & dataFinal$is_ds_con == 1,]) > 0) {
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 1 & dataFinal$is_ds_inf == 1 & dataFinal$is_ds_con == 1,]$keepFinal <- 1 # Else, Keep if adjusted and match common studygroup (w/i analytic dataset) informant and construct
  } 
  else if(nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 1 & dataFinal$is_ds_inf == 1 & dataFinal$is_ds_con == 0,]) > 0){
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 1 & dataFinal$is_ds_inf == 1 & dataFinal$is_ds_con == 0,]$keepFinal <- 1 # Else, Keep if adjusted and match common studygroup (w/i analytic dataset) informant 
  } 
  else if (nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 0 & dataFinal$is_d_inf == 1 & dataFinal$is_d_con == 1,]) > 0) {
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 0 & dataFinal$is_d_inf == 1 & dataFinal$is_d_con == 1,]$keepFinal <- 1 # Else, Keep if unadjusted and match common analytic dataset informant and construct
  } 
  else if(nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 0 & dataFinal$is_d_inf == 1 & dataFinal$is_d_con == 0,]) > 0){
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 0 & dataFinal$is_d_inf == 1 & dataFinal$is_d_con == 0,]$keepFinal <- 1 # Else, Keep if unadjusted and match common analytic dataset informant 
  } 
  else if (nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 0 & dataFinal$is_ds_inf == 1 & dataFinal$is_ds_con == 1,]) > 0) {
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 0 & dataFinal$is_ds_inf == 1 & dataFinal$is_ds_con == 1,]$keepFinal <- 1 # Else, Keep if unadjusted and match common studygroup (w/i analytic dataset) informant and construct
  } 
  else if(nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 0 & dataFinal$is_ds_inf == 1 & dataFinal$is_ds_con == 0,]) > 0){
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & dataFinal$prev_strat == escount$prev_strat[i] & dataFinal$studygroup == escount$studygroup[i] & dataFinal$is_single == 0 & dataFinal$is_adj == 0 & dataFinal$is_ds_inf == 1 & dataFinal$is_ds_con == 0,]$keepFinal <- 1 # Else, Keep if unadjusted and match common studygroup (w/i analytic dataset) informant 
  } 
}


##### KEEP THE EFFECT SIZES FLAGGED AS KEEPFINAL=1
dataFinal <- dataFinal[dataFinal$keepFinal==1,]

##### Find weighted mean effect sizes and variances for studies with > 1 ES
agg_wtmeans <-dataFinal %>%
  group_by(studygroup, dv2_new, prev_strat) %>%
  summarise(wmean = weighted.mean(es_preadj, eff_es3_tukey), #Change back to tukey 2 when done with winsorizing experiment! (6/23)
            wvar = weighted.mean(esvar, eff_es3_tukey))
agg_wtmeans <- as.data.frame(agg_wtmeans)
dataFinal <- merge(dataFinal, agg_wtmeans, by=c("studygroup","dv2_new","prev_strat"))
# Below is really an overly complated way to select one of the rows - it doesn't need to be random because they should all be the same. But it works.
dataFinal$rand_wt <- sample(nrow(dataFinal), size=nrow(dataFinal), replace=FALSE) # Add a column of random numbers to the dataset. replace=FALSE, so no duplicate values.
agg_wtmeans <- aggregate(dataFinal$rand_wt, by=list(dataFinal$studygroup, dataFinal$dv2_new, dataFinal$prev_strat), min)# Smallest randomly generated number within each studygroup+dv+ps combo.
colnames(agg_wtmeans)[c(1:4)] <- c("studygroup","dv2_new", "prev_strat","min_rand_wt")
dataFinal <- merge(dataFinal, agg_wtmeans, by=c("studygroup","dv2_new","prev_strat")) # Merge the smallest randomly generated number (by studygroup+dv+ps) back into the dataset
dataFinal$select_wt <- ifelse(dataFinal$rand_wt == dataFinal$min_rand_wt, 1, 0) # We'll keep the 1 in our analysis and drop the 0
dataFinal <- subset(dataFinal, select = -c(rand_wt, min_rand_wt))
dataFinal <- dataFinal[dataFinal$select_wt == 1,]
rm(agg_wtmeans)

# Final variable names: effect size = wmean, variance = wvar

save(dataFinal,file="R_ANALYTIC_DATASETv2.Rda")
# This file can be used later. Or just move on to the analysis now!

dataexp<-subset(dataFinal, select=c(studygroup, esid, universalid, varno, dv5, prev_strat, dv2_new, decade, 
                                    h2, h8h9, h33, h35_d1, h14i, routine_eval_dev, h25alt_d1, h22i, 
                                    h21_d6, h26h85, h27h86, h37r_d1, h38i, h68i, h73i, h70i, h40_d1, h41_d1, 
                                    h57_d1, h59alt_d1, h64i, s5, s8i, SESFinal, s2i, es3, h68_count, h73_count, 
                                    h70r_count, h42, ageavg, agerange, safe_level, sc100_d101, sc100_d102, sc100_d103,
                                    sc100_d104,sc100_d105,sc100_d106,sc100_d107,sc100_d108,sc100_d109,sc100_d110,
                                    sc100_d111,sc100_d112,sc100_d113,sc100_d114,sc100_d115,sc100_d116,sc100_d117,sc100_d118,
                                    sc100_d119,sc100_d120,sc100_d121,sc100_d122,sc100_d123,sc100_d124,sc100_d125,
                                    sc100_d126,sc100_d127,sc100_d128,sc100_d129,sc100_d130,sc100_d131,sc100_d132,sc100_d133,sc100_d134,
                                    sc100_d135,sc100_d137,sc100_d139))
write.csv(dataexp, file="datafinal.csv")


####################OLD CODE GRAVEYARD############################################################
