Meta-Analysis: Exploring Treatment Effectiveness
================

The first half of this R code was used to prepare a data file for analysis in a large meta-analysis project. Some preliminary analyses follow in the second half.

Data Processing
===============

Data cleaning was performed in SPSS and the cleaned dataset was imported into R. Information about the file and processes were generally stored outside of the R code (e.g, SPSS data dictionaries, Word documents). At the time (circa 2013), we were beginning the transition to R to take advantage of the `metafor` package.

There are things that I would approach differently were I writing this code in 2019, including:

1.  Less reliance on `for` loops.
2.  Increase the use of functions, as a means of self-documentation and increasing readability and maintainability. (e.g. Even with comments, the effect size calculation portions of this program are very hard to read and could be very difficult to debug. This seems like a good place to separate the calculation-intensive code from the program-flow related code)
3.  Overall, write the code in a way that allows it to be run all at once in a reproducible manner. There are portions of the code that must be run twice, others that should be run once, some that should be skipped, etc. Some of this is due to the nature of the stage of the project (exploratory, still trying to decide what models to run, etc.), but the majority should be completely avoided.

### A note on the dataset:

Each row of the dataset contained information about one effect size. Information about study, participant, and outcome measure characteristics had previously been merged (in SPSS) into the effect-size level dataset.

Brief data dictionaries will be displayed before code chunks as new variables are introduced.

Import data
-----------

``` r
# Initialize environment -----

# Load packages
library(foreign)
library(metafor)
library(dplyr)
library(robustbase)

# Disable scientific notation in output
options(scipen = 999)

# Set random number seed for reproducibility
set.seed(62784)

# Import SPSS data files -----

# Pretest data
pre  <- read.spss("Pretest.sav",use.value.labels=F,to.data.frame=T)

# Posttest data
post <- read.spss("Posttest.sav",use.value.labels = F,to.data.frame = T)
```

Load custom winsorizing functions
---------------------------------

``` r
# Functions -----

# inner_tukey -----
# This function will repace outlier values with Tukey's fence, 
# where Tukey's fence is calculated using 1.5*IQR
inner_tukey <- function(x, data){
  arguments <- as.list(match.call()) # Get user arguments as a list
  x = eval(arguments$x, data) # Recreates "data$x" from argument list
  # Print five number summary to console
  print(round(quantile(x)[1:5], digits=4))
  # Get low fence (Q1 - 1.5*IQR)
  low_tukey <- quantile(x, names=FALSE)[2]-(1.5*(quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2])) 
  # Get high fence (Q3 + 1.5*IQR)
  high_tukey <- quantile(x, names=FALSE)[4]+(1.5*(quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2])) 
  # Flag outliers based on high and low fences (yes if outlier)
  high_outlier <- ifelse(x >= high_tukey, "yes", "no")
  low_outlier <- ifelse(x <= low_tukey, "yes", "no") 
  if(all(x==x[1])){ # If all values in vector are identical ...
    # Return the original vector.
    es_tukey <- x
  } else {# Else...
    # Keep non-outlier values values (temporary NA for outliers)
    es_tukey <- ifelse(high_outlier=="no" & low_outlier=="no", x, NA)
    # If high outlier, substitute value with high fence
    es_tukey <- ifelse(high_outlier=="yes", high_tukey, es_tukey) 
    # If low outier, substitute value with low fence
    es_tukey <- ifelse(low_outlier=="yes", low_tukey, es_tukey) 
  }
  # Return the winsorized vector
  return(es_tukey) # Return the effect sizes
}

# outer_tukey -----
# This function will repace outlier values with Tukey's fence, 
# where Tukey's fence is calculated using 3*IQR
outer_tukey <- function(x, data){
  arguments <- as.list(match.call())
  x = eval(arguments$x, data) 
  # Print five number summary to console
  print(round(quantile(x)[1:5], digits=4))
  # Get low fence (Q1 - 3*IQR)
  low_tukey <- quantile(x, names=FALSE)[2]-(3*(quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2]))
  # Get high fence (Q3 + 3*IQR)
  high_tukey <- quantile(x, names=FALSE)[4]+(3*(quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2]))
  # Flag outliers based on high and low fences (yes if outlier)
  high_outlier <- ifelse(x >= high_tukey, "yes", "no")
  low_outlier <- ifelse(x <= low_tukey, "yes", "no")
  if(all(x==x[1])){ # If all values in vector are identical ...
    # Return the original vector.
    es_tukey <- x
  } else { # Else...
    # Keep non-outlier values values (temporary NA for outliers)
    es_tukey <- ifelse(high_outlier=="no" & low_outlier=="no", x, NA) 
    # If high outlier, substitute value with high fence
    es_tukey <- ifelse(high_outlier=="yes", high_tukey, es_tukey) 
    # If low outier, substitute value with low fence
    es_tukey <- ifelse(low_outlier=="yes", low_tukey, es_tukey)
  }
  # Return the winsorized vector
  return(es_tukey) 
}

# outer_tukey_val -----
# This function will repace outlier values with the most extreme non-outlier, 
# where Tukey's fence is calculated using 3*IQR
outer_tukey_val <- function(x, data){
  arguments <- as.list(match.call())
  x = eval(arguments$x, data) 
  # Print five number summary to console
  print(round(quantile(x)[1:5], digits=4))
  # Get low fence (Q1 - 3*IQR)
  low_tukey <- quantile(x, names=FALSE)[2]-(3*(quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2]))
  # Get high fence (Q3 + 3*IQR)
  high_tukey <- quantile(x, names=FALSE)[4]+(3*(quantile(x, names=FALSE)[4]-quantile(x, names=FALSE)[2])) 
  # Flag outliers based on high and low fences (yes if outlier)
  high_outlier <- ifelse(x >= high_tukey, "yes", "no")
  low_outlier <- ifelse(x <= low_tukey, "yes", "no") 
  if(all(x==x[1])){ # If all values in vector are identical ...
    # Return the original vector.
    es_tukey <- x
  } else { # Else...
    # Keep non-outlier values values (temporary NA for outliers)
    es_tukey <- ifelse(high_outlier=="no" & low_outlier=="no", x, NA) 
    # If high outlier, substitute value with highest non-outlier in vector
    es_tukey <- ifelse(high_outlier=="yes", max(es_tukey, na.rm=TRUE), es_tukey)
    # If low outlier, substitute value with lowest non-outlier in vector
    es_tukey <- ifelse(low_outlier=="yes", min(es_tukey, na.rm=TRUE), es_tukey)
  }
  # Return the winsorized vector
  return(es_tukey)
}

# inner_tukey_skew -----
# This function will repace outlier values with the skewness-adjusted fence, 
# using a coefficient of 1.5
inner_tukey_skew <- function(x, data){
  arguments <- as.list(match.call()) 
  x = eval(arguments$x, data)
  # Get skewness-adjusted fences 
  low_tukey <- adjboxStats(x, coef=1.5)$fence[1]  # low fence 
  high_tukey <-  adjboxStats(x, coef=1.5)$fence[2] # high fence
  # Flag outliers based on high and low fences (yes if outlier)
  high_outlier <- ifelse(x >= high_tukey, "yes", "no") 
  low_outlier <- ifelse(x <= low_tukey, "yes", "no") 
  if(all(x==x[1])){ # If all values in vector are identical ...
    # Return the original vector.
    es_tukey <- x
  } else { # Else...
    # Keep non-outlier values (temporary NA for outliers)
    es_tukey <- ifelse(high_outlier=="no" & low_outlier=="no", x, NA) 
    # If high outlier, substitute value with high fence
    es_tukey <- ifelse(high_outlier=="yes", high_tukey, es_tukey)
    # If low outier, substitute value with low fence
    es_tukey <- ifelse(low_outlier=="yes", low_tukey, es_tukey)
  }
  # Return the winsorized vector
  return(es_tukey)
}

# outer_tukey_skew -----
# This function will repace outlier values with the skewness-adjusted fence, 
# using a coefficient of 3
outer_tukey_skew <- function(x, data){
  arguments <- as.list(match.call()) 
  x = eval(arguments$x, data) 
  # Get skewness-adjusted fences
  low_tukey <- adjboxStats(x, coef=3)$fence[1]  
  high_tukey <-  adjboxStats(x, coef=3)$fence[2] 
  # Flag outliers based on high and low fences (yes if outlier)
  high_outlier <- ifelse(x > high_tukey, "yes", "no") 
  low_outlier <- ifelse(x < low_tukey, "yes", "no")
  if(all(x==x[1])){ # If all values in vector are identical ...
    # Return the original vector.
    es_tukey <- x
  } else { # Else...
    # Keep non-outlier values (temporary NA for outliers)
    es_tukey <- ifelse(high_outlier=="no" & low_outlier=="no", x, NA)
    # If high outlier, substitute value with high fence
    es_tukey <- ifelse(high_outlier=="yes", high_tukey, es_tukey)
    # If low outier, substitute value with low fence
    es_tukey <- ifelse(low_outlier=="yes", low_tukey, es_tukey)
  }
  # Return the winsorized vector
  return(es_tukey)
}
```

Analytic dataset key
--------------------

In the context of this specific project, the dataset is better thought of as a collection of smaller analytic datasets (with some rows possibly belonging to no analytic datasets). Ultimately, analyses will be performed separately on each analytic dataset.

Analytic datasets were defined by unique combinations of the following variables:

-   `dv2_new` Outcome domain (e.g. school performance, self-esteem)
-   `prev_strat` Treatment program prevention strategy (primary prevention or secondary prevention)

``` r
# DATA AND MODEL INFORMATION DATAFRAME -----
# This creates a data frame that contains information about the dataset.
# DO NOT CHANGE unless to add/remove analyses or change dv2_new/prev_strat coding.
key <- data.frame(model=rep("",38), dv=rep("", 38), ps=rep("", 38), es=rep("", 38), 
                  stringsAsFactors=FALSE) 
n <- 0
# Columns:
# Model Output Name (semi-descriptive), dv2_new, 
# prev_strat (1=universal, 2=indicated), effect size type (SMD, OR)
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
```

Begin pseudo-loop
-----------------

The following code was run twice - once on the pretest dataset and once on the posttest dataset - and handles samples size and effect size corrections and conversions. This "pseudo-loop" will end once it's time to merge the corrected pretest effect sizes into the posttest dataset, at which point the posttest effect sizes can be adjusted based on their pretest values.

``` r
# Perform code below twice - once with the pre-test and once with the post-test

dataset <- post

dataset <- pre
```

Impute cluster information for cluster-randomized trials with missing cluster data
----------------------------------------------------------------------------------

Some studies in the dataset were cluster randomized, but not all cluster randomized studies reported the number of clusters present in the analytic sample. In these cases, the median number of clusters for the analytic dataset was used, where the median was computed separately for each analytic dataset. The median was also calculated separately for the intervention and comparison groups.

<table style="width:58%;">
<colgroup>
<col width="19%" />
<col width="38%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>h10</td>
<td>number of tx clusters</td>
</tr>
<tr class="even">
<td>h12</td>
<td>number of ct clusters</td>
</tr>
<tr class="odd">
<td>h8</td>
<td>group assignment method (1: individual assignment; 2,3: cluster assignment)</td>
</tr>
</tbody>
</table>

``` r
# Calculate the number of clusters in each group -----
# Do this when: 
# > 1) # of clusters is  missing, and 
# > 2) the study used cluster-level assignment

# New variables:
dataset$cluster_err_tx <- 0 # 1 if cluster randomized and # of clusters is missing
dataset$cluster_err_ct <- 0
dataset$h10_new <- NA # Number of clusters, including "imputed" data (from medians)
dataset$h12_new <- NA

# Flag effect sizes from cluster-randomized studies with missing 'number of clusters'  
# > Note: 999,9999,888,and 8888 all mean missing/invalid
dataset<- within(dataset,{ 
  cluster_err_tx[( h10 %in% c(999,9999,888,8888) | is.na(h10) ) & h8 %in% c(2,3)] <- 1 
  cluster_err_ct[( h12 %in% c(999,9999,888,8888) | is.na(h12) ) & h8 %in% c(2,3)] <- 1
})

# Find the median number of clusters for treatment groups within each analytic dataset
agg_cluster_sz_tx <-dataset %>%
  # Keep cases from cluster randomized trials that do not have a missing # clusters
  filter(h8 %in% c(2,3) & cluster_err_tx==0) %>% 
  # Group by outcome domain and prevention strategy (analytic dataset)
  group_by(dv2_new, prev_strat) %>%
  # For each outcome domain and prevention strategy combo, 
  # return the median treatment # of clusters
  summarise(med_num_cluster_tx = round(median(h10), digits=0))
agg_cluster_sz_tx <- as.data.frame(agg_cluster_sz_tx)

# Same as the above code, but for the comparison group
agg_cluster_sz_ct <-dataset %>%
  filter(h8 %in% c(2,3) & cluster_err_ct==0) %>%
  group_by(dv2_new, prev_strat) %>%
  summarise(med_num_cluster_ct = round(median(h12), digits=0))
agg_cluster_sz_ct <- as.data.frame(agg_cluster_sz_ct)

# Merge the median # of clusters back into the original dataset
# > Merging by outcome domain and prevention strategy (analytic dataset)
# > all.x = TRUE: keep rows (in orig dataset) with no match
dataset <- merge(dataset, agg_cluster_sz_tx, by=c("dv2_new","prev_strat"), all.x = TRUE) 
dataset <- merge(dataset, agg_cluster_sz_ct, by=c("dv2_new","prev_strat"), all.x = TRUE) 

# Remove objects containing the summary data (medians)
rm(agg_cluster_sz_tx)
rm(agg_cluster_sz_ct)

# If the # of clusters was missing, use the median *for that analytic dataset*,
# else, keep old data
dataset <- within(dataset,{
  h10_new <- ifelse(cluster_err_tx == 1, med_num_cluster_tx, h10)
  h12_new <- ifelse(cluster_err_ct == 1, med_num_cluster_ct, h12)
})
```

Cluster correct sample sizes
----------------------------

| Variable | Description       |
|----------|-------------------|
| es1      | tx sample size    |
| es2      | ct sample size    |
| es3      | total sample size |

``` r
# Cluster correct the sample sizes -----

# New variables:
dataset$rho <- NA # ICC
dataset$eff_es3 <- NA # Cluster corrected total n

dataset <- within(dataset,{
  # If DV is an achievement outcome, use ICC of 0.2, 
  # else use 0.1 (e.g., behavioral and attitudinal outcomes)
  rho <- ifelse(dv2_new == 21, 0.2, 0.1)
  # Calculate avg number of units per cluster, for tx and ct groups
  # > avg units = sample size / # clusters
  avg_units_tx_cluster <- es1/h10_new
  avg_units_ct_cluster <- es2/h12_new
  # Take the inverse of the avg number of units per cluster
  tx_cluster_inv <- 1/avg_units_tx_cluster
  ct_cluster_inv <- 1/avg_units_ct_cluster
  # Calculate the weighted average of the inverse
  cluster_inv_wtavg <- ( tx_cluster_inv * h10_new + ct_cluster_inv * h12_new ) / ( h10_new + h12_new )
  # Take the inverse of theat weighted average
  # This creates the updated average number of units per cluster
  cluster_size_synth <- 1/cluster_inv_wtavg
  # Multiply by the total number of clusters to get the updated sample size
  total_sample_size <- cluster_size_synth * ( h10_new + h12_new )
  # If individual assignment, keep the old sample size
  # Else, compute the cluster corrected sample size
  eff_es3 <- ifelse(h8==1, es3,  
                    total_sample_size / ( 1 + (cluster_size_synth-1) * rho)) 
})
```

Winsorize cluster-corrected sample size outliers
------------------------------------------------

We discussed many ways to handle outlier winsorization. For sample size outliers, we decided to bring outliers to the inner fence while adjusting for skewness.

``` r
# Winsorize the cluster-corrected sample size outliers -----

# Create the variables we will use:
dataset$eff_es3_tukey <- NA # Winsorized cluster corrected total n

# Perform winsorization separately within each analytic dataset
# > Each row of the 'key' object provides information about the analytic datasets of interest 
# > that can be used in filtering (outcome (dv) and prevention strategy (ps))

for(i in 1:nrow(key)){ 
  dataset[dataset$dv2_new==key$dv[i] & dataset$prev_strat==key$ps[i], ]$eff_es3_tukey <- inner_tukey_skew(
    eff_es3, dataset[dataset$dv2_new==key$dv[i] & dataset$prev_strat==key$ps[i], ]
    )
}
```

Calculate intervention and comparison sample sizes from the winsorized total sample size
----------------------------------------------------------------------------------------

``` r
# Recalculate group sample sizes using cluster corrected and winsorized total n  -----
# > These will be calculated based on the original assignment ratios

# New variables:
dataset$eff_es1_tukey <- NA # Winsorized cluster corrected treatment n
dataset$eff_es2_tukey <- NA # Winsorized cluster corrected comparison n

# Create winsorized sample size for each group
# >> These group n's are calculated with the cluster adjusted N
# >> (original group n / original total n) * winsorized total N
dataset <- within(dataset,{
  eff_es1_tukey <- (es1/es3)*eff_es3_tukey
  eff_es2_tukey <- (es2/es3)*eff_es3_tukey 
})
```

Calculate cell frequencies for 2x2 contingency tables using the winsorized sample size
--------------------------------------------------------------------------------------

In this dataset, some 2x2 tables contained count data and others contained proportion data.

| Group | OBSERVED | NOT OBS. |
|-------|----------|----------|
| TX    | A        | B        |
| CT    | C        | D        |

<table style="width:58%;">
<colgroup>
<col width="19%" />
<col width="38%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>xa</td>
<td>observed, tx (n)</td>
</tr>
<tr class="even">
<td>xb</td>
<td>not observed, tx (n)</td>
</tr>
<tr class="odd">
<td>xc</td>
<td>observed, ct (n)</td>
</tr>
<tr class="even">
<td>xd</td>
<td>not observed, ct (n)</td>
</tr>
<tr class="odd">
<td>xa_pr</td>
<td>observed, tx (proportion)</td>
</tr>
<tr class="even">
<td>xb_pr</td>
<td>not observed, tx (proportion)</td>
</tr>
<tr class="odd">
<td>xc_pr</td>
<td>observed, ct (proportion)</td>
</tr>
<tr class="even">
<td>xd_pr</td>
<td>not observed, ct (proportion)</td>
</tr>
<tr class="odd">
<td>es4</td>
<td>Continuous or dichotomous ES? (1: Standardized Mean Difference (continuous); 2: Odds Ratio (dichotomous))</td>
</tr>
<tr class="even">
<td>es23</td>
<td>Data used to calculate ES (2: 2x2 table, frequencies; 3: 2x2 table, proportions/percents)</td>
</tr>
</tbody>
</table>

``` r
# For Odds Ratios with 2x2 frequency tables, re-compute the cell frequencies -----
# > This is similar to what was done above to recreate group sample sizes

# Step 1: Fill in incomplete 2x2 tables -----
# > We're purposefully using the old sample sizes here to complete the *original* tables
dataset <- within(dataset,{
  # If xa is missing and xb is not, xa = ORIGINAL tx group n - xb, etc.
  xa <- ifelse(is.na(xa) & !is.na(xb), es1-xb, xa) 
  xb <- ifelse(is.na(xb) & !is.na(xa), es1-xa, xb)
  xc <- ifelse(is.na(xc) & !is.na(xd), es2-xd, xc)
  xd <- ifelse(is.na(xd) & !is.na(xc), es2-xc, xd)
  # If xa proportion is missing and xb proportion is not, xa_pr = 1 - xb_pr, etc.
  xa_pr <- ifelse(is.na(xa_pr) & !is.na(xb_pr), 1-xb_pr, xa_pr) 
  xb_pr <- ifelse(is.na(xb_pr) & !is.na(xa_pr), 1-xa_pr, xb_pr)
  xc_pr <- ifelse(is.na(xc_pr) & !is.na(xd_pr), 1-xd_pr, xc_pr)
  xd_pr <- ifelse(is.na(xd_pr) & !is.na(xc_pr), 1-xc_pr, xd_pr)
})

# Step 2: Calculate percent cells  -----
# > Will be needed to redistribute the updated sample size among cells.
# > Convert any proportions to percentages.
# > Saving them as percentages temporarily differentiates them from the proportions that
# >> were originally stored in the variables. They will be converted to equivalent
# >> scales later.

# Convert frequency ORs to percents
# >> Use ratio of cell n to group n to get percents
dataset[dataset$es23==2,] <- within(dataset[dataset$es23==2,],{
  xa_pr <- xa/(xa+xb)*100
  xb_pr <- xb/(xa+xb)*100
  xc_pr <- xc/(xc+xd)*100
  xd_pr <- xd/(xc+xd)*100
})

# Handle proportion cells
# >> Convert proprotions to percents

# New variables:
dataset$convertORPR <- NA # Flags an OR 2x2 table that was originally proportional
# (i.e., not converted to percents from frequencies)

dataset[dataset$es23==2 | dataset$es23==3,] <- within(
  dataset[dataset$es23==2 | dataset$es23==3,],{
    # Flag contingency tables coded as proportions 
    # >> (We'll ignore the percents we just calculated using frequency tables)
    convertORPR <- ifelse(xa_pr + xb_pr + xc_pr + xd_pr < 2.5, 1, 0)
    # >> 2.5 allows for some padding for rounding (technically none should sum to >2)
    # Convert flagged entries to percents
    xa_pr <- ifelse(convertORPR == 1, xa_pr*100, xa_pr)
    xb_pr <- ifelse(convertORPR == 1, xb_pr*100, xb_pr)
    xc_pr <- ifelse(convertORPR == 1, xc_pr*100, xc_pr)
    xd_pr <- ifelse(convertORPR == 1, xd_pr*100, xd_pr)
  })

# Step 3: Calculate the 2x2 table cell counts for all odds ratios 
# >> Use the percentages and winsorized cluster corrected sample sizes

# New variables:
dataset$xa_eff <- NA # N for 2x2 table cell A, winsorized and cluster corrected
dataset$xb_eff <- NA
dataset$xc_eff <- NA
dataset$xd_eff <- NA

dataset[dataset$es23==2 | dataset$es23==3,] <- within(
  dataset[dataset$es23==2 | dataset$es23==3,],{ 
    xa_eff <- (xa_pr/100) * eff_es1_tukey
    xb_eff <- (xb_pr/100) * eff_es1_tukey
    xc_eff <- (xc_pr/100) * eff_es2_tukey
    xd_eff <- (xd_pr/100) * eff_es2_tukey
  })
```

Recompute odds ratios using the new cell frequencies
----------------------------------------------------

<table style="width:58%;">
<colgroup>
<col width="19%" />
<col width="38%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>es17</td>
<td>Which group performed &quot;better&quot;? (Treatment, comparison, or eqivalent) [Note: This variable is used to ensure that all effect sizes have consistent directionality - i.e., positive effect sizes always indicate that the intervention group performed better, regardless of whether higher values on the measure indicated a better or worse response]</td>
</tr>
<tr class="even">
<td>es4</td>
<td>Continuous or dichotomous ES? (1: Standardized Mean Difference (continuous); 2: Odds Ratio (dichotomous))</td>
</tr>
<tr class="odd">
<td>es23</td>
<td>Data used to calculate ES (2: 2x2 table, frequencies; 3: 2x2 table, proportions/percents; 11: Hand calculated by coder; 17: Author calculated)</td>
</tr>
</tbody>
</table>

``` r
# Step 4: Recompute the odds ratios

# New variables:
dataset$es81_eff <- NA # Odds ratio calculated from corrected Ns
# This is a temporary variable - it does not yet account for directionality
# It also does not include author- or coder- calculated odds ratios

dataset <- within(dataset,{
  es81_eff <- # OR uncorrected for favored group AND not including author calculated ORs
    ifelse(es23 %in% c(2,3) & xa_eff>0 & xb_eff>0 & xc_eff>0 & xd_eff>0, 
           (xa_eff*xd_eff)/(xb_eff*xc_eff),
    ifelse(es23 %in% c(2,3) & xa_eff==0 | xb_eff==0 | xc_eff==0 | xd_eff==0, 
           ((xa_eff+.5)*(xd_eff+.5))/((xb_eff+.5)*(xc_eff+.5)),
    NA))
})

# Step 5: Recalculate the OR based on favored group

# New variables:
dataset$es81_new <- NA # Odds ratios recomputed and corrected for outcome directionality

dataset <- within(dataset, {
  es81_new <- 
    ifelse (es23==17 & es4==2 & es17 %in% c(1,2), 
            # If author- or hand-calculated OR, just use the original odds ratio
            es81, 
    ifelse (es23 %in% c(2, 3, 11, 17) & es17==1 & es81_eff > 1, 
            es81_eff,
    ifelse (es23 %in% c(2, 3, 11, 17) & es17==1 & es81_eff < 1, 
            exp(abs(log(es81_eff))),
    ifelse (es23 %in% c(2, 3, 11, 17) & es17==2 & es81_eff < 1, 
            es81_eff,
    ifelse (es23 %in% c(2, 3, 11, 17) & es17==2 & es81_eff > 1, 
            exp(-1*(log( es81_eff ))),
    # OR is 1 if both groups are favored equally or if the OR is 1.
    ifelse (es23 %in% c(2, 3, 11, 17) & es17==3, 1, 
    ifelse (es23 %in% c(2, 3, 11, 17) & es81_eff==1, 1, 
    NA)))))))
})
```

Recompute SMDs using the new sample sizes
-----------------------------------------

| Variable | Description |
|----------|-------------|
| \[TODO\] |             |

``` r
# Step 6: Recompute SMDs and change them based on favored group

# New variables:
dataset$es21_new <- NA # SMD recomputed and corrected for outcome directionality

dataset <- within(dataset, {
  es21_new <- 
    # for author-reported SMDs*/ es23=17
    ifelse(es4==1 & es23==17, 
           es21,
    # for SMDs from means and standard deviations
    ifelse(es4==1 & eff_es3_tukey>0 & 
             !is.na(es9) & !is.na(es10) & !is.na(es12) & !is.na(es13) & 
             es17==1,
           abs(es11/sqrt(((es12*es12)*(eff_es1_tukey-1)+(es13*es13)*(eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    ifelse(eff_es3_tukey>0 & 
             !is.na(es9) & !is.na(es10) & !is.na(es12) & !is.na(es13) & 
             es17==2,
           0-abs(es11/sqrt(((es12*es12)*(eff_es1_tukey-1)+(es13*es13)*(eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    # for SMDs from means and variances
    ifelse(es4==1 & eff_es3_tukey>0 & 
             !is.na(es9) & !is.na(es10) & !is.na(es40) & !is.na(es41) 
           & es17==1,
           abs(es11/sqrt(((es40)*(eff_es1_tukey-1)+(es41)*(eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    ifelse(eff_es3_tukey>0 & 
             !is.na(es9) & !is.na(es10) & !is.na(es40) & !is.na(es41) & 
             es17==2,
           0-abs(es11/sqrt(((es40)*(eff_es1_tukey-1)+(es41)*(eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    # for SMDs from means and standard errors
    ifelse(es4==1 & eff_es3_tukey>0 & 
             !is.na(es9) & !is.na(es10) & !is.na(es42) & !is.na(es43) & 
             es17==1,
           abs(es11/sqrt((((es42*(sqrt(eff_es1_tukey)))*(es42*(sqrt(eff_es1_tukey))))*
                            (eff_es1_tukey-1)+((es43*(sqrt(eff_es2_tukey)))*(es43*(sqrt(eff_es2_tukey))))*
                            (eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    ifelse(eff_es3_tukey>0 & 
             !is.na(es9) & !is.na(es10) & !is.na(es42) & !is.na(es43) & 
             es17==2,
           0-abs(es11/sqrt((((es42*(sqrt(eff_es1_tukey)))*(es42*(sqrt(eff_es1_tukey))))*
                              (eff_es1_tukey-1)+((es43*(sqrt(eff_es2_tukey)))*(es43*(sqrt(eff_es2_tukey))))*
                              (eff_es2_tukey-1))/((eff_es1_tukey+eff_es2_tukey)-2))),
    # for SMDs from t-tests (es23==5)
    ifelse(es4==1 & eff_es1_tukey>1 & eff_es2_tukey>1 & es23==5 & !is.na(es5) & es17==1,
           abs(es5*sqrt((eff_es1_tukey+eff_es2_tukey)/(eff_es1_tukey*eff_es2_tukey))),
    ifelse(es4==1 & eff_es1_tukey>1 & eff_es2_tukey>1 & es23==5 & !is.na(es5) & es17==2,
           0-abs(es5*sqrt((eff_es1_tukey+eff_es2_tukey)/(eff_es1_tukey*eff_es2_tukey))),
    # for SMDs from F-tests (es23==6)
    ifelse(es4==1 & eff_es1_tukey>1 & eff_es2_tukey>1 & es23==6 & !is.na(es7) & es17==1,
           abs((sqrt(es7))*sqrt((eff_es1_tukey+eff_es2_tukey)/(eff_es1_tukey*eff_es2_tukey))),
    ifelse(es4==1 & eff_es1_tukey>1 & eff_es2_tukey>1 & es23==6 & !is.na(es7) & es17==2,
           0-abs((sqrt(es7))*sqrt((eff_es1_tukey+eff_es2_tukey)/(eff_es1_tukey*eff_es2_tukey))),
    # for SMDs from chi-squares and no marginals or cell info (es23==9)
    ifelse(es4==1 & !is.na(es6) & es23==9 & es17==1, abs(2*sqrt(es6/(eff_es3_tukey-es6))),
    ifelse(es4==1 & !is.na(es6) & es23==9 & es17==2, 0-abs(2*sqrt(es6/(eff_es3_tukey-es6))),
    # for hand computed effect sizes
    ifelse(!is.na(es70) & es17==1, abs(es70), 
    ifelse(!is.na(es70) & es17==2, 0-abs(es70), 
    # for zero effect sizes
    ifelse(es4==1 & es17==3, 0,NA))))))))))))))))
})
```

Compute Hedges' G, Variance, and Inverse Variance
-------------------------------------------------

``` r
# Compute Hedges' G and Variance using new SMDs

# New variables:
dataset$adjes <- NA # Hedges' G
dataset$V <- NA # Variance of Hedges' G
dataset$W <- NA # Inverse of the variance 

dataset[dataset$es4==1,] <- within(dataset[dataset$es4==1,], {
  adjes <- (1-(3/((4*(eff_es3_tukey))-9)))*es21_new
  V <- ((eff_es1_tukey+eff_es2_tukey)/(eff_es1_tukey*eff_es2_tukey))+((adjes**2)/(2*(eff_es1_tukey+eff_es2_tukey)))
  W <- 1/V
})
```

Compute the Log Odds Ratio and Variance
---------------------------------------

``` r
# Compute LogOR and Variance using new ORs

# Prep:
# Adjust cases in which the OR (es81_new) is 0 to be a very small number, so we dont' get -Inf logs. 
# (This code will show an error message if there are no matches - this warning can be ignored)

dataset[dataset$es81_new==0 & dataset$es4==2,] <- within(dataset[dataset$es81_new==0 & dataset$es4==2,], {
  es81_new <- .00001
})

# New variables
dataset$LOR <- NA # Log odds ratio
dataset$seLOR <- NA # SE of the log OR
dataset$varLOR <- NA # Variance of the log OR

# Compute LOR 
# > Add a small number (0.5) to cells with 0 frequencies when computing the variance
dataset[dataset$es4==2,] <- within(dataset[dataset$es4==2,], {
  LOR <- log(es81_new)
  seLOR <- ifelse(xa_eff==0 | xb_eff==0 | xc_eff==0 | xd_eff==0, 
                  sqrt( (1/(xa_eff+.5))+(1/(xb_eff+.5))+(1/(xc_eff+.5))+(1/(xd_eff+.5)) ),
                  sqrt( (1/xa_eff)+(1/xb_eff)+(1/xc_eff)+(1/xd_eff) ) 
                  )
  varLOR <- seLOR * seLOR
})
```

Convert Effect Sizes
--------------------

This converts all Hedges' G values to log odds ratios, and vice versa (along with their associated variances).

This is done because each analytic dataset contains both continuous and dichotomous effect sizes. In order for them all be analysed together, they need to be converted to the same scale.

``` r
# Convert ES types (G->LOR & LOR->G) for all cases -----

# Prep:
# If ES is dichotomous and there is a LOR effect size, set adjes to NA
dataset[dataset$es4==2 & !is.na(dataset$LOR),]$adjes <- NA 
# If ES is continuous and there is an adjes effect size, set LOR to NA
dataset[dataset$es4==1 & !is.na(dataset$adjes),]$LOR <- NA 
# (This is to clean up any stray calculations that may have occurred in previous steps)

# Remove cases where the effect sizes couldn't be computed 
# (i.e., remove dichotomous ES with missing LOR, and continuous ES with missing adjes)
dataset<- dataset[(dataset$es4==2 & !is.na(dataset$LOR)) | (dataset$es4==1 & !is.na(dataset$adjes)),]
# This effects 0 SMDs and 2 ORs

# New Variables:
# (These are all used to hold temporary values)
dataset$LORconvert <- NA # Log OR form SMD 
dataset$LORconvertV <- NA # Log OR Variance from SMD
dataset$SMDconvert <- NA # SMD from Log OR (Not yet converted to Hedges G)
dataset$SMDconvertV <- NA # SMD Variance from Log OR variance

dataset <- within(dataset,{
  # Continuous to dichotomous
  LORconvert <- es21_new*(pi/sqrt(3)) # Log odds ratio calculated from SMD
  LORconvertV <- V*((pi*pi)/3) # LOR Variance convertd from Hedges' G Variance
  # Dichotomous to continuous
  SMDconvert <- LOR*(sqrt(3)/pi) # Calculate SMD from Log OR
  SMDconvertV <- varLOR*(3/(pi*pi)) # Calculate SMD Variance from Log OR variance
  # Copy the "converted" effect sizes into the effect size variables
  # (Now the "original" and "converted" versions are all in the same variables)
  varLOR <- ifelse(is.na(LOR), LORconvertV, varLOR)
  LOR <- ifelse(is.na(LOR), LORconvert, LOR)
  V <- ifelse(is.na(adjes), SMDconvertV, V)
  adjes <- ifelse(is.na(adjes), 
                  (1-(3/((4*(eff_es3_tukey))-9)))*SMDconvert, adjes) # Converting to G on the fly
})
```

Determine majority OR/SMD within each analytic dataset
======================================================

If the majority of effect sizes within an analytic dataset were originally continuous, Hedges G will be used in the final analysis. If the majority of were dichotomous, the log of the odds ratio will be used.

``` r
# Determine effect size type majority within each analytic dataset -----

# RUN ON POST TEST ONLY ##################################################
# FIGURE OUT ANALYTIC ES TYPE. PRETEST WILL DEFAULT TO THE POST-TEST DECISION
# Standardize EStype use within analytic datasets

# Convert effect size type variable into two dummy variables
dataset$isSMD <- ifelse(dataset$es4==1, 1, 0)
dataset$isLOR <- ifelse(dataset$es4==2, 1, 0)

# Get counts of effect size types within each analytic dataset
escount <- dataset %>% 
  group_by(dv2_new, prev_strat) %>% 
  summarize(NumSMD=sum(isSMD), NumOR=sum(isLOR))
escount <- as.data.frame(escount)

# Determine which effect size type will be considered primary
# 1 = use G, 2 = use LOR
escount$analyticestype <- ifelse(escount$NumSMD>=escount$NumOR, 1, 2)

# Drop columns we don't want to merge
# > We only need merge variables and "analyticestype"
escount <- escount[-c(3,4)] 

# POSTTEST ONLY OVER #####################################################

# Merge the summary dataframe (created using posttest) into data set
dataset <- merge(dataset, escount, 
                 by=c("dv2_new","prev_strat"), all.x = TRUE)

# New Variables:
dataset$es <- NA    # Either G or LogOR, depending on majority es type
dataset$esvar <- NA # Same as above, for variance
dataset$ORtoSMD <- NA # Flag if dichotomous ES converted to continuous
dataset$SMDtoOR <- NA # Flag if continuous ES converted to dichotomous
dataset$ESsame <- NA  # Flag if ES type did not change

dataset <- within(dataset,{
  es <- ifelse(analyticestype==1, adjes, LOR)
  esvar <- ifelse(analyticestype==1, V, varLOR)
  ORtoSMD <- ifelse(isLOR==1 & analyticestype==1, 1, 0)
  SMDtoOR <- ifelse(isSMD==1 & analyticestype==2, 1, 0)
  ESsame <- ifelse( (isSMD==1 & analyticestype==1) | (isLOR==1 & analyticestype==2), 1, 0)
})

# Drop cases with missing variances 
# (If they could have been recovered, it would've been during conversion)
dataset <- dataset[!is.na(dataset$esvar),]
# (These are all odds ratios - 10 cases in posttest dataset)
```

Winsorize effect sizes
----------------------

``` r
# Winsorize the Effect Size Outliers within analytic datasets -----
# > Winsorize to the outer fence, don't account for skewness

# New variables:
dataset$estukey1 <- NA # Winsorized effect size

for(i in 1:nrow(key)){
  dataset[dataset$dv2_new==key$dv[i] & dataset$prev_strat==key$ps[i], ]$estukey1 <- 
    outer_tukey(es, dataset[dataset$dv2_new==key$dv[i] & dataset$prev_strat==key$ps[i], ])
}
```

... End pseudo-loop
-------------------

And so ends the portion of the code that was run twice (once on the pretest and once on the posttest).

``` r
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
```

Merge pre-test and post-test datasets
-------------------------------------

``` r
# Merge pre-test effect sizes into post-test dataset -----

# Keep only the effect size and variables needed for the merge
preFinalMerge <- subset(preFINAL, 
                        select = c(studyid, substudy, source, breakid, subgrp, typevar, varno, 
                                   estukey1))
# Add pre suffix to es variable name 
# > This is the only pretest variable that will be appended to the posttest dataset
colnames(preFinalMerge)[8] <- "estukey1_pre" 

dataFinal <- merge(postFINAL, preFinalMerge, 
                   by=c("studyid", "substudy", "source", "breakid", "subgrp", "typevar", "varno"), 
                   all.x = TRUE) 
```

Adjust posttest effect sizes using pretest effect size
------------------------------------------------------

<table style="width:58%;">
<colgroup>
<col width="19%" />
<col width="38%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>es50</td>
<td>Was this posttest effect size adjusted based on pretest or demographics? (2 and 4 indicate adjustement )</td>
</tr>
</tbody>
</table>

``` r
# Adjust post-test scores using pre-test scores (when appropriate) ------

# New Variables:
dataFinal$es_preadj <- NA # Pretest-adjusted effect size

dataFinal$es_preadj <- ifelse(is.na(dataFinal$estukey1_pre) | dataFinal$es50 %in% c(2, 4), 
                              # If there is no pretest OR ES already adjusted, use the posttest
                              dataFinal$estukey1, 
                              # Else, subtract the pretest from the posttest
                              dataFinal$estukey1 - dataFinal$estukey1_pre) 

# Move on to the analysis (or complete optional filtering)
# In either case: dataFinal is your final dataset from this stage
```

Optional: Filter the dataset before analysis
--------------------------------------------

Each analytic dataset has many effect sizes, and there may be variability in the effect sizes that we we don't want - e.g., are self-reported, and parent- or teacher-reported outcomes really comparable?

This block of code creates a filtered version of the dataset with the goal of decreasing this kind of measurement-related error. The filtering occurrs within each study, within each analytic dataset.

``` r
# Optional: Filter dataset to decrease unwanted measurement-related heterogeneity

# Function to find the mode:
Mode <- function(x) {
  unique_ls <- unique(x) # Get list of unique values
  unique_ls[
    which.max # Get location of most common value
    (tabulate (match(x, unique_ls))) # Get number of appearances of each value
    ]
}

# Create finer-tuned construct variable
# > This contains both the the macro (dv2_new) and micro (dv1) construct 
# > E.g., Internalizing problems is a macro-level construct, anxiety is a micro-level construct
# > (Macro constructs were initially created by grouping together similar micro constructs)
# > Note that micro construct IDs are NOT unique ACROSS macro constructs
dataFinal$macro_micro <- dataFinal$dv1 + dataFinal$dv2/100 #  Max is 99, so can divide by 100.

# Get summary information:
# Number of effect sizes within each study within each dataset
escount <- dataFinal %>% 
  group_by(dv2_new, prev_strat, studygroup) %>% 
  summarize(ds_num=length(esid))
escount <- as.data.frame(escount)
# Each row of 'escount' contains:
# 1) dv2_new, prev_strat, study id, and # of ES in that combo

# New Variables:
escount$d_inf <- NA # Most common informant within each study within each dataset
escount$d_con <- NA # Most common construct (micro+macro) within each dataset
escount$ds_inf <- NA # Most common informant within each study within each dataset
escount$ds_con <- NA # Most common construct wiithin each study within each dataset

# Runs once for each study within each dataset
for(i in 1:nrow(escount)){ 
  # Get the most common informant and construct across all effect sizes within each analytic dataset
  # >> Informant
  escount[i,5] <- Mode(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                                   dataFinal$prev_strat==escount$prev_strat[i],]$dv5)
  # >> Construct (micro/macro combo)
  escount[i,6] <- Mode(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                                   dataFinal$prev_strat==escount$prev_strat[i],]$macro_micro)
  # Get the most common informant and construct across all effect sizes within 
  #each study within each analytic dataset
  # >> Informant
  escount[i,7] <- Mode(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                                   dataFinal$prev_strat==escount$prev_strat[i] & 
                                   dataFinal$studygroup==escount$studygroup[i],]$dv5) 
  # >> Construct  (micro/macro combo)
  escount[i,8] <- Mode(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                                   dataFinal$prev_strat==escount$prev_strat[i] & 
                                   dataFinal$studygroup==escount$studygroup[i],]$macro_micro)
}

# Merge the modes into the main dataset by analytic dataset (dv2+prev) and study
dataFinal <- merge(dataFinal, escount, by=c("dv2_new", "prev_strat", "studygroup"))

# Filter dataset to keep "similar" effect sizes -----

# New Variables:
# Is there only 1 ES in the study in the dataset? (1 if YES)
dataFinal$is_single <- ifelse(dataFinal$ds_num==1, 1, 0) 
# Is this ES adjusted for pretest values? (1 if YES)
dataFinal$is_adj <- ifelse(!is.na(dataFinal$estukey1_pre) | dataFinal$es50 %in% c(2, 4), 
                           1, 0)
# Does the informant match the most popular for the analytic dataset? (1 if YES)
dataFinal$is_d_inf <- ifelse(dataFinal$dv5 == dataFinal$d_inf, 1, 0) 
 # Does the construct match the most popular for the analytic dataset? (1 if YES)
dataFinal$is_d_con <- ifelse(dataFinal$macro_micro == dataFinal$d_con, 1, 0)
# Does the informant match the most popular for the analytic dataset and study? (1 if YES)
dataFinal$is_ds_inf <- ifelse(dataFinal$dv5 == dataFinal$ds_inf, 1, 0) 
# Does the construct match the most popular for the analytic dataset and study? (1 if YES)
dataFinal$is_ds_con <- ifelse(dataFinal$macro_micro == dataFinal$ds_con, 1, 0)   
# keepFinal will be set to 1 if the ES is to be kept
dataFinal$keepFinal <- 0

# Runs once for each study within each dataset
for(i in 1:nrow(escount)){
  # If the only ES is a singleton, keep it
  if (nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                     dataFinal$prev_strat == escount$prev_strat[i] & 
                     dataFinal$studygroup == escount$studygroup[i] & 
                     dataFinal$is_single == 1,]) > 0) {
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                dataFinal$prev_strat == escount$prev_strat[i] & 
                dataFinal$studygroup == escount$studygroup[i] & 
                dataFinal$is_single == 1,]$keepFinal <- 1 
  } 
  # Else, if any ES are adjusted and match the most
  # common analytic dataset informant and micro construct, keep them
  else if (nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                          dataFinal$prev_strat == escount$prev_strat[i] & 
                          dataFinal$studygroup == escount$studygroup[i] & 
                          dataFinal$is_single == 0 & 
                          dataFinal$is_adj == 1 & 
                          dataFinal$is_d_inf == 1 & 
                          dataFinal$is_d_con == 1,]) > 0) {
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                dataFinal$prev_strat == escount$prev_strat[i] & 
                dataFinal$studygroup == escount$studygroup[i] & 
                dataFinal$is_single == 0 & 
                dataFinal$is_adj == 1 & 
                dataFinal$is_d_inf == 1 & 
                dataFinal$is_d_con == 1,]$keepFinal <- 1 
  } 
  # Else, if any ES are adjusted and match the most
  # common analytic dataset informant, keep them
  else if(nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                         dataFinal$prev_strat == escount$prev_strat[i] & 
                         dataFinal$studygroup == escount$studygroup[i] & 
                         dataFinal$is_single == 0 & 
                         dataFinal$is_adj == 1 & 
                         dataFinal$is_d_inf == 1 & 
                         dataFinal$is_d_con == 0,]) > 0){
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                dataFinal$prev_strat == escount$prev_strat[i] & 
                dataFinal$studygroup == escount$studygroup[i] & 
                dataFinal$is_single == 0 & 
                dataFinal$is_adj == 1 & 
                dataFinal$is_d_inf == 1 & 
                dataFinal$is_d_con == 0,]$keepFinal <- 1 
  } 
  # Else, if any ES are adjusted and match the most
  # common study-level (within each analytic dataset) informant and micro construct, keep them
  else if (nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                          dataFinal$prev_strat == escount$prev_strat[i] & 
                          dataFinal$studygroup == escount$studygroup[i] & 
                          dataFinal$is_single == 0 & 
                          dataFinal$is_adj == 1 & 
                          dataFinal$is_ds_inf == 1 & 
                          dataFinal$is_ds_con == 1,]) > 0) {
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                dataFinal$prev_strat == escount$prev_strat[i] & 
                dataFinal$studygroup == escount$studygroup[i] & 
                dataFinal$is_single == 0 & 
                dataFinal$is_adj == 1 & 
                dataFinal$is_ds_inf == 1 & 
                dataFinal$is_ds_con == 1,]$keepFinal <- 1 
  } 
  # Else, if any ES are adjusted and match the
  # common study-level (within each analytic dataset) informant, keep them
  else if(nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                         dataFinal$prev_strat == escount$prev_strat[i] & 
                         dataFinal$studygroup == escount$studygroup[i] & 
                         dataFinal$is_single == 0 & 
                         dataFinal$is_adj == 1 & 
                         dataFinal$is_ds_inf == 1 & 
                         dataFinal$is_ds_con == 0,]) > 0){
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                dataFinal$prev_strat == escount$prev_strat[i] & 
                dataFinal$studygroup == escount$studygroup[i] & 
                dataFinal$is_single == 0 & 
                dataFinal$is_adj == 1 & 
                dataFinal$is_ds_inf == 1 & 
                dataFinal$is_ds_con == 0,]$keepFinal <- 1 
  } 
  # Else, if any ES are unadjusted and match the
  # common analytic dataset informant and micro construct, keep them
  else if (nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                          dataFinal$prev_strat == escount$prev_strat[i] & 
                          dataFinal$studygroup == escount$studygroup[i] & 
                          dataFinal$is_single == 0 & 
                          dataFinal$is_adj == 0 & 
                          dataFinal$is_d_inf == 1 & 
                          dataFinal$is_d_con == 1,]) > 0) {
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                dataFinal$prev_strat == escount$prev_strat[i] & 
                dataFinal$studygroup == escount$studygroup[i] & 
                dataFinal$is_single == 0 & 
                dataFinal$is_adj == 0 & 
                dataFinal$is_d_inf == 1 & 
                dataFinal$is_d_con == 1,]$keepFinal <- 1 #
  } 
  # Else, if any ES are unadjusted and match the
  # common analytic dataset informant, keep them
  else if(nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                         dataFinal$prev_strat == escount$prev_strat[i] & 
                         dataFinal$studygroup == escount$studygroup[i] & 
                         dataFinal$is_single == 0 & 
                         dataFinal$is_adj == 0 & 
                         dataFinal$is_d_inf == 1 & 
                         dataFinal$is_d_con == 0,]) > 0){
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                dataFinal$prev_strat == escount$prev_strat[i] & 
                dataFinal$studygroup == escount$studygroup[i] & 
                dataFinal$is_single == 0 & 
                dataFinal$is_adj == 0 & 
                dataFinal$is_d_inf == 1 & 
                dataFinal$is_d_con == 0,]$keepFinal <- 1  
  } 
  # Else, if any ES are unadjusted and match the
  # common study-level (within each analytic dataset) informant and micro construct, keep them
  else if (nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                          dataFinal$prev_strat == escount$prev_strat[i] & 
                          dataFinal$studygroup == escount$studygroup[i] & 
                          dataFinal$is_single == 0 & 
                          dataFinal$is_adj == 0 & 
                          dataFinal$is_ds_inf == 1 & 
                          dataFinal$is_ds_con == 1,]) > 0) {
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                dataFinal$prev_strat == escount$prev_strat[i] & 
                dataFinal$studygroup == escount$studygroup[i] & 
                dataFinal$is_single == 0 & 
                dataFinal$is_adj == 0 & 
                dataFinal$is_ds_inf == 1 & 
                dataFinal$is_ds_con == 1,]$keepFinal <- 1 
  } 
    # Else, if any ES are unadjusted and match the
  # common study-level (within each analytic dataset) informant, keep them
  else if(nrow(dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                         dataFinal$prev_strat == escount$prev_strat[i] & 
                         dataFinal$studygroup == escount$studygroup[i] & 
                         dataFinal$is_single == 0 & 
                         dataFinal$is_adj == 0 & 
                         dataFinal$is_ds_inf == 1 & 
                         dataFinal$is_ds_con == 0,]) > 0){
    dataFinal[dataFinal$dv2_new==escount$dv2_new[i] & 
                dataFinal$prev_strat == escount$prev_strat[i] & 
                dataFinal$studygroup == escount$studygroup[i] & 
                dataFinal$is_single == 0 & 
                dataFinal$is_adj == 0 & 
                dataFinal$is_ds_inf == 1 & 
                dataFinal$is_ds_con == 0,]$keepFinal <- 1 
  } 
}


# KEEP THE EFFECT SIZES FLAGGED AS KEEPFINAL=1
dataFinal <- dataFinal[dataFinal$keepFinal==1,]
```

Modeling
========

This code is from early in the modeling process. It contains multiple methods of running meta-analytic models, and none of the models contain covariates.

``` r
data <- dataFinal
```

Analysis Options
----------------

This section of the code contains multiple chunks that were designed be run on an as-needed basis.

Some of our preliminary analyses will inlude all effect sizes. However, others will only include one effect size per study within each analytic dataset.

There are two ways to achieve this:

1.  Within each analytic dataset, randomly select one effect size from each study

``` r
# OPTION 1: RANDOMLY SELECTED EFFECT SIZE

# Run this block if you are running single-level analyses in which
# the effect size is randomly selected

# Add a column of random numbers to the dataset without replacement
data$rand_rnd <- sample(nrow(data), size=nrow(data), replace=FALSE)  

# Find the smallest randomly generated number within each study within each analytic dataset
agg_rnd <- aggregate(data$rand_rnd, by=list(data$studygroup, data$dv2_new, data$prev_strat), min) 
colnames(agg_rnd)[c(1:4)] <- c("studygroup","dv2_new", "prev_strat","min_rand_rnd")

# Merge that information back into the main dataset
data <- merge(data, agg_rnd, by=c("studygroup","dv2_new","prev_strat"))

# Flag rows in which the random number matches the merged minimum random number
# >> select_rnd == 1. This is how we'll identify them later
data$select_rnd <- ifelse(data$rand_rnd == data$min_rand_rnd, 1, 0) 
data <- subset(data, select = -c(rand_rnd, min_rand_rnd))

rm(agg_rnd)
```

1.  Within each analytic dataset, compute the weighted mean and weighted variance for each study.

``` r
# OPTIONAL: WEIGHTED MEAN EFFECT SIZE AND WEIGHTED MEAN VARIANCE DATA FRAME

# Run this block if you are running single level analyses in which 
# the effect size and variance are weighted means of the study group

# Compute the weighted mean effect size and variance within each 
# study within each analytic dataset
agg_wtmeans <-data %>%
  group_by(studygroup, dv2_new, prev_strat) %>%
  summarise(wmean = weighted.mean(es_preadj, eff_es3_tukey),
            wvar = weighted.mean(esvar, eff_es3_tukey))
agg_wtmeans <- as.data.frame(agg_wtmeans)

# Merge the weighted means into the main dataset
data <- merge(data, agg_wtmeans, by=c("studygroup","dv2_new","prev_strat"))


# Select one of the rows for use in the analysis (same method as above)

# Add a column of random numbers to the dataset without replacement
data$rand_wt <- sample(nrow(data), size=nrow(data), replace=FALSE)

# Find the smallest randomly generated number within each study within each analytic dataset
agg_wtmeans <- aggregate(data$rand_wt, by=list(data$studygroup, data$dv2_new, data$prev_strat), min)
colnames(agg_wtmeans)[c(1:4)] <- c("studygroup","dv2_new", "prev_strat","min_rand_wt")

# Merge that information back into the main dataset
data <- merge(data, agg_wtmeans, by=c("studygroup","dv2_new","prev_strat"))

# Flag rows in which the random number matches the merged minimum random number
# >> select_rnd == 1. This is how we'll identify them later
data$select_wt <- ifelse(data$rand_wt == data$min_rand_wt, 1, 0)
data <- subset(data, select = -c(rand_wt, min_rand_wt))

rm(agg_wtmeans)
```

Alternatively, this chunk of code is designed to be used when conducting analyses on all effect sizes (multiple per study). It flags studies that only have one effect size so they can be removed from the analysis (if desired).

``` r
# OPTIONAL: NO STUDY GROUP SINGLETONS DATA FRAME

# Run this block if you are running multi-level analyses and wish to 
# drop study groups with only one effect size (within each analytic dataset)

data$rand_wt <- sample(nrow(data), size=nrow(data), replace=FALSE)
# We know rand_wt is a variable with 0 missing values. Get its length.
agg_singletons <- aggregate(data$rand_wt, 
                            by=list(data$studygroup, data$dv2_new, data$prev_strat), length)
colnames(agg_singletons)[c(1:4)] <- c("studygroup","dv2_new", "prev_strat","EScount")

# single = 1 if it's a study singleton, 0 if not. We will drop the 1, keep the 0.
agg_singletons$single <- ifelse(agg_singletons$EScount==1, 1, 0) 

# Merge that information back into the main dataset
data <- merge(data, agg_singletons, by=c("studygroup","dv2_new","prev_strat"))
data <- subset(data, select = -c(rand_wt, EScount))

rm(agg_singletons)
```

Run models
----------

As mentioned above, these are all early models with no moderators. All possible analyses are included in this chunk of code. They were designed to be run piecemeal as desired.

``` r
# ANALYSES -----
# The following blocks of code perform different types of analyses.
# Do NOT run all blocks - only run the block for the analysis you want to perform.

# If multiple analyses are to be run - run one analysis, then export it (see next section of code) 
# >> then run the next analysis, then export it, etc.

# Each analysis will be preceeded by a brief description. 
# Be sure to read the rma.mv/uni.mv syntax closely, 
# as you may wish to make small changes (e.g., to the method argument).

# NO MODERATORS ----------

# ANALYSIS: Multilevel analysis ----
# >> full dataset
# >> Using Hedges G and log OR
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=es_preadj, V=esvar, 
           intercept=TRUE, 
           data=data, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
    ))
  # Print model summary to console
  print(summary(get(key$model[i]))) 
}

### ANALYSIS: Multilevel analysis ----
# > using weighted mean (wmean) and variance (wvar) from espicker syntax
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=wmean, V=wvar, 
           intercept=TRUE, 
           data=data, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
  ))
  # Print model summary to console
  print(summary(get(key$model[i])))
}

# ANALYSIS: Multilevel analysis ----
# >> Using Hedges G and log OR
# >> weights=1
data$wone <- 1
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=es_preadj, V=esvar, 
           W=wone, 
           intercept=TRUE, 
           data=data, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i], 
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
  ))
  # Print model summary to console
  print(summary(get(key$model[i]))
}

# ANALYSIS: Multilevel analysis ----
# >> Using non-weighted Hedges G and log OR
# no singletons (single==0)
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.mv(yi=es_preadj, V=esvar,  
           intercept=TRUE, 
           data=data, 
           subset = dv2_new==key$dv[i] & prev_strat==key$ps[i] & single==0, 
           sparse = TRUE, 
           random = list(~1|studygroup,~1|esid), 
           method = "REML")
  ))
  # Print model summary to console
  print(summary(get(key$model[i])))
}

# ANALYSIS: Single level analysis ----
# Weighted mean effect size and variance (wmean, wvar) 
# One randomly selected effect size per study
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.uni(yi=wmean, vi=wvar, 
            intercept=TRUE, 
            data=data,
            subset=dv2_new==key$dv[i] & prev_strat==key$ps[i] & select_wt==1,
            method="REML")
  ))
  # Print model summary to console
  print(summary(get(key$model[i])))
}

# ANALYSIS: Single level analysis ----
# Using Hedges G and log OR
# One randomly selected effect size per study
for (i in 1:nrow(key)){
  assign(eval(key$model[i]),eval(
    rma.uni(yi=es_preadj, vi=esvar, 
            intercept=TRUE, 
            data=data,
            subset=dv2_new==key$dv[i] & prev_strat==key$ps[i] & select_rnd==1,
            method="REML")
  ))
  # Print model summary to console
  print(summary(get(key$model[i])))
}
```

Export model information
------------------------

``` r
# EXPORT ----- 

# After running each analysis, export the output.

output <-data.frame()

# Add a brief description of the output to be used in the .csv filename
# >> Date information will be appended automatically
fn <- "MODEL"

# Run one of the two code blocks that follow -----
# >> The first works on rma.mv output (Run if Multi Level). 
# >> The second works on rma.uni output (Run if Single Level)

# 1. Run if Multi-Level -----
for(i in 1:nrow(key)){ # As usual, we're iterating across analytic datasets
  # Get information from the rma.mv  object:
  temp_ci <- confint(eval(as.name(key$model[i])))
  # Model Name
  output[i,1] <- key$model[i]
  # Study Level Tau^2
  output[i,2] <- round(eval(as.name(key$model[i]))[[8]][1], digits=4) 
  # Study Level Sqrt(Tau^2)
  output[i,3] <- round(sqrt(eval(as.name(key$model[i]))[[8]][1]), digits=4) 
  # Study tau^2, CI lower bound
  output[i,4] <- round(temp_ci[[1]][[1]][3], digits=4) 
  # Study tau^2, CI upper bound
  output[i,5] <- round(temp_ci[[1]][[1]][5], digits=4) 
  # Study Levels
  output[i,6] <- eval(as.name(key$model[i]))$s.nlevels[1] 
  # Effect Size Level Tau^2
  output[i,7] <- round(eval(as.name(key$model[i]))[[8]][2], digits=4)
  # Effect Size Level Sqrt(Tau^2)
  output[i,8] <- round(sqrt(eval(as.name(key$model[i]))[[8]][2]), digits=4) 
  # Effect Size tau^2, CI lower bound
  output[i,9] <- round(temp_ci[[2]][[1]][3], digits=4) 
  # Effect Size tau^2, CI upper bound
  output[i,10] <- round(temp_ci[[2]][[1]][5], digits=4) 
  # Effect Size Levels
  output[i,11] <- eval(as.name(key$model[i]))$s.nlevels[2] 
  # Total variance
  output[i,12] <- round(eval(as.name(key$model[i]))[[8]][1] + 
                          eval(as.name(key$model[i]))[[8]][2], 
                        digits=4) 
  # % Variance at study level
  output[i,13] <- round(eval(as.name(key$model[i]))[[8]][1] / 
                          (eval(as.name(key$model[i]))[[8]][1] +
                             eval(as.name(key$model[i]))[[8]][2]) * 100, 
                        digits=4) 
  # % Variance at effect size level
  output[i,14] <- round(eval(as.name(key$model[i]))[[8]][2] / 
                          (eval(as.name(key$model[i]))[[8]][1] + 
                             eval(as.name(key$model[i]))[[8]][2])*100, 
                        digits=4) 
  # Q
  output[i,15] <- round(eval(as.name(key$model[i]))[[20]], digits=4)
  # Q df (Effect Size levels - 1)
  output[i,16] <- eval(as.name(key$model[i]))[[14]]-1 
  # Model effect size mean
  output[i,17] <- ifelse(key$es[i]=="SMD", 
                         # IF analyzing SMDs can use as-is:
                              round(eval(as.name(key$model[i]))[[1]][1], digits=4), 
                         # IF ORs, use exp(mean)
                              exp(round(eval(as.name(key$model[i]))[[1]][1], digits=4)))
  # Standard Error
  output[i,18] <- round(eval(as.name(key$model[i]))[[2]][1], digits=4) 
  # Model name/label provided by the key
  output[i,19] <- toString(eval(as.name(key$model[1]))$call)
}

# Attach column names to the output
colnames(output) <- c("Model", 
                          "SG Tau^2", "SG Sqrt(Tau^2)","SG CI lower", "SG CI upper", "SG k", 
                          "ES Tau^2", "ES Sqrt(Tau^2)", "ES CI lower", "ES CI upper","ES n", 
                          "Total Variance", "% at SG level", "% at ES level",
                          "Q", "Q df", "Mean ES", "SE", "Model")

# Save the output to a .csv
write.csv(output, file= paste0("RMA-MV_",fn,"_",format(Sys.time(), "%Y%m%d_%H%M%S_"),".csv") )

# Run if Single Level -----
for(i in 1:nrow(key)){ # As usual, we're iterating across analytic datasets
  # Get information from the rma.uni  object:
  temp_ci <- confint(eval(as.name(key$model[i])))
  # Model Name
  output[i,1] <- key$model[i]
  # Tau^2
  output[i,2] <- round(eval(as.name(key$model[i]))[[8]], digits=4) 
  # Sqrt(Tau^2)
  output[i,3] <- round(sqrt(eval(as.name(key$model[i]))[[8]]), digits=4) 
  # tau^2, CI lower bound
  output[i,4] <- round(temp_ci[[1]][5], digits=4) 
  # tau^2, CI upper bound
  output[i,5] <- round(temp_ci[[1]][9], digits=4) 
  # Study Group Levels (k)
  output[i,6] <- eval(as.name(key$model[i]))[[11]]  
  # Q
  output[i,7] <- round(eval(as.name(key$model[i]))[[18]], digits=4) 
  # Q df (Effect Size levels - 1)
  output[i,8] <- eval(as.name(key$model[i]))[[11]]-1 
  # Model effect size mean
  output[i,9] <- ifelse(i<35, 
                        # IF analyzing SMDs can use as-is:
                        round(eval(as.name(key$model[i]))[[1]][1], digits=4),
                        # IF ORs, use exp(mean)
                        exp(round(eval(as.name(key$model[i]))[[1]][1], digits=4)))
  # Model effect size standard error 
  output[i,10] <- round(eval(as.name(key$model[i]))[[2]], digits=4) 
  output[i,11] <- toString(eval(as.name(key$model[1]))$call)
}

# Attach column names to the output
colnames(output) <- c("Model", "Tau^2", "Sqrt(Tau^2)","CI lower", "CI upper", 
                      "SG k", "Q", "Q df", "Mean ES", "SE","Model")

# Save the output to a .csv
write.csv(output, file= paste0("RMA-UNI_",fn,"_",format(Sys.time(), "%Y%m%d_%H%M%S_"),".csv") )
```
