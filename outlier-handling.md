Outlier Handling Methods
================

Background and Description
--------------------------

Over the course of the project we took different approaches to outlier handling:

1.  Replace outliers with "narrow" Tukey fences (Q<sub>1</sub>-1.5\*(IQR) and Q<sub>3</sub>+1.5\*(IQR))
2.  Replace outliers with "wide" Tukey fences (Q<sub>1</sub>-3\*(IQR) and Q<sub>3</sub>+3\*(IQR))
3.  Replace outliers with the most extreme non-outlier values present in the vector, defined by "narrow" Tukey fences
4.  Replace outliers with the most extreme non-outlier values present in the vector, defined by "wide" Tukey fences

In addition, an alternative to Tukey's fences was identified. This method adjusts for skewed distributions and is implemented in the `robustbase` package (see also: Hubert and Vandervieren (2004)). Combined with the above four approaches, that leaves eight approaches to outlier handling.

The function below combines all approaches into one function.

Usage
-----

`wt_winsor(x, coef = 1.5, to = "fence", skew = FALSE)`

Arguments
---------

<table style="width:88%;">
<colgroup>
<col width="19%" />
<col width="68%" />
</colgroup>
<thead>
<tr class="header">
<th>Argument</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>x</code></td>
<td>a numeric vector</td>
</tr>
<tr class="even">
<td><code>coef</code></td>
<td>number to be multiplied by the interquartile range</td>
</tr>
<tr class="odd">
<td><code>to</code></td>
<td>identifies the point to which outliers should be adjusted. If &quot;fence&quot;, replaces outliers with the outlier cutoff. If &quot;value&quot;, replaces outliers with the most extreme non-oulier present in the dataset</td>
</tr>
<tr class="even">
<td><code>skew</code></td>
<td>If <code>FALSE</code>, uses Tukey's fence values. If <code>TRUE</code> uses values from <code>robustbase</code> as outlined in Hubert and Vandervieren (2004)</td>
</tr>
</tbody>
</table>

This function ignores `NA` values. If you do not wish to run it on a vector that includes missing data, that check must be manually performed.

This function returns a numeric vector.

Function
========

``` r
wt_winsor <- function(x, coef=1.5, to="fence", skew = FALSE){
  
  # Stop function if inappropriate argument value
  if (coef < 0 | !is.numeric(coef)) 
    stop("'coef' must be a positive number")
  
  if (!is.numeric(x))
    stop("'x' must be a numeric vector")
  
  if(!to %in% c("fence", "value"))
    stop("Invalid argument: 'to' must be 'fence' or 'value'")
  
  if(!is.logical(skew))
    stop("'skew' must specify a logical value")
  
  # If want to adjust for skewness, need to use robustbase package
  # if(skew == TRUE) 
    # [--- SEE NOTE 1 BELOW ---]
    # requireNamespace("robustbase") ????
    # require(robustbase) ???
  
  if(skew == FALSE){
    # If not adjusting for skewness,
    # Get 1st and 3rd quartiles and calculate the interquartile range 
    q1 <- quantile(x, na.rm=TRUE)[[2]] # Quartile 1 (Q1), 25% 
    q3 <- quantile(x, na.rm=TRUE)[[4]] # Quartile 3 (Q3), 75%
    iqr <- q3 - q1 
    # Calculate Tukey fences
    low_tukey  <- q1 - (coef * iqr) 
    high_tukey <- q3 + (coef * iqr)
  } else if(skew == TRUE){
    # If adjusting for skewness,
    # Use fences provided by adjboxStats
    low_tukey  <- robustbase::adjboxStats(x, coef=coef)$fence[1] 
    high_tukey <- robustbase::adjboxStats(x, coef=coef)$fence[2]
  }
  
  # Flag high and low outliers (TRUE if outlier)
  high_outlier <- ifelse(x > high_tukey, TRUE, FALSE)
  low_outlier <- ifelse(x < low_tukey, TRUE, FALSE) 
  
  # Create vector of corrected values
  if(all(x==x[min(which(!is.na(x)))], na.rm=TRUE)){
    # If all values of x are equal, return the original vector
    # (Compares vector against the first non-NA value in the vector)
    # [--- SEE NOTE 2 ---]
    es_tukey <- x
  } else if (to == "fence") { # If using fence
    # Replace value with high/low fence for high/low outliers, respectively 
    # Keep original values for non-outliers
    es_tukey <- ifelse(low_outlier, low_tukey,
                       ifelse(high_outlier, high_tukey, x))
  } else if (to == "value") { # If using most extreme non-outlier
    # Keep original values for non-outliers
    es_tukey <- ifelse(!high_outlier & !low_outlier, x, NA)
    # If is outlier, replace with most extreme non-outlier
    es_tukey <- ifelse(high_outlier, max(es_tukey, na.rm=TRUE), 
                       ifelse(low_outlier, min(es_tukey, na.rm=TRUE),
                              es_tukey))
  }
  
  # Show a warning and continue if coefficient is outside of traditional values
  if (coef > 3 | coef < 1.5)
    warning("Output calcuated using a 'coef' outside of expected range (1.5 to 3)")
  
  #return(cbind(x, low_tukey, high_tukey, high_outlier, low_outlier, es_tukey)) # remove line after testing complete
  return(es_tukey)
}
```

Notes:
------

1.  Need to read up on best practices for handling non-base packages in functions. Starting points:
    -   <https://stackoverflow.com/questions/23232791/is-it-a-good-practice-to-call-functions-in-a-package-via>
    -   <https://stackoverflow.com/questions/15258398/whats-the-impact-of-requiring-a-package-inside-a-function-if-the-package-is-alre>
    -   <https://stackoverflow.com/questions/46270860/importing-a-library-inside-a-function>
    -   <https://www.r-bloggers.com/difference-between-library-and-require-in-r/>
2.  It would be more efficient to move this chunk of code near the top of the function. Also, surely the logic can be less obtuse.

This function is still under development.

Example
-------

Create a small dataset that should contain outliers:

``` r
set.seed(1234)
rand_val <- rnorm(n=10, mean=20, sd=5)
ind <- which(rand_val %in% sample(rand_val, 2))
rand_val[ind]<-rand_val[ind]+100
dat <- data.frame(value = rand_val)
```

Function output:

``` r
dat$value_winsorized <- wt_winsor(dat$value)
dat
```

    ##        value value_winsorized
    ## 1   13.96467         13.96467
    ## 2   21.38715         21.38715
    ## 3  125.42221         30.37655
    ## 4  108.27151         30.37655
    ## 5   22.14562         22.14562
    ## 6   22.53028         22.53028
    ## 7   17.12630         17.12630
    ## 8   17.26684         17.26684
    ## 9   17.17774         17.17774
    ## 10  15.54981         15.54981
