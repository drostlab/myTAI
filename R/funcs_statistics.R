#' @title Function to compute the TAI or TDI destruction score of the global hourglass pattern based on a pre-computed TAI or TDI vector.
#' @description
#' 
#'         This function reduces the destruction of an hourglass shaped pattern to a single score value.
#'         
#'         Based on a given \code{\link{TAI}} or \code{\link{TDI}} pattern the given vector is being
#'         divided into three developmental modules: early, mid, and late. The corrisponding \code{\link{TAI}} or \code{\link{TDI}}
#'         values in each developmental module are accumulated using the \emph{scoringMethod} argument ("max-min" or "mean-mean").
#'         
#'         In more detail:
#'                 
#'         (1) for a given \code{\link{TAI}} or \code{\link{TDI}} vector \emph{tai_profile} or \emph{tdi_profile}, we classify each value of \emph{tai_profile} or \emph{tdi_profile} into its corresponding developmental module early, mid, or late.
#'         
#'         (2) accumulate the \emph{tai_profile} or \emph{tdi_profile} values in each developmental module using the arithmetic mean (\code{\link{mean}}) in case scoringMethod = "mean-mean", or accumulate the \emph{tai_profile} or \emph{tdi_profile} values in each developmental module using \code{\link{max}} for the early and late module and \code{\link{min}} for the mid module in case scoringMethod = "max-min".
#'         
#'         (3) then reduce the three values for each developmental module by computing the difference between: early - mid, and late - mid.
#'         
#'         (4) the two difference values are referred to as a_early and a_late. 
#v         
#'         
#'         Each developmental module now has an accumulated representation value which is being reduced to one value using the
#'         \emph{method} argument ("max", "min", or "mean"). 
#'         
#v         The "max", "min", or "mean" parameters refer to the following reduction procedure:
#'                 
#'                 Given the two accumulated values for each hourglass module: a_early and a_late,
#'         we reduce the two given values by:
#'                 
#'                 \emph{"max"}
#'         
#'         \eqn{S = max{a_early,a_late}}
#'         
#'         \emph{"min"}
#'         
#'         \eqn{S = min{a_early,a_late}}
#'         
#'         \emph{"mean"}
#'         
#'         \eqn{S = mean{a_early,a_late}}
#'         
#'         All together this results in a global score \emph{S}.
#'         This global score \emph{S} is being returned by this function \code{\link{gpScore}}.
#'         
#' @param age_vals a numeric vector containing \code{\link{TAI}} or \code{\link{TDI}} values for each developmental stage s.
#' @param early a numeric vector including the numeric stage values that correspond to the early phase of development.
#' @param mid a numeric vector including the numeric stage values that correspond to the middle phase of development.
#' @param late a numeric vector including the numeric stage values that correspond to the late phase of development.
#' @param method to determine the two value reduction value, resulting in the global score S: "max", or "min", or "mean".
#' @param scoringMethod method to determine the module accumulation value: "max-min" or "mean-mean".
#' @details
#' 
#' The gpScore is a heuristic score enabling to  construct a test statistic to determine  
#' the significance of a present (phylotranscriptomic) hourglass pattern.
#' 
#' @return a numeric value representing the hourglass destruction score.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{ReductiveHourglassTest}}, \code{\link{TAI}}, \code{\link{TDI}}
#' @examples \dontrun{
#' 
#'  # read standard phylotranscriptomics data
#'  data(PhyloExpressionSetExample)
#'  data(DivergenceExpressionSetExample)
#'
#'  # example PhyloExpressionSet:
#'
#'  # compute the TAI profile
#'  TAIs <- TAI(PhyloExpressionSetExample)
#'
#'  # compute the global hourglass destruction score for the TAIs profile using reduction method: mean(mean-mean)
#'  gScore <- gpScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,method = "mean",scoringMethod = "mean-mean")
#'
#'
#'  # example DivergenceExpressionSet:
#'
#'  # compute the TDI profile
#'  TDIs <- TDI(DivergenceExpressionSetExample)
#'
#'  # compute the global hourglass destruction score for the TDIs profile using reduction method: mean(mean-mean)
#'  gScore <- gpScore(age_vals = TDIs,early = 1:2,mid = 3:5,late = 6:7,method = "mean",scoringMethod = "mean-mean")
#' }
#' @export
gpScore <- function(age_vals,early,mid,late,method,scoringMethod)
{
  
  Score.Early <- vector(mode = "numeric", length = 1)
  Score.Late <- vector(mode = "numeric", length = 1)
  age_valsEarly <- vector(mode = "numeric",length = length(early))
  age_valsMid <- vector(mode = "numeric",length = length(mid))
  age_valsLate <- vector(mode = "numeric",length = length(late))
  

  age_valsEarly <- age_vals[early]
  age_valsMid <- age_vals[mid]
  age_valsLate <- age_vals[late]
  
  
  
  if(scoringMethod == "max-min"){
      Score.Early <- max(age_valsEarly) - min(age_valsMid)
      Score.Late <- max(age_valsLate) - min(age_valsMid) 
  }
  
  if(scoringMethod == "mean-mean"){
      Score.Early <- mean(age_valsEarly) - mean(age_valsMid)
      Score.Late <- mean(age_valsLate) - mean(age_valsMid)    
  }

  if(method == "max")
     return(max(Score.Early, Score.Late))
  if(method == "min")
     return(min(Score.Early, Score.Late))
  if(method == "mean")
     return(mean(Score.Early, Score.Late))


}


#' @title Function to perform the Flat Line Test.
#' @description This function computes N phylotranscriptomics profiles of a randomly sampled (row-permutations: random permutations of phylostratum or divergence-stratum assignments) 
#' PhyloExpressionSets or DivergenceExpressionSets.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param permutations a numeric value specifying the number of permutations to be performed for the \emph{FlatLineTest}.
#' @param plotHistogram a boolean value specifying whether a detailed statistical analysis concerning the goodness of fit shall be performed.
#' @param parallel a boolean value specifying whether goodness of fit computations shall be performed in parallel using the \emph{doParallel} package.
#' @param runs specify the number of runs to be performed for goodness of fit computations. In most cases runs = 100 is a reasonable choice.
#' @details Internally the function performs N phylotranscriptomics pattern computations (\code{\link{TAI}} or \code{\link{TDI}}) based on sampled PhyloExpressionSets or DivergenceExpressionSets (see \code{\link{bootMatrix}}). The test statistics is being developed as follows:
#'
#' The variance \emph{V_pattern} of the S phylotranscriptomics values defines the test statistic for the \code{\link{FlatLineTest}}. The basic assumption is, that the variance of a flat line should be equivalent to zero for a perfect flat line. Any deviation from a flat line can be measured with a variance value > 0. 
#'
#' To determine the null distribution of \emph{V_p}, all PS or DS values within each developmental stage s are randomly permuted, S surrogate phylotranscriptomics values are computed from this permuted dataset, and a surrogate value of \emph{V_p} is computed from these S phylotranscriptomics values. This permutation process is repeated N times, yielding a histogram of \emph{V_p}. 
#'
#' After applying a \emph{Lilliefors Kolmogorov-Smirnov Test for gamma distribution}, \emph{V_p} is approximated by a \emph{gamma distribution}. The two parameters of the \emph{gamma distribution} are estimated by the function \emph{fitdistrplus::fitdist} from the \emph{fitdistrplus} package by \emph{moment matching estimation}. The fitted \emph{gamma distribution} is considered the null distribution of \emph{V_pattern}, and the p-value of the observed value of \emph{V_p} is computed from this null distribution.
#'
#' In case the parameter \emph{plotHistogram = TRUE}, a multi-plot is generated showing:
#'        
#' (1) A Cullen and Frey skewness-kurtosis plot generated by fitdistrplus::descdist().
#'
#' (2) A histogram of V_p combined with the density plot using the Maximum Likelihood estimated parameters returned by the MASS::fitdistr() function using a gamma distribution.
#'
#' (3) A plot showing the p-values for N independent runs to verify that a specific p-value is biased by a specific permutation order.
#'
#' The \emph{goodness of fit} for the random vector \emph{V_p} is quantified statistically by an adapted Lilliefors (Kolmogorov-Smirnov) test for gamma distributions.
#' @return a list object containing the list elements:
#' \item{p.value}{the p-value quantifying the statistical significance (deviation from a flat line) of the given phylotranscriptomics pattern.}
#' \item{std.dev}{the standard deviation of the N sampled phylotranscriptomics patterns for each developmental stage S.}
#' @references Quint M et al. (2012). "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#'
#' M. L. Delignette-Muller, R. Pouillot, J.-B. Denis and C. Dutang (2014), fitdistrplus: help to fit of a parametric distribution to non-censored or censored data.
#'
#' Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-159.
#'
#' Evans M, Hastings N and Peacock B (2000) Statistical distributions. John Wiley and Sons Inc.
#'
#' Sokal RR and Rohlf FJ (1995) Biometry. W.H. Freeman and Company, USA, pp. 111-115.
#'
#' Juergen Gross and bug fixes by Uwe Ligges (2012). nortest: Tests for Normality. R package version 1.0-2. http://CRAN.R-project.org/package=nortest
#'
#' Dallal, G.E. and Wilkinson, L. (1986): An analytic approximation to the distribution of Lilliefors test for normality. The American Statistician, 40, 294–296.
#'
#' Stephens, M.A. (1974): EDF statistics for goodness of fit and some comparisons. Journal of the American Statistical Association, 69, 730–737.
#'
#' http://stackoverflow.com/questions/4290081/fitting-data-to-distributions?rq=1
#'
#' http://stats.stackexchange.com/questions/45033/can-i-use-kolmogorov-smirnov-test-and-estimate-distribution-parameters
#' 
#' http://cran.r-project.org/doc/contrib/Ricci-distributions-en.pdf
#'
#' http://cran.r-project.org/doc/contrib/Ricci-distributions-en.pdf
#'
#' @author Hajk-Georg Drost
#' @note In case there are extreme outlier expression values stored in the dataset (PhyloExpressionSet or DivergenceExpressionSet),
#'  the internal MASS::fitdistr function that is based on the \code{\link{bootMatrix}} output might return a warning:
#'   "In densfun(x, parm[1], parm[2], ...) : NaNs were produced" which indicates that permutation results caused by extreme outlier expression values 
#'   cannot be fitted accordingly. This warning will not be printed out when the corresponding outlier values are extracted from the dataset.
#' @seealso \code{\link{TAI}}, \code{\link{TDI}}, \code{\link{PlotPattern}}, \code{\link{bootMatrix}}
#' @examples \dontrun{
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#'
#' # example PhyloExpressionSet using 1000 permutations
#' FlatLineTest(PhyloExpressionSetExample, permutations = 1000, plotHistogram = FALSE)
#'
#' # example DivergenceExpressionSet using 1000 permutations
#' FlatLineTest(DivergenceExpressionSetExample, permutations = 1000, plotHistogram = FALSE)
#'
#' # perform the Kolmogorov-Smirnov-Test and plot the resulting histogram for a PhyloExpressionSet
#' FlatLineTest(PhyloExpressionSetExample, permutations = 1000, plotHistogram = TRUE, runs = 100)
#'
#' # in case you have a multi-core processor and the doParallel package installed on you system,
#' # you can run the tests shown above in parallel 
#' # perform the Kolmogorov-Smirnov-Test and plot the resulting histogram for a PhyloExpressionSet
#' # FlatLineTest(PhyloExpressionSetExample, permutations = 1000, plotHistogram = TRUE, parallel = TRUE, runs = 100)
#' 
#'
#' # Example: finding outlier expression levels that badly influence the permutation test statistic
#' # this plot visualizes the distribution and variance of permuted TAI values for each developmental stage
#' # In case there are outlier expression levels in a specific developmental stage,
#' # the corresponding boxplot of this stage also shows an unusual fluctuation of TAI values
#' # compared to other stages
#' boxplot(bootMatrix(PhyloExpressionSetExample) , ylab = "TAI")
#' 
#' # analogous for TDI permutation results
#' boxplot(bootMatrix(DivergenceExpressionSetExample) , ylab = "TDI")
#' 
#' } 
#' @export
FlatLineTest <- function(ExpressionSet,permutations=1000,plotHistogram=FALSE,parallel=FALSE,runs=10)
{
  
  is.ExpressionSet(ExpressionSet)
  
  if((plotHistogram == TRUE) & is.null(runs))
    stop("Please specify the number of runs to be performed for the goodness of fit computations.")
  
  nCols <- dim(ExpressionSet)[2]
  resMatrix <- matrix(NA_real_, permutations,(nCols-2))
  var_values <- vector(mode = "numeric", length = permutations)
  sd_values <- vector(mode = "numeric",length = nCols-2)
  #random_mean_age <- vector(mode = "numeric", length = permutations)
  age.real <- vector(mode = "numeric",length = nCols-2)
  age.real <- TAI(ExpressionSet)
  ### compute the real variance of TAIs of the observed TAI/TDI-Hourglass pattern
  real.var <- var(age.real)
  ### sample only the phylostrata (row-permutations) 
  
  resMatrix <- bootMatrix(ExpressionSet, permutations)
  var_values <- apply(resMatrix,1,var)
  #random_mean_age <- apply(resMatrix,2,mean)
  ### estimate the maximum likelihood parameter (shape,rate) 
  ### of the gamma distributed variance values
  #gamma_MLE <- MASS::fitdistr(var_values,"gamma")
  gamma_MME <- fitdistrplus::fitdist(var_values,"gamma",method = "mme")
  ### estimate shape:
  #shape <- est.shape(var_values)
  shape <- gamma_MME$estimate[1]
  ### estimate the rate:
  #rate <- est.rate(var_values)
  rate <- gamma_MME$estimate[2]
  
  if(plotHistogram == TRUE){
    ### Perform a Kolmogorov–Smirnov test to quantify a
    ### Gamma-distribution as Null-Distribution
    #KS.Test.Output <- ks.test(var_values,"pgamma",shape = shape,rate = rate)
    # plot histogram of standard deviations
    gammaDensity <- function(x){
      
       return(dgamma(x = x,shape = shape,rate = rate))
       
    }
    
    par(mfrow = c(1,3))
    # plot a Cullen and Frey graph
    fitdistrplus::descdist(var_values, boot = permutations)
    # plot the histogram and the fitted curve
    curve(gammaDensity,xlim = c(min(var_values),max(c(var_values,real.var))),col = "steelblue",lwd=5,xlab = "Variances",ylab="Frequency", main = paste0("permutations = ",permutations))
    histogram <- hist(var_values,prob = TRUE,add = TRUE, breaks = permutations / (0.01 * permutations))
    rug(var_values)
    abline(v = real.var, lty = 1, lwd = 4, col = "darkred")
    
    p.vals_vec <- vector(mode = "numeric", length = runs)
    
    if(parallel == TRUE){
      ### Parallellizing the sampling process using the 'doMC' and 'parallel' package
      ### register all given cores for parallelization
      ### detectCores(all.tests = TRUE, logical = FALSE) returns the number of cores available on a multi-core machine
      cores <- parallel::makeForkCluster(parallel::detectCores(all.tests = FALSE, logical = FALSE))
      doParallel::registerDoParallel(cores)
      
      ### Perform the sampling process in parallel
      p.vals_vec <- as.vector(foreach::foreach(i = 1:runs,.combine="c") %dopar% {FlatLineTest(ExpressionSet)$p.value})
      
      ### close the cluster connection
      ### The is important to be able to re-run the function N times
      ### without getting cluster connection problems
      stopCluster(cores)
    }
    
    if(parallel == FALSE){
      # sequential computations of p-values 
       if(runs >= 10){
            # initializing the progress bar
            progressBar <- txtProgressBar(min = 1,max = runs,style = 3)
    
       }

      for(i in 1:runs){
        p.vals_vec[i] <- FlatLineTest(ExpressionSet)$p.value

        if(runs >= 10){
             # printing out the progress
             setTxtProgressBar(progressBar,i)
        }
      }
    }
    
    plot(p.vals_vec,type = "l" , lwd = 6, ylim = c(0,1), col = "darkblue", xlab = "Runs", ylab = "p-value", main = paste0("runs = ",runs))
    abline(h = 0.05, lty = 2, lwd = 3, col = "darkred")
    
  }
  
   pval <- pgamma(real.var,shape = shape,rate = rate,lower.tail = FALSE)
   sd_values <- apply(resMatrix,2,sd)

   return(list(p.value = pval,std.dev = sd_values))
}


ReductiveHourglassTest <- function(ExpressionSet,modules = NULL,permutations=1000, lillie.test = FALSE, plotHistogram=FALSE,parallel=FALSE,runs=10)
{
  
  
  is.ExpressionSet(ExpressionSet)
    
  if(length(modules) != 3)
    stop("Please specify three modules: early, mid, and late to perform the ReductiveHourglassTest.")
  
  if(length(unlist(modules)) != (dim(ExpressionSet)[2] - 2))
    stop("The number of stages classified into the three modules does not match the total number of stages stored in the given ExpressionSet.")
  
  if(is.null(modules))
    stop("Please specify the three modules: early, mid, and late using the argument 'module = list(early = ..., mid = ..., late = ...)'.")
  
  if(any(table(unlist(modules)) > 1))
    stop("Intersecting modules are not defined for the ReductiveHourglassTest.")
  
  nCols <- dim(ExpressionSet)[2]
  score_vector <- vector(mode = "numeric",length = permutations)
  resMatrix <- matrix(NA_real_, permutations,(nCols-2))
  real_age <- vector(mode = "numeric",length = nCols-2)
  real_age <- TAI(ExpressionSet)
  ### compute the real reductive hourglass scores of the observed phylotranscriptomics pattern
  real_score <- gpScore(real_age,early = modules[[1]],mid = modules[[2]],late = modules[[3]],
                        method = "min",scoringMethod = "mean-mean")
  
  ### compute the bootstrap matrix 
  resMatrix <- bootMatrix(ExpressionSet, permutations)
  
  ### compute the global phylotranscriptomics destruction scores foe each sampled age vector
  score_vector <- apply(resMatrix, 1 ,gpScore,early = modules[[1]],mid = modules[[2]],late = modules[[3]],method = "min",scoringMethod = "mean-mean")
  
  
  
  # parameter estimators using MASS::fitdistr
  param <- fitdistrplus::fitdist(score_vector,"norm", method = "mme")
  mu <- param$estimate[1]
  sigma <- param$estimate[2]
  #mu <- mean(score_vector)
  #sigma <- sd(score_vector)
  
  
    
  if(plotHistogram == TRUE){
     # plot histogram of scores
     normDensity <- function(x){
      
       return(dnorm(x,mu,sigma))
    
     }
    
    if(lillie.test == TRUE)
       par(mfrow = c(2,2))
    if(lillie.test == FALSE)
       par(mfrow = c(1,2))
    
    fitdistrplus::descdist(score_vector, boot = permutations)
    curve(normDensity,xlim = c(min(score_vector),max(score_vector)),col = "steelblue",lwd = 5,xlab = "Scores",ylab="Frequency")
    hist(score_vector,prob = TRUE,add = TRUE, breaks = permutations / (0.01 * permutations))
    rug(score_vector)
    legend("topleft", legend = "A", bty = "n")
    
    p.vals_vec <- vector(mode = "numeric", length = runs)
    lillie_vec <- vector(mode = "logical", length = runs)
    rht <- vector(mode = "list", length = 3)


    if(parallel == TRUE){

      # parallellizing the sampling process using the 'doMC' and 'parallel' package
      # register all given cores for parallelization
      # detectCores(all.tests = TRUE, logical = FALSE) returns the number of cores available on a multi-core machine
      cores <- parallel::makeForkCluster(parallel::detectCores(all.tests = FALSE, logical = FALSE))
      doParallel::registerDoParallel(cores)
      
      # perform the sampling process in parallel
      parallel_results <- as.data.frame(foreach::foreach(i = 1:runs,.combine = "rbind") %dopar% {
        
        if(lillie.test == TRUE)
            ReductiveHourglassTest(ExpressionSet = ExpressionSet,permutations = permutations,lillie.test = TRUE, plotHistogram = FALSE,modules = list(early = modules[[1]],mid = modules[[2]],late = modules[[3]]))[c(1,3)]
        if(lillie.test == FALSE)
          ReductiveHourglassTest(ExpressionSet = ExpressionSet,permutations = permutations,lillie.test = FALSE, plotHistogram = FALSE,modules = list(early = modules[[1]],mid = modules[[2]],late = modules[[3]]))[c(1,3)]
        
        })
    
      # close the cluster connection
      # The is important to be able to re-run the function N times
      # without getting cluster connection problems
      stopCluster(cores)
     
     #print(parallel_results$p.value[1])
     #print(parallel_results$lillie.test[1])
     p.vals_vec <- parallel_results$p.value[1]
     
     if(lillie.test == TRUE)
        lillie_vec <- parallel_results$lillie.test[[1]]

    }
    
    if(parallel == FALSE){

      # sequential computations of p-values 
       if(runs >= 10){
            # initializing the progress bar
            progressBar <- txtProgressBar(min = 1,max = runs,style = 3)
    
       }

        for(i in 1:runs){
            if(lillie.test == TRUE)
                rht <- ReductiveHourglassTest(ExpressionSet = ExpressionSet,permutations = permutations,lillie.test = TRUE, plotHistogram = FALSE,modules = list(early = modules[[1]],mid = modules[[2]],late = modules[[3]]),runs=NULL)
            if(lillie.test == FALSE)
                rht <- ReductiveHourglassTest(ExpressionSet = ExpressionSet,permutations = permutations,lillie.test = FALSE, plotHistogram = FALSE,modules = list(early = modules[[1]],mid = modules[[2]],late = modules[[3]]),runs=NULL)
            
            p.vals_vec[i] <- rht$p.value
            if(lillie.test == TRUE)
                lillie_vec[i] <- rht$lillie.test

            if(runs >= 10){
                 # printing out the progress
                 setTxtProgressBar(progressBar,i)
            }
        }
    }

    plot(p.vals_vec,type = "l" , lwd = 6, ylim = c(0,1), col = "darkblue", xlab = "Runs", ylab = "p-value")
    abline(h = 0.05, lty = 2, lwd = 3)
    legend("topleft", legend = "B", bty = "n")
    
    if(lillie.test == TRUE){
        tbl <- table(factor(lillie_vec, levels = c("FALSE","TRUE")))
        barplot(tbl/sum(tbl) , beside = TRUE, names.arg = c("FALSE", "TRUE"), ylab = "relative frequency", main = paste0("runs = ",runs))
        legend("topleft", legend = "C", bty = "n")
    }
  }
    
    
    #if(real_score >= 0)
    pval <- pnorm(real_score,mean = mu,sd = sigma,lower.tail = FALSE)
    
    #if(real_score < 0)
    #pval <- pnorm(real_score,mean=mu,sd=sigma,lower.tail=TRUE)
    ### computing the standard deviation of the sampled TAI values for each stage separately
    sd_vals <- vector(mode = "numeric",length = nCols-2)
    sd_vals <- apply(resMatrix,2,sd)
    
    if(lillie.test == TRUE){
        # perform Lilliefors K-S-Test
        lillie_p.val <- nortest::lillie.test(score_vector)$p.value
        # does the Lilliefors test pass the criterion
        lillie_bool <- (lillie_p.val > 0.05)
  
        if((lillie_p.val < 0.05) & (plotHistogram == FALSE)){
            warning("Lilliefors (Kolmogorov-Smirnov) test for normality did not pass the p > 0.05 criterion!")
        }
    }
   
    if(lillie.test == TRUE)
        return(list(p.value = pval,std.dev = sd_vals,lillie.test = lillie_bool))
    if(lillie.test == FALSE)
        return(list(p.value = pval,std.dev = sd_vals,lillie.test = NA))
}


#' @title Function to compute the statistical significance of each replicate combination.
#' @description In case a PhyloExpressionSet or DivergenceExpressionSet stores replicates for each
#' developmental stage or experiment, this function allows to compute the p-values quantifying
#' the statistical significance of the underlying pattern for all combinations of replicates.
#' 
#' The intention of this analysis is to validate that there exists no sequence of replicates 
#' (for all possible combination of replicates) that results in a non-significant pattern,
#' when the initial pattern with combined replicates was shown to be significant.
#' 
#' A small Example: 
#' 
#'      
#' Assume PhyloExpressionSet stores 2 developmental stages with 3 replicates measured for each stage.
#' The 6 replicates in total are denoted as: 1.1, 1.2, 1.3, 2.1, 2.2, 2.3. Now the function computes the
#' statistical significance of each pattern derived by the corresponding combination of replicates, e.g.
#'
#' 1.1 + 2.1 -> p-value for combination 1
#'
#' 1.1 + 2.2 -> p-value for combination 2
#'
#' 1.1 + 2.3 -> p-value for combination 3
#'
#' 1.2 + 2.1 -> p-value for combination 4
#'
#' 1.2 + 2.2 -> p-value for combination 5
#'
#' 1.2 + 2.3 -> p-value for combination 6
#'
#' 1.3 + 2.1 -> p-value for combination 7
#'
#' 1.3 + 2.2 -> p-value for combination 8
#'
#' 1.3 + 2.3 -> p-value for combination 9
#' 
#' This procedure yields 9 p-values for the \eqn{2^3} (\eqn{#stages^#replicates}) replicate combinations.
#' 
#' Note, that in case you have a large amount of stages/experiments and a large amount of replicates
#' the computation time will increase by \eqn{#stages^#replicates}. For 11 stages and 4 replicates, 4^11 = 4194304 p-values have to be computed. Each p-value computation itself is based on a permutation test running with 1000 or more permutations. Be aware that this might take some time.
#'
#' The p-value vector returned by this function can then be used to plot the p-values to see
#' whether an critical value \eqn{\alpha} is exeeded or not (e.g. \eqn{\alpha = 0.05}).
#' 
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param replicates a numeric vector storing the number of replicates within each developmental stage or experiment.
#' In case replicate stores only one value, then the function assumes that each developmental stage or experiment
#' stores the same number of replicates.
#' @param TestStatistic a string defining the type of test statistics to be used to quantify the statistical significance the present phylotranscriptomics pattern.
#' Possible values can be: \code{TestStatistic} = "FlatLineTest" : Statistical test for the deviation from a flat line.
#' \code{TestStatistic} = "ReductiveHourglassTest" : Statistical test for the existence of a hourglass shape (high-low-high pattern).
#' @param permutations a numeric value specifying the number of permutations to be performed for the \code{\link{FlatLineTest}} or \code{\link{ReductiveHourglassTest}}.
#' @param parallel a boolean value specifying whether parallel processing (multicore processing) shall be performed.
#' @details The function receives a standard PhyloExpressionSet or DivergenceExpressionSet object and a vector storing the number of replicates present in each stage or experiment. Based on these arguments the function computes all possible replicate combinations using the \code{\link{expand.grid}} function and performs a permutation test (either a \code{\link{FlatLineTest}} or \code{\link{ReductiveHourglassTest}}) for each replicate combination. The \emph{permutation} parameter of this function specifies the number of permutations that shall be performed for each permutation test. When all p-values are computed, a numeric vector storing the corresponding p-values for each replicate combination is returned.  
#' 
#' In other words, for each replicate combination present in the PhyloExpressionSet or DivergenceExpressionSet object, the TAI or TDI pattern of the corresponding replicate combination is tested for its statistical significance based on the underlying test statistic.
#' 
#' This function is also able to perform all computations in parallel using multicore processing. The underlying statistical tests are written in C++ and optimized for fast computations.
#' 
#' @return a numeric vector storing the p-values returned by the underlying test statistic for all possible replicate combinations.
#' @references Drost et al. 2014, Active maintenance of phylotranscriptomic hourglass patterns in animal and plant embryogenesis. MBE. 2014
#' @author Hajk-Georg Drost
#' @seealso \code{\link{expand.grid}}, \code{\link{FlatLineTest}}, \code{\link{ReductiveHourglassTest}}
#' @examples \dontrun{
#' 
#' # load a standard PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'
#' # we assume that the PhyloExpressionSetExample consists of 3 developmental stages 
#' # and 2 replicates for stage 1, 3 replicates for stage 2, and 2 replicates for stage 3
#' p.vector <- combinatorialSignificance(PhyloExpressionSetExample, replicates = c(2,3,2), TestStatistic = "FlatLineTest", permutations = 1000, parallel = FALSE)
#'
#'
#' # now we assume that the PhyloExpressionSetExample consists of 3 developmental stages 
#' # and 2 replicates for each stage
#' # in this case typing replicates = 2 is enough to allow the function to assume that
#' # each stage has 2 replicates
# here we also compute the p-values using multicore processing
#' p.vector <- combinatorialSignificance(PhyloExpressionSetExample[ , 1:8], replicates = 2, TestStatistic = "FlatLineTest", permutations = 1000, parallel = TRUE)
#'
#'
#' }
#' @export
combinatorialSignificance <- function(ExpressionSet,replicates,TestStatistic = "FlatLineTest", permutations = 1000, parallel = FALSE)
{
  
  is.ExpressionSet(ExpressionSet)
  
  if(!is.element(TestStatistic, c("FlatLineTest","ReductiveHourglassTest"))){
    stop("Please enter a correct string for the test statistic: 'FlatLineTest' or 'ReductiveHourglassTest'.")
  }
  
  ncols <- dim(ExpressionSet)[2]
  
  if(length(replicates) == 1){
  
     if((ncols - 2) %% replicates != 0)
        stop("The number of stages and the number of replicates do not match.")
  
     nStages <- (ncols - 2) / replicates
     replicateName.List <- lapply(lapply(1:nStages,rep,times = replicates),paste0,paste0(".",1:replicates))
     stageNames <- as.vector(unlist(replicateName.List))
  
     colnames(ExpressionSet)[3:ncols] <- stageNames
  
     # compute all possible combinations
     combinatorialMatrix <- expand.grid(replicateName.List, stringsAsFactors = FALSE)
  }
  
  if(length(replicates) > 1){
    
    nStages <- length(replicates)
    
    if(sum(replicates) != (ncols - 2))
      stop("The number of stages and the number of replicates do not match.")
    
    f <- function(x){ 
      unlist(lapply(lapply(x,rep,times = replicates[x]),paste0,paste0(".",1:replicates[x])))
    }
    
    
    replicateName.List <- lapply(1:nStages,f)
    stageNames <- as.vector(unlist(replicateName.List))
    colnames(ExpressionSet)[3:ncols] <- stageNames
    
    combinatorialMatrix <- expand.grid(replicateName.List, stringsAsFactors = FALSE)
    
  }
  
  nCombinations <- dim(combinatorialMatrix)[1]
  p.vals <- vector(mode = "numeric",length = nCombinations)
  first_cols_names <- as.character(colnames(ExpressionSet)[1:2])
  
  if(parallel == TRUE){
    # parallellizing the sampling process using the 'doMC' and 'parallel' package
    # register all given cores for parallelization
    # detectCores(all.tests = TRUE, logical = FALSE) returns the number of cores available on a multi-core machine
    cores <- makeForkCluster(detectCores(all.tests = FALSE, logical = FALSE))
    registerDoParallel(cores)
    
    # perform the sampling process in parallel
    p.vals <- as.vector(foreach(i = 1:nCombinations,.combine = "c") %dopar% {
      
      FlatLineTest(as.data.frame(ExpressionSet[c(first_cols_names,as.character(combinatorialMatrix[i , ]))]), permutations = permutations)$p.value
      
    })
    
    # close the cluster connection
    # The is important to be able to re-run the function N times
    # without getting cluster connection problems
    stopCluster(cores)
  }
  
  
  if(parallel == FALSE){
    # sequential computations of p-values 
    if(nCombinations > 10){
      # initializing the progress bar
      progressBar <- txtProgressBar(min = 1,max = nCombinations,style = 3)
      
    }
    
    for(i in 1:nCombinations){
      
      p.vals[i] <- FlatLineTest(as.data.frame(ExpressionSet[c(first_cols_names,as.character(combinatorialMatrix[i , ]))]), permutations = permutations)$p.value
      
      if(nCombinations > 10){
        # printing out the progress
        setTxtProgressBar(progressBar,i)
      }
      
    }
  }
  
  return(p.vals)
}

