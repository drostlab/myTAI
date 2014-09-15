
### This function is used within the statitical test function: ReductiveHourglassTest
### to perform faster computations on an existing sampled TXI matrix
gpScore <- function(age_vals,early,mid,late,method,scoringMethod)
{
  
  Score.Early <- 1
  Score.Late <- 1
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


### Statistical Test for the deviation of a given
### Hourglass-Pattern from a flat line
FlatLineTest <- function(ExpressionSet,permutations=1000,plotHistogram=FALSE,parallel=FALSE,runs=10)
{
  
  # require(MASS, quietly = TRUE)

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


# compute all combinations of replicates of an expressionset
# and compute the significance of a non-flat line TAI/TDI battern for
# each replicate combination

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

