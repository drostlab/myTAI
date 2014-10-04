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

