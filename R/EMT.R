.packageName <- "EMT"



"multinomial.test" <- 
function(observed, prob, useChisq = FALSE, MonteCarlo = FALSE, ntrial = 1e6, atOnce = 1e6) 
{
    if(!is.vector(observed, mode = "numeric")) stop(" Observations have to be stored in a vector, e.g. 'observed <- c(5,2,1)'")
    if(!is.vector(prob, mode = "numeric")) stop(" Probabilities have to be stored in a vector, e.g. 'prob <- c(0.25, 0.5, 0.25)'")
    if(round(sum(prob), digits = 3) != 1) stop("Wrong input: sum of probabilities must not deviate from 1.")
    if(length(observed) != length(prob)) stop(" Observations and probabilities must have equal dimensions.")

    size = sum(observed)
    groups = length(observed)
    numEvents = choose(size + groups - 1, groups - 1)  
    
    cat(paste("\n The model includes", numEvents, "different events.\n\n"))

    if ( MonteCarlo == FALSE ) {
        if (useChisq == FALSE) {
            res <- ExactMultinomialTest(observed, prob, size, groups, numEvents)
        } else {
            res <- ExactMultinomialTestChisquare(observed, prob, size, groups, numEvents)
        }
    } else {
        if(ntrial < 10*numEvents) 
        cat(" The chosen number of trials is rather low, you might consider increasing the parameter 'ntrials'.\n\n")
        flush.console()
        
        if (useChisq == FALSE) {
            res <- MonteCarloMultinomialTest(observed, prob, size, groups, numEvents, ntrial, atOnce)
        } else {
            res <- MonteCarloMultinomialTestChisquare(observed, prob, size, groups, numEvents, ntrial, atOnce)
        }
    }
    invisible(res)
} 




"plotMultinom" <- 
  function(listMultinom) 
  {
    if(!is.list(listMultinom)) stop(" First argument must be a list (output of function 'multinomial.test').")
    if(is.null(listMultinom$allProb) | is.null(listMultinom$criticalValue)) {
      stop(" A barplot of probabilities cannot be shown for this dataset.\n\n")
    }  
    
    epsilon = sqrt(.Machine$double.eps) 
    h <- listMultinom$allProb
    cols = rep("blue", length(h))
    cols[h - listMultinom$criticalValue <= epsilon] = "red"
    
    barplot(h, main = "Probability vs. Events", xlab = "", ylab = "", space = 1, font.main = 1, las = 2, col = cols)
    mtext("Events (sorted)", side = 1)
    mtext(paste("p.value =", signif(listMultinom$p.value,4)), side = 3, col = "blue", cex = 0.9)
    
    if (sum(grep("Carlo", listMultinom$id))) mtext(paste(" Trials: ", listMultinom$ntrial), side = 4, col = "blue", cex = 0.9)  
    invisible(listMultinom)  			
  }




"ExactMultinomialTest" <- 
function(observed, prob, size, groups, numEvents) 
{
    epsilon = sqrt(.Machine$double.eps) 
  
    pObs <- dmultinom(observed, size = size, prob) 	
    eventMat <- findVectors(groups, size)    		
    if( nrow(eventMat) != numEvents ) stop("Wrong number of events calculated. \n This is a bug.")

    eventProb <- apply(eventMat, 1, function(x) dmultinom(x, size = size, prob = prob))
    if(round(sum(eventProb), digits = 3) != 1) stop("Incorrect values for probabilities. \n This is a bug.")
    p.value <- sum(eventProb[eventProb <= pObs + epsilon])
    
    head <- paste("\n Exact Multinomial Test\n\n")
    tab <- as.data.frame(cbind(numEvents, signif(pObs, digits = 6), signif(p.value, digits = 6)))
    colnames(tab) <- c("   Events","   pObs","   p.value")
    cat(head); print(tab, row.names = FALSE)

    invisible(list(id = "Exact Multinomial Test", size = size, groups = groups, numEvents = numEvents, 
         stat = "lowP", allProb = sort(eventProb, decreasing = TRUE), criticalValue = pObs,
         ntrial = NULL, p.value = p.value))
} 




"ExactMultinomialTestChisquare" <- 
  function(observed, prob, size, groups, numEvents) 
  {
    epsilon <- sqrt(.Machine$double.eps)
    
    expectedFreq  <- size * prob                     	
    chi2Obs <- chisqStat(observed, expectedFreq)  
    pObs <- dmultinom(observed, size = size, prob) 	
    
    eventMat <- findVectors(groups, size) 			
    if( nrow(eventMat) != numEvents ) stop("Incorrect number of events calculated. \n This is a bug.")
    
    eventProb <- apply(eventMat, 1, function(x) dmultinom(x, size = size, prob = prob)) 
    eventChi2 <- apply(eventMat, 1, function(x) chisqStat(x, expectedFreq)) 
    
    if(round(sum(eventProb), digits = 3) != 1) stop("Incorrect values for probabilities. \n This is a bug.")
    eventPandChi2 <- cbind(eventProb, eventChi2)
    
    p.value <- sum(eventPandChi2[eventPandChi2[,2] + epsilon >= chi2Obs, 1])
    
    criticalValue <- eventPandChi2[eventPandChi2[,2] == chi2Obs, 1][1]; attr(criticalValue, "names") <- NULL
    
    head <- paste("\n Exact Multinomial Test, chiSquared\n\n")
    tab <- as.data.frame(cbind(numEvents, signif(chi2Obs, digits = 6), signif(p.value, digits = 6)))
    colnames(tab) <- c("   Events","   chi2Obs","   p.value")
    cat(head); print(tab, row.names = FALSE)
    
    invisible(list(id = "Exact Multinomial Test, chiSquared", size = size, groups = groups, numEvents = numEvents,
                   stat = "highChisq", allProb = sort(eventProb, decreasing = TRUE), criticalValue = criticalValue, 
                   ntrial = NULL, p.value = p.value))
  } 




"MonteCarloMultinomialTest" <- 
  function(observed, prob, size, groups, numEvents, ntrial, atOnce) 
  {
    IDofObs <- paste0(observed, collapse = "")
    sumLowFreq = 0
    loops = floor(ntrial/atOnce)
    if(ntrial%%atOnce > 0) loops = loops + 1    
    
    for (i in 1:loops) {
      res <- rmultinom(n = atOnce, size = size, prob = prob) 			
      eventID <- apply(res, 2, function(x) paste0(x, collapse = ""))    		
      frequencyTable <- table(eventID)    							    				
      if (!(IDofObs %in% rownames(frequencyTable))) {
        freqObs = 0   										  
      } else {
        freqObs <- frequencyTable[rownames(frequencyTable) == IDofObs][[1]]  	
      } 
      sumLowFreq  <- sumLowFreq + sum(frequencyTable[as.vector(frequencyTable) <= freqObs])
      cat(" Number of withdrawals accomplished: ", prettyNum(i*atOnce, scientific = FALSE, big.mark = ","),"\n")
      flush.console() 
    }
    cat("\n Number of withdrawals less frequent than the observation: ", sumLowFreq,"\n\n")
    
    p.value = (sumLowFreq + 1)/(ntrial + 1)   
    
    if(numEvents > 100 | loops > 1) {
      allProb = NULL 
      criticalValue = NULL
    } else {
      allProb = sort(as.vector(frequencyTable), decreasing = TRUE)/ntrial  
      criticalValue = frequencyTable[IDofObs]/ntrial; attr(criticalValue, "names") = NULL  
    }
    
    head <- paste("Monte Carlo Multinomial Test, distance measure: frequency\n\n")
    tab <- as.data.frame(cbind(numEvents, signif(freqObs/ntrial, digits = 6), signif(p.value, digits = 6)))
    colnames(tab) <- c("   Events","   fObs","   p.value")
    cat(head); print(tab,row.names = FALSE)
    
    invisible(list(id = "Monte Carlo Multinomial Test", size = size, groups = groups, numEvents = numEvents,
                   stat = "lowF", allProb = allProb, criticalValue = criticalValue, ntrial = ntrial, p.value = p.value))
  } 




"MonteCarloMultinomialTestChisquare" <- 
  function(observed, prob, size, groups, numEvents, ntrial, atOnce) 
  { 
    epsilon = sqrt(.Machine$double.eps)   
    expectedFreq  = size * prob                     			
    chi2Obs = chisqStat(observed, expectedFreq)
    bigChis = 0
    loops = floor(ntrial/atOnce)
    if(ntrial%%atOnce > 0) loops = loops + 1    
    
    for (i in 1:loops) {
      res <- rmultinom(n = atOnce, size = size, prob = prob) 		
      chi2all <- apply(res, 2, function(x) chisqStat(x, expectedFreq))	
      bigChis <- bigChis + sum(chi2all >= chi2Obs - epsilon) 
      cat(" Number of withdrawals accomplished: ", prettyNum(i*atOnce, scientific = FALSE, big.mark = ","),"\n")
      flush.console()
    }
    cat("\n Number of withdrawals with chiSquared values bigger than the observation: ", bigChis,"\n\n")
    
    p.value = (bigChis + 1)/(ntrial + 1) 
    
    if(numEvents > 100 | loops > 1) {
      allProb = NULL  
      criticalValue = NULL
    } else {
      frequencyTable <- table(chi2all)
      allProb = sort(as.vector(frequencyTable)/ntrial, decreasing = TRUE) 
      criticalValue <- frequencyTable[as.character(chi2Obs)]/ntrial; attr(criticalValue, "names") = NULL
    }

    head <- paste("\n Monte Carlo Multinomial Test, distance measure: chiSquared\n\n")
    tab <- as.data.frame(cbind(numEvents, signif(chi2Obs, digits = 6), signif(p.value, digits = 6)))
    colnames(tab) <- c("   Events","   chi2Obs","   p.value")
    cat(head); print(tab,row.names = FALSE)
    
    invisible(list(id = "Monte Carlo Multinomial Test, chiSquared", size = size, groups = groups, numEvents = numEvents,
                   stat = "chi2", allProb = allProb, criticalValue = criticalValue, ntrial = ntrial, p.value = p.value))			   
  } 



"chisqStat" <- 
function(observed, expected) 
{
  chisq = sum((observed - expected)^2/expected)
  invisible(chisq)   	
}



"findVectors" <- 			
function(groups, size) 
{  
    if (groups == 1) {
        mat = size    
    } else { 
        mat <- matrix(rep(0, groups - 1), nrow = 1)  
        for (i in 1:size) {
            mat <- rbind(mat, findVectors(groups - 1, i))  
        } 
        mat <- cbind(mat, size - rowSums(mat))   
    } 
    invisible(mat)
}



## 2010; uwe.menzel@math.uu.se   uwemenzel@gmail.com
## 2021; uwe.menzel@matstat.org  uwemenzel@gmail.com
## 2023; uwe.menzel@matstat.org  uwemenzel@gmail.com



