.packageName <- "EMT"



"multinomial.test" <- 
function(observed, prob, useChisq = FALSE, MonteCarlo = FALSE, ntrial = 1e6) 
{
    if(!is.vector(observed, mode = "numeric")) stop(" Observations have to be stored in a vector, e.g. 'observed <- c(5,2,1)'")
    if(!is.vector(prob, mode = "numeric")) stop(" Probabilities have to be stored in a vector, e.g. 'prob <- c(0.25, 0.5, 0.25)'")
    if(round(sum(prob), digits = 1) != 1) stop("Wrong input: sum of probabilities must not deviate from 1.")
    if(length(observed) != length(prob)) stop(" Observations and probabilities must have equal dimensions.")

    size = sum(observed)
    groups = length(observed)
    numEvents = choose(size + groups - 1, groups - 1)  
    
    cat(paste("\n The model includes", numEvents, "different events.\n\n"))
    if(ntrial < 10*numEvents) cat(" The chosen number of trials is rather low, should be at least 10 times the numver of events.\n\n")

    if ( MonteCarlo == FALSE ) {
        if (useChisq == FALSE) {
            res <- ExactMultinomialTest(observed, prob, size, groups, numEvents)
        } else {
            res <- ExactMultinomialTestChisquare(observed, prob, size, groups, numEvents)
        }
    } else {
        if ( ntrial < numEvents ) {
            cat(" \n WARNING: Number of withdrawels is lower than the number of possible outcomes. 
                This might yield unreliable results!\n\n")}
            flush.console()
        
        if (useChisq == FALSE) {
            res <- MonteCarloMultinomialTest(observed, prob, size, groups, numEvents, ntrial)
        } else {
            res <- MonteCarloMultinomialTestChisquare(observed, prob, size, groups, numEvents, ntrial)
        }
    }
    invisible(res)
} 



"plotMultinom" <- 
function(listMultinom) 
{
    if(!is.list(listMultinom)) stop(" First argument must be a list (output of function 'multinomial.test').")
    if(is.null(listMultinom$allProb)) stop(" A barplot of probabilities cannot be shown for Monte Carlo methods.\n\n")
    if(listMultinom$numEvents > 100) stop(" A barplot is not made when the number of events is higher than 100.\n\n")
        
    h <- listMultinom$allProb
    cols = rep("blue", length(h))
    cols[h <= listMultinom$criticalValue] = "red"
        
    barplot(h, main = "Probability vs. Events", xlab = "", ylab = "", space = 1, font.main = 1, las = 2, col = cols)
    mtext("Events (sorted)", side = 1)
    mtext(paste("p.value =",listMultinom$p.value), side = 3, col = "blue", cex = 0.9)
    
    if (sum(grep("Carlo", listMultinom$id))) mtext(paste(" Trials: ", listMultinom$ntrial), side = 4, col = "blue", cex = 0.9)  
    invisible(listMultinom)  			
}





"ExactMultinomialTest" <- 
function(observed, prob, size, groups, numEvents) 
{
    pObs = dmultinom(observed, size = size, prob) 	
    eventMat <- findVectors(groups, size)    		
    if( nrow(eventMat) != numEvents ) stop("Wrong number of events calculated. \n This is probably a bug.")

    eventProb <- apply(eventMat, 1, function(x) dmultinom(x, size = size, prob = prob)) 
    eventProb[abs(eventProb - pObs) < .Machine$double.eps^0.5] <- pObs
    p.value = sum(eventProb[eventProb <= pObs])

    if(round(sum(eventProb), digits = 2) != 1) stop("Wrong values for probabilities. \n This is probably a bug.")

    head <- paste("\n Exact Multinomial Test\n\n")
    tab <- as.data.frame(cbind(numEvents, round(pObs, digits = 4), round(p.value, digits = 4)))
    colnames(tab) <- c("   Events","   pObs","   p.value")
    cat(head); print(tab, row.names = FALSE)

    invisible(list(id = "Exact Multinomial Test", size = size, groups = groups, numEvents = numEvents, 
         stat = "lowP", allProb = sort(eventProb, decreasing = TRUE), criticalValue = pObs,
         ntrial = NULL, p.value = round(p.value, digits = 4)))
} 




"ExactMultinomialTestChisquare" <- 
function(observed, prob, size, groups, numEvents) 
{
    expectedFreq  = size * prob                     	
    chi2Obs = chisqStat(observed, expectedFreq)  
    pObs = dmultinom(observed, size = size, prob) 	
    
    eventMat <- findVectors(groups, size) 			
    if( nrow(eventMat) != numEvents ) stop("Wrong number of events calculated. \n This is probably a bug.")
    eventProb <- apply(eventMat, 1, function(x) dmultinom(x, size = size, prob = prob)) 
    eventChi2 <- apply(eventMat, 1, function(x) chisqStat(x, expectedFreq))   			
    eventProb[abs(eventProb - pObs)    < .Machine$double.eps^0.5] <- pObs
    eventChi2[abs(eventChi2 - chi2Obs) < .Machine$double.eps^0.5] <- chi2Obs
    eventPandChi2 <- cbind(eventProb, eventChi2)

    if(round(sum(eventProb), digits = 2) != 1) stop("Wrong values for probabilities. \n This is probably a bug.")
    p.value <- sum((eventPandChi2[eventPandChi2[,2] >= chi2Obs,])[,1])
    cV = max(eventPandChi2[eventPandChi2[,2] >= chi2Obs,][,1])  

    head <- paste("\n Exact Multinomial Test, Chisquare\n\n")
    tab <- as.data.frame(cbind(numEvents, round(chi2Obs, digits = 4), round(p.value, digits = 4)))
    colnames(tab) <- c("   Events","   chi2Obs","   p.value")
    cat(head); print(tab, row.names = FALSE)

    invisible(list(id = "Exact Multinomial Test, Chisquare", size = size, groups = groups, numEvents = numEvents,
         stat = "highChisq", allProb = sort(eventProb, decreasing = TRUE), criticalValue = cV, 
         ntrial = NULL, p.value = round(p.value, digits = 4)))
} 


"MonteCarloMultinomialTest" <- 
function(observed, prob, size, groups, numEvents, ntrial) 
{
    pObs = dmultinom(observed, size = size, prob)
    lowPs = 0
    for (i in 1:ntrial) {
      res <- rmultinom(n = 1, size = size, prob = prob)[,1]   
      pDraw = dmultinom(res, size = size, prob)  
      if(abs(pDraw - pObs) < .Machine$double.eps^0.5) pDraw = pObs
      if(pDraw <= pObs) lowPs <- lowPs + 1
      if(i%%100000 == 0) cat(" Number of withdrawals: ", prettyNum(i, scientific = FALSE, big.mark = ","),"\n")
    }
    p.value = (lowPs + 1)/(ntrial + 1)  
    
    head <- paste("\n Monte Carlo Multinomial Test\n\n")
    tab <- as.data.frame(cbind(numEvents, round(pObs, digits = 4), round(p.value, digits = 4)))
    colnames(tab) <- c("   Events","   pObs","   p.value")
    cat(head); print(tab, row.names = FALSE)
    
    invisible(list(id = "Monte Carlo Multinomial Test", size = size, groups = groups, numEvents = numEvents,
                   stat = "lowF", allProb = NULL, criticalValue = pObs,
                   ntrial = ntrial, p.value = round(p.value, digits = 4)))    
} 



"MonteCarloMultinomialTestChisquare" <- 
function(observed, prob, size, groups, numEvents, ntrial) 
{ 
    expectedFreq  = size * prob                    			
    chi2Obs = chisqStat(observed, expectedFreq)
    
    bigChis = 0
    for (i in 1:ntrial) {
      res <- rmultinom(n = 1, size = size, prob = prob)[,1]   
      chi2Draw = chisqStat(res, expectedFreq)  
      if(abs(chi2Draw - chi2Obs) < .Machine$double.eps^0.5) chi2Draw = chi2Obs
      if(chi2Draw >= chi2Obs) bigChis <- bigChis + 1
      if(i%%100000 == 0) cat(" Number of withdrawals: ", prettyNum(i, scientific = FALSE, big.mark = ","),"\n")
    } 
    p.value = (bigChis + 1)/(ntrial + 1)
    
    head <- paste("\n Monte Carlo Multinomial Test, Chisquare\n\n")
    tab <- as.data.frame(cbind(numEvents, round(chi2Obs, digits = 4), round(p.value, digits = 4)))
    colnames(tab) <- c("   Events","   chi2Obs","   p.value")
    cat(head); print(tab, row.names = FALSE)
    
    invisible(list(id = "Monte Carlo Multinomial Test, Chisquare", size = size, groups = groups, numEvents = numEvents, 
                   stat = "highChi2", allProb = NULL, criticalValue = chi2Obs,
                   ntrial = ntrial, p.value = round(p.value, digits = 4)))			   
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
## 2021 uwe.menzel@matstat.org  uwemenzel@gmail.com

