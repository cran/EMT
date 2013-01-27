\name{plotMultinom}
\alias{plotMultinom}
\title{ Plot the Probability distribution fot the Exact Multinomial Test}
\description{
  This function takes the results of \code{multinomial.test} as input and plots the calculated probability
  distribution.  
}
\usage{
plotMultinom(listMultinom, showmax = 50)
}
\arguments{
  \item{listMultinom}{ a list created by running the function \code{multinomial.test}. }
  \item{showmax}{ maximum number of bars to show in the plot (to avoid long tails). }
}
\details{
The function \code{multinomial.test} creates an output list that is used as input in \code{plotMultinom} to depict some results.
If the default approach was used, the figure shows the exact probabilities of all possible outcomes of the experiment.
If the Monte Carlo approach was used, the relative frequencies of the outcomes are shown as occurred during the simulated withdrawals.
The probabilities/relative frequencies are shown in descending order from the left to the right.
}
\value{
The first argument (\emph{listMultinom}) is returned without modification. 
} 
\author{ Uwe Menzel   <uwemenzel@gmail.com> }
\note{ 
For better visibility, the parameter \code{showmax} excludes very long right tails from the plot.
However, the default value of \code{showmax} should be incremented to get a significant plot if the number of possible outcomes is big.
}
\seealso{  
The Multinomial Distribution: \code{\link{multinomial.test}}
}
\examples{

## Load the EMT package:
library(EMT)


## input and calculation of p-values:
observed <- c(5,2,1)   		
prob <- c(0.25, 0.5, 0.25) 	
out <- multinomial.test(observed, prob)        


## Plot the probability distribution:
plotMultinom(out)
plotMultinom(out, showmax = 30)  # suppress part of the tail in the plot 

}
\keyword{ htest }
\keyword{ multivariate }



