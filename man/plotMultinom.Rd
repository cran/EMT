\name{plotMultinom}
\alias{plotMultinom}
\title{ Plot the Probability distribution fot the Exact Multinomial Test}
\description{
  This function takes the results of \code{multinomial.test} as input and plots the calculated probability
  distribution.  
}
\usage{
plotMultinom(listMultinom)
}
\arguments{
  \item{listMultinom}{ a list created by running the function \code{multinomial.test}. }
}
\details{
The function \code{plotMultinom} displays a barplot of the probabilities for the individual events.
The probabilities are shown in descending order from the left to the right.
Events contributing to the p-value are marked red.
Plots are only made if the number of different events is lower than or equal to 100 and for low number of trials in Monte Carlo simulations.
}
\value{
The first argument (\emph{listMultinom}) is returned without modification. 
} 
\author{ Uwe Menzel   <uwemenzel@gmail.com> }

\seealso{  
The Multinomial Distribution: \code{\link{multinomial.test}}
}
\examples{

## Load the EMT package:
library(EMT)

## input and calculation of p-values:
observed <- c(5,2,1)   		
prob <- c(0.25, 0.5, 0.25) 
out <- multinomial.test(observed, prob)     # p.value = 0.0767

## Plot the probability distribution:
plotMultinom(out)

}
\keyword{ htest }
\keyword{ multivariate }



