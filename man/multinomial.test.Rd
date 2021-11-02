\name{multinomial.test}
\alias{multinomial.test}
\title{ Exact Multinomial Test: Goodness-of-Fit Test for Discrete Multivariate Data}
\description{
  Goodness-of-fit tests for discrete multivariate data. 
  It is tested if a given observation is likely to have occurred under the assumption of an ab-initio model. 
  Monte Carlo methods are provided to make the function capable of solving high-dimensional problems.  
}
\usage{
multinomial.test(observed, prob, useChisq = FALSE,
               MonteCarlo = FALSE, ntrial = 1e6)
}
\arguments{
  \item{observed}{ vector describing the observation: contains the \emph{observed numbers} of items in each category. }
  \item{prob}{ vector describing the model: contains the \emph{hypothetical probabilities} corresponding to each category. }
  \item{useChisq}{ if \code{TRUE}, Pearson's chisquare is used as a distance measure between observed and expected frequencies. }
  \item{MonteCarlo}{ if \code{TRUE}, a Monte Carlo approach is used. }
  \item{ntrial}{ number of simulated samples in the Monte Carlo approach. }
}
\details{
The Exact Multinomial Test is a Goodness-of-fit test for discrete multivariate data. 
It is tested if a given observation is likely to have occurred under the assumption of an ab-initio model. 
In the experimental setup belonging to the test, \emph{n} items fall into \emph{k} categories with certain probabilities
(sample size \emph{n} with \emph{k} categories). 
The \bold{observation}, described by the vector \code{observed}, indicates how many items have been observed in each category.  
The \bold{model}, determined by the vector \code{prob}, assigns to each category the hypothetical probability that an item falls into it. 
Now, if the observation is unlikely to have occurred under the assumption of the model, it is advisible to
regard the model as \emph{not} valid. The p-value estimates how likely the observation is, given the model.   
In particular, low p-values suggest that the model is \emph{not} valid.
The \bold{default approach} used by \code{multinomial.test} obtains the p-values by 
calculating the exact probabilities of \emph{all} possible outcomes given \code{n} and \code{k},
using the multinomial probability distribution function \code{\link{dmultinom}} provided by R. 
Then, by default, the p-value is obtained by summing the probabilities of all outcomes which are less likely
than the observed outcome (or equally likely as the observed outcome), i.e. by summing all \eqn{p(i) <= p(observed)}
(distance measure based on probabilities).
Alternatively, the p-value can be obtained by summing the probabilities of all outcomes connected with a chisquare no smaller than
the chisquare connected with the actual observation (distance measure based on chisquare). 
The latter is triggered by setting \code{useChisq = TRUE}.
Having a sample of size \emph{n} in an experiment with \emph{k} categories, the number of distinct
possible outcomes is the binomial coefficient \code{choose(n+k-1,k-1)}. This number grows rapidly with increasing parameters \emph{n} and
\emph{k}. If the parameters grow too big, numerical calculation might fail because of time or
memory limitations. 
In this case, usage of a \bold{Monte Carlo approach} provided by \code{multinomial.test} is suggested.
A Monte Carlo approach, activated by setting \code{MonteCarlo = TRUE},
simulates withdrawal of \emph{ntrial} samples of size \emph{n} from the hypothetical distribution specified by the vector \code{prob}. 
The default value for \emph{ntrial} is \code{100000} but might be incremented for big \emph{n} and \emph{k}. 
The advantage of the Monte Carlo approach is that memory requirements and running time are essentially determined by \emph{ntrial} 
but not by \emph{n} or \emph{k}. 
By default, the p-value is then obtained by summing the relative frequencies of occurrence of unusual outcomes, i.e. of
outcomes occurring less frequently than the observed one (or equally frequent as the observed one).    
Alternatively, as above, Pearson's chisquare can be used as a distance measure by setting \code{useChisq = TRUE}.
}
\value{
  \item{id}{textual description of the method used.}
  \item{size}{sample size \emph{n}, equals the sum of the components of the vector \code{observed}.}
  \item{groups}{number of categories \emph{k} in the experiment, equals the number of components of the vector \code{observed}.} 
  \item{numEvents}{number of different events for the model considered.}  
  \item{stat}{textual description of the distance measure used.}
  \item{allProb}{vector containing the probabilities (rel. frequencies for the Monte Carlo approach) of all possible outcomes (might be huge for big \emph{n} and \emph{k}).}
  \item{criticalValue}{the critical value of the hypothesis test.}
  \item{ntrial}{number of trials if the Monte Carlo approach was used, \code{NULL} otherwise.}
  \item{p.value}{the calculated p-value rounded to four significant digits.}  
} 
\references{
  H. Bayo Lawal (2003)
  \emph{Categorical data analysis with SAS and SPSS applications}, Volume 1, Chapter 3
  ISBN: 978-0-8058-4605-8

  Read, T. R. C. and Cressie, N. A. C. (1988). 
  \emph{Goodness-of-fit statistics for discrete multivariate data.} 
  Springer, New York.
}
\author{ Uwe Menzel <uwemenzel@gmail.com>}
\note{ For two categories (\eqn{k = 2}), the test is called Exact Binomial Test.
} 
\seealso{  
The Multinomial Distribution: \code{\link{dmultinom}}
}
\examples{

## Load the EMT package:
library(EMT)

## Input data for a three-dimensional case:
observed <- c(5,2,1)   		 
prob <- c(0.25, 0.5, 0.25) 	 

## Calculate p-value using default options:
out <- multinomial.test(observed, prob)        
# p.value = 0.0767

## Plot the probabilities for each event:
plotMultinom(out)

## Calculate p-value for the same input using Pearson's chisquare:
out <- multinomial.test(observed, prob, useChisq = TRUE)  
# p.value = 0.0596 ; not the same!

## Test the hypothesis that all sides of a dice have the same probabilities:
prob <- rep(1/6, 6)	                 
observed <- c(4, 5, 2, 7, 0, 1)			
out <- multinomial.test(observed, prob) 
# p.value = 0.0357 -> better get another dice!

# the same problem using a Monte Carlo approach:
\dontrun{
  out <- multinomial.test(observed, prob, MonteCarlo = TRUE, ntrial = 5e+6)   
}

}
\keyword{ htest }
\keyword{ multivariate }



