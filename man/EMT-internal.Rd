\name{EMT-internal}
\alias{ExactMultinomialTest}
\alias{ExactMultinomialTestChisquare}
\alias{MonteCarloMultinomialTest}
\alias{MonteCarloMultinomialTestChisquare}
\alias{chisqStat}
\alias{findVectors}
\title{ Internal functions for the EMT package }
\description{
  Internal functions for the EMT package
}
\usage{
ExactMultinomialTest(observed, prob, size, groups, numEvents)
ExactMultinomialTestChisquare(observed, prob, size, groups, numEvents)
MonteCarloMultinomialTest(observed, prob, size, groups, 
                         numEvents, ntrial)
MonteCarloMultinomialTestChisquare(observed, prob, size, groups,
                         numEvents, ntrial)
chisqStat(observed,expected)
findVectors(groups,size)
}
\arguments{
  \item{observed}{ vector describing the observation: contains the \emph{observed numbers} of items in each category. }
  \item{prob}{ vector describing the model: contains the \emph{hypothetical probabilities} corresponding to each category. }
  \item{expected}{vector containing the expected numbers of items in each category under the assumption that the model is valid.}
  \item{size}{ sample size, sum of the components of the vector \code{observed}. }
  \item{groups}{ number of categories in the experiment. }
  \item{numEvents}{ number of possible outcomes of the experiment. }
  \item{ntrial}{ number of simulated samples in the Monte Carlo approach. }
}
\details{
These functions are not intended to be called by the user.
}



