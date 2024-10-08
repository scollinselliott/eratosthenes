% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eratosthenes.R
\name{seq_adj}
\alias{seq_adj}
\alias{seq_adj.character}
\title{Adjust Sequence to Target}
\usage{
seq_adj(input, target)

\method{seq_adj}{character}(input, target)
}
\arguments{
\item{input}{A vector of elements in a sequence.}

\item{target}{A vector of elements in a sequence, containing at least three of the same elements as \code{input}.}
}
\value{
A vector of the adjusted sequence.
}
\description{
Given an "input" sequence of elements and another "target" seqeunce that contains fewer elements in a different order, shift the order of the input sequence to match that of the target, keeping all other elements as proximate to one another as possible. This adjusted ranking is accomplished using piecewise linear interpolation between joint elements ranks. That is, joint rankings are plotted, with input rankings along the \eqn{x} axis and target rankings on the \eqn{y} axis. Remaining rankings in the input sequence are assigned a ranking of \eqn{y} based on the piecewise linear function between joint rankings. If the rank order of elements in the target are identical to those in the input, the result is identical to the input. A minimum number of three joint elements in both the input and target are required.
}
\examples{
x <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J") # the input sequence
y <- c("D", "A", "J") # the target sequence

seq_adj(x, y)

}
