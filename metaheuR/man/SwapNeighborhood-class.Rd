% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/SwapNeighborhood.R
\docType{class}
\name{SwapNeighborhood-class}
\alias{SwapNeighborhood-class}
\title{An S4 class to represent swap neighborhood for permutations}
\description{
An S4 class to represent swap neighborhood for permutations
}
\section{Slots}{

\describe{
\item{\code{base}}{An object of class \code{\linkS4class{Permutation}}, which will be the solution whose neighborhood we will explore}

\item{\code{position.list}}{A permutation indicating the order in which the exploration will be carried out. Each time a neighbohr is created, this list is used to get the first position of the swap movement (the other position will be the one inmediately after that in the base permutation)}

\item{\code{random}}{A logical value indicating whether the exploration is at random or not}

\item{\code{id}}{Numeric value indicating the current position (in the \code{position.list}) that will be used to generate a new neighbor}
}}

