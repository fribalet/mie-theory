% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mietheory.R
\name{ALegendr}
\alias{ALegendr}
\title{Calculate associated Legendre polynomial functions}
\usage{
ALegendr(ang, nmax)
}
\arguments{
\item{ang}{Numeric vector of scattering angles (in radians).}

\item{nmax}{Integer specifying the maximum order of the functions.}
}
\value{
A list containing:
 \itemize{
   \item `p`: A matrix of dimensions `nmax` by `length(ang)`, where each column 
         represents the `p` function values for a specific angle.
   \item `t`: A matrix of dimensions `nmax` by `length(ang)`, where each column 
         represents the `t` function values for a specific angle.
}
}
\description{
This function calculates the required functions of the associated Legendre 
polynomials, which are needed to determine the scattering angle dependence 
of scattered intensity. It computes the functions from order 1 to `nmax` 
for the given scattering angle(s) (`ang`) using upwards recursion.
}
\examples{
ALegendr(ang = seq(0, pi, length.out = 10), nmax = 5)
}
