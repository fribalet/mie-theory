#' Calculate Ricatti-Bessel functions of the first kind
#'
#' This function calculates the Ricatti-Bessel functions of the first kind 
#' from order 1 to `nmax` for the given radius parameter(s) `r`, using 
#' downwards recursion.
#'
#' @param r Numeric vector of radius parameters.
#' @param nmax Integer specifying the maximum order of the Ricatti-Bessel function.
#'
#' @return A matrix of dimensions `nmax` by `length(r)`, where each column 
#'   represents the Ricatti-Bessel functions of the first kind for a 
#'   specific radius parameter.
#'
#' @seealso \code{\link{RB2}} for calculating Ricatti-Bessel functions of the second kind.
#' @examples
#' RB1(r = 1:10, nmax = 5)

RB1 <- function(r, nmax) {
  # Ensure r is a column vector
  r <- as.vector(r)
  
  # Compute stopping index nst
  nst <- ceiling(nmax + sqrt(101 + max(r)))
  
  # Initialize phi matrix
  phi <- matrix(0, nrow = nst, ncol = length(r))
  
  # Set initial condition
  phi[nst - 1, ] <- 1e-10
  
  # Compute phi using backward recurrence (descending loop)
  if (nst > 2) {
    for (n in seq(nst - 2, 1, by = -1)) {  # Count down
      phi[n, ] <- ((2 * n + 3) * phi[n + 1, ] / r) - phi[n + 2, ]
    }
  }
  
  # Compute phi0, ensuring no division by zero
  r_safe <- ifelse(r == 0, 1e-10, r)
  phi0 <- (3 * phi[1, ] / r_safe) - phi[2, ]
  phi0 <- sin(r) / phi0
  
  # Extract first nmax rows and multiply by phi0
  phi <- phi[1:nmax, , drop = FALSE] * matrix(rep(phi0, each = nmax), nrow = nmax)
  
  return(phi)
}

#' Calculate Ricatti-Bessel functions of the second kind
#'
#' This function calculates the Ricatti-Bessel functions of the second kind 
#' from order 1 to `nmax` for the given radius parameter(s) `r`, using 
#' upwards recursion.
#'
#' @param r Numeric vector of radius parameters.
#' @param nmax Integer specifying the maximum order of the Ricatti-Bessel function.
#'
#' @return A matrix of dimensions `nmax` by `length(r)`, where each column 
#'   represents the Ricatti-Bessel functions of the second kind for a 
#'   specific radius parameter.
#'
#' @seealso \code{\link{RB1}} for calculating Ricatti-Bessel functions of the first kind.
#'
#' @examples
#' RB2(r = 1:10, nmax = 5)
RB2 <- function(r, nmax) {
  # Ensure r is a vector
  r <- as.vector(r)
  
  # Handle division by zero (avoid r = 0)
  r_safe <- ifelse(r == 0, 1e-10, r)
  
  # Initialize zeta matrix
  zeta <- matrix(0, nrow = nmax, ncol = length(r))
  
  # Compute first two rows of zeta
  zeta[1, ] <- -cos(r) / r_safe - sin(r)
  zeta[2, ] <- 3 * zeta[1, ] / r_safe + cos(r)
  
  # Compute zeta using recurrence relation
  if (nmax > 2) {
    for (n in 3:nmax) {
      zeta[n, ] <- ((2 * n - 1) * zeta[n - 1, ] / r_safe) - zeta[n - 2, ]
    }
  }
  
  return(zeta)
}

#' Calculate scattering coefficients for Mie calculations
#'
#' This function calculates the scattering coefficients required for Mie 
#' calculations (e.g., scattering intensity). It uses the Ricatti-Bessel 
#' functions to compute the coefficients `a` and `b` from order 1 to `nmax` 
#' for the given refractive index (`r`) and particle size parameter(s) (`x`).
#'
#' @param r Numeric vector of refractive indices. Can be a scalar.
#' @param x Numeric vector of particle size parameters. Can be a scalar.
#' @param nmax Integer specifying the maximum order of the 
#'   Ricatti-Bessel function.
#'
#' @return A list containing: 
#' \itemize{
#'  \item `a`: A complex matrix of dimensions `nmax` by `length(r)` (or `length(x)` 
#'          if `x` is a vector), where each column represents the `a` 
#'          scattering coefficients for a specific particle parameter (pair).
#'  \item `b`: A complex matrix of dimensions `nmax` by `length(r)` (or `length(x)` 
#'          if `x` is a vector), where each column represents the `b` 
#'          scattering coefficients for a specific particle parameter (pair).
#' }
#' @details 
#' If both `r` and `x` are vectors, they must have the same length to form 
#' particle parameter pairs. If either `r` or `x` is a scalar, it will be 
#' expanded to match the length of the other vector.
#'
#' @examples
#' ScatCoef(r = 1.5, x = 1:10, nmax = 5)
ScatCoef <- function(r, x, nmax) {

  r <- as.vector(r)
  x <- as.vector(x)
  
  # Expand scalar x to match r's size if needed
  if (length(x) == 1) {
    x <- rep(x, length(r))
  }
  
  # Ensure r and x have the same dimensions
  if (length(r) > 1 & length(x) != length(r)) {
    stop("Dimensions of x & m must be the same or scalar")
  }
  
  # Create matrix N for indices 1 to nmax
  N <- matrix(rep(1:nmax, length(x)), nrow = nmax, ncol = length(x))
  
  # Compute Ricatti-Bessel functions
  phi <- RB1(x, nmax)             # First kind
  phir <- RB1(r * x, nmax)        # First kind with complex argument
  zeta <- RB2(x, nmax)            # Second kind
  xi <- phi - 1i * zeta           # Complex xi function
  
  # Compute phi, phir, and zeta for n-1
  phin_1 <- rbind(sin(x), phi[1:(nmax - 1), ])
  phirn_1 <- rbind(sin(r * x), phir[1:(nmax - 1), ])
  zetan_1 <- rbind(-cos(x), zeta[1:(nmax - 1), ])
  
  # Ensure r and x are properly expanded for element-wise division
  if (length(r) > 1) {
    r <- matrix(rep(r, each = nmax), nrow = nmax, ncol = length(r))
  }
  if (length(x) > 1) {
    x <- matrix(rep(x, each = nmax), nrow = nmax, ncol = length(x))
  }
  
  # Compute derivatives with element-wise division
  dphi <- phin_1 - (N * phi) / x
  dphir <- phirn_1 - (N * phir) / (r * x)
  dzeta <- zetan_1 - (N * zeta) / x
  dxi <- dphi - 1i * dzeta
  
  # Compute scattering coefficients
  denom_a <- r * phir * dxi - xi * dphir
  denom_b <- phir * dxi - r * xi * dphir
  
  # Handle potential division by zero
  denom_a[denom_a == 0] <- NaN
  denom_b[denom_b == 0] <- NaN
  
  a <- (r * phir * dphi - phi * dphir) / denom_a
  b <- (phir * dphi - r * phi * dphir) / denom_b
  
  return(list(a = a, b = b))
}


#' Calculate associated Legendre polynomial functions
#'
#' This function calculates the required functions of the associated Legendre 
#' polynomials, which are needed to determine the scattering angle dependence 
#' of scattered intensity. It computes the functions from order 1 to `nmax` 
#' for the given scattering angle(s) (`ang`) using upwards recursion.
#'
#' @param ang Numeric vector of scattering angles (in radians).
#' @param nmax Integer specifying the maximum order of the functions.
#'
#' @return A list containing:
#'  \itemize{
#'    \item `p`: A matrix of dimensions `nmax` by `length(ang)`, where each column 
#'          represents the `p` function values for a specific angle.
#'    \item `t`: A matrix of dimensions `nmax` by `length(ang)`, where each column 
#'          represents the `t` function values for a specific angle.
#' }
#' @examples
#' ALegendr(ang = seq(0, pi, length.out = 10), nmax = 5)
ALegendr <- function(ang, nmax) {
  # Initialize matrices p and t
  p <- matrix(NA, nrow = nmax, ncol = length(ang))
  t <- matrix(NA, nrow = nmax, ncol = length(ang))
  
  p[1,] <- rep(1, length(ang)) # Equivalent to ones(1, size(ang, 2)) in Matlab
  t[1,] <- cos(ang)
  p[2,] <- 3 * cos(ang)
  t[2,] <- 2 * cos(ang) * p[2,] - 3
  
  for (n in 3:nmax) {
    p[n,] <- ((2 * n - 1) * cos(ang) * p[n-1,] - n * p[n-2,]) / (n - 1)
    t[n,] <- n * cos(ang) * p[n,] - (n + 1) * p[n-1,]
  }
  
  return(list(p = p, t = t))
}


#' Calculate scattered light intensity for a sphere
#'
#' This function calculates the scattered light intensity for a sphere, 
#' given its size parameter (`x`), refractive index relative to the medium (`m`), 
#' and scattering angle (`ang`). The incident light is assumed to be polarized.
#'
#' @param m Numeric vector of refractive index ratios. Can be a scalar.
#' @param x Numeric vector of size parameters. Can be a scalar.
#' @param ang Numeric vector of scattering angles (in radians).
#'
#' @return A list containing:
#' \itemize{
#'  \item `I1`: A matrix representing the parallel component of the scattered 
#'          light intensity. If `m` and/or `x` are vectors, `I1` will be a 3D 
#'          array with dimensions corresponding to `length(ang)`, `length(m)`, 
#'          and `length(x)`.
#'  \item `I2`: A matrix representing the perpendicular component of the scattered 
#'          light intensity.  If `m` and/or `x` are vectors, `I2` will be a 3D 
#'          array with dimensions corresponding to `length(ang)`, `length(m)`, 
#'          and `length(x)`.
#' }
#' @details
#' If `m` and/or `x` are vectors, the output matrices `I1` and `I2` will be 3D 
#' arrays, with the `m` and/or `x` parameters spanning the third dimension.
#'
#' @examples
#' Intensity(m = 1.5, x = 1:10, ang = seq(0, pi, length.out = 20))
Intensity <- function(m, x, ang) {
  
  # Ensure x and m are vectors of the same length
  if (length(x) == 1) {
    x <- rep(x, length(r))
  }
  if (length(m) == 1) {
    m <- rep(m, length(x))
  }
  
  # Compute the necessary number of terms
  nc <- ceiling(max(x) + 4.05 * (max(x)^(1/3)) + 2)
  n <- 1:nc
  E <- outer((2 * n + 1) / (n * (n + 1)), rep(1, length(x)))
  
  # Compute Legendre polynomials
  legendre_result <- ALegendr(ang, nc)
  p <- legendre_result$p
  t <- legendre_result$t
  
  old_warn <- options("warn")
  options(warn = -1)  # Suppress warnings
  
  # Compute scattering coefficients
  scat_coef_result <- ScatCoef(m, x, nc)
  a <- scat_coef_result$a
  b <- scat_coef_result$b
  
  # Check for invalid (NaN) results due to excessive terms for small particles
  invalid <- which(apply(is.na(rbind(a, b)), 2, any))
  
  while (length(invalid) > 0) {
    
    # Set invalid coefficients to zero
    a[, invalid] <- 0
    b[, invalid] <- 0
    
    # Recalculate `nc2` based on the invalid values
    nc2 <- ceiling(max(x[invalid]) + 4.05 * (max(x[invalid])^(1/3)) + 2)
    
    scat_coefs_new <- ScatCoef(m[invalid], x[invalid], nc2)
    a[1:nc2, invalid] <- scat_coefs_new$a
    b[1:nc2, invalid] <- scat_coefs_new$b
    
    
    # Re-evaluate invalid entries
    invalid <- which(apply(is.na(rbind(a, b)), 2, any))
    
    # Remove cases where x or m is zero
    if (length(x) >= max(invalid)) {
      invalid <- invalid[x[invalid] != 0]
    } else if (x == 0) {
      invalid <- integer(0)
    }
    
    if (length(m) >= max(invalid, 0)) {
      invalid <- invalid[m[invalid] != 0]
    } else if (m == 0) {
      invalid <- integer(0)
    }
  }
  
  options(old_warn)  # Restore previous warning settings
  
  # Compute intensity
  a <- a * E
  b <- b * E
  
  S1 <- t(a) %*% p + t(b) %*% t
  S2 <- t(a) %*% t + t(b) %*% p
  
  I1 <- Re(S1 * Conj(S1))  # Parallel component
  I2 <- Re(S2 * Conj(S2))  # Perpendicular component
  
  return(list(I1 = I1, I2 = I2))
}


#' Calculate Mie scattering intensity for a sphere at finite angles
#'
#' This function calculates the Mie scattering intensity for a sphere with 
#' given parameters across a range of particle diameters. It returns both the 
#' parallel and perpendicular components of the scattered light.
#'
#' @param np Refractive index of the particle.
#' @param nm Refractive index of the medium.
#' @param wl Wavelength of the incident light in nanometers.
#' @param na Numerical aperture of the objective lens.
#' @param dmin Minimum diameter of the particle in microns (must be > 0).
#' @param dmax Maximum diameter of the particle in microns.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'  \item `Diameter`: Particle diameter in microns.
#'  \item `AvgIpara`: Average parallel component of scattered intensity.
#'  \item `AvgIperp`: Average perpendicular component of scattered intensity.
#'  }
#' @export
#' @examples
#' 
#' # Refractive indices
#' np <- 1.6 # Beads
#'# n_particle <- 1.36 # lower limit for Phytoplankton
#'# n_particle <- 1.38 # average for Phytoplankton
#'# n_particle <- 1.41 # upper limit for Phytoplankton
#' nm <-   1.3371 # Water
#'
#' # Laser wavelength (in nanometers)
#' wl <- 457 # (SeaFlow)       
#' # wl <- 488 # (InFlux)       
#'
#' # Numerical aperture of the objective
#' na <- 0.55 # (SeaFlow)
#' # na <- 0.42  # (InFlux)
#'
#' # Particle size range 
#' dmin <- 0.025   # Minimum radius (microns)
#' dmax <- 50      # Maximum radius (microns)
#'
#' df <- FiniteAngleCalc(np, nm, lam, na, dmin, dmax)
#'
#' ggplot2::ggplot(df, ggplot2::aes(x = Diameter, y = AvgIpara)) +
#'   ggplot2::geom_line() +
#'   ggplot2::scale_y_log10() +
#'   ggplot2::labs(y = 'Relative Intensity', x = 'Diameter (µm)') +
#'   ggplot2::theme_bw() +
#'   ggplot2::theme(legend.position = "bottom") +
#'   ggplot2::ggtitle('Intensity vs Diameter')
#'   
#' # Save output
#' # readr::write_csv(df, paste0('mie_',round(m,3),'n.csv'))
FiniteAngleCalc <- function(np, nm, wl, na, dmin, dmax) {
  
  # Calculate parameters
  m <- np / nm  # Refractive index ratio
  lam <- wl * 1e-9  # Convert wavelength to microns
  k <- 2 * pi / lam   # wavenumber
  r <- 10^(seq(log10(dmin/2), log10(dmax/2), length.out = 2000)) * 1e-6  # Particle radius
  x <- k * r  # Size parameter
  max_angle <- -asin(na/nm)* 180/pi # Calculate scattering angles
  angle <- seq(max_angle, -5.7, by = 0.1)  # Scattering angles (-5.7 to account for the scatter bar)
  ang <- angle * pi / 180  # Convert to radians
  
  # Calculate intensities
  intensity_results <- Intensity(m, x, ang)
  
  # Extract I1 and I2 (fix variable names)
  Ipara <- intensity_results$I1  # Parallel component
  Iperp <- intensity_results$I2  # Perpendicular component
  
  # Average the angular intensities
  avgIpara <- rowSums(Ipara, na.rm = TRUE) / length(angle)
  avgIperp <- rowSums(Iperp, na.rm = TRUE) / length(angle)
  
  # Create a data frame
  df <- dplyr::tibble(Diameter = 2 * r * 1e6, AvgIpara = avgIpara, AvgIperp = avgIperp)
  return(df)
}