# 'RB1 is a function used to calculate the Ricatti-Bessel functions of the first 
#' and second # 'kinds respectively. The function provides the radial component 
#' of the scattering coefficients.
# 'RB1(r, nmax) calculates the Ricatti-Bessel function of the first
# 'kind from order 1 to nmax, for radius parameter(s) r, by downwards
# 'recursion. The returns are nmax-by-nr matrices, where nr is the dimension
# 'of the vector r, each column representing the results for one radius.
#' 
#' RB2(r, nmax) similarly calculates the Ricatti-Bessel function by upwards recursion.
#' 
#' @param r radius.
#' @param nmax maximum order of the Ricatti-Bessel function.
#' @return nmax-by-nr matrices, where nr is the dimension
# 'of the vector r, each column representing the results for one radius.
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

# 'RB2 is a function used to calculate the Ricatti-Bessel functions of the first 
#' and second # 'kinds respectively. The function provides the radial component 
#' of the scattering coefficients.
# 'RB2(r, nmax) calculates the Ricatti-Bessel function of the first
# 'kind from order 1 to nmax, for radius parameter(s) r, by upwards
# 'recursion. The returns are nmax-by-nr matrices, where nr is the dimension
# 'of the vector r, each column representing the results for one radius.
#' 
#' RB1(r, nmax) similarly calculates the Ricatti-Bessel function by downwards recursion.

#' @param r radius.
#' @param nmax maximum order of the Ricatti-Bessel function.
#' @return nmax-by-nr matrices, where nr is the dimension
# 'of the vector r, each column representing the results for one radius.
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

#' ScatCoef is a function used to calculate the scattering coefficients required
#' by Mie calculations (e.g. scattering intensity)
#' This function requires the Ricatti-Bessel functions. 
#' It calculates the scattering coefficients from order 1 to nmax for refractive
#' index and particle size parameter(s) m and x. The function returns a and b 
#' are complex nmax-by-nx matrices, where nx is the dimension of the m and/or x, 
#' with each column representing one particle parameter(pair). 
#' m, x or both can be scalar. If both are vectors, they must be of the same 
#' dimension to form particle parameter pairs.
#' @param r radius.
#' @param x refractive index.
#' @param nmax maximum order of the Ricatti-Bessel function.
#' @return returns a list of a and b, that are nmax-by-nr matrices, where nr is the dimension
# 'of the vector r, each column representing the results for one radius.
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


#' The ALegendr function returns the required functions of the associated Legendre
#' polynomials. This function is required to provide the scattering angle 
#' dependency of the scattered intensity.
#' The function calculates the required functions from order 1 to nmax, 
#' for angle(s) q, by upwards recursion. The returns p and t are nmax-by-nq
#' matrices, where nq is the dimension of the vector q, each column representing 
#' the results for one angle.
#' @param ang scattering angle.
#' @param nmax maximum order of the Ricatti-Bessel function.
#' @return returns a list of p and t, that are nmax-by-nq matrices, where nq is 
#' the dimension of the vector q, each column representing  the results for one angle.
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


#' The Intensity function returns Stoke's vector of scattered light.
#' the function calculates  the scattered Light for a sphere, size x, 
#' refractive index relative to medium m, at angle ang.
#' The incident light is polarized. If m and/or x are vectors, return will 3D,
#' with m and/or x parameters spanning third dimension.
#' @param m the Refractive index ratio of the sphere to the medium.
#' @param x The Size parameter
#' @param ang scattering angle (in radians).
#' @return returns a list of p and t, that are nmax-by-nq matrices, where nq is 
#' the dimension of the vector q, each column representing  the results for one angle.
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


#' The function FiniteAngleCalc calculates Mie scattering intensity for a sphere
#' of given parameters at the required angle.
#' 
#' The function returns both the parallel and perpendicular components of the 
#' scattered light across the range of particle diameter.
#'
#'@param np The refractive index of the particle.
#'@param nm The refractive index of the medium.
#'@param wl The wavelength of the incident light (in nanometer).
#'@param na The numerical aperture of the objective lens.
#'@param dmin The minimum diameter of the particle (in microns). Needs to be larger than 0
#'@param dmax The maximum diameter of the particle (in microns).
#' @return returns a table of the parallel and perpendicular components of the 
#' scattered light across the range of particle diameter.
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