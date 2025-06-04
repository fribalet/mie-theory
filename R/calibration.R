# Load required libraries
library(DEoptim)
library(tidyverse)
library(purrr)

# Source Mie theory calculations
source("mietheory.R")

#-----------------
# DATA PREPARATION 
#-----------------

# Load and prepare ALL calibration data (beads + cultures)
# Include both FSC and SSC data with separate numerical apertures
summary_aggregate <- read_csv("summary_aggregate.csv") %>%
  filter((Experiment == "Beads" & pop == "beads") | 
           (Experiment == "Cultures" & pop != "beads")) %>%
  filter(channel == "FSC" | channel == "SSC") %>%
  filter(diameter > 0.2) %>% 
  mutate(
    na_fsc = case_when(
      grepl("SeaFlow", Instrument) ~ 0.55,
      grepl("Influx", Instrument) ~ 0.42,
      TRUE ~ NA),
    na_ssc = case_when(
      grepl("Influx", Instrument) ~ 0.6,
      grepl("Attune", Instrument) ~ 1.2, 
      grepl("Cytoflex", Instrument) ~ 1.3, 
      TRUE ~ NA)) %>%
  # Add refractive index (np) based on sample type
  # https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/product/documents/817/733/lb11pis.pdf
  # Refractive Index for Polystyrene-based Latex Particles: 1.602 at 486 nm / 1.5905 at 589 nm
  mutate(np = case_when(
    Experiment == "Beads" ~ 1.603,
    Experiment == "Cultures" ~ 1.38,
    TRUE ~ NA
  )) %>%
  filter(!is.na(geom), geom > 0)

# Plot data
p1 <- summary_aggregate %>% 
  mutate(Instrument = factor(Instrument, 
                             levels = c("Influx 1", "Influx 2", "Influx 3",
                                        "Biorad", "Cytek", "Cytoflex", 
                                        "Attune 1", "Attune 2", "Calibur", "Accuri"))) %>%
  ggplot(aes(diameter, geom, color = Experiment)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_pointrange(aes(ymin = geom / sd, 
                      ymax = geom * sd)) +
  geom_pointrange(aes(xmin = diameter - sd_diameter, 
                      xmax = diameter + sd_diameter)) +  
  scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  labs(x = "Diameter (μm)", 
       y = "Scatter (normalized to 1 μm beads)") +
  facet_grid(channel ~ Instrument, scales = "free_y") +
  coord_flip() +
  my_theme

print(p1)




#------------------------------------
# MIE THEORY PARAMETERS
#------------------------------------

# Physical constants
nm <- 1.3371   # Water refractive index
wl <- 488      # Laser wavelength (nm)
dmin <- 0.3    # Minimum diameter (μm)
dmax <- 10     # Maximum diameter (μm)

#---------------------------------
# MODIFIED OPTIMIZATION FUNCTION
#---------------------------------

#' Calculate weighted sum of squared errors for combined bead + culture data
#' Handles both FSC and SSC predictions with separate numerical apertures
#' 
#' @param data Combined data frame with columns: diameter, geom, channel, na, np
#' @param params Vector [c_fsc, b_fsc, c_ssc, b_ssc, na_fsc, na_ssc] - optimization parameters
#' @param optimize_na_fsc Logical - whether to optimize FSC numerical aperture
#' @param optimize_na_ssc Logical - whether to optimize SSC numerical aperture
#' @return Weighted mean squared error
sigma.combined <- function(data, params, optimize_na_fsc = FALSE, optimize_na_ssc = FALSE) {
  
  c_fsc <- params[1]   # FSC scaling constant
  b_fsc <- params[2]   # FSC power law exponent  
  c_ssc <- params[3]   # SSC scaling constant
  b_ssc <- params[4]   # SSC power law exponent
  
  # Determine NA values
  if (optimize_na_fsc) {
    na_fsc_opt <- params[5]
  } else {
    na_fsc_opt <- NA
  }
  
  if (optimize_na_ssc) {
    na_ssc_opt <- ifelse(optimize_na_fsc, params[6], params[5])
  } else {
    na_ssc_opt <- NA
  }
  
  # Calculate Mie predictions for each particle type and scatter type
  mie_results <- data %>%
    group_by(np, channel) %>%
    do({
      current_np <- unique(.$np)
      current_channel <- unique(.$channel)
      
      # Determine which NA to use
      if (current_channel == "FSC") {
        # Use optimized NA if available, otherwise use data NA
        if (optimize_na_fsc) {
          current_na <- na_fsc_opt
        } else {
          current_na <- unique(.$na_fsc[!is.na(.$na_fsc)])
          if (length(current_na) == 0) return(tibble())
          current_na <- current_na[1]
        }
        
        mie_pred <- FiniteAngleCalc(current_np, nm, wl, current_na, .$diameter, fsc = TRUE)
        mutate(., 
               predicted_intensity = mie_pred$AvgIpara,
               predicted_scatter = (predicted_intensity/c_fsc)^b_fsc)
        
      } else { # SSC
        # Use optimized NA if available, otherwise use data NA
        if (optimize_na_ssc) {
          current_na <- na_ssc_opt
        } else {
          current_na <- unique(.$na_ssc[!is.na(.$na_ssc)])
          if (length(current_na) == 0) return(tibble())
          current_na <- current_na[1]
        }
        
        mie_pred <- FiniteAngleCalc(current_np, nm, wl, current_na, .$diameter, fsc = FALSE)
        mutate(., 
               predicted_intensity = mie_pred$AvgIpara,
               predicted_scatter = (predicted_intensity/c_ssc)^b_ssc)
      }
    }) %>%
    ungroup()
  
  # Skip if no valid predictions were generated
  if (nrow(mie_results) == 0) {
    return(Inf)
  }
  
  ###############################################################
  # Weight measurements by biological relevance and sample type #
  ###############################################################
  
  weights <- case_when(
    # Higher weight for beads (calibration standards)
    mie_results$Experiment == "Beads" & 
      mie_results$diameter >= 0.7 & mie_results$diameter <= 5 ~ 1,
    # Extra weight for biologically relevant sizes in cultures
    mie_results$Experiment == "Cultures" & 
      mie_results$diameter >= 0.7 & mie_results$diameter <= 5 ~ 1,
    TRUE ~ 5)
  
  # Calculate weighted mean squared error on log scale
  sigma <- weighted.mean(
    (log(mie_results$geom) - log(mie_results$predicted_scatter))^2, 
    weights, 
    na.rm = TRUE
  )
  
  return(sigma)
}

#-------------------------------------------------
# MODIFIED OPTIMIZATION FUNCTION
#-------------------------------------------------

optimize_combined <- function(data) {
  
  ####################
  # Parameter bounds #
  ####################
  # Check what numerical apertures are available
  available_na_fsc <- unique(data$na_fsc[data$channel == "FSC" & !is.na(data$na_fsc)])
  available_na_ssc <- unique(data$na_ssc[data$channel == "SSC" & !is.na(data$na_ssc)])
  
  # Determine if we need to optimize NA values
  optimize_na_fsc <- length(available_na_fsc) == 0
  optimize_na_ssc <- length(available_na_ssc) == 0
  
  cat("  Available NA for FSC:", ifelse(length(available_na_fsc) > 0, available_na_fsc, "None - will optimize"), "\n")
  cat("  Available NA for SSC:", ifelse(length(available_na_ssc) > 0, available_na_ssc, "None - will optimize"), "\n")
  
  # Set up bounds for scaling parameters
  c_fsc_lower <- 1; c_fsc_upper <- 5000
  b_fsc_lower <- 0.1; b_fsc_upper <- 1.5
  c_ssc_lower <- 0.1; c_ssc_upper <- 100
  b_ssc_lower <- 0.1; b_ssc_upper <- 1.5
  
  # Set up bounds for NA parameters if needed
  lower_bounds <- c(c_fsc_lower, b_fsc_lower, c_ssc_lower, b_ssc_lower)
  upper_bounds <- c(c_fsc_upper, b_fsc_upper, c_ssc_upper, b_ssc_upper)
  
  if (optimize_na_fsc) {
    # FSC NA bounds
    na_center <- 0.3
    na_lower <- round(max(0.2, na_center / 2), 2)
    na_upper <- round(min(nm, na_center * 2), 2)
    lower_bounds <- c(lower_bounds, na_lower)
    upper_bounds <- c(upper_bounds, na_upper)
    cat("  FSC NA optimization range:", na_lower, "to", na_upper, "\n")
  }
  
  if (optimize_na_ssc) {
    # SSC NA bounds
    na_center <- 0.9
    na_lower <- round(max(0.2, na_center / 2), 2)
    na_upper <- round(min(0.95 * nm, na_center * 2), 2)
    lower_bounds <- c(lower_bounds, na_lower)
    upper_bounds <- c(upper_bounds, na_upper)
    cat("  SSC NA optimization range:", na_lower, "to", na_upper, "\n")
  }
  
  # Create objective function with closure
  objective_function <- function(params) {
    sigma.combined(data, params, optimize_na_fsc, optimize_na_ssc)
  }
  
  cat("  Starting optimization...\n")
  
  optimal_result <- tryCatch({
    DEoptim(
      fn = objective_function,
      lower = lower_bounds, 
      upper = upper_bounds,
      control = DEoptim.control(
        itermax = 1000,      
        reltol = 1e-4, 
        steptol = 50,
        trace = 50
      )
    )
  }, error = function(e) {
    cat("    Combined optimization failed:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(optimal_result)) {
    # Return results along with optimization flags
    result <- list(
      params = optimal_result$optim$bestmem,
      optimize_na_fsc = optimize_na_fsc,
      optimize_na_ssc = optimize_na_ssc,
      available_na_fsc = available_na_fsc,
      available_na_ssc = available_na_ssc
    )
    return(result)
  }
  return(NULL)
}

#--------------------------
# MAIN OPTIMIZATION WORKFLOW
#--------------------------

all_instruments <- unique(summary_aggregate$Instrument)

# Run combined optimization for each instrument
results <- all_instruments %>%
  set_names() %>%
  map_dfr(function(instrument_name) {
    
    cat("Processing instrument:", instrument_name, "\n")
    cat(paste(rep("-", nchar(instrument_name) + 21), collapse = ""), "\n")
    
    # Get all data for this instrument (beads + cultures, FSC + SSC)
    instrument_data <- summary_aggregate %>% 
      filter(Instrument == instrument_name) %>%
      arrange(Experiment, channel, diameter)
    
    cat("  Bead points:", sum(instrument_data$Experiment == "Beads"), "\n")
    cat("  Culture points:", sum(instrument_data$Experiment == "Cultures"), "\n")
    
    if (nrow(instrument_data) < 5) {  # Need more points for 5-parameter fit
      cat("  WARNING: Insufficient data points\n\n")
      return(tibble())
    }
    
    # Optimize using combined data
    optimization_result <- optimize_combined(instrument_data)
    
    if (is.null(optimization_result)) {
      cat("  ERROR: Combined optimization failed\n\n")
      return(tibble())
    }
    
    # Extract results
    params <- optimization_result$params
    optimize_na_fsc <- optimization_result$optimize_na_fsc
    optimize_na_ssc <- optimization_result$optimize_na_ssc
    available_na_fsc <- optimization_result$available_na_fsc
    available_na_ssc <- optimization_result$available_na_ssc
    
    cat("  Combined optimization results:\n")
    cat("    FSC scaling factor (c):", round(params[1], 2), "\n")
    cat("    FSC power exponent (b):", round(params[2], 3), "\n") 
    cat("    SSC scaling factor (c):", round(params[3], 2), "\n")
    cat("    SSC power exponent (b):", round(params[4], 3), "\n")
    
    # Determine final NA values
    if (optimize_na_fsc) {
      na_fsc_final <- params[5]
      cat("    FSC numerical aperture (optimized):", round(na_fsc_final, 4), "\n")
    } else {
      na_fsc_final <- available_na_fsc[1]
      cat("    FSC numerical aperture (fixed):", round(na_fsc_final, 4), "\n")
    }
    
    if (optimize_na_ssc) {
      na_ssc_final <- ifelse(optimize_na_fsc, params[6], params[5])
      cat("    SSC numerical aperture (optimized):", round(na_ssc_final, 4), "\n")
    } else {
      na_ssc_final <- available_na_ssc[1]
      cat("    SSC numerical aperture (fixed):", round(na_ssc_final, 4), "\n")
    }
    
    # Generate theoretical curves
    diameter_range <- 10^(seq(log10(dmin), log10(dmax), length.out = 200))
    
    # Create all combinations of particle types and scatter types
    combinations <- expand_grid(
      np = c(1.603, 1.38),
      particle_type = c("Beads", "Cultures"),
      channel = c("FSC", "SSC")
    ) %>%
      filter((np == 1.603 & particle_type == "Beads") | 
               (np == 1.38 & particle_type == "Cultures"))
    
    # Generate curves for each combination
    curves <- map_dfr(1:nrow(combinations), function(i) {
      combo <- combinations[i, ]
      
      # Determine which NA and parameters to use
      if (combo$channel == "FSC") {
        current_na <- na_fsc_final
        is_fsc <- TRUE
      } else { # SSC
        current_na <- na_ssc_final
        is_fsc <- FALSE
      }
      
      # Skip if no valid NA
      if (is.na(current_na)) {
        return(tibble())
      }
      
      # Calculate Mie predictions
      mie_curve <- FiniteAngleCalc(combo$np, nm, wl, current_na, diameter_range, fsc = is_fsc) %>%
        mutate(
          Instrument = instrument_name,
          particle_type = combo$particle_type,
          channel = combo$channel,
          np = combo$np,
          c_fsc = params[1], b_fsc = params[2],
          c_ssc = params[3], b_ssc = params[4],
          na_fsc = na_fsc_final,
          na_ssc = na_ssc_final,
          calibrated_Scatter = if (is_fsc) {
            (AvgIpara/params[1])^params[2]
          } else {
            (AvgIpara/params[3])^params[4]
          }
        )
      
      return(mie_curve)
    })
    
    return(curves)
  })


write_csv(results, "mie_calibration.csv")

#-----------------
# PLOTTING RESULTS  
#-----------------
results <- read_csv("mie_calibration.csv") %>%
  mutate(Instrument = factor(Instrument, 
                             levels = c("Influx 1", "Influx 2", "Influx 3",
                                        "Biorad", "Cytek", "Cytoflex", 
                                        "Attune 1", "Attune 2", "Calibur", "Accuri")))

# Summary of parameters
results %>% 
  group_by(Instrument) %>%
  reframe(c_fsc = mean(c_fsc), 
          b_fsc = mean(b_fsc),
          c_ssc = mean(c_ssc),
          b_ssc = mean(b_ssc),
          na_fsc = mean(na_fsc, na.rm = TRUE),
          na_ssc = mean(na_ssc, na.rm = TRUE))


# Plot combined results with data overlay for both FSC and SSC
p2 <- summary_aggregate %>% 
  mutate(Instrument = factor(Instrument, 
                             levels = c("Influx 1", "Influx 2", "Influx 3",
                                        "Biorad", "Cytek", "Cytoflex", 
                                        "Attune 1", "Attune 2", "Calibur", "Accuri"))) %>%
  ggplot(aes(diameter, geom, color = Experiment)) +
  geom_line(data = results,
            aes(x = Diameter, y = calibrated_Scatter,
                color = particle_type), size = 1) +
  geom_point(alpha = 0.5, size = 2) +
  geom_linerange(aes(ymin = geom / sd, 
                     ymax = geom * sd)) +
  geom_linerange(aes(xmin = diameter - sd_diameter, 
                     xmax = diameter + sd_diameter)) +  
  scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  labs(x = "Diameter (μm)", 
       y = "Light Scatter Intensity (normalized to 1 μm beads)") +
  facet_wrap(channel ~ Instrument) +
  coord_flip() +
  my_theme 

print(p2)

#---------------
# Export results
#---------------
ggsave("calibration_plot.png", p2, width = 15, height = 10, dpi = 300)

