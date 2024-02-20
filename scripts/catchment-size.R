################################################
# Sampling Strategies: Catchment Size Analysis #
################################################

#==========#
# PREAMBLE #
#==========#

# Libraries
library(tidyverse)
library(testit)

# Import plot formatting
source("scripts/aux_format-plots.R")

# Fixed parameters
f <- 1e-6 # Pathogen reads from infected individual
z <- 0 # Non-human fraction
r_tot <- 1e9 # Total reads

#=========================================#
# PART 1: VALIDATION OF MEAN AND VARIANCE #
#=========================================#

# Pick simple fixed values for n & p
N <- 1000000
P <- 0.01
n_iter <- 1e6

# Calculate number of infecteds
K <- rbinom(n_iter, N, P)

# Calculate lambda values
L <- r_tot * f * K / (N+z)

# Calculate number of pathogen reads
R <- rpois(n_iter, L)

# Predicted properties
E_R_exp <- r_tot * f * P
E_R_sq_exp <- E_R_exp^2 + E_R_exp + r_tot^2*f^2*P*(1-P)/N
Var_R_exp <- E_R_sq_exp - E_R_exp^2

# Observed properties
E_R_obs <- mean(R)
E_R_sq_obs <- mean(R^2)
Var_R_obs <- var(R)

cat("E(R):", E_R_exp, "vs", E_R_obs, "\n")
cat("E(R^2):", E_R_sq_exp, "vs", E_R_sq_obs, "\n")
cat("Var(R):", Var_R_exp, "vs", Var_R_obs, "\n")

#===========================================#
# PART 2: PROBABILITY OF DETECTION OVERVIEW #
#===========================================#

#---------------------
# Variable parameters
#---------------------

n <- round(10^seq(0,4,0.25)) # Catchment size
p <- signif(10^seq(-4.5,-0.5,0.5), 2) # Prevalence in population
T <- round(10^seq(0,2,0.5)) # Read thresholds for detection

#---------------------
# Auxiliary functions
#---------------------

compute_infection_table <- function(catchment_size, # Catchment size
                                    prevalence, # Prevalence in population
                                    pfrac_indiv, # Fraction of pathogen reads from infected individual
                                    nfrac_nonhuman, # Fraction of nonhuman material
                                    nreads_tot){ # Total reads
  # Calculate probability distribution across number of infected individuals
  tab_inf <- expand_grid(catchment_size = catchment_size,
                         n_infected = 0:max(catchment_size),
                         prevalence = prevalence,
                         nfrac_nonhuman = nfrac_nonhuman,
                         pfrac_indiv = pfrac_indiv,
                         nreads_tot = nreads_tot) %>%
    filter(n_infected <= catchment_size) %>%
    mutate(pfrac_sample = pfrac_indiv * n_infected / (catchment_size + nfrac_nonhuman), # Fraction of nucleic acid in sample originating from pathogen
           lambda = nreads_tot * pfrac_sample, # Rate of pathogen reads from sequencing Poisson process
           p_k = dbinom(n_infected, catchment_size, prevalence)) # Probability of observing k infected individuals given catchment & prevalence
  return(tab_inf)
}

compute_threshold_table <- function(tab_infected,
                                    threshold){
  # Calculate probability of detection given values for threshold and K
  tab_thresh <- expand_grid(threshold = threshold,
                            j = 0:max(T)) %>% # Intermediate value for calculating Poisson CDF
    filter(j <= threshold) %>%
    full_join(tab_infected, by = character()) %>%
    mutate(cdf_val = lambda^j/factorial(j)) # Contribution to CDF of Poisson process
  tab_thresh_summ <- tab_thresh %>%
    group_by(threshold, catchment_size, n_infected, prevalence, nfrac_nonhuman,
             pfrac_indiv, nreads_tot, pfrac_sample, lambda, p_k) %>%
    summarise(cdf_sum = sum(cdf_val), .groups = "drop") %>% # Summation part of CDF
    mutate(cdf = cdf_sum * exp(-lambda), # Overall CDF of Poisson process
           p_ident_k = 1-cdf) # Probability of identification given k infected individuals
  return(tab_thresh_summ)
}

compute_detection_probability <- function(tab_threshold){
  # Compute overall probability of identification
  tab_threshold %>% 
    group_by(threshold, catchment_size, prevalence, nfrac_nonhuman, pfrac_indiv,
             nreads_tot) %>%
    summarise(p_ident = sum(p_ident_k * p_k), .groups = "drop")
}

#-----------------------
# Compute probabilities
#-----------------------

tab_infected <- compute_infection_table(n,p,f,z,r_tot)
tab_threshold <- compute_threshold_table(tab_infected, T)
tab_detection <- compute_detection_probability(tab_threshold)

#-----------------------
# Predict asymptotes
#-----------------------

tab_asymptotes <- tab_detection %>% group_by(threshold, prevalence, pfrac_indiv, nreads_tot) %>%
  summarize %>% mutate(asymptote = 1-ppois(threshold, nreads_tot * prevalence * pfrac_indiv))


#--------------
# Plot results
#--------------

# Define aesthetics
g_aes_calc <- aes(x=catchment_size, y=p_ident,
                  colour = threshold, group = threshold)

# Generate plot
tab_plot <- tab_detection %>% 
  mutate(threshold = fct_inorder(as.character(threshold)),
         label = paste0("Prevalence: 10^", round(log10(prevalence), 1)))
tab_plot_asymptotes <- tab_asymptotes %>%
  mutate(threshold = fct_inorder(as.character(threshold)),
         label = paste0("Prevalence: 10^", round(log10(prevalence), 1)))
g_calc <- ggplot(tab_plot, g_aes_calc) + 
  geom_hline(data=tab_plot_asymptotes, aes(yintercept=asymptote, color=threshold), linetype="dashed", alpha=0.5) +
  geom_line(size = 1) +
  scale_x_log10(name = "Catchment size") +
  scale_y_continuous(name = "Probability of detection", expand = c(0,0),
                     breaks = seq(0,1,0.2)) +
  scale_color_manual(values = unname(palette_primary), name = "Threshold (# reads)") +
  facet_wrap(~label, labeller = label_parsed) +
  theme_base

#==========================================#
# PART 3: PROBABILITY OF DETECTION ZOOM-IN #
#==========================================#

#---------------------
# Variable parameters
#---------------------

n2 <- round(10^seq(0,5,0.25)) # Catchment size
p2 <- 0.01 # Fix prevalence
T2 <- c(1, 5, 10, 12, 15, 20, 30, 50, 100) # Thresholds

#-----------------------
# Compute probabilities
#-----------------------

tab_infected_2 <- compute_infection_table(n2,p2,f,z,r_tot)
tab_threshold_2 <- compute_threshold_table(tab_infected_2, T2)
tab_detection_2 <- compute_detection_probability(tab_threshold_2)
tab_asymptotes_2 <- tab_detection_2 %>% group_by(threshold, prevalence, pfrac_indiv, nreads_tot) %>%
  summarize %>% mutate(asymptote = 1-ppois(threshold, nreads_tot * prevalence * pfrac_indiv))

#--------------
# Plot results
#--------------

# Generate plot
tab_plot_2 <- tab_detection_2 %>%
  arrange(threshold) %>%
  mutate(threshold = fct_inorder(as.character(threshold)),
         label = paste0("Prevalence: 10^", round(log10(prevalence), 1)))
tab_plot_asymptotes_2 <- tab_asymptotes_2 %>%
  mutate(threshold = fct_inorder(as.character(threshold)),
         label = paste0("Prevalence: 10^", round(log10(prevalence), 1)))

g_calc_2 <- ggplot(tab_plot_2, g_aes_calc) + 
  geom_hline(data=tab_plot_asymptotes_2, aes(yintercept=asymptote, color=threshold), linetype="dashed", alpha=0.5) +
  geom_line(size = 1) +
  scale_x_log10(name = "Catchment size") +
  scale_y_continuous(name = "Probability of detection", expand = c(0,0),
                     breaks = seq(0,1,0.2)) +
  scale_color_manual(values = unname(palette_secondary), name = "Threshold (# reads)") +
  facet_wrap(~label, labeller = label_parsed) +
  theme_base
g_calc_2

#========================#
# PART 4: SIGNAL DENSITY #
#========================#

#---------------------
# Set parameters
#---------------------

n3 <- c(30, 100, 1000, 10000)
p3 <- 0.01
R  <- 0:1000

#-----------------------------
# Calculate poisson densities
#-----------------------------

# Number of infected individuals
tab_infected_3 <- compute_infection_table(n3,p3,f,z,r_tot) %>%
  filter(p_k > 1e-100) %>%
  full_join(tibble(nreads_path = R), by = character()) %>%
  mutate(p_r_k = dpois(nreads_path, lambda))

# Probability of pathogen read counts
tab_reads_3 <- tab_infected_3 %>%
  group_by(catchment_size, prevalence, nfrac_nonhuman, pfrac_indiv, nreads_tot,
           nreads_path) %>%
  summarise(p_r = sum(p_k * p_r_k), .groups = "drop")

#-----------------------------
# Plot
#-----------------------------

tab_plot_3 <- tab_reads_3 %>% arrange(catchment_size) %>%
  mutate(catchment_size = fct_inorder(as.character(catchment_size)))

g_aes_reads <- aes(x=nreads_path, y=p_r, color = catchment_size,
                   group = catchment_size)
g_reads <- ggplot(tab_plot_3, g_aes_reads) + geom_line(size = 1) +
  geom_vline(xintercept = c(5,20), colour = "red", linetype = "dashed", size = 0.7) +
  scale_x_log10(name = "# pathogen reads") +
  scale_y_continuous(name = "Probability", lim = c(0,1), expand = c(0,0)) +
  scale_color_manual(values = unname(palette_primary), name = "Catchment size") +
  coord_cartesian(ylim = c(0,0.15)) +
  theme_base

#========================#
# PART 5: EXPORT RESULTS #
#========================#

save_fig("figures/p_detection_broad.png", g_calc, width = 6, aspect = 1)
save_fig("figures/p_detection_narrow.png", g_calc_2, width = 4, aspect = 1)
save_fig("figures/p_detection_narrow_2.png", g_calc_2 + theme(aspect.ratio = 1/2), width = 6, aspect = 0.8)
save_fig("figures/density.png", g_reads, width = 4, aspect = 1)
