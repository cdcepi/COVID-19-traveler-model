rm(list = ls())

require(deSolve)
require(tidyverse)

################################
# Within host viral dynamics ODE model
# adapted from https://www.medrxiv.org/content/10.1101/2020.08.07.20169920v2
# and https://github.com/ashish2goyal/SARS_CoV_2_Super_Spreader_Event
################################

SARS_COV_model <- function(t, x, params) {
  with(as.list(x), {
    ddt_S <- -beta * V * S
    ddt_I <- beta * V * S - delta * I^k * I - m * E^r * I / (E^r + E50^r)
    ddt_V <- p * I - c * V

    ddt_M1 <- w * I * M1 - q * M1
    ddt_M2 <- q * (M1 - M2)
    ddt_E <- q * M2 - de * E

    der <- c(ddt_S, ddt_I, ddt_V, ddt_E, ddt_M1, ddt_M2)

    list(der)
  })
}

################################
# 1. Read parameters for within host model
################################

Parameters <- read.csv("Goyal dynamics/SimulatedParameters_no_added_heterogentity.txt", header = TRUE)

alpha <- c(2)
AUC_factor <- c(4) ## lambda is 10^AUC_factor

n_simulation <- 10000 ### We simulate each parameter combination for 1000 transmitters

tzero <- 1
dT <- 0.1 # Sampling frequency in days ; this should be less than 1 -- noted only rho will change with change in dT as rho_final = rho_paper/dT
times <- seq(tzero, tzero + 30, by = dT)
Vdata <- as.data.frame(times)
ProbV_data <- as.data.frame(times)

sample_idx <- sample(1:n_simulation, replace = FALSE) 

for (i in sample_idx) { ## running the loop for 1000 transmitters

  ## parameters for the within host model
  beta <- Parameters$beta[i]
  delta <- Parameters$delta[i]
  k <- Parameters$k[i]
  p <- Parameters$p[i]
  m <- Parameters$m[i]
  w <- Parameters$w[i]
  E50 <- Parameters$E50[i]
  r <- Parameters$r[i]
  q <- Parameters$q[i]
  de <- Parameters$de[i]
  c <- Parameters$c[i]


  # Initial conditions for the within-host model
  S_0 <- 1e7
  I_0 <- 1
  V_0 <- p * I_0 / c
  #V_0 <- params$p * I_0 / params$c
  E_0 <- 0
  M1_0 <- 1
  M2_0 <- 0

  ## The incubation period (time from exposure to symtptoms) is gamma distributed -- an incubation period with a mean of 6.4 and SD of 2.3 days [4], and a mean of 4.8 and a SD 2.6 days [7]  --- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/
  mean_incubation_period <- 5.2 # mean from https://www.acpjournals.org/doi/10.7326/M20-0504 and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/
  std_incubation_period <- 2.8 # days --  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/
  incubation_period_infector <- pmax(0, rgamma(
    n = 1, shape = 3.45,
    scale = (1 / 0.66)
  )) # rate=0.66, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/

  incubation_period_infectee <- pmax(0, rgamma(
    n = 1, shape = 3.45,
    scale = (1 / 0.66)
  )) # rate=0.66, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/


  ##### Now we simulate viral dynamics for period after viral load kicks off with one infected cell
  init.x <- c(S = S_0, I = I_0, V = V_0, E = E_0, M1 = M1_0, M2 = M2_0)

  viral_load <- as.data.frame(lsodar(init.x, times, SARS_COV_model, parms=c()))$V
  viral_load[times[1:length(viral_load)] > tzero + 20] <- 0 ## assumption that viral loads do not persists beyond 20 days
  if (length(viral_load) < length(times)) viral_load <- c(viral_load, rep(0, length(times) - length(viral_load)))

  prob_shedding <- (pmax(0, (viral_load)))^alpha / ((10^AUC_factor)^alpha + (pmax(0, (viral_load)))^alpha) ### defining infectiousness according to viral loads
  prob_shedding[is.na(prob_shedding)] <- 0

  Vdata <- cbind(Vdata, viral_load)
  ProbV_data <- cbind(ProbV_data, prob_shedding)
}

### plot some simulations to check
plot(Vdata$times, log10(Vdata[, 2]), type = "n")
for (i in 2:30) lines(Vdata$times, log10(Vdata[, i]), col=adjustcolor('black', 0.1))

plot(ProbV_data$times, ProbV_data[, 2], type = "n")
for (i in 2:30) lines(ProbV_data$times, ProbV_data[, i], col=adjustcolor('black', 0.1))

#######################################################################
### convert to density function
prob_infectious_long <- as_tibble(ProbV_data, .name_repair = 'unique') %>%
  pivot_longer(starts_with('prob_shedding'),
    names_to = "sim", values_to = "p_infectious") %>%
  mutate(p_infectious = ifelse(is.na(p_infectious), 0, p_infectious))

d_infectious_Goyal_unnormalized <- 
  approxfun(density(filter(prob_infectious_long, p_infectious >= 0.5)$times, 
  from=-100, to=100))

inf_Goyal_norm <- sum(d_infectious_Goyal_unnormalized(-100:100))
  
d_infectious_Goyal <- function(t) d_infectious_Goyal_unnormalized(t) / inf_Goyal_norm

sum(d_infectious_Goyal(-10:30))

save(d_infectious_Goyal, d_infectious_Goyal_unnormalized, inf_Goyal_norm,
  file="Goyal dynamics/Goyal_infectiousness.RData")

### plot to check
dpi <- seq(0, 21, by=0.1)
plot(dpi, d_infectious_Goyal(dpi), type='l')


