### code adapted from https://github.com/cmmid/travel_screening_strategies

rm(list=ls())

require(tidyverse)

source('Clifford dynamics/GammaParamsFromQuantiles.R')
source('https://raw.githubusercontent.com/cmmid/travel_screening_strategies/568cd63ca152e5221f0d2e5ab34686c463fd337c/he.R')
source('https://raw.githubusercontent.com/cmmid/travel_screening_strategies/568cd63ca152e5221f0d2e5ab34686c463fd337c/wolfel.R')

n_travellers <- 10000

gamma2mv <- function(shape, rate=NULL, scale=NULL){
  if (is.null(rate)){
    rate <- 1/scale
  }

  list(mean = shape/rate,
       var  = shape/rate^2)
}

mv2gamma <- function(mean, var){
  list(shape = mean^2/var,
       rate  = mean/var,
       scale = var/mean)
}

time_to_event <- function(n, mean, var){
  if (var > 0){
    parms <- mv2gamma(mean, var)
    return(rgamma(n, shape = parms$shape, rate = parms$rate))
  } else{
    return(rep(mean, n))
  }
}

pathogen <- list(
  symptomatic = 
    # review paper Byrne et al. (2020) https://doi.org/10.1101/2020.04.25.20079889
    # define T1 as infection to beginning of presymptomatic infectious period
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      gamma.parms.from.quantiles(q = c(5.1, 11.5),
                                 p = c(0.5, 0.975)) %>%
        {list(shape = .$shape, scale = .$scale)}  %>% 
        {gamma2mv(.$shape, scale = .$scale)} %>% 
        set_names(., c("mu_inc", "sigma_inc")),
      
      # Li et al https://www.nejm.org/doi/full/10.1056/nejmoa2001316
      {c(9.1, 14.7)} %>% 
        set_names(., c("mu_inf", "sigma_inf"))),
  
  asymptomatic = 
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      gamma.parms.from.quantiles(q = c(5.1, 11.5),
                                 p = c(0.5,0.975)) %>%
        {list(shape = .$shape, scale = .$scale)}  %>% 
        {gamma2mv(.$shape, scale = .$scale)} %>% 
        set_names(., c("mu_inc", "sigma_inc")),
      # https://doi.org/10.1101/2020.04.25.20079889
      list(
        mu_inf    =  6,
        sigma_inf = 12))) %>%
  map(~data.frame(.x), .id = "type")

incubation_times <- crossing(idx  = 1:n_travellers,
                               type = c("symptomatic",
                                        "asymptomatic") %>%
                                 factor(x = .,
                                        levels = .,
                                        ordered = T)) %>%
    split(.$type) %>%
    map2_df(.x = .,
            .y = pathogen,
            ~mutate(.x,
                    exp_to_onset   = time_to_event(n = n(),
                                                   mean = .y$mu_inc, 
                                                   var  = .y$sigma_inc),
                    onset_to_recov = time_to_event(n = n(),
                                                   mean = .y$mu_inf, 
                                                   var  = .y$sigma_inf))) 
  
# infectious period from time of onset to no longer infectious
gamma.parms.from.quantiles(q = c(5.1, 11.5),
                           p = c(0.5, 0.975)) %>%
  {list(shape = .$shape, scale = .$scale)} -> inc_parms

incubation_times %<>% 
  mutate(u = runif(n = nrow(.), 0.01, 0.99)) %>%
  mutate(inf_from_onset = 
           approx(x    = wolfel_pred$y, 
                  y    = wolfel_pred$day, 
                  xout = u)$y,
         pre_symp_lead  = 
           approx(x    = HE$p,
                  y    = HE$delay,
                  xout = pmin(1 - 1e-5,
                              pmax(1e-5,
                                   pgamma(q = exp_to_onset,
                                          shape = inc_parms$shape,
                                          scale = inc_parms$scale))))$y
  )

incubation_times <- incubation_times %>% 
  mutate(onset     = exp_to_onset,
         inf_start = onset - pre_symp_lead,
         inf_end   = ifelse(type == "asymptomatic",
                            exp_to_onset + onset_to_recov,
                            exp_to_onset + inf_from_onset),
         symp_end  = ifelse(type == "asymptomatic",
                            onset, # but really never matters because asymptomatics are never symptomatic!
                            exp_to_onset + onset_to_recov),
         inf_dur   = inf_end - inf_start,
         symp_dur  = symp_end - onset)

times <- seq(0, 35, by=0.1)
infectious_times <- tibble()
for (i in 1:length(times)) {
  infectious_times <- bind_rows(infectious_times, 
    incubation_times %>%
      mutate(time = times[i]) %>%
      filter(inf_start <= time & time <= inf_end)
    )
}
hist(infectious_times$time)

d_infectious_Clifford_unnormalized <- approxfun(density(infectious_times$time, 
  from=-100, to=100))

inf_Clifford_norm <- sum(d_infectious_Clifford_unnormalized(-100:100))

d_infectious_Clifford <- function(t) d_infectious_Clifford_unnormalized(t) / inf_Clifford_norm

sum(d_infectious_Clifford(-10:30))

save(d_infectious_Clifford, d_infectious_Clifford_unnormalized, inf_Clifford_norm,
  file="Clifford dynamics/Clifford_infectiousness.RData")

dpi <- seq(0, 35, by=0.1)
plot(dpi, d_infectious_Clifford(dpi), type='l')


