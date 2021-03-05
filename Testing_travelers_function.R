require(testthat)

generate_p_symp_default <- function(prop_asx = 0.3, meanlog = 1.6, sdlog = 0.5) {
  function(t) (1 - prop_asx) * plnorm(t, meanlog, sdlog)
}

generate_p_test_pos_ag <- function(d_infectious, max_sensitivity = 0.80) {
  function(t) {
    inf_max <- optimize(d_infectious, c(0, 10), maximum=T)$objective
    max_sensitivity * d_infectious(t)/inf_max
  }
}
  
travel_trans_reduction <- function(
    test_time=NA,
    infection_window = c(0, NA), # relative to travel
    transmission_window = c(0, 28), # relative to travel
    quarantine_length = NA,  # relative to travel
    quarantine_effectiveness = 1,
    d_infectious = function(t) dgamma(t, shape=7.2, scale=0.80),
    p_symp = NA,
    p_test_pos = NA,
    additional_tests = NA
    ) {
  if ((!is.na(test_time[1]) | !is.na(additional_tests[1])) & !is.function(p_test_pos)) {
    stop("must specify p_test_pos()")
  }
  symptom_monitoring <- is.function(p_symp)
  infectiousness_fxn <- 
    if (is.na(infection_window[2])) {
      function(t, tau) d_infectious(t - tau)
    } else {
      function(t, tau) {
        dunif(tau, infection_window[1], infection_window[2]) * d_infectious(t - tau)
      }
    }
  p_not_detected <-
    if (symptom_monitoring & is.na(test_time[1]) & is.na(quarantine_length)) {
      function(t, tau) 1 - p_symp(t - tau)
    } else if (!symptom_monitoring & is.na(test_time[1]) & !is.na(quarantine_length)) {
      function(t, tau) {
        (1 - quarantine_effectiveness * (as.numeric(0 < t & t < quarantine_length)))
      }
    } else if (!symptom_monitoring & !is.na(test_time[1]) & is.na(quarantine_length)) {
      function(t, tau) {
        (1 - (as.numeric(t > test_time) * p_test_pos(test_time - tau)))
      }
    } else if (symptom_monitoring & !is.na(test_time[1]) & is.na(quarantine_length)) {
      function(t, tau) {
        (1 - p_symp(t - tau)) *
        (1 - (as.numeric(t > test_time) * p_test_pos(test_time - tau)))
      }
    } else if (symptom_monitoring & is.na(test_time[1]) & !is.na(quarantine_length)) {
      function(t, tau) {
        (1 - p_symp(t - tau)) *
        (1 - quarantine_effectiveness * (as.numeric(0 < t & t < quarantine_length)))
      }
    } else if (!symptom_monitoring & !is.na(test_time[1]) & !is.na(quarantine_length)) {
      function(t, tau) {
        (1 - (as.numeric(t > test_time) * p_test_pos(test_time - tau))) * 
        (1 - quarantine_effectiveness * (as.numeric(0 < t & t < quarantine_length)))
      }
    } else if (symptom_monitoring & !is.na(test_time[1]) & !is.na(quarantine_length)) {
      function(t, tau) {
        (1 - p_symp(t - tau)) *
        (1 - (as.numeric(t > test_time) * p_test_pos(test_time - tau))) *
        (1 - quarantine_effectiveness * (as.numeric(0 < t & t < quarantine_length)))
      }
    } else {
      function(t, tau) { 1 }
    }
  p_not_detected_multitest <- function(t, tau) {
    p_not_detected(t, tau) * ifelse(is.na(additional_tests[1]), 1, 
        prod(1 - (as.numeric(t > additional_tests) * p_test_pos(additional_tests - tau))))
  }
  p_not_detected_multitest <- Vectorize(p_not_detected_multitest, 't')

  if (is.na(infection_window[2])) {
    risk_prevented <- integrate(function(t) infectiousness_fxn(t, infection_window[1]) * 
            (1 - p_not_detected_multitest(t, infection_window[1])),
          lower=transmission_window[1], upper=transmission_window[2], abs.tol=0.001)$value
    total_risk <- 
        integrate(function(t) infectiousness_fxn(t, infection_window[1]),
          lower=transmission_window[1], upper=transmission_window[2], abs.tol=0.001)$value
  } else {
    risk_prevented <- integrate(function(tau) {
        sapply(tau, function(tau) {
          integrate(function(t) infectiousness_fxn(t, tau) * 
              (1 - p_not_detected_multitest(t, tau)),
            lower=transmission_window[1], upper=transmission_window[2], abs.tol=0.001)$value
        })
      }, lower=infection_window[1], upper=infection_window[2], abs.tol=0.001)$value
    total_risk <- integrate(function(tau) {
      sapply(tau, function(tau) {
        integrate(function(t) infectiousness_fxn(t, tau),
          lower=transmission_window[1], upper=transmission_window[2], abs.tol=0.001)$value
      })
    }, lower=infection_window[1], upper=infection_window[2], abs.tol=0.001)$value
  }
  return(risk_prevented / total_risk)
}

travel_trans_reduction_vec <- 
  Vectorize(travel_trans_reduction, vectorize.args = 'test_time')


### function tests
p_test_pos_test <- generate_p_test_pos_ag(function(t) dgamma(t, shape=7.2, scale=0.80), 0.99)

test_that("removing symptomatics works", {
  expect_equivalent(round(travel_trans_reduction(
    p_symp = generate_p_symp_default()), 3), 0.389)
})

test_that("removing test positive works", {
  expect_equivalent(round(travel_trans_reduction(3, 
    p_test_pos = p_test_pos_test), 3), 0.471)
})

test_that("removing test positive & symptomatics works", {
  expect_equivalent(round(travel_trans_reduction(3, 
    p_symp = generate_p_symp_default(prop_asx = 0.3, meanlog = 1.6, sdlog = 0.5), 
    p_test_pos = p_test_pos_test), 3), 0.664)
})

test_that("vectorization works", {
  expect_equivalent(round(travel_trans_reduction_vec(c(0, 5, 10), 
    p_symp = generate_p_symp_default(prop_asx = 0.3, meanlog = 1.6, sdlog = 0.5), 
    p_test_pos = p_test_pos_test), 3),
    c(0.389, 0.679, 0.391))
})

test_that("quarantine + removing test positive & symptomatics works", {
  expect_equivalent(round(travel_trans_reduction_vec(c(4, 8), 
    p_symp = generate_p_symp_default(prop_asx = 0.3, meanlog = 1.6, sdlog = 0.5), 
    p_test_pos = p_test_pos_test, 
    quarantine_length = 7), 3), c(0.986, 0.921))
})

test_that("removing test positive & symptomatics works + 2nd test", {
  expect_equivalent(round(travel_trans_reduction_vec(c(4, 8), 
    p_symp = generate_p_symp_default(prop_asx = 0.3, meanlog = 1.6, sdlog = 0.5), 
    p_test_pos = p_test_pos_test, 
    additional_tests = 0), 3), c(0.759, 0.413))
})

test_that("multiple tests", {
  expect_equivalent(round(travel_trans_reduction_vec(c(3, 5), 
    p_symp = generate_p_symp_default(), 
    p_test_pos = p_test_pos_test, 
    additional_tests = c(0, 5)), 3), c(0.807, 0.682))
})

test_that("multiple tests v2", {
  expect_equivalent(round(travel_trans_reduction(
    p_symp = generate_p_symp_default(), 
    p_test_pos = p_test_pos_test, 
    additional_tests = c(0, 3, 5)), 3), 0.807)
})

