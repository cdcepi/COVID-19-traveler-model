rm(list = ls())

source('Testing_travelers_function.R')

load("Clifford dynamics/Clifford_RTPCR.RData") # p_test_pos_Clifford()
test_options <- c('RTPCR', 'Ag')

load("Goyal dynamics/Goyal_infectiousness.RData")
load("Clifford dynamics/Clifford_infectiousness.RData")
d_inf_Gamma <- function(t) dgamma(t, shape=7.2, scale=0.80)

times <- seq(0, 30, by=0.1)

png('estimates/Fig1.png', 6, 7, units='i', res=256)
# A
par(mfrow=c(2, 1), mar=c(3, 3, 1, 0.5), mgp=c(1.5, 0.3, 0), tck=-0.02)
plot(times, d_inf_Gamma(times), type='l', lwd=2, xlim=c(0, 21), col='black', 
  xlab="Time since infection (days)", 
  ylab="Relative infectiousness", bty='n', axes=F)
lines(times, d_infectious_Goyal(times), lwd=2, col='darkorange')
lines(times, d_infectious_Clifford(times), lwd=2, col='blue')
axis(1)

legend('topright', legend=c('Gamma', 
  'Goyal et al.', 'Clifford et al.'), 
  lty=c(1, 1, 1), lwd=2, col=c('black', 'darkorange', 'blue'), bty='n')
mtext("A", 3, 0, adj=-0.1, font=2)

# B
plot(times, p_test_pos_Clifford(times), type='l', lwd=2, xlim=c(0, 21), col='darkgreen', 
  xlab="Time since infection (days)", 
  ylab="Probability of positive test (%)", bty='n', ylim=c(0, 1))
lines(times, generate_p_test_pos_ag(d_inf_Gamma, 0.80)(times), lwd=2, lty=3, col='black')
lines(times, generate_p_test_pos_ag(d_infectious_Goyal, 0.80)(times), lwd=2, lty=3, col='darkorange')
lines(times, generate_p_test_pos_ag(d_infectious_Clifford, 0.80)(times), lwd=2, lty=3, col='blue')
legend('topright', legend=c('RT-PCR', 'Ag 80% sens. - Gamma', 
  'Ag 80% sens. - Goyal et al.', 'Ag 80% sens. - Clifford et al.'), 
  lty=c(1, 3, 3, 3), lwd=2,  col=c('darkgreen', 'black', 'darkorange', 'blue'), bty='n')
mtext("B", 3, 0, adj=-0.1, font=2)
dev.off()

