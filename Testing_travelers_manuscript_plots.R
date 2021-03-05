rm(list = ls())

require(tidyverse)

source('Testing_travelers_function.R')

load("Clifford dynamics/Clifford_RTPCR.RData") # p_test_pos_Clifford()

test_options <- c('RTPCR', 'Ag')

load("Goyal dynamics/Goyal_infectiousness.RData")
load("Clifford dynamics/Clifford_infectiousness.RData")

d_inf_options <- list(
  gamma = function(t) dgamma(t, shape=7.2, scale=0.80), 
  Goyal = d_infectious_Goyal,
  Clifford = d_infectious_Clifford
)

#Manuscript Results - Reducing transmission risk after a known exposure - Fig 2
protocols2 <- tibble()
for (inf_i in 1:length(d_inf_options)) {
  for (this_test in test_options) {
    for (this_q in c(0, 7, 14)) {
      for (this_adherence in c(1)) {
        for (this_symp_monitor in c(T, F)) {
          for (this_test_time in c(NA, 1:7)) {
            protocols2 <- bind_rows(protocols2,
               tibble(
                 q = this_q,
                 q_effectiveness = this_adherence,
                 symp_monitor = this_symp_monitor,
                 test_time = this_test_time,
                 test = this_test,
                 inf = names(d_inf_options)[inf_i],
                 reduction = travel_trans_reduction(
                   test_time = this_test_time,
                   quarantine_length = this_q,
                   quarantine_effectiveness = this_adherence,
                   d_infectious = d_inf_options[[inf_i]],
                   p_symp = ifelse(this_symp_monitor, generate_p_symp_default(), NA), 
                   p_test_pos = ifelse(this_test == 'RTPCR', p_test_pos_Clifford,
                     generate_p_test_pos_ag(d_inf_options[[inf_i]], 0.80)),
                   infection_window = c(0, NA))
               )
            )
          }
        }
      }
    }
  }
}

#save(protocols2, file="estimates/fig2protocols.rdata")
load("estimates/fig2protocols.rdata")

protocol_table2 <- protocols2 %>%
  mutate(any_test = ifelse(is.na(test_time), F, T)) %>%
  group_by(q, q_effectiveness, test_time, symp_monitor) %>%
  summarize(
    median_reduction = median(reduction),
    mean_reduction = mean(reduction),
    min_reduction = min(reduction),
    max_reduction = max(reduction))

### Figure 2a
#symp alone, quar alone, testing alone, symp + testing 
quar_alone7 = subset(protocol_table2, q == 7 & is.na(test_time) & symp_monitor==F)
quar_alone14 = subset(protocol_table2, q == 14 & is.na(test_time) & symp_monitor==F)
symp_alone = subset(protocol_table2, q == 0 & is.na(test_time) & symp_monitor==T)
test_alone = subset(protocol_table2, q == 0 & test_time %in% 1:7 & symp_monitor==F)
test_symp = subset(protocol_table2, q == 0 & test_time %in% 1:7 & symp_monitor==T)
test_quar = subset(protocol_table2, q == 7 & test_time %in% 1:7 & symp_monitor==F)
test_quar_symp = subset(protocol_table2, q == 7 & test_time %in% 1:7 & symp_monitor==T)

ds = rbind(symp_alone, quar_alone7, quar_alone14, test_alone, test_symp, test_quar, test_quar_symp)

bards = cbind(c(ds$median_reduction[1], ds$median_reduction[2], ds$median_reduction[3], 
  rep(NA, 4)), matrix(c(rep(0, 21), ds$median_reduction[4:31]), ncol=7, byrow=T))
minds = cbind(c(ds$min_reduction[1], ds$min_reduction[2], ds$min_reduction[3], 
  rep(NA, 4)), matrix(c(rep(0, 21), ds$min_reduction[4:31]), ncol=7, byrow=T))
maxds =cbind(c(ds$max_reduction[1], ds$max_reduction[2], ds$max_reduction[3], 
  rep(NA, 4)), matrix(c(rep(0, 21), ds$max_reduction[4:31]), ncol=7, byrow=T))
rownames(bards) <- rownames(minds) <- rownames(maxds) <- 
  c("Symp. alone", "Quar. alone7", "Quar. alone14", "Test. alone", "Test/symp. screen", "Test/quar.", "Test/quar/symp")
colnames(bards) <- colnames(minds) <- colnames(maxds) <- 
  c("No test", "1", "2", "3", "4", "5","6", "7")

png('estimates/Fig2.png', 6, 4, units='i', res=256)
par(mar=c(3, 3, 0.5, 0.5), mgp=c(1.5, 0.2, 0), tck=-0.02)
fig2 = barplot(bards*100, beside=T, ylim=c(-0.5, 125), 
    col=c("gold", "tomato3", adjustcolor("tomato3", 0.5), "blue", "seagreen", "purple4", "maroon"), 
    xlab="Time from infection to test (days)",
    ylab="Reduction in transmission risk (%)", cex.names=1, axes=F)
abline(h=0, lwd=2, col="black")
arrows(x0=fig2, y0=maxds*100, y1=minds*100, lwd=1.5, angle = 90, code=3, length = 0.03)
axis(2, at=seq(0, 100, by=20))
legend(x=0, y=128, c("Symptom monitoring", "7d quarantine", "14d quarantine", "",
    "Test", "Test + symptom monitoring", "Test + 7d quarantine", "Test + symp. monitoring + 7d quarantine"), 
  border=c(rep('black', 3), 'white', rep('black', 4),
    "blue", "seagreen", "purple4", "maroon"),
  fill=c("gold", "tomato3", adjustcolor("tomato3", 0.5), adjustcolor('white', 0), 
      "blue", "seagreen", "purple4", "maroon"), 
    cex=0.75, bty='n', horiz=F, ncol=2)
dev.off()


#Manuscript Results - Transmission risk during travel - Fig 3
protocols3a <- tibble()
for (inf_i in 1:length(d_inf_options)) {
  for (this_test in test_options) {
    for (this_symp_monitor in c(T, F)) {
      for (this_test_time in c(NA, -7:-1)) {
        protocols3a <- bind_rows(protocols3a,
           tibble(
             symp_monitor = this_symp_monitor,
             test_time = this_test_time,
             test = this_test,
             inf = names(d_inf_options)[inf_i],
             reduction = travel_trans_reduction(
               test_time = this_test_time,
               infection_window = c(-8, -1),
               transmission_window = c(-1, 0),
               d_infectious = d_inf_options[[inf_i]],
               p_symp = ifelse(this_symp_monitor, generate_p_symp_default(), NA),
               p_test_pos = ifelse(this_test == 'RTPCR', p_test_pos_Clifford,
                 generate_p_test_pos_ag(d_inf_options[[inf_i]], 0.80))
             )
           )
        )
      }
    }
  }
}


#Fig 3b
protocols3b <- tibble()
for (inf_i in 1:length(d_inf_options)) {
  for (this_test in "Ag") {
    for (this_max_sens in c(0.8, 0.95)) {
      for (this_symp_monitor in c(F)) {
        for (this_test_time in seq(-7, -1, by=0.2)) {
          protocols3b <- bind_rows(protocols3b,
            tibble(
              symp_monitor = this_symp_monitor,
              test_time = this_test_time,
              test = this_test,
              test_sens = this_max_sens,
              inf = names(d_inf_options)[inf_i],
              reduction = travel_trans_reduction(
                test_time = this_test_time,
                infection_window = c(-8, -1),
                transmission_window = c(-1, 0),
                d_infectious = d_inf_options[[inf_i]],
                p_symp = ifelse(this_symp_monitor, generate_p_symp_default(), NA),
                p_test_pos = ifelse(this_test == 'RTPCR', p_test_pos_Clifford,
                  generate_p_test_pos_ag(d_inf_options[[inf_i]], this_max_sens))
              )
            )
          )
        }
      }
    }
  }
}

#save(protocols3a, protocols3b, file="estimates/fig3protocols.rdata")
load("estimates/fig3protocols.rdata")

protocol_table3a <- protocols3a %>%
  mutate(any_test = ifelse(is.na(test_time), F, T)) %>%
  group_by(test_time, symp_monitor) %>%
  summarize(
    median_reduction = median(reduction),
    mean_reduction = mean(reduction),
    min_reduction = min(reduction),
    max_reduction = max(reduction))

#symp alone, quar alone, testing alone, symp + testing 
symp_alone = subset(protocol_table3a, is.na(test_time) & symp_monitor == T)
test_alone = subset(protocol_table3a, test_time %in% -7:-1 & symp_monitor == F)
test_symp = subset(protocol_table3a, test_time %in% -7:-1  & symp_monitor == T)

ds = rbind(symp_alone, test_alone, test_symp)

bards3a = cbind(c(ds$median_reduction[1], rep(NA, 2)), matrix(c(rep(0, 7), ds$median_reduction[2:15]), ncol=7, byrow=T))
minds3a = cbind(c(ds$min_reduction[1], rep(NA, 2)), matrix(c(rep(0, 7), ds$min_reduction[2:15]), ncol=7, byrow=T))
maxds3a = cbind(c(ds$max_reduction[1], rep(NA, 2)), matrix(c(rep(0, 7), ds$max_reduction[2:15]), ncol=7, byrow=T))
rownames(bards3a) <- rownames(minds3a) <- rownames(maxds3a) <- 
  c("Symp. alone", "Test. alone", "Test/symp. screen")
colnames(bards3a) <- colnames(minds3a) <- colnames(maxds3a) <- 
  c("No test", "-6", "-5", "-4", "-3", "-2", "-1", "0")

protocol_table3b <- protocols3b %>%
  group_by(test_time, test_sens) %>%
  summarize(
    median_reduction = median(reduction),
    mean_reduction = mean(reduction),
    min_reduction = min(reduction),
    max_reduction = max(reduction))

test_alone95 = subset(protocol_table3b, test_sens == 0.95)
test_alone80 = subset(protocol_table3b, test_sens == 0.80)

ds = rbind(test_alone95, test_alone80)
test_times <- unique(protocol_table3b$test_time)

bards = matrix(ds$median_reduction, ncol=length(test_times), byrow=T)
minds = matrix(ds$min_reduction, ncol=length(test_times), byrow=T)
maxds = matrix(ds$max_reduction, ncol=length(test_times), byrow=T)
rownames(bards) <- rownames(minds) <- rownames(maxds) <- c("Test (95% sensitivity)", "Test (80% sensitivity)")
colnames(bards) <- colnames(minds) <- colnames(maxds) <- test_times

png('estimates/Fig3.png', 6, 7, units='i', res=256)
# 3A
par(mfrow=c(2, 1), mar=c(3, 3, 1, 0.5), mgp=c(1.5, 0.3, 0), tck=-0.02)
fig3a = barplot(bards3a*100, beside=T, ylim=c(-0.5, 105), 
  col=c("gold", "blue", "seagreen"), 
  xlab="Time from test to departure (days)",
  ylab="Reduction in transmission risk (%)", cex.names=1)
abline(h=0, lwd=2, col="black")
arrows(x0=fig3a, y0=maxds3a*100, y1=minds3a*100, lwd=1.5, angle = 90, code=3, length = 0.05)
legend("topleft", c("Symp. monitoring", "Test alone", "Test + symp. monitoring"), 
  fill=c("gold", "blue", "seagreen"), cex=1, bty='n')
mtext("A", 3, 0, adj=-0.1, font=2)
# 3B
par(mar=c(3, 3, 0.5, 0.5), mgp=c(1.5, 0.5, 0))
plot(test_times + 1, 100 * bards[1,], type="l", xlim=c(-6, 0.25), ylim=c(0, 100),
     lwd=2, col="darkorange", xlab="Time from test to departure (days)", 
     ylab="Reduction in transmission risk (%)", bty='n')
lines(test_times + 1, 100 * bards[2,], lwd=2, col='blue')
polygon(c(test_times + 1, rev(test_times + 1)), c(100*minds[1,], rev(100*maxds[1,])), 
  col=adjustcolor('darkorange', 0.75), border=adjustcolor('darkorange', 0.75))
polygon(c(test_times + 1, rev(test_times + 1)), c(100*minds[2,], rev(100*maxds[2,])), 
  col=adjustcolor('blue', 0.5), border=adjustcolor('blue', 0.5))
legend("topleft",
       c("95% sensitivity antigen test",  "80% sensitivity antigen test"),
       fill=c("darkorange", "blue"), bty="n", cex=1)
mtext("B", 3, 0, adj=-0.1, font=2)
dev.off()

#Manuscript - transmission risk after travel - Fig 4
protocols4 <- tibble()
for (inf_i in 1:length(d_inf_options)) {
  for (this_test in test_options) {
    for (this_q in c(0, 7)) {
      for (this_symp_monitor in c(T, F)) {
        for (this_test_time in c(NA, -1:7)) {
          protocols4 <- bind_rows(protocols4,
            tibble(
              q = this_q,
              symp_monitor = this_symp_monitor,
              test_time = this_test_time,
              test = this_test,
              inf = names(d_inf_options)[inf_i],
              reduction = travel_trans_reduction(
                test_time = this_test_time,
                quarantine_length = this_q,
                infection_window = c(-7, 0),
                transmission_window = c(0, 28),
                d_infectious = d_inf_options[[inf_i]],
                p_symp = ifelse(this_symp_monitor, generate_p_symp_default(), NA),
                p_test_pos = ifelse(this_test == 'RTPCR', p_test_pos_Clifford,
                  generate_p_test_pos_ag(d_inf_options[[inf_i]], 0.80))
                )
              )
            )
        }
      }
    }
  }
}

#save(protocols4, file="estimates/fig4protocols.rdata")
load("estimates/fig4protocols.rdata")

protocol_table4 <- protocols4 %>%
  mutate(any_test = ifelse(is.na(test_time), F, T)) %>%
  group_by(test_time, symp_monitor, q) %>%
  summarize(
    median_reduction = median(reduction),
    mean_reduction = mean(reduction),
    min_reduction = min(reduction),
    max_reduction = max(reduction))


#Figure 4
symp_alone = subset(protocol_table4, q == 0 & is.na(test_time) & symp_monitor == T)
quar_alone7 = subset(protocol_table4, q == 7 & is.na(test_time) & symp_monitor == F)
quar_symp7 = subset(protocol_table4, q == 7 & is.na(test_time) & symp_monitor == T)
test_alone = subset(protocol_table4, q == 0 & test_time %in% -1:7 & symp_monitor == F)
test_symp = subset(protocol_table4, q == 0 & test_time %in% -1:7 & symp_monitor == T)
test_quar7 = subset(protocol_table4, q == 7 & test_time %in% -1:7 & symp_monitor == F)
test_symp_quar7 = subset(protocol_table4, q == 7 & test_time %in% -1:7 & symp_monitor == T)

ds = rbind(symp_alone, quar_alone7, quar_symp7, test_alone, test_symp, test_quar7, test_symp_quar7)

bards = cbind(c(ds$median_reduction[1:3], rep(NA, 4)), 
  matrix(c(rep(0, 27), ds$median_reduction[4:39]), ncol=9, byrow=T))

minds = cbind(c(ds$min_reduction[1], ds$min_reduction[2], ds$min_reduction[3], rep(NA, 4)), matrix(c(rep(0, 27), ds$min_reduction[4:39]), ncol=9, byrow=T))
maxds =cbind(c(ds$max_reduction[1], ds$max_reduction[2], ds$max_reduction[3], rep(NA, 4)), matrix(c(rep(0, 27), ds$max_reduction[4:39]), ncol=9, byrow=T))
rownames(bards) <- rownames(minds) <- rownames(maxds) <- 
  c("Symp. alone", "Quar. alone7", "Quar. symp7", "Test. alone", 
    "Test/symp. screen", "Test/quar.", "Test/symp/quar")
colnames(bards) <- colnames(minds) <- colnames(maxds) <- 
  c("No test", "-1", "0", "1", "2", "3", "4", "5", "6", "7")


png('estimates/Fig4.png', 6, 4, units='i', res=256)
par(mar=c(3, 3, 0.5, 0.5), mgp=c(1.5, 0.5, 0))
fig4 = barplot(bards*100, beside=T, ylim=c(-0.5, 125), 
        col=c("gold", "tomato3", "darkorange", "blue", "seagreen", "purple4", "maroon"), 
        xlab="Time from arrival to test (days)",
        ylab="Reduction in transmission risk (%)", cex.names=1, axes=F)
abline(h=0, lwd=2, col="black")
arrows(x0=fig4, y0=maxds*100, y1=minds*100, lwd=1.5, angle = 90, code=3, length = 0.02)
axis(2, at=seq(0, 100, by=20))
legend(x=0, y=128, c("Symptom monitoring", "7d quarantine", "7d quarantine + symptom monitoring", "",
    "Test", "Test + symptom monitoring", "Test + 7d quarantine", "Test + symp. monitoring + 7d quarantine"), 
  border=c(rep('black', 3), adjustcolor('white', 0), rep('black', 4)), 
  fill=c("gold", "tomato3", "darkorange", adjustcolor('white', 0),
      "blue", "seagreen", "purple4", "maroon"), 
  cex=0.75, bty='n', horiz=F, ncol=2)
dev.off()

#Manuscript - transmission risk after travel, varying quarantine adherence - Fig 5
protocols5 <- tibble()
for (inf_i in 1:length(d_inf_options)) {
  for (this_test in test_options) {
    for (this_q in c(7, 10, 14)) {
      for (this_adherence in c(1, 0.5)) {
        for (this_symp_monitor in c(T)) {
          for (this_test_time in c(NA, -1:7)) { 
            protocols5 <- bind_rows(protocols5,
               tibble(
                 q = this_q,
                 q_effectiveness = this_adherence,
                 symp_monitor = this_symp_monitor,
                 test_time = this_test_time,
                 test = this_test,
                 inf = names(d_inf_options)[inf_i],
                 reduction = travel_trans_reduction(
                   test_time = this_test_time,
                   quarantine_length = this_q,
                   quarantine_effectiveness = this_adherence,
                   infection_window = c(-7, 0),
                   transmission_window = c(0, 28),
                   d_infectious = d_inf_options[[inf_i]],
                   p_symp = ifelse(this_symp_monitor, generate_p_symp_default(), NA),
                   p_test_pos = ifelse(this_test == 'RTPCR', p_test_pos_Clifford,
                     generate_p_test_pos_ag(d_inf_options[[inf_i]], 0.80))
                 )
              )
            )
          }
        }
      }
    }
  }
}

# save(protocols5, file="estimates/fig5protocols.rdata")
load("estimates/fig5protocols.rdata")

protocol_table5 <- protocols5 %>%
  filter(test != 'default_95') %>%
  mutate(any_test = ifelse(is.na(test_time), F, T)) %>%
  group_by(test_time, symp_monitor, q, q_effectiveness) %>%
  summarize(
    median_reduction = median(reduction),
    mean_reduction = mean(reduction),
    min_reduction = min(reduction),
    max_reduction = max(reduction)) %>%
  ungroup()

ds = filter(protocol_table5, symp_monitor == T) %>%
  ungroup() %>%
  mutate(
    no_test = is.na(test_time)) %>% 
  arrange(desc(no_test), test_time, desc(q_effectiveness), q)

# bards = cbind(c(ds$median_reduction[1], ds$median_reduction[2], ds$median_reduction[3], 
#   rep(NA, 8)), matrix(c(rep(0, 12), ds$median_reduction[4:27]), ncol=5, byrow=T))
bards = matrix(ds$median_reduction, ncol=10, byrow=F)
minds = matrix(ds$min_reduction, ncol=10, byrow=F)
maxds = matrix(ds$max_reduction, ncol=10, byrow=F)
rownames(bards) <- rownames(minds) <- rownames(maxds) <- 
  c("7 day quar. 100% adh.", "10 day quar. 100% adh.", "14 day quar. 100% adh.", 
    "7 day quar. 50% adh.", "10 day quar. 50% adh.", "14 day quar. 50% adh.")
colnames(bards) <- colnames(minds) <- colnames(maxds) <- 
  c("No test", "-1", "0", "1", "2", "3", "4", "5", "6", "7")

png('estimates/Fig5.png', 6, 4, units='i', res=256)
par(mar=c(3, 3, 0.5, 0.5), mgp=c(1, 0.1, 0), tck=-0.01)
fig5 = barplot(bards*100, beside=T, ylim=c(-0.1, 115), axes=F, 
    col=c("navy", "mediumorchid4", "darkseagreen4", 
      adjustcolor("navy", 0.5), adjustcolor("mediumorchid4", 0.5), adjustcolor("darkseagreen4", 0.5)), 
    xlab="Time from arrival to test (days)",
    ylab="Reduction in transmission risk (%)", cex.names=0.6)
abline(h=0, lwd=2, col="black")
arrows(x0=fig5, y0=maxds*100, y1=minds*100, lwd=1.5, angle = 90, code=3, length = 0.02)
axis(2)
legend('topleft', c("7 day quar. 100% adh.", "10 day quar. 100% adh.", "14 day quar. 100% adh.", 
    "7 day quar. 50% adh.", "10 day quar. 50% adh.", "14 day quar. 50% adh."), 
    fill=c("navy", "mediumorchid4", "darkseagreen4", adjustcolor("navy", 0.5), 
      adjustcolor("mediumorchid4", 0.5), adjustcolor("darkseagreen4", 0.5))[c(1, 4, 2, 5, 3, 6)], 
    cex=0.75, bty='n', horiz=F, ncol=3)
dev.off()

