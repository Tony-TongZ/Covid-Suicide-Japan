############################################################
## Decision on Time Lags
############################################################
Sys.setenv(
  OMP_NUM_THREADS="1",
  OPENBLAS_NUM_THREADS="1",
  MKL_NUM_THREADS="1",
  VECLIB_MAXIMUM_THREADS="1"
)
set.seed(123)

############################################################
## 0. Packages
############################################################
library(splines)
library(brms)
library(cmdstanr)
cmdstanr::check_cmdstan_toolchain()
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(rlang)

############################################################
## 1. suicide vs unemployment & depression
############################################################

# standardization
excess_depression_s<-as.vector(scale(excess_depression, center = TRUE, scale = TRUE))
rate_of_change_s<-as.vector(scale(rate_of_change, center = TRUE, scale = TRUE))

###### Spline-Lag Basis Construction
L <- 7
K <- 4
make_spline_lag <- function(P, L = 7, K = 4) {
  n <- length(P)
  laggedP <- embed(P, L)
  lag_index <- 2:L
  B <- bs(lag_index, df = K, intercept = TRUE)
  S <- laggedP[, lag_index] %*% B
  return(S)
}

exun_df <- cbind(exun_male, exun_female)
SX_list <- lapply(1:ncol(exun_df), function(j) make_spline_lag(exun_df[, j], L = 7, K = 4))
excess_depression_s_2 <- cbind(excess_depression_s, excess_depression_s)
SZ_list <- lapply(1:ncol(excess_depression_s_2), function(j) make_spline_lag(excess_depression_s_2[, j], L = 7, K = 4))

###### Data Preparation with Spline-Lag Covariates
pt <- pt_monthly_means
sui_male_obs2   <- sui_male_obs[7:24, ]
sui_male_pre2   <- sui_male_pre[7:24, ]
sui_female_obs2 <- sui_female_obs[7:24, ]
sui_female_pre2 <- sui_female_pre[7:24, ]

exsui_male_new <- data.frame(
  t = rep(1:18, 6),
  Yobs = c(sui_male_obs2$au20, sui_male_obs2$a20, sui_male_obs2$a30,
           sui_male_obs2$a40,  sui_male_obs2$a50, sui_male_obs2$a60),
  Ypre = c(sui_male_pre2$au20, sui_male_pre2$a20, sui_male_pre2$a30,
           sui_male_pre2$a40,  sui_male_pre2$a50, sui_male_pre2$a60),
  age = c(rep("u20", 18), rep("20", 18), rep("30", 18),
          rep("40", 18), rep("50", 18), rep("60", 18)),
  gender = rep("Male", 18 * 6)
)

exsui_female_new <- data.frame(
  t = rep(1:18, 6),
  Yobs = c(sui_female_obs2$au20, sui_female_obs2$a20, sui_female_obs2$a30,
           sui_female_obs2$a40,  sui_female_obs2$a50, sui_female_obs2$a60),
  Ypre = c(sui_female_pre2$au20, sui_female_pre2$a20, sui_female_pre2$a30,
           sui_female_pre2$a40,  sui_female_pre2$a50, sui_female_pre2$a60),
  age = c(rep("u20", 18), rep("20", 18), rep("30", 18),
          rep("40", 18), rep("50", 18), rep("60", 18)),
  gender = rep("Female", 18 * 6)
)

exun_male_new <- data.frame(
  t = rep(1:18, 6),
  X = c(exun_male$a10[7:24], exun_male$a20[7:24], exun_male$a30[7:24],
        exun_male$a40[7:24], exun_male$a50[7:24], exun_male$a60[7:24]),
  age = c(rep("u20", 18), rep("20", 18), rep("30", 18),
          rep("40", 18), rep("50", 18), rep("60", 18)),
  gender = rep("Male", 18 * 6)
)

exun_female_new <- data.frame(
  t = rep(1:18, 6),
  X = c(exun_female$a10[7:24], exun_female$a20[7:24], exun_female$a30[7:24],
        exun_female$a40[7:24], exun_female$a50[7:24], exun_female$a60[7:24]),
  age = c(rep("u20", 18), rep("20", 18), rep("30", 18),
          rep("40", 18), rep("50", 18), rep("60", 18)),
  gender = rep("Female", 18 * 6)
)

data1 <- data.frame(
  t = c(exun_male_new$t, exun_female_new$t),
  X = c(exun_male_new$X, exun_female_new$X),
  Yobs = c(exsui_male_new$Yobs, exsui_female_new$Yobs),
  Ypre = c(exsui_male_new$Ypre, exsui_female_new$Ypre),
  Z = rep(excess_depression_s[1:18], 12),
  V = rep(rate_of_change_s[1:18], 12),
  age = c(exun_male_new$age, exun_female_new$age),
  gender = c(exun_male_new$gender, exun_female_new$gender),
  SX_1 = c(SX_list[[1]][, 1], SX_list[[2]][, 1], SX_list[[3]][, 1], SX_list[[4]][, 1], SX_list[[5]][, 1], SX_list[[6]][, 1],
           SX_list[[7]][, 1], SX_list[[8]][, 1], SX_list[[9]][, 1], SX_list[[10]][, 1], SX_list[[11]][, 1], SX_list[[12]][, 1]),
  SX_2 = c(SX_list[[1]][, 2], SX_list[[2]][, 2], SX_list[[3]][, 2], SX_list[[4]][, 2], SX_list[[5]][, 2], SX_list[[6]][, 2],
           SX_list[[7]][, 2], SX_list[[8]][, 2], SX_list[[9]][, 2], SX_list[[10]][, 2], SX_list[[11]][, 2], SX_list[[12]][, 2]),
  SX_3 = c(SX_list[[1]][, 3], SX_list[[2]][, 3], SX_list[[3]][, 3], SX_list[[4]][, 3], SX_list[[5]][, 3], SX_list[[6]][, 3],
           SX_list[[7]][, 3], SX_list[[8]][, 3], SX_list[[9]][, 3], SX_list[[10]][, 3], SX_list[[11]][, 3], SX_list[[12]][, 3]),
  SX_4 = c(SX_list[[1]][, 4], SX_list[[2]][, 4], SX_list[[3]][, 4], SX_list[[4]][, 4], SX_list[[5]][, 4], SX_list[[6]][, 4],
           SX_list[[7]][, 4], SX_list[[8]][, 4], SX_list[[9]][, 4], SX_list[[10]][, 4], SX_list[[11]][, 4], SX_list[[12]][, 4]),
  SZ_1 = rep(SZ_list[[1]][, 1], 12),
  SZ_2 = rep(SZ_list[[1]][, 2], 12),
  SZ_3 = rep(SZ_list[[1]][, 3], 12),
  SZ_4 = rep(SZ_list[[1]][, 4], 12)
)

data1_new <- data.frame(
  t = c(exun_male_new$t, exun_female_new$t, rep(NA, 6 * 12)),
  Yobs = c(exsui_male_new$Yobs, exsui_female_new$Yobs, rep(NA, 6 * 12)),
  Ypre = c(exsui_male_new$Ypre, exsui_female_new$Ypre, rep(NA, 6 * 12)),
  Z = rep(excess_depression_s[1:24], 12),
  V = rep(rate_of_change_s[1:24], 12),
  age = c(exun_male_new$age, exun_female_new$age, rep(NA, 6 * 12)),
  gender = c(exun_male_new$gender, exun_female_new$gender, rep(NA, 6 * 12)),
  SX_1 = c(SX_list[[1]][, 1], SX_list[[2]][, 1], SX_list[[3]][, 1], SX_list[[4]][, 1], SX_list[[5]][, 1], SX_list[[6]][, 1],
           SX_list[[7]][, 1], SX_list[[8]][, 1], SX_list[[9]][, 1], SX_list[[10]][, 1], SX_list[[11]][, 1], SX_list[[12]][, 1], rep(NA, 6 * 12)),
  SX_2 = c(SX_list[[1]][, 2], SX_list[[2]][, 2], SX_list[[3]][, 2], SX_list[[4]][, 2], SX_list[[5]][, 2], SX_list[[6]][, 2],
           SX_list[[7]][, 2], SX_list[[8]][, 2], SX_list[[9]][, 2], SX_list[[10]][, 2], SX_list[[11]][, 2], SX_list[[12]][, 2], rep(NA, 6 * 12)),
  SX_3 = c(SX_list[[1]][, 3], SX_list[[2]][, 3], SX_list[[3]][, 3], SX_list[[4]][, 3], SX_list[[5]][, 3], SX_list[[6]][, 3],
           SX_list[[7]][, 3], SX_list[[8]][, 3], SX_list[[9]][, 3], SX_list[[10]][, 3], SX_list[[11]][, 3], SX_list[[12]][, 3], rep(NA, 6 * 12)),
  SX_4 = c(SX_list[[1]][, 4], SX_list[[2]][, 4], SX_list[[3]][, 4], SX_list[[4]][, 4], SX_list[[5]][, 4], SX_list[[6]][, 4],
           SX_list[[7]][, 4], SX_list[[8]][, 4], SX_list[[9]][, 4], SX_list[[10]][, 4], SX_list[[11]][, 4], SX_list[[12]][, 4], rep(NA, 6 * 12)),
  SZ_1 = c(rep(SZ_list[[1]][, 1], 12), rep(NA, 6 * 12)),
  SZ_2 = c(rep(SZ_list[[1]][, 2], 12), rep(NA, 6 * 12)),
  SZ_3 = c(rep(SZ_list[[1]][, 3], 12), rep(NA, 6 * 12)),
  SZ_4 = c(rep(SZ_list[[1]][, 4], 12), rep(NA, 6 * 12)),
  X = c(exun_male$a10, exun_male$a20, exun_male$a30, exun_male$a40, exun_male$a50, exun_male$a60,
        exun_female$a10, exun_female$a20, exun_female$a30, exun_female$a40, exun_female$a50, exun_female$a60),
  P = c(pt$V1, pt$V2, pt$V3,
        pt$V4, pt$V5, pt$V6,
        pt$V8, pt$V9, pt$V10,
        pt$V11, pt$V12, pt$V13),
  gender_X = c(rep("Male", 24 * 6), rep("Female", 24 * 6)),
  age_X = c(rep("u20", 24), rep("20", 24), rep("30", 24),
            rep("40", 24), rep("50", 24), rep("60", 24),
            rep("u20", 24), rep("20", 24), rep("30", 24),
            rep("40", 24), rep("50", 24), rep("60", 24))
)

###### Joint Bayesian Model: Spline-Lag Distributed Effects
fit_joint1 <- brm(
  bf(X ~ P + (P | age_X) + (P | gender_X),family = gaussian()) +
  bf(Z ~ V,family = gaussian()) +
  bf(Yobs ~ Ypre + SX_1+SX_2+SX_3+SX_4 + (SX_1+SX_2+SX_3+SX_4 | age) + (SX_1+SX_2+SX_3+SX_4 | gender) + SZ_1+SZ_2+SZ_3+SZ_4,family = gaussian()) +
  set_rescor(FALSE),
  data = data1_new,
  chains = 4, cores = 4, threads = threading(1),
  seed = 123, iter = 2000,
  control = list(adapt_delta = 0.9999, max_treedepth = 20),
  backend = "cmdstanr"
)

### Lag Effect of X on Y: Fixed + Group-Specific Effects
posterior <- fixef(fit_joint1)

sx_coef <- c(
  posterior["Yobs_SX_1", "Estimate"],
  posterior["Yobs_SX_2", "Estimate"],
  posterior["Yobs_SX_3", "Estimate"],
  posterior["Yobs_SX_4", "Estimate"]
)

ranef <- ranef(fit_joint1)

sx_ranef_age <- sapply(1:K, function(k) {
  ranef$age[, , paste0("Yobs_SX_", k)][, "Estimate"]
})

sx_ranef_gender <- sapply(1:K, function(k) {
  ranef$gender[, , paste0("Yobs_SX_", k)][, "Estimate"]
})

ages <- rownames(sx_ranef_age)
genders <- rownames(sx_ranef_gender)

lag_list <- list()
lag_index <- 2:L
B <- bs(lag_index, df = K, intercept = TRUE)

for (a in ages) {
  for (g in genders) {
    coef_group <- sx_coef + sx_ranef_age[a, ] + sx_ranef_gender[g, ]
    effect <- as.numeric(B %*% coef_group)
    
    lag_list[[paste(a, g, sep = "_")]] <- data.frame(
      lag = 2:L,
      effect = effect,
      age = a,
      gender = g
    )
  }
}

lag_df <- do.call(rbind, lag_list)

lag_df$age <- factor(
  lag_df$age,
  levels = unique(lag_df$age),
  labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "<20")
)

lag_df$age <- factor(
  lag_df$age,
  levels = c("<20", "20-29", "30-39", "40-49", "50-59", "60-69")
)

ggplot(lag_df, aes(x = lag, y = effect, color = age)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(size = 1.2) +
  facet_wrap(~factor(gender, levels = c("Male", "Female")), scales = "free_y") +
  geom_point(size = 2) +
  labs(
    x = "Time lag (months)",
    y = "Estimated lag effect of X on Y",
    color = "Age (years)"
  ) +
  scale_x_continuous(
    breaks = unique(lag_df$lag),
    labels = 1:6
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_blank(),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    strip.text = element_text(size = 30, face = "bold")
  )
# ggsave("results/figures/Figure_XonY.jpg",width=10,height=8,dpi=300)

### Lag Effect of Z on Y: Fixed Effects
sz_coef <- c(
  posterior["Yobs_SZ_1", "Estimate"],
  posterior["Yobs_SZ_2", "Estimate"],
  posterior["Yobs_SZ_3", "Estimate"],
  posterior["Yobs_SZ_4", "Estimate"]
)

lag_index <- 2:L
B <- bs(lag_index, df = K, intercept = TRUE)

effect_Z <- as.numeric(B %*% sz_coef)

lag_df_Z <- data.frame(
  lag = lag_index,
  effect = effect_Z
)

ggplot(lag_df_Z, aes(x = lag, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(size = 1.2, color = "blue") +
  geom_point(size = 2, color = "blue") +
  labs(
    x = "Time lag (months)",
    y = "Estimated lag effect of Z on Y"
  ) +
  scale_x_continuous(
    breaks = unique(lag_df$lag),
    labels = 1:6
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_blank(),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24)
  )
# ggsave("results/figures/Figure_ZonY.jpg",width=10,height=8,dpi=300)

############################################################
## Lag Effect of Z on Y for Age ≥70
############################################################

# Spline-lag covariates for age ≥70
SZ_70 <- make_spline_lag(excess_depression_s, L = 7, K = 4)

# Excess suicides among age groups <70
Yij <- cbind(
  sui_male_obs2 - sui_male_pre2,
  sui_female_obs2 - sui_female_pre2
)[,-c(7,14)]

# Data for age ≥70
data2 <- data.frame(
  t = c(1:18, 1:18),
  Yobs = c(sui_male_obs2$aa70, sui_female_obs2$aa70),
  Ypre = c(sui_male_pre2$aa70, sui_female_pre2$aa70),
  extra = c(
    rowSums(sweep(Yij, 2, C[7,-c(7,14)], "*")),
    rowSums(sweep(Yij, 2, C[14,-c(7,14)], "*"))
  ),
  gender = c(rep("Male",18), rep("Female",18)),
  SZ_1 = rep(SZ_70[,1], 2),
  SZ_2 = rep(SZ_70[,2], 2),
  SZ_3 = rep(SZ_70[,3], 2),
  SZ_4 = rep(SZ_70[,4], 2)
)

# Model for age ≥70
fit_joint2 <- brm(
  bf(
    Yobs ~ Ypre +
      extra + (extra | gender) +
      SZ_1 + SZ_2 + SZ_3 + SZ_4,
    family = gaussian()
  ),
  data = data2,
  chains = 4, cores = 4, threads = threading(1),
  seed = 123, iter = 2000,
  control = list(adapt_delta = 0.9999, max_treedepth = 20),
  backend = "cmdstanr"
)

# Lag effect of Z on Y for age ≥70
posterior_70 <- fixef(fit_joint2)

sz_coef_70 <- c(
  posterior_70["SZ_1","Estimate"],
  posterior_70["SZ_2","Estimate"],
  posterior_70["SZ_3","Estimate"],
  posterior_70["SZ_4","Estimate"]
)

lag_index <- 2:L
B <- bs(lag_index, df = K, intercept = TRUE)

effect_Z_70 <- as.numeric(B %*% sz_coef_70)

lag_df_Z_70 <- data.frame(
  lag = lag_index,
  effect = effect_Z_70
)

max_lag_Z_70 <- lag_df_Z_70 %>%
  slice_max(
    order_by = effect,
    n = 1,
    with_ties = FALSE
  )

ggplot(lag_df_Z_70, aes(x = lag, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(size = 1.2, color = "blue") +
  geom_point(size = 2, color = "blue") +
  labs(
    x = "Time lag (months)",
    y = "Estimated lag effect of Z on Y"
  ) +
  scale_x_continuous(
    breaks = unique(lag_df_Z_70$lag),
    labels = 1:6
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_blank(),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24)
  )

# ggsave("results/figures/Figure_ZonY_70.jpg",width=10,height=8,dpi=300)

############################################################
## 2. unemployment vs control effort
############################################################
pt_df <- cbind(
  pt_monthly_means[1:6],
  pt_monthly_means[8:13]
)

unemp_df <- as.data.frame(
  cbind(
    exun_male$a10, exun_male$a20, exun_male$a30,
    exun_male$a40, exun_male$a50, exun_male$a60,
    exun_female$a10, exun_female$a20, exun_female$a30,
    exun_female$a40, exun_female$a50, exun_female$a60
  )
)

colnames(pt_df) <- colnames(unemp_df) <- c(
  "Male <20", "Male 20-29", "Male 30-39",
  "Male 40-49", "Male 50-59", "Male 60-69",
  "Female <20", "Female 20-29", "Female 30-39",
  "Female 40-49", "Female 50-59", "Female 60-69"
)

ccf_min <- c()
ccf_max <- c()

for (i in seq_len(ncol(pt_df))) {
  ok <- stats::complete.cases(pt_df[[i]], unemp_df[[i]])
  x1 <- as.numeric(pt_df[[i]][ok])
  y1 <- as.numeric(unemp_df[[i]][ok])
  
  ccf_res <- stats::ccf(
    x1, y1,
    lag.max = NULL,
    type = "correlation",
    plot = FALSE,
    na.action = na.omit
  )
  
  ccf_min <- c(ccf_min, min(ccf_res$acf))
  ccf_max <- c(ccf_max, max(ccf_res$acf))
}

y_limits <- c(min(ccf_min) - 0.05, max(ccf_max) + 0.05)

plot_ccf_pair <- function(x, y, var_name, lag_max = NULL, type = "correlation") {
  ok <- stats::complete.cases(x, y)
  x1 <- as.numeric(x[ok])
  y1 <- as.numeric(y[ok])
  n_eff <- length(x1)
  
  ccf_res <- stats::ccf(
    x1, y1,
    lag.max = lag_max,
    type = type,
    plot = FALSE,
    na.action = na.omit
  )
  
  df <- data.frame(
    lag = as.numeric(ccf_res$lag),
    ccf = as.numeric(ccf_res$acf)
  )
  
  ci <- 2 / sqrt(n_eff)
  
  ggplot(df, aes(x = lag, y = ccf)) +
    geom_col(width = 0.9) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(ci, -ci), linetype = "dashed") +
    labs(title = var_name) +
    scale_y_continuous(limits = y_limits) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 20),
      axis.title = element_blank()
    )
}

stopifnot(ncol(pt_df) == ncol(unemp_df))

var_names <- colnames(pt_df) %||% paste0("V", seq_len(ncol(pt_df)))

plot_list <- vector("list", ncol(pt_df))

for (i in seq_len(ncol(pt_df))) {
  plot_list[[i]] <- plot_ccf_pair(
    x = pt_df[[i]],
    y = unemp_df[[i]],
    var_name = var_names[i],
    lag_max = NULL,
    type = "correlation"
  )
}

big_plot <- wrap_plots(plot_list, nrow = 2, ncol = 6) &
  theme(plot.margin = margin(t = 10, r = 10, b = 15, l = 15))

final_plot <- ggdraw(big_plot) +
  draw_label("Time lag (months)", x = 0.5, y = -0.02, vjust = -1, size = 28) +
  draw_label("Cross-correlation", x = 0, y = 0.5, angle = 90, vjust = 1.5, size = 28)
# ggsave("results/figures/Figure_PonX.jpg",final_plot,width=24,height=8,dpi=300)

############################################################
## 3. depression vs growth rate
############################################################
ccf_result<-ccf(rate_of_change_s,excess_depression_s,lag.max=NULL,type=c("correlation","covariance"),plot=TRUE) 
df_ccf <- data.frame(
  Lag = ccf_result$lag,
  Correlation = ccf_result$acf
)
N <- min(length(rate_of_change_s), length(excess_depression_s))
ci <- 1.96 / sqrt(N)
max_idx <- which.max(abs(df_ccf$Correlation))
max_lag <- df_ccf$Lag[max_idx]
max_corr <- df_ccf$Correlation[max_idx]
ggplot(df_ccf, aes(x = Lag, y = Correlation)) +
  geom_col(width = 0.9) +       
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(ci, -ci), linetype = "dashed") +
  theme_classic() +                                
  theme(
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24)
  ) +
  labs(
    x = "Time lag (months)",
    y = "Cross-correlation"
  )
# ggsave("results/figures/Figure_VonZ.jpg",width=10,height=8,dpi=300)

############################################################
# End of script
############################################################