############################################################
## Random Effects Model Fitting and Post-processing
## Suicide and Socioeconomic Covariate Analysis
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
library(brms)
library(cmdstanr)
cmdstanr::check_cmdstan_toolchain()
library(ggplot2)
library(scales)
library(cowplot)
library(patchwork)

############################################################
## 1. Model construction
############################################################

###### age<70
pt <- pt_monthly_means
sui_male_obs1    <- sui_male_obs[7:24, ]
sui_male_pre1    <- sui_male_pre[7:24, ]
sui_female_obs1  <- sui_female_obs[7:24, ]
sui_female_pre1  <- sui_female_pre[7:24, ]

# Construct observed / predicted suicides (male / female)
exsui_male_new <- data.frame(
  t = rep(1:18, 6),
  Yobs = c(sui_male_obs1$au20, sui_male_obs1$a20, sui_male_obs1$a30,
           sui_male_obs1$a40,  sui_male_obs1$a50, sui_male_obs1$a60),
  Ypre = c(sui_male_pre1$au20, sui_male_pre1$a20, sui_male_pre1$a30,
           sui_male_pre1$a40,  sui_male_pre1$a50, sui_male_pre1$a60),
  age = c(rep("u20", 18), rep("20", 18), rep("30", 18),
          rep("40", 18), rep("50", 18), rep("60", 18)),
  gender = rep("Male", 18 * 6)
)

exsui_female_new <- data.frame(
  t = rep(1:18, 6),
  Yobs = c(sui_female_obs1$au20, sui_female_obs1$a20, sui_female_obs1$a30,
           sui_female_obs1$a40,  sui_female_obs1$a50, sui_female_obs1$a60),
  Ypre = c(sui_female_pre1$au20, sui_female_pre1$a20, sui_female_pre1$a30,
           sui_female_pre1$a40,  sui_female_pre1$a50, sui_female_pre1$a60),
  age = c(rep("u20", 18), rep("20", 18), rep("30", 18),
          rep("40", 18), rep("50", 18), rep("60", 18)),
  gender = rep("Female", 18 * 6)
)

# Construct unemployment covariate (male / female)
exun_male_new <- data.frame(
  t = rep(1:18, 6),
  X = c(exun_male$a10[6:23], exun_male$a20[6:23], exun_male$a30[6:23],
        exun_male$a40[6:23], exun_male$a50[6:23], exun_male$a60[6:23]),
  age = c(rep("u20", 18), rep("20", 18), rep("30", 18),
          rep("40", 18), rep("50", 18), rep("60", 18)),
  gender = rep("Male", 18 * 6)
)

exun_female_new <- data.frame(
  t = rep(1:18, 6),
  X = c(exun_female$a10[4:21], exun_female$a20[4:21], exun_female$a30[4:21],
        exun_female$a40[4:21], exun_female$a50[4:21], exun_female$a60[4:21]),
  age = c(rep("u20", 18), rep("20", 18), rep("30", 18),
          rep("40", 18), rep("50", 18), rep("60", 18)),
  gender = rep("Female", 18 * 6)
)

# standardization
excess_depression_s<-as.vector(scale(excess_depression, center = TRUE, scale = TRUE))
rate_of_change_s<-as.vector(scale(rate_of_change, center = TRUE, scale = TRUE))

# Main dataset
data1 <- data.frame(
  t = c(exun_male_new$t, exun_female_new$t),
  X = c(exun_male_new$X, exun_female_new$X),
  Yobs = c(exsui_male_new$Yobs, exsui_female_new$Yobs),
  Ypre = c(exsui_male_new$Ypre, exsui_female_new$Ypre),
  Z = rep(excess_depression_s[1:18], 12),
  age = c(exun_male_new$age, exun_female_new$age),
  gender = c(exun_male_new$gender, exun_female_new$gender)
)

# Expanded dataset
data1_new <- data.frame(
  t = c(exun_male_new$t, exun_female_new$t, rep(NA, 6 * 12)),
  Yobs = c(exsui_male_new$Yobs, exsui_female_new$Yobs, rep(NA, 6 * 12)),
  Ypre = c(exsui_male_new$Ypre, exsui_female_new$Ypre, rep(NA, 6 * 12)),
  X = c(exun_male_new$X, exun_female_new$X, rep(NA, 6 * 12)),
  Z = c(rep(excess_depression_s[1:18], 12), rep(NA, 6 * 12)),
  age = c(exun_male_new$age, exun_female_new$age, rep(NA, 6 * 12)),
  gender = c(exun_male_new$gender, exun_female_new$gender, rep(NA, 6 * 12)),
  Zall = rep(excess_depression_s[1:24], 12),
  Xall = c(exun_male$a10, exun_male$a20, exun_male$a30, exun_male$a40, exun_male$a50, exun_male$a60,
           exun_female$a10, exun_female$a20, exun_female$a30, exun_female$a40, exun_female$a50, exun_female$a60),
  P = c(pt_monthly_means$V1, pt_monthly_means$V2, pt_monthly_means$V3,
        pt_monthly_means$V4, pt_monthly_means$V5, pt_monthly_means$V6,
        pt_monthly_means$V8, pt_monthly_means$V9, pt_monthly_means$V10,
        pt_monthly_means$V11, pt_monthly_means$V12, pt_monthly_means$V13),
  V = rep(rate_of_change_s[1:24], 12),
  genderall = c(rep("Male", 24 * 6), rep("Female", 24 * 6)),
  ageall = c(rep("u20", 24), rep("20", 24), rep("30", 24),
             rep("40", 24), rep("50", 24), rep("60", 24),
             rep("u20", 24), rep("20", 24), rep("30", 24),
             rep("40", 24), rep("50", 24), rep("60", 24))
)

# model
fit_joint1 <- brm(
  bf(Xall ~ P + (P | ageall) + (P | genderall),family = gaussian()) +
  bf(Zall ~ V,family = gaussian()) +
  bf(Yobs ~ Ypre + X + Z + (X | age) + (X | gender),family = gaussian()) +
  set_rescor(FALSE),
  data = data1_new,
  chains = 4, cores = 4, threads = threading(1),
  seed = 123, iter = 2000,
  control = list(adapt_delta = 0.9999, max_treedepth = 20),
  backend = "cmdstanr"
)

y_pred1 <- posterior_predict(fit_joint1, newdata = data1, resp = "Yobs")
pred_mean1 <- apply(y_pred1, 2, mean)
pred_low1  <- apply(y_pred1, 2, quantile, probs = 0.025)
pred_high1 <- apply(y_pred1, 2, quantile, probs = 0.975)

fixef(fit_joint1)["Yobs_Z", ]  # Z → Y
fixef(fit_joint1)["Zall_V", ]  # V → Z

###### age>70
Yij<-cbind(sui_male_obs1-sui_male_pre1,sui_female_obs1-sui_female_pre1)[,-c(7,14)]
data2<-data.frame(t=c(1:18,1:18),
                  Yobs=c(sui_male_obs1$aa70,sui_female_obs1$aa70),
                  Ypre=c(sui_male_pre1$aa70,sui_female_pre1$aa70),
                  extra=c(rowSums(sweep(Yij, 2, C[7,-c(7,14)], "*")),rowSums(sweep(Yij, 2, C[14,-c(7,14)], "*"))),
                  Z=rep(excess_depression_s[6:23],2),
                  gender=c(rep("Male",18),rep("Female",18)))

# model
fit_joint2 <- brm(
  bf(Yobs ~ Ypre + extra + (extra | gender) + Z,family = gaussian()),
  data = data2,
  chains = 4, cores = 4, threads = threading(1),
  seed = 123, iter = 2000,
  control = list(adapt_delta = 0.9999, max_treedepth = 20),
  backend = "cmdstanr"
)

y_pred2 <- posterior_predict(fit_joint2, newdata = data2, resp = "Yobs")
pred_mean2 <- apply(y_pred2, 2, mean)
pred_low2  <- apply(y_pred2, 2, quantile, probs = 0.025)
pred_high2 <- apply(y_pred2, 2, quantile, probs = 0.975)

############################################################
## 2. Posterior Prediction with Group-Specific Effects
############################################################
X_post <- posterior_epred(
  fit_joint1,
  resp = "Xall",
  newdata = data1_new,
  re_formula = ~(P | ageall) + (P | genderall)
)

X_df <- data1_new %>%
  mutate(
    X_mean = apply(X_post, 2, mean),
    X_low  = apply(X_post, 2, quantile, probs = 0.025),
    X_high = apply(X_post, 2, quantile, probs = 0.975)
  )

X_df$ageall <- factor(
  X_df$ageall,
  levels = c("u20", "20", "30", "40", "50", "60"),
  labels = c(
    "Age < 20 years", "Age 20–29 years", "Age 30–39 years",
    "Age 40–49 years", "Age 50–59 years", "Age 60–69 years"
  )
)

ymin_X <- min(X_df$X_low, na.rm = TRUE)
ymax_X <- max(X_df$X_high, na.rm = TRUE)

age_levels_X <- levels(X_df$ageall)
plot_list_X <- map(age_levels_X, function(age_group) {
  df_sub <- X_df %>% filter(ageall == age_group)
  
  ggplot(df_sub, aes(x = P, y = X_mean, color = genderall, fill = genderall)) +
    geom_line(aes(group = interaction(ageall, genderall)), size = 1) +
    geom_ribbon(aes(ymin = X_low, ymax = X_high), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("Male" = "#1F77B4", "Female" = "#FF7F0E")) +
    scale_fill_manual(values = c("Male" = "#1F77B4", "Female" = "#FF7F0E")) +
    scale_x_continuous(
      breaks = c(0, 0.2, 0.4, 0.6, 0.8)
    ) +
    expand_limits(x = 0.8) +
    coord_cartesian(ylim = c(ymin_X, ymax_X)) +
    ggtitle(age_group) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 20),
      strip.text = element_text(size = 20),
      axis.ticks.length = unit(0.25, "cm"),
      plot.title = element_text(size = 24, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      legend.position = "none"
    )
})

Y_df_input <- data1_new %>% filter(!is.na(X) & !is.na(Yobs) & !is.na(age) & !is.na(gender))
X_seq <- seq(min(Y_df_input$X), max(Y_df_input$X), length.out = 100)

age_levels_Y <- unique(Y_df_input$age)
gender_levels_Y <- unique(Y_df_input$gender)

newdata_pred <- expand.grid(
  X = X_seq,
  age = age_levels_Y,
  gender = gender_levels_Y,
  Ypre = mean(Y_df_input$Ypre, na.rm = TRUE),
  Z = mean(Y_df_input$Z, na.rm = TRUE)
)

Y_post_smooth <- posterior_epred(
  fit_joint1,
  resp = "Yobs",
  newdata = newdata_pred,
  re_formula = ~(X | age) + (X | gender)
)

Y_df_smooth <- newdata_pred %>%
  mutate(
    Y_mean = apply(Y_post_smooth, 2, mean),
    Y_low  = apply(Y_post_smooth, 2, quantile, probs = 0.025),
    Y_high = apply(Y_post_smooth, 2, quantile, probs = 0.975)
  )

Y_df_smooth$age <- factor(
  Y_df_smooth$age,
  levels = c("u20", "20", "30", "40", "50", "60"),
  labels = c(
    "Age < 20 years", "Age 20–29 years", "Age 30–39 years",
    "Age 40–49 years", "Age 50–59 years", "Age 60–69 years"
  )
)

ymin_Y <- min(Y_df_smooth$Y_low - Y_df_smooth$Ypre, na.rm = TRUE)
ymax_Y <- max(Y_df_smooth$Y_high - Y_df_smooth$Ypre, na.rm = TRUE)

age_levels_Y <- levels(Y_df_smooth$age)
plot_list_Y <- map(age_levels_Y, function(age_group) {
  df_sub <- Y_df_smooth %>% filter(age == age_group)
  
  ggplot(df_sub, aes(x = X, y = Y_mean - Ypre, color = gender, fill = gender)) +
    geom_line(aes(group = interaction(age, gender)), size = 1) +
    geom_ribbon(aes(ymin = Y_low - Ypre, ymax = Y_high - Ypre), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("Male" = "#1F77B4", "Female" = "#FF7F0E")) +
    scale_fill_manual(values = c("Male" = "#1F77B4", "Female" = "#FF7F0E")) +
    scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8, 10, 12)) +
    coord_cartesian(ylim = c(ymin_Y, ymax_Y)) +
    ggtitle(age_group) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 20),
      strip.text = element_text(size = 20),
      axis.ticks.length = unit(0.25, "cm"),
      plot.title = element_text(size = 24, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      legend.position = "none"
    )
})

group_names <- c(
  "Male <20", "Male 20-29", "Male 30-39", "Male 40-49", "Male 50-59", "Male 60-69", "Male ≥70",
  "Female <20", "Female 20-29", "Female 30-39", "Female 40-49", "Female 50-59", "Female 60-69", "Female ≥70"
)

pred_mean <- c(pred_mean1[1:108], pred_mean2[1:18], pred_mean1[109:216], pred_mean2[19:36])
pred_low  <- c(pred_low1[1:108],  pred_low2[1:18],  pred_low1[109:216],  pred_low2[19:36])
pred_high <- c(pred_high1[1:108], pred_high2[1:18], pred_high1[109:216], pred_high2[19:36])

plot_data <- data.frame(
  obs_id = 1:length(pred_mean),
  group  = rep(1:14, each = 18),
  pred_mean = pred_mean,
  pred_low  = pred_low,
  pred_high = pred_high,
  actual = c(data1$Yobs[1:108], data2$Yobs[1:18], data1$Yobs[109:216], data2$Yobs[19:36])
)

df_mean <- plot_data %>%
  dplyr::select(group, pred_mean) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(row = dplyr::row_number()) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = pred_mean
  ) %>%
  dplyr::select(-row)

df_low <- plot_data %>%
  dplyr::select(group, pred_low) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(row = dplyr::row_number()) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = pred_low
  ) %>%
  dplyr::select(-row)

df_upper <- plot_data %>%
  dplyr::select(group, pred_high) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(row = dplyr::row_number()) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = pred_high
  ) %>%
  dplyr::select(-row)

colnames(df_mean) <- colnames(df_low) <- colnames(df_upper) <-
  colnames(suibase) <- colnames(suibase_upper) <- colnames(suibase_lower) <-
  colnames(suiobs) <- group_names

base  <- rbind(suibase, df_mean)
upper <- rbind(suibase_upper, df_upper)
lower <- rbind(suibase_lower, df_low)

n_per_group <- 84
plot_data2 <- do.call(rbind, lapply(1:14, function(g) {
  data.frame(
    group = g,
    group_name = group_names[g],
    obs_id = 1:n_per_group,
    mean = base[, g],
    low = lower[, g],
    upper = upper[, g],
    obs = suiobs[, g]
  )
}))

ymin2 <- min(plot_data2$low, plot_data2$obs, na.rm = TRUE)
ymax2 <- max(plot_data2$upper, plot_data2$obs, na.rm = TRUE)

create_unified_theme <- function() {
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_blank(),
      plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
      plot.margin = margin(t = 10, r = 10, b = 15, l = 15),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank()
    )
}

plot_list2 <- lapply(1:14, function(g) {
  df <- subset(plot_data2, group == g)
  
  ggplot(df, aes(x = obs_id)) +
    geom_ribbon(aes(ymin = low, ymax = upper), alpha = 0.1,
                fill = "pink", color = "red", linetype = "dotted") +
    geom_line(aes(y = mean), color = "red", size = 1) +
    geom_vline(xintercept = 66, linetype = "dashed") +
    geom_point(aes(y = obs), shape = 16, fill = "black", size = 1) +
    scale_x_continuous(breaks = c(1, 12, 24, 36, 48, 60, 72, 84)) +
    scale_y_continuous(limits = c(ymin2, ymax2)) +
    labs(title = group_names[g]) +
    create_unified_theme()
})

big_plot2 <- wrap_plots(plot_list2, nrow = 2, ncol = 7) &
  theme(plot.margin = margin(t = 10, r = 10, b = 15, l = 15))

legend_plot <- get_legend(
  plot_list_X[[1]] +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 28),
      legend.text = element_text(size = 24)
    )
)

legend_grob <- as_grob(legend_plot)

legend_centered <- plot_grid(
  NULL, legend_grob, NULL,
  ncol = 1,
  rel_heights = c(0.25, 0.5, 0.25)
)

plot_list_block <- plot_grid(
  wrap_plots(plot_list_X, nrow = 1, ncol = 6),
  legend_centered,
  ncol = 2,
  rel_widths = c(6, 1)
)

plot_list_Y_block <- plot_grid(
  wrap_plots(plot_list_Y, nrow = 1, ncol = 6),
  plot_grid(NULL, NULL, ncol = 1),
  ncol = 2,
  rel_widths = c(6, 1)
)

big_block_with_xlab <- plot_grid(big_plot2, NULL, ncol = 1, rel_heights = c(1, 0.05))
plot_list_block_with_xlab <- plot_grid(plot_list_block, NULL, ncol = 1, rel_heights = c(1, 0.05))
plot_list_Y_block_with_xlab <- plot_grid(plot_list_Y_block, NULL, ncol = 1, rel_heights = c(1, 0.05))

final_big <- plot_grid(
  big_block_with_xlab,
  plot_list_block_with_xlab,
  plot_list_Y_block_with_xlab,
  ncol = 1,
  rel_heights = c(1, 0.5, 0.5)
)

final_big <- ggdraw(final_big, xlim = c(-0.025, 1), ylim = c(-0.01, 1)) +
  draw_label("Months from January 2015 to December 2021", x = 0.5, y = 0.53, size = 28) +
  draw_label("Suicides", x = -0.005, y = 0.75, angle = 90, size = 28) +
  draw_label("Required effort for control", x = 0.5, y = 0.26, size = 28) +
  draw_label("Excess unemployed \npersons (×10,000)", x = -0.01, y = 0.38, angle = 90, size = 28) +
  draw_label("Excess unemployed persons (×10,000)", x = 0.5, y = 0.01, size = 28) +
  draw_label("Excess suicides", x = -0.005, y = 0.13, angle = 90, size = 28)
# ggsave("results/figures/Figure_PosteriorPrediction.jpg",final_big,width=28,height=16,dpi=300)

############################################################
## 3. Senario Loop
############################################################
seeds <- 1:5
n_seed <- length(seeds)

estsui_PV_upper_all    <- array(NA_real_, dim = c(n_seed, 61))
estsui_PV_lower_all    <- array(NA_real_, dim = c(n_seed, 61))
estsui_PV_df_upper_all <- array(NA_real_, dim = c(n_seed, 61, 14))
estsui_PV_df_lower_all <- array(NA_real_, dim = c(n_seed, 61, 14))

for (s in seq_along(seeds)) {

  set.seed(seeds[s])
  
  # upper
  estsui_PV_upper    <- rep(0, 61)
  estsui_PV_df_upper <- matrix(0, nrow = 61, ncol = 14)

  for (i in 1:61) {
    pt <- pt_monthly_CI[[i]]$upper
    pt <- rbind(rep(0, 14), pt)
    pt <- as.data.frame(pt)

    rate_of_change_s <- as.vector(scale(rate_of_change_m_lower[, i], center = TRUE, scale = TRUE))

    data <- data1_new
    data$P <- c(pt$V1, pt$V2, pt$V3, pt$V4, pt$V5, pt$V6,
                pt$V8, pt$V9, pt$V10, pt$V11, pt$V12, pt$V13)
    data$V <- rep(rate_of_change_s[1:24], 12)

    x_pred <- posterior_linpred(
      fit_joint1, newdata = data, resp = "Xall",
      transform = TRUE, allow_new_levels = TRUE
    )
    pred_x <- apply(x_pred, 2, mean)

    z_pred <- posterior_linpred(
      fit_joint1, newdata = data, resp = "Zall",
      transform = TRUE, allow_new_levels = TRUE
    )
    pred_z <- apply(z_pred, 2, mean)

    mat_X <- matrix(pred_x, nrow = 24, ncol = 12)
    mat_Z <- matrix(pred_z, nrow = 24, ncol = 12)

    data <- data1
    data$X <- c(mat_X[6:23, 1], mat_X[6:23, 2], mat_X[6:23, 3], mat_X[6:23, 4],
                mat_X[6:23, 5], mat_X[6:23, 6], mat_X[4:21, 7], mat_X[4:21, 8],
                mat_X[4:21, 9], mat_X[4:21, 10], mat_X[4:21, 11], mat_X[4:21, 12])
    data$Z <- c(mat_Z[1:18, 1], mat_Z[1:18, 2], mat_Z[1:18, 3], mat_Z[1:18, 4],
                mat_Z[1:18, 5], mat_Z[1:18, 6], mat_Z[1:18, 7], mat_Z[1:18, 8],
                mat_Z[1:18, 9], mat_Z[1:18, 10], mat_Z[1:18, 11], mat_Z[1:18, 12])

    y_pred1 <- posterior_predict(fit_joint1, newdata = data, resp = "Yobs")
    pred_mean1 <- apply(y_pred1, 2, mean)
    pred_low1  <- apply(y_pred1, 2, quantile, probs = 0.025)
    pred_high1 <- apply(y_pred1, 2, quantile, probs = 0.975)

    sui_obs <- matrix(pred_mean1, nrow = 18, byrow = FALSE)
    sui_pre <- cbind(sui_male_pre1[, -7], sui_female_pre1[-7])
    Yij <- sui_obs - sui_pre

    data <- data.frame(
      t      = data2$t,
      Yobs   = data2$Yobs,
      Ypre   = data2$Ypre,
      extra  = c(
        rowSums(sweep(Yij, 2, C[7,  -c(7, 14)],  "*")),
        rowSums(sweep(Yij, 2, C[14, -c(7, 14)], "*"))
      ),
      Z      = rep(mat_Z[6:23, 1], 2),
      gender = c(rep("Male", 18), rep("Female", 18))
    )

    y_pred2 <- posterior_predict(fit_joint2, newdata = data, resp = "Yobs")
    pred_mean2 <- apply(y_pred2, 2, mean)
    pred_low2  <- apply(y_pred2, 2, quantile, probs = 0.025)
    pred_high2 <- apply(y_pred2, 2, quantile, probs = 0.975)

    pred_mean_all <- c(
      pred_high1[1:108], pred_high2[1:18],
      pred_high1[109:216], pred_high2[19:36]
    ) - c(
      data1$Ypre[1:108], data2$Ypre[1:18],
      data1$Ypre[109:216], data2$Ypre[19:36]
    )

    estsui_PV_df_upper[i, ] <- tapply(pred_mean_all, ceiling(seq_along(pred_mean_all) / 18), sum)
    estsui_PV_upper[i]      <- sum(pred_mean_all)
  }
  
  # lower
  estsui_PV_lower    <- rep(0, 61)
  estsui_PV_df_lower <- matrix(0, nrow = 61, ncol = 14)

  for (i in 1:61) {
    pt <- pt_monthly_CI[[i]]$lower
    pt <- rbind(rep(0, 14), pt)
    pt <- as.data.frame(pt)

    rate_of_change_s <- as.vector(scale(rate_of_change_m_upper[, i], center = TRUE, scale = TRUE))
    
    data <- data1_new
    data$P <- c(pt$V1, pt$V2, pt$V3, pt$V4, pt$V5, pt$V6,
                pt$V8, pt$V9, pt$V10, pt$V11, pt$V12, pt$V13)
    data$V <- rep(rate_of_change_s[1:24], 12)

    x_pred <- posterior_linpred(
      fit_joint1, newdata = data, resp = "Xall",
      transform = TRUE, allow_new_levels = TRUE
    )
    pred_x <- apply(x_pred, 2, mean)

    z_pred <- posterior_linpred(
      fit_joint1, newdata = data, resp = "Zall",
      transform = TRUE, allow_new_levels = TRUE
    )
    pred_z <- apply(z_pred, 2, mean)

    mat_X <- matrix(pred_x, nrow = 24, ncol = 12)
    mat_Z <- matrix(pred_z, nrow = 24, ncol = 12)

    data <- data1
    data$X <- c(mat_X[6:23, 1], mat_X[6:23, 2], mat_X[6:23, 3], mat_X[6:23, 4],
                mat_X[6:23, 5], mat_X[6:23, 6], mat_X[4:21, 7], mat_X[4:21, 8],
                mat_X[4:21, 9], mat_X[4:21, 10], mat_X[4:21, 11], mat_X[4:21, 12])
    data$Z <- c(mat_Z[1:18, 1], mat_Z[1:18, 2], mat_Z[1:18, 3], mat_Z[1:18, 4],
                mat_Z[1:18, 5], mat_Z[1:18, 6], mat_Z[1:18, 7], mat_Z[1:18, 8],
                mat_Z[1:18, 9], mat_Z[1:18, 10], mat_Z[1:18, 11], mat_Z[1:18, 12])

    y_pred1 <- posterior_predict(fit_joint1, newdata = data, resp = "Yobs")
    pred_mean1 <- apply(y_pred1, 2, mean)
    pred_low1  <- apply(y_pred1, 2, quantile, probs = 0.025)
    pred_high1 <- apply(y_pred1, 2, quantile, probs = 0.975)

    sui_obs <- matrix(pred_mean1, nrow = 18, byrow = FALSE)
    sui_pre <- cbind(sui_male_pre1[, -7], sui_female_pre1[-7])
    Yij <- sui_obs - sui_pre

    data <- data.frame(
      t      = data2$t,
      Yobs   = data2$Yobs,
      Ypre   = data2$Ypre,
      extra  = c(
        rowSums(sweep(Yij, 2, C[7,  -c(7, 14)],  "*")),
        rowSums(sweep(Yij, 2, C[14, -c(7, 14)], "*"))
      ),
      Z      = rep(mat_Z[6:23, 1], 2),
      gender = c(rep("Male", 18), rep("Female", 18))
    )

    y_pred2 <- posterior_predict(fit_joint2, newdata = data, resp = "Yobs")
    pred_mean2 <- apply(y_pred2, 2, mean)
    pred_low2  <- apply(y_pred2, 2, quantile, probs = 0.025)
    pred_high2 <- apply(y_pred2, 2, quantile, probs = 0.975)

    pred_mean_all <- c(
      pred_low1[1:108], pred_low2[1:18],
      pred_low1[109:216], pred_low2[19:36]
    ) - c(
      data1$Ypre[1:108], data2$Ypre[1:18],
      data1$Ypre[109:216], data2$Ypre[19:36]
    )

    estsui_PV_df_lower[i, ] <- tapply(pred_mean_all, ceiling(seq_along(pred_mean_all) / 18), sum)
    estsui_PV_lower[i]      <- sum(pred_mean_all)
  }

  estsui_PV_upper_all[s, ] <- estsui_PV_upper
  estsui_PV_lower_all[s, ] <- estsui_PV_lower
  estsui_PV_df_upper_all[s, , ] <- estsui_PV_df_upper
  estsui_PV_df_lower_all[s, , ] <- estsui_PV_df_lower
}

### seed-averaged results
estsui_PV_upper <- apply(estsui_PV_upper_all, 2, mean, na.rm = TRUE)
estsui_PV_lower <- apply(estsui_PV_lower_all, 2, mean, na.rm = TRUE)

estsui_PV_df_upper <- apply(estsui_PV_df_upper_all, c(2, 3), mean, na.rm = TRUE)
estsui_PV_df_lower <- apply(estsui_PV_df_lower_all, c(2, 3), mean, na.rm = TRUE)

estsui_PV_mean    <- (estsui_PV_lower + estsui_PV_upper) / 2
estsui_PV_df_mean <- (estsui_PV_df_lower + estsui_PV_df_upper) / 2
#
estsui_PV_lower[estsui_PV_lower < 0] <- 0

### plot
df <- data.frame(
  l_value = 0:60,
  est = rev(estsui_PV_mean),
  lower = rev(estsui_PV_lower),
  upper = rev(estsui_PV_upper)
)
ggplot(df, aes(x = l_value, y = est)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray30", alpha = 0.3) +  
  geom_line(color = "black", size = 1) +  
  geom_line(aes(y = lower), color = "black", linetype = "dashed", size = 0.8) + 
  geom_line(aes(y = upper), color = "black", linetype = "dashed", size = 0.8) + 
  coord_cartesian(ylim = c(0, max(df$upper))) +
  scale_y_continuous(labels = label_comma()) +
  scale_x_continuous(
    breaks = c(0, 10, 20, 30, 40, 50, 60)
  ) +
  labs(
    x = "Additional relative control effort (ρ)",
    y = "Estimated excess suicides"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  )
# ggsave("results/figures/Figure_Excesssuicides.jpg",width=10,height=8,dpi=300)

############################################################
## 4. Total Excess Deaths (COVID-19 + Excess Suicides)
############################################################
###### Total death
all<-estsui_PV_mean+estdeathall
all_lower<-estsui_PV_upper+estdeathall_lower
all_upper<-estsui_PV_lower+estdeathall_upper
min<-min(all,all_lower,all_upper)
max<-max(all,all_lower,all_upper)
df <- data.frame(
  l_value = 0:60,
  est = rev(all),
  lower = rev(all_lower),
  upper = rev(all_upper)
)
ymin <- min-10
ymax <- max
ggplot(df, aes(x = l_value, y = est)) +
  geom_line(color = "black", size = 1) +  
  geom_line(aes(y = lower), color = "blue", linetype = "dashed", size = 0.8) + 
  geom_line(aes(y = upper), color = "red", linetype = "dashed", size = 0.8) + 
  annotate("point", x = 100-(which(all==min(all))+39), y = min(all), color = "black", size = 2.5) +
  annotate("point", x = 100-(which(all_lower==min(all_lower))+39), y = min(all_lower), color = "blue", size = 2.5) +
  annotate("point", x = 100-(which(all_upper>=1)[1]-1+39), y = min(all_upper), color = "red", size = 2.5) +
  scale_y_continuous(limits = c(ymin, ymax), labels = label_comma()) +
  scale_x_continuous(
    breaks = c(0, 10, 20, 30, 40, 50, 60)
  ) +
  labs(
    x = "Additional relative control effort (ρ)",
    y = "Estimated excess deaths"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24),
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.margin = margin(8, 5.5, 5.5, 5.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  )
# ggsave("results/figures/Figure_death.jpg",width=10,height=8,dpi=300)

###### YLL
YLL <- read.csv(
  "data/YLL.csv",
  header = TRUE,
  row.names = 1
)
YLL_male<-YLL$Male
YLL_female<-YLL$Female

YLL_covid_male<-c(mean(YLL_male[1:4]),mean(YLL_male[5:6]),mean(YLL_male[7:8]),mean(YLL_male[9:10]),mean(YLL_male[11:12]),mean(YLL_male[13:14]),mean(YLL_male[15:21]))
YLL_covid_female<-c(mean(YLL_female[1:4]),mean(YLL_female[5:6]),mean(YLL_female[7:8]),mean(YLL_female[9:10]),mean(YLL_female[11:12]),mean(YLL_female[13:14]),mean(YLL_female[15:21]))
YLL_covid<-c(YLL_covid_male,YLL_covid_female)

YLL_sui_male<-c(mean(YLL_male[3:4]),mean(YLL_male[5:6]),mean(YLL_male[7:8]),mean(YLL_male[9:10]),mean(YLL_male[11:12]),mean(YLL_male[13:14]),mean(YLL_male[15:21]))
YLL_sui_female<-c(mean(YLL_female[3:4]),mean(YLL_female[5:6]),mean(YLL_female[7:8]),mean(YLL_female[9:10]),mean(YLL_female[11:12]),mean(YLL_female[13:14]),mean(YLL_female[15:21]))
YLL_sui<-c(YLL_sui_male,YLL_sui_female)

### YLL Calculation
le_covid<-rep(0,61)
for(i in 1:61){
  le_covid[i]<-sum(estdeath_df[i,]*YLL_covid)
}
le_covid_lower<-rep(0,61)
for(i in 1:61){
  le_covid_lower[i]<-sum(estdeath_df_lower[i,]*YLL_covid)
}
le_covid_upper<-rep(0,61)
for(i in 1:61){
  le_covid_upper[i]<-sum(estdeath_df_upper[i,]*YLL_covid)
}
##
le_sui<-rep(0,61)
for(i in 1:61){
  le_sui[i]<-sum(estsui_PV_df_mean[i,]*YLL_sui)
}
le_sui_lower<-rep(0,61)
for(i in 1:61){
  le_sui_lower[i]<-sum(estsui_PV_df_lower[i,]*YLL_sui)
}
le_sui_upper<-rep(0,61)
for(i in 1:61){
  le_sui_upper[i]<-sum(estsui_PV_df_upper[i,]*YLL_sui)
}
#
le_sui_lower[le_sui_lower<0]<-0

### Total YLL
all<-le_covid+le_sui
all_lower<-le_covid_lower+le_sui_upper
all_upper<-le_covid_upper+le_sui_lower
min<-min(all,all_lower,all_upper)
max<-max(all,all_lower,all_upper)
df <- data.frame(
  l_value = 0:60,
  est = rev(all),
  lower = rev(all_lower),
  upper = rev(all_upper)
)
ymin <- min-1000
ymax <- max
ggplot(df, aes(x = l_value, y = est)) +
  geom_line(color = "black", size = 1) +  
  geom_line(aes(y = lower), color = "blue", linetype = "dashed", size = 0.8) + 
  geom_line(aes(y = upper), color = "red", linetype = "dashed", size = 0.8) + 
  annotate("point", x = 100-(which(all==min(all))+39), y = min(all), color = "black", size = 2.5) +
  annotate("point", x = 100-(which(all_lower==min(all_lower))+39), y = min(all_lower), color = "blue", size = 2.5) +
  annotate("point", x = 100-(which(all_upper>0)[1]-1+39), y = min(all_upper), color = "red", size = 2.5) +
  scale_y_continuous(limits = c(ymin, ymax), labels = label_comma()) +
  scale_x_continuous(
    breaks = c(0, 10, 20, 30, 40, 50, 60)
  ) +
  labs(
    x = "Additional relative control effort (ρ)",
    y = "Estimated YLL"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24),
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.margin = margin(8, 5.5, 5.5, 5.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  )
# ggsave("results/figures/Figure_YLL.jpg",width=10,height=8,dpi=300)

############################################################
# End of script
############################################################