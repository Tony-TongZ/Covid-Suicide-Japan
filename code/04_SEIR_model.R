############################################################
## SEIR Model Fitting and Post-processing
## COVID-19 Epidemic Analysis
############################################################
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
set.seed(123)

############################################################
## 0. Packages
############################################################
library(dplyr)
library(tidyr)
library(lubridate)
library(deSolve)
library(nloptr)
library(numDeriv)
library(Matrix)
library(MASS)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyverse)
library(scales)

############################################################
## 1. Total daily infections
############################################################
InfData <- read.csv(
  "data/newly_confirmed_cases_daily_MHLW.csv",
  header = TRUE
)

InfDatabefore <- InfData[
  which(InfData$Date == "2020/1/16") : which(InfData$Date == "2020/9/1"),
  1:2
]

InfDataafter <- InfData[
  which(InfData$Date == "2020/9/2") : which(InfData$Date == "2021/12/31"),
  1:2
]

############################################################
## 2. Daily age–sex distribution (early period)
############################################################
InfData1 <- read.csv(
  "data/case2020_detail_prefectures_reports.csv",
  header = TRUE
)

InfData1 <- InfData1 %>%
  filter(!is.na(AgeSex), !is.na(Sex)) %>%
  mutate(pref = 1) %>%
  arrange(Date) %>%
  group_by(Date, AgeSex, Sex) %>%
  summarise(case = sum(pref), .groups = "drop") %>%
  mutate(Date = as.Date(Date))

## Helper: aggregate age groups
aggregate_age <- function(df, ages) {
  df %>%
    filter(AgeSex %in% ages) %>%
    group_by(Date) %>%
    summarise(case = sum(case), .groups = "drop")
}

## Male
male <- InfData1 %>% filter(Sex == "M")
male_groups <- list(
  maleu20 = c(0, 10),
  male20  = 20,
  male30  = 30,
  male40  = 40,
  male50  = 50,
  male60  = 60,
  malea70 = c(70, 80, 90, 100)
)
male_df <- lapply(male_groups, aggregate_age, df = male)

## Female
female <- InfData1 %>% filter(Sex == "F")
female_groups <- list(
  femaleu20 = c(0, 10),
  female20  = 20,
  female30  = 30,
  female40  = 40,
  female50  = 50,
  female60  = 60,
  femalea70 = c(70, 80, 90, 100)
)
female_df <- lapply(female_groups, aggregate_age, df = female)

## Fill missing dates
all_dates <- data.frame(
  Date = seq(as.Date("2020-01-15"), as.Date("2020-09-01"), by = "day")
)

fill_zero <- function(df) {
  out <- merge(all_dates, df, by = "Date", all.x = TRUE)
  out$case[is.na(out$case)] <- 0
  out
}

male_df   <- lapply(male_df, fill_zero)
female_df <- lapply(female_df, fill_zero)

############################################################
## 3. Daily age–sex counts (early period)
############################################################
InfData1_mf <- data.frame(
  date = male_df$maleu20$Date,
  maleu20 = male_df$maleu20$case,
  male20  = male_df$male20$case,
  male30  = male_df$male30$case,
  male40  = male_df$male40$case,
  male50  = male_df$male50$case,
  male60  = male_df$male60$case,
  malea70 = male_df$malea70$case,
  femaleu20 = female_df$femaleu20$case,
  female20  = female_df$female20$case,
  female30  = female_df$female30$case,
  female40  = female_df$female40$case,
  female50  = female_df$female50$case,
  female60  = female_df$female60$case,
  femalea70 = female_df$femalea70$case
)

InfData1_prop <- apply(InfData1_mf[,-1], 1, function(x) {
  if (sum(x) == 0) rep(0, length(x)) else x / sum(x)
}) |> t() |> as.data.frame()

df1 <- InfData1_prop[-1, ] * InfDatabefore$ALL
df1 <- cbind(date = InfData1_mf$date[-1], df1)

############################################################
## 4. Weekly age–sex distribution (later period)
############################################################
InfData2_raw <- read.csv(
  "data/newly_confirmed_cases_detail_weekly_MHLW.csv",
  header = TRUE
)[2:71, 1:21]

colnames(InfData2_raw) <- c(
  "Week",
  "Male Under 10","Male 10s","Male 20s","Male 30s","Male 40s","Male 50s","Male 60s","Male 70s","Male 80s","Male Over 90",
  "Female Under 10","Female 10s","Female 20s","Female 30s","Female 40s","Female 50s","Female 60s","Female 70s","Female 80s","Female Over 90"
)

InfData2_raw[as.matrix(InfData2_raw) == "*"] <- 0

InfData2 <- data.frame(
  week      = InfData2_raw$Week,
  maleu20   = as.numeric(InfData2_raw$`Male Under 10`) + as.numeric(InfData2_raw$`Male 10s`),
  male20    = as.numeric(InfData2_raw$`Male 20s`),
  male30    = as.numeric(InfData2_raw$`Male 30s`),
  male40    = as.numeric(InfData2_raw$`Male 40s`),
  male50    = as.numeric(InfData2_raw$`Male 50s`),
  male60    = as.numeric(InfData2_raw$`Male 60s`),
  malea70   = as.numeric(InfData2_raw$`Male 70s`) +
    as.numeric(InfData2_raw$`Male 80s`) +
    as.numeric(InfData2_raw$`Male Over 90`),
  femaleu20 = as.numeric(InfData2_raw$`Female Under 10`) + as.numeric(InfData2_raw$`Female 10s`),
  female20  = as.numeric(InfData2_raw$`Female 20s`),
  female30  = as.numeric(InfData2_raw$`Female 30s`),
  female40  = as.numeric(InfData2_raw$`Female 40s`),
  female50  = as.numeric(InfData2_raw$`Female 50s`),
  female60  = as.numeric(InfData2_raw$`Female 60s`),
  femalea70 = as.numeric(InfData2_raw$`Female 70s`) +
    as.numeric(InfData2_raw$`Female 80s`) +
    as.numeric(InfData2_raw$`Female Over 90`)
)

InfData2_prop <- as.data.frame(
  t(apply(InfData2[, -1], 1, function(x) x / sum(x)))
)

InfData2_prop <- cbind(InfData2$week, InfData2_prop)

InfData2_prop <- InfData2_prop %>%
  separate(col = 1, into = c("StartDate", "EndDate"), sep = "~") %>%
  mutate(
    StartDate = as.Date(trimws(StartDate), format = "%Y/%m/%d"),
    EndDate   = as.Date(trimws(EndDate),   format = "%Y/%m/%d")
  ) %>%
  rowwise() %>%
  mutate(Date = list(seq(StartDate, EndDate, by = "day"))) %>%
  ungroup() %>%
  dplyr::select(-StartDate, -EndDate) %>% 
  unnest(cols = c(Date)) %>%
  relocate(Date)

df2 <- InfData2_prop[1:486, -1] * InfDataafter$ALL
df2 <- cbind(InfData2_prop$Date[1:486], df2)
colnames(df2)[1] <- "date"
############################################################
## 5. Aggregate full daily age–sex infections
############################################################
colnames(df1) <- colnames(df2)
df <- rbind(df1, df2)
df <- df[17:716, ]
rownames(df) <- NULL
############################################################
## 6. Define epidemic waves
############################################################
wave_periods <- list(
  wave1 = c(as.Date("2020-02-01"), as.Date("2020-07-09")),
  wave2 = c(as.Date("2020-06-30"), as.Date("2020-11-06")),
  wave3 = c(as.Date("2020-10-28"), as.Date("2021-03-06")),
  wave4 = c(as.Date("2021-02-25"), as.Date("2021-07-14")),
  wave5 = c(as.Date("2021-06-15"), as.Date("2021-12-31"))
)

extract_wave <- function(df, start_date, end_date) {
  df[df$date >= start_date & df$date <= end_date, -1]
}

wave1df <- extract_wave(df, wave_periods$wave1[1], wave_periods$wave1[2])
wave2df <- extract_wave(df, wave_periods$wave2[1], wave_periods$wave2[2])
wave3df <- extract_wave(df, wave_periods$wave3[1], wave_periods$wave3[2])
wave4df <- extract_wave(df, wave_periods$wave4[1], wave_periods$wave4[2])
wave5df <- extract_wave(df, wave_periods$wave5[1], wave_periods$wave5[2])

############################################################
## 7. Population structure
############################################################
population <- data.frame(
  age = c("4","9","14","19","24","29","34","39","44","49",
          "54","59","64","69","74","79","84","89","94","99","100"),
  male = c(2413,2588,2736,2949,3303,3238,3404,3779,4313,4954,
           4333,3887,3675,4073,4256,3201,2230,1317,504,97,10),
  female = c(2295,2468,2605,2807,3093,3036,3257,3677,4207,4848,
             4286,3896,3771,4334,4761,3962,3151,2386,1313,419,67)
)

population_group <- data.frame(
  age = c("a20","20","30","40","50","60","a70"),
  male = 1000 * c(
    sum(population$male[1:4]),
    sum(population$male[5:6]),
    sum(population$male[7:8]),
    sum(population$male[9:10]),
    sum(population$male[11:12]),
    sum(population$male[13:14]),
    sum(population$male[15:21])
  ),
  female = 1000 * c(
    sum(population$female[1:4]),
    sum(population$female[5:6]),
    sum(population$female[7:8]),
    sum(population$female[9:10]),
    sum(population$female[11:12]),
    sum(population$female[13:14]),
    sum(population$female[15:21])
  )
)

popN <- c(population_group$male, population_group$female)
popsum <- sum(popN)
popN_mat <- outer(popN, popN, "/")

############################################################
## 8. Contact matrix construction
############################################################
C <- matrix(
  c(
    0.8,0.6,0.1,0,0.3,0.8,0.3,0.1,0.2,0.1,
    0.6,2.7,0.3,0.1,0.3,2.1,0.8,0.1,0.2,0.2,
    0.1,0.3,1.5,0.3,0.2,0.8,0.9,0.4,0.2,0.1,
    0,0.1,0.3,1.8,0.4,0.3,0.5,0.6,0.1,0.1,
    0.3,0.3,0.2,0.4,2.5,1.2,1.3,1.5,1,0.7,
    0.8,2.1,0.8,0.3,1.2,2.6,1.6,1,1.6,0.8,
    0.3,0.8,0.9,0.5,1.3,1.6,3.3,1.6,1.3,1.3,
    0.1,0.1,0.4,0.6,1.5,1,1.6,2.5,1.2,0.9,
    0.2,0.2,0.2,0.1,1,1.6,1.3,1.2,2.9,1.4,
    0.1,0.2,0.1,0.1,0.7,0.8,1.3,0.9,1.4,2.1
  ),
  nrow = 10, byrow = TRUE
)

rownames(C) <- colnames(C) <- c(
  "0-4","5-9","10-14","15-19",
  "20-29","30-39","40-49","50-59","60-69","70-"
)

group_mapping <- list(
  "0-4"   = 4,
  "5-9"   = 9,
  "10-14" = 14,
  "15-19" = 19,
  "20-29" = c(24,29),
  "30-39" = c(34,39),
  "40-49" = c(44,49),
  "50-59" = c(54,59),
  "60-69" = c(64,69),
  "70-"   = c(74,79,84,89,94,99,100)
)

group_pop <- sapply(group_mapping, function(ages) {
  rows <- which(population$age %in% ages)
  sum(population$male[rows] + population$female[rows])
})

merge_mapping <- list(
  "0-19"  = c("0-4","5-9","10-14","15-19"),
  "20-29" = "20-29",
  "30-39" = "30-39",
  "40-49" = "40-49",
  "50-59" = "50-59",
  "60-69" = "60-69",
  "70+"   = "70-"
)

new_group_pop <- sapply(merge_mapping, function(groups) {
  sum(group_pop[groups])
})

new_C <- matrix(0, nrow = 7, ncol = 7)
rownames(new_C) <- colnames(new_C) <- names(new_group_pop)

for (i in 1:7) {
  for (j in 1:7) {
    numer <- 0
    denom <- 0
    for (old_i in merge_mapping[[i]]) {
      for (old_j in merge_mapping[[j]]) {
        numer <- numer + group_pop[old_i] * C[old_i, old_j]
      }
      denom <- denom + group_pop[old_i]
    }
    new_C[i, j] <- numer / denom
  }
}

# Expand to sex-specific 14 × 14 matrix and normalize
C_star <- rbind(
  cbind(new_C, new_C),
  cbind(new_C, new_C)
)

rownames(C_star) <- colnames(C_star) <- c(
  paste0(names(new_group_pop), "_male"),
  paste0(names(new_group_pop), "_female")
)

C <- sweep(C_star, 1, rowSums(C_star), FUN = "/")

############################################################
## 9. SEIR model
############################################################
###### Wave-specific initial values (sigma, q, alpha)

## Overdispersion parameter (sigma)
sigma <- read.csv(
  "data/sigma_initial.csv",
  header = TRUE
)[, 1]

sigma1 <- sigma[1]
sigma2 <- sigma[2]
sigma3 <- sigma[3]
sigma4 <- sigma[4]
sigma5 <- sigma[5]

## Control parameters (q)
q <- read.csv(
  "data/q_initial.csv",
  header = TRUE
)[, 2]

initq1 <- q[1:15]
initq2 <- q[16:27]
initq3 <- q[28:39]
initq4 <- q[40:50]
initq5 <- q[51:70]

## Group-specific exposure parameters (alpha)
alpha <- read.csv(
  "data/alpha_initial.csv",
  header = TRUE
)[, 2]

initalpha1 <- alpha[1:14]
initalpha2 <- alpha[15:28]
initalpha3 <- alpha[29:42]
initalpha4 <- alpha[43:56]
initalpha5 <- alpha[57:70]

###### model construction and estimation

## R0
calculate_R0 <- function(alpha, q_vec) {
  D_E <- 4.5
  D_I <- 1.5
  ratio <- 3  
  alpha_mat <- matrix(alpha, nrow = length(alpha), ncol = length(alpha), byrow = FALSE)
  R0_values <- numeric(length(q_vec))
  for (k in 1:length(q_vec)) {
    beta_E_mat <- q_vec[k] * alpha_mat * C
    beta_I_mat <- ratio * beta_E_mat
    K <- beta_E_mat * popN_mat * D_E +
      beta_I_mat * popN_mat * D_I
    R0_values[k] <- max(Re(eigen(K)$values))
  }
  return(list(R0_max = max(R0_values), R0_values = R0_values))
}

### wave1
seir<-function(time,y,params){

  q <- if (time <= 10) params$q1 else
    if (time <= 20) params$q2 else
      if (time <= 30) params$q3 else
        if (time <= 40) params$q4 else
          if (time <= 50) params$q5 else
            if (time <= 60) params$q6 else
              if (time <= 70) params$q7 else
                if (time <= 80) params$q8 else
                  if (time <= 90) params$q9 else
                    if (time <= 100) params$q10 else
                      if (time <= 110) params$q11 else
                        if (time <= 120) params$q12 else
                          if (time <= 130) params$q13 else
                            if (time <= 140) params$q14 else
                              if (time <= 150) params$q15 else
                                params$q16
  beta_mat <- q * C

  epsilon=params[["epsilon"]]
  gamma=params[["gamma"]]
  alpha1<-params[["alpha1"]]
  alpha2<-params[["alpha2"]]
  alpha3<-params[["alpha3"]]
  alpha4<-params[["alpha4"]]
  alpha5<-params[["alpha5"]]
  alpha6<-params[["alpha6"]]
  alpha7<-params[["alpha7"]]
  alpha8<-params[["alpha8"]]
  alpha9<-params[["alpha9"]]
  alpha10<-params[["alpha10"]]
  alpha11<-params[["alpha11"]]
  alpha12<-params[["alpha12"]]
  alpha13<-params[["alpha13"]]
  alpha14<-params[["alpha14"]]

  alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14)

  S<-y[1:n_group]
  E<-y[(n_group+1):(2*n_group)]
  I<-y[(2*n_group+1):(3*n_group)]
  R<-y[(3*n_group+1):(4*n_group)]

  dS<-numeric(n_group)
  dE<-numeric(n_group)
  dI<-numeric(n_group)
  dR<-numeric(n_group)

  for(i in 1:14){
    infection_force <- sum(beta_mat[i,] * S[i] * alpha[i] * (3*I + E)/popN)
    dS[i] <- -infection_force
    dE[i] <- infection_force - epsilon * E[i]
    dI[i] <- epsilon * E[i] - gamma * I[i]
    dR[i] <- gamma * I[i]
  }

  solution=c(dS,dE,dI,dR)
  return(list(solution))
}
nonlinear_constraint <- function(par) {
  alpha <- par[17:30]
  q_vec <- par[1:16]
  sigma <- par[31]
  result <- calculate_R0(alpha, q_vec)
  R0_max <- result$R0_max
  constraint1 <- R0_max - 2.5
  constraint2 <- sigma - (sigma1+0.01)
  return(c(constraint1,constraint2))
}
n_group<-14
init=c(population_group$male,population_group$female,rep(0.01,n_group),rep(0,n_group),rep(0,n_group))
times=seq(1,nrow(wave1df),by=1)
epsilon=1/4.5
f<-function(par){

  paras<-list(q1=par[1],q2=par[2],q3=par[3],q4=par[4],q5=par[5],q6=par[6],q7=par[7],q8=par[8],q9=par[9],q10=par[10],
              q11=par[11],q12=par[12],q13=par[13],q14=par[14],q15=par[15],q16=par[16],
              alpha1=par[17],alpha2=par[18],alpha3=par[19],alpha4=par[20],alpha5=par[21],alpha6=par[22],alpha7=par[23],
              alpha8=par[24],alpha9=par[25],alpha10=par[26],alpha11=par[27],alpha12=par[28],alpha13=par[29],alpha14=par[30],
              epsilon=1/4.5,gamma=1/1.5)
  sigma <- par[31]
  sigma_vec <- rep(sigma , 14)

  out <- tryCatch(ode(y=init, times=times, func=seir, parms=paras, method="lsoda", rtol=1e-8, atol=1e-8), error=function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) return(1e30)
  out=as.data.frame(out)
  out<-out[,-1]
  estE<-as.matrix(out[,(n_group+1):(2*n_group)])
  lambda=epsilon*estE

  Iobs <- as.matrix(wave1df)
  lnL <- numeric(14)
  for(i in 1:14){
    lambda_i <- pmax(lambda[, i], 1e-10)
    obs_i <- Iobs[, i]
    sigma_i <- sigma_vec[i]

    lnL[i] <- sum(
      lgamma(obs_i + 1/sigma_i) - lgamma(1/sigma_i) +
        obs_i * log(lambda_i) + obs_i * log(sigma_i) -
        (obs_i + 1/sigma_i) * log(1 + lambda_i * sigma_i)
    )
  }

  obj <- -sum(lnL)
  if (!is.finite(obj)) obj <- 1e30
  return(obj)
}
intipar<-c(initq1,0.02,initalpha1,sigma1)
model<-nloptr(
  x0 = intipar,
  eval_f = f,
  eval_g_ineq = nonlinear_constraint,
  opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-6, "maxeval" = 1000)
)
estpar<-model$solution
allpara<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
              q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],q12=estpar[12],q13=estpar[13],q14=estpar[14],
              q15=estpar[15],q16=estpar[16],
              alpha1=estpar[17],alpha2=estpar[18],alpha3=estpar[19],alpha4=estpar[20],alpha5=estpar[21],alpha6=estpar[22],alpha7=estpar[23],
              alpha8=estpar[24],alpha9=estpar[25],alpha10=estpar[26],alpha11=estpar[27],alpha12=estpar[28],alpha13=estpar[29],alpha14=estpar[30],
              epsilon=1/4.5,gamma=1/1.5)
out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
out=as.data.frame(out)
out<-out[,-1]
estdailyI<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
alpha<-estpar[17:30]
q_vec<-estpar[1:15]
result <- calculate_R0(alpha, q_vec)
R0_values <- result$R0_values
#
q1<-q_vec
alpha1<-alpha
R01<-R0_values
estdailyI1<-estdailyI[1:150,]
estsigma1<-estpar[31]
initnew <- as.numeric(out[150,])
##CI
get_inv <- function(mat) {
  cov_mat <- ginv(mat)
  eig_vals <- eigen(cov_mat, symmetric = TRUE)$values
  is_posdef <- all(eig_vals > 0)
  if (!is_posdef) {
    message("Matrix is not positive definite; applying nearPD for correction.")
    cov_mat_pd <- nearPD(cov_mat)$mat
    cov_mat_pd <- as.matrix(cov_mat_pd)
  } else {
    cov_mat_pd <- cov_mat
  }
  return(cov_mat_pd)
}
hess_mat<-hessian(func=f,x=estpar)
cov_mat<-get_inv(hess_mat)
boot_pars <- mvrnorm(1000,mu=estpar,Sigma=cov_mat)
q_boot <- boot_pars[, 1:16]
alpha_boot <- boot_pars[, 17:30]
sigma_boot <- boot_pars[, 31]
R0_values<-matrix(0,nrow=1000,ncol=length(R0_values))
for(i in 1:1000){
  result <- calculate_R0(alpha_boot[i,], q_boot[i,1:ncol(R0_values)])
  R0_values[i,]<-result$R0_values
}
estdailyI_boot<-list()
for(i in 1:1000){
  estpar<-c(q_boot[i,],alpha_boot[i,])
  allpara<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
                q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],q12=estpar[12],q13=estpar[13],q14=estpar[14],
                q15=estpar[15],q16=estpar[16],
                alpha1=estpar[17],alpha2=estpar[18],alpha3=estpar[19],alpha4=estpar[20],alpha5=estpar[21],alpha6=estpar[22],alpha7=estpar[23],
                alpha8=estpar[24],alpha9=estpar[25],alpha10=estpar[26],alpha11=estpar[27],alpha12=estpar[28],alpha13=estpar[29],alpha14=estpar[30],
                epsilon=1/4.5,gamma=1/1.5)
  out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
  out=as.data.frame(out)
  out<-out[,-1]
  estdailyI_boot[[i]]<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
}
estdailyI_boot_nb <- vector("list", length(estdailyI_boot))
for(i in seq_along(estdailyI_boot)) {
  lambda_mat <- estdailyI_boot[[i]]
  sigma_vec <- rep(sigma_boot[i],14)
  obs_mat <- matrix(NA, nrow = nrow(lambda_mat), ncol = ncol(lambda_mat))
  for(g in 1:ncol(lambda_mat)) {
    mu_vec <- lambda_mat[, g]
    size_param <- 1 / sigma_vec[g]
    mu_vec <- pmax(mu_vec, 1e-6)
    size_param <- max(size_param, 1e-6)
    obs_mat[, g] <- rnbinom(length(mu_vec), size = size_param, mu = mu_vec)
  }
  estdailyI_boot_nb[[i]] <- obs_mat
}
array_boot <- array(unlist(estdailyI_boot_nb),
                    dim = c(nrow(estdailyI_boot_nb[[1]]),
                            ncol(estdailyI_boot_nb[[1]]),
                            length(estdailyI_boot_nb)))
lower_CI <- matrix(NA, nrow = dim(array_boot)[1], ncol = dim(array_boot)[2])
upper_CI <- matrix(NA, nrow = dim(array_boot)[1], ncol = dim(array_boot)[2])
for (t in 1:dim(array_boot)[1]) {
  for (g in 1:dim(array_boot)[2]) {
    vals <- array_boot[t, g, ]
    lower_CI[t, g] <- quantile(vals, 0.025)
    upper_CI[t, g] <- quantile(vals, 0.975)
  }
}
#
q_boot1<-q_boot[,1:15]
alpha_boot1<-alpha_boot
sigma_boot1<-sigma_boot
q_CI1 <- apply(q_boot, 2, quantile, probs = c(0.025, 0.975))[,1:15]
alpha_CI1 <- apply(alpha_boot, 2, quantile, probs = c(0.025, 0.975))
R0_CI1 <- apply(R0_values, 2, quantile, probs = c(0.025, 0.975))
estsigma_CI1 <- quantile(sigma_boot, probs = c(0.025, 0.975))
estdailyI1_lower_CI<-lower_CI[1:150,]
estdailyI1_upper_CI<-upper_CI[1:150,]

###wave2
seir<-function(time,y,params){

  q <- if (time <= 10) params$q1 else
    if (time <= 20) params$q2 else
      if (time <= 30) params$q3 else
        if (time <= 40) params$q4 else
          if (time <= 50) params$q5 else
            if (time <= 60) params$q6 else
              if (time <= 70) params$q7 else
                if (time <= 80) params$q8 else
                  if (time <= 90) params$q9 else
                    if (time <= 100) params$q10 else
                      if (time <= 110) params$q11 else
                        if (time <= 120) params$q12 else
                          params$q13
  beta_mat <- q * C

  epsilon=params[["epsilon"]]
  gamma=params[["gamma"]]
  alpha1<-params[["alpha1"]]
  alpha2<-params[["alpha2"]]
  alpha3<-params[["alpha3"]]
  alpha4<-params[["alpha4"]]
  alpha5<-params[["alpha5"]]
  alpha6<-params[["alpha6"]]
  alpha7<-params[["alpha7"]]
  alpha8<-params[["alpha8"]]
  alpha9<-params[["alpha9"]]
  alpha10<-params[["alpha10"]]
  alpha11<-params[["alpha11"]]
  alpha12<-params[["alpha12"]]
  alpha13<-params[["alpha13"]]
  alpha14<-params[["alpha14"]]

  alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14)

  S<-y[1:n_group]
  E<-y[(n_group+1):(2*n_group)]
  I<-y[(2*n_group+1):(3*n_group)]
  R<-y[(3*n_group+1):(4*n_group)]

  dS<-numeric(n_group)
  dE<-numeric(n_group)
  dI<-numeric(n_group)
  dR<-numeric(n_group)

  for(i in 1:14){
    infection_force <- sum(beta_mat[i,] * S[i] * alpha[i] * (3*I + E)/popN)
    dS[i] <- -infection_force
    dE[i] <- infection_force - epsilon * E[i]
    dI[i] <- epsilon * E[i] - gamma * I[i]
    dR[i] <- gamma * I[i]
  }

  solution=c(dS,dE,dI,dR)
  return(list(solution))
}
nonlinear_constraint <- function(par) {
  alpha <- par[14:27]
  q_vec <- par[1:13]
  sigma <- par[28]
  result <- calculate_R0(alpha, q_vec)
  R0_max <- result$R0_max
  constraint1 <- R0_max - 2.5
  constraint2 <- sigma - (sigma2+0.01)
  return(c(constraint1,constraint2))
}
init=initnew
times=seq(1,nrow(wave2df),by=1)
f<-function(par){

  paras<-list(q1=par[1],q2=par[2],q3=par[3],q4=par[4],q5=par[5],q6=par[6],q7=par[7],q8=par[8],q9=par[9],q10=par[10],
              q11=par[11],q12=par[12],q13=par[13],
              alpha1=par[14],alpha2=par[15],alpha3=par[16],alpha4=par[17],alpha5=par[18],alpha6=par[19],alpha7=par[20],
              alpha8=par[21],alpha9=par[22],alpha10=par[23],alpha11=par[24],alpha12=par[25],alpha13=par[26],alpha14=par[27],
              epsilon=1/4.5,gamma=1/1.5)
  sigma <-par[28]
  sigma_vec <- rep(sigma , 14)

  out <- tryCatch(ode(y=init, times=times, func=seir, parms=paras, method="lsoda", rtol=1e-8, atol=1e-8), error=function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) return(1e30)
  out=as.data.frame(out)
  out<-out[,-1]
  estE<-as.matrix(out[,(n_group+1):(2*n_group)])
  lambda=epsilon*estE

  Iobs <- as.matrix(wave2df)
  lnL <- numeric(14)
  for(i in 1:14){
    lambda_i <- pmax(lambda[, i], 1e-10)
    obs_i <- Iobs[, i]
    sigma_i <- sigma_vec[i]

    lnL[i] <- sum(
      lgamma(obs_i + 1/sigma_i) - lgamma(1/sigma_i) +
        obs_i * log(lambda_i) + obs_i * log(sigma_i) -
        (obs_i + 1/sigma_i) * log(1 + lambda_i * sigma_i)
    )
  }

  obj <- -sum(lnL)
  if (!is.finite(obj)) obj <- 1e30
  return(obj)
}
intipar<-c(initq2,0.02,initalpha2,sigma2)
model<-nloptr(
  x0 = intipar,
  eval_f = f,
  eval_g_ineq = nonlinear_constraint,
  opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-6, "maxeval" = 1000)
)
estpar<-model$solution
allpara<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
              q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],q12=estpar[12],q13=estpar[13],
              alpha1=estpar[14],alpha2=estpar[15],alpha3=estpar[16],alpha4=estpar[17],alpha5=estpar[18],alpha6=estpar[19],alpha7=estpar[20],
              alpha8=estpar[21],alpha9=estpar[22],alpha10=estpar[23],alpha11=estpar[24],alpha12=estpar[25],alpha13=estpar[26],alpha14=estpar[27],
              epsilon=1/4.5,gamma=1/1.5)
out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
out=as.data.frame(out)
out<-out[,-1]
estdailyI<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
alpha<-estpar[14:27]
q_vec<-estpar[1:12]
result <- calculate_R0(alpha, q_vec)
R0_values <- result$R0_values
#
q2<-q_vec
alpha2<-alpha
R02<-R0_values
estdailyI2<-estdailyI[1:120,]
estsigma2<-estpar[28]
initnew <- as.numeric(out[120,])
##CI
hess_mat<-hessian(func=f,x=estpar)
cov_mat<-get_inv(hess_mat)
boot_pars <- mvrnorm(1000,mu=estpar,Sigma=cov_mat)
q_boot <- boot_pars[, 1:13]
alpha_boot <- boot_pars[, 14:27]
sigma_boot <- boot_pars[, 28]
R0_values<-matrix(0,nrow=1000,ncol=length(R0_values))
for(i in 1:1000){
  result <- calculate_R0(alpha_boot[i,], q_boot[i,1:ncol(R0_values)])
  R0_values[i,]<-result$R0_values
}
estdailyI_boot<-list()
for(i in 1:1000){
  estpar<-c(q_boot[i,],alpha_boot[i,])
  allpara<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
                q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],q12=estpar[12],q13=estpar[13],
                alpha1=estpar[14],alpha2=estpar[15],alpha3=estpar[16],alpha4=estpar[17],alpha5=estpar[18],alpha6=estpar[19],alpha7=estpar[20],
                alpha8=estpar[21],alpha9=estpar[22],alpha10=estpar[23],alpha11=estpar[24],alpha12=estpar[25],alpha13=estpar[26],alpha14=estpar[27],
                epsilon=1/4.5,gamma=1/1.5)
  out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
  out=as.data.frame(out)
  out<-out[,-1]
  estdailyI_boot[[i]]<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
}
estdailyI_boot_nb <- vector("list", length(estdailyI_boot))
for(i in seq_along(estdailyI_boot)) {
  lambda_mat <- estdailyI_boot[[i]]
  sigma_vec <- rep(sigma_boot[i],14)
  obs_mat <- matrix(NA, nrow = nrow(lambda_mat), ncol = ncol(lambda_mat))
  for(g in 1:ncol(lambda_mat)) {
    mu_vec <- lambda_mat[, g]
    size_param <- 1 / sigma_vec[g]
    mu_vec <- pmax(mu_vec, 1e-6)
    size_param <- max(size_param, 1e-6)
    obs_mat[, g] <- rnbinom(length(mu_vec), size = size_param, mu = mu_vec)
  }
  estdailyI_boot_nb[[i]] <- obs_mat
}
array_boot <- array(unlist(estdailyI_boot_nb),
                    dim = c(nrow(estdailyI_boot_nb[[1]]),
                            ncol(estdailyI_boot_nb[[1]]),
                            length(estdailyI_boot_nb)))
lower_CI <- matrix(NA, nrow = dim(array_boot)[1], ncol = dim(array_boot)[2])
upper_CI <- matrix(NA, nrow = dim(array_boot)[1], ncol = dim(array_boot)[2])
for (t in 1:dim(array_boot)[1]) {
  for (g in 1:dim(array_boot)[2]) {
    vals <- array_boot[t, g, ]
    lower_CI[t, g] <- quantile(vals, 0.025)
    upper_CI[t, g] <- quantile(vals, 0.975)
  }
}
#
q_boot2<-q_boot[,1:12]
alpha_boot2<-alpha_boot
sigma_boot2<-sigma_boot
q_CI2 <- apply(q_boot, 2, quantile, probs = c(0.025, 0.975))[,1:12]
alpha_CI2 <- apply(alpha_boot, 2, quantile, probs = c(0.025, 0.975))
R0_CI2 <- apply(R0_values, 2, quantile, probs = c(0.025, 0.975))
estsigma_CI2 <- quantile(sigma_boot, probs = c(0.025, 0.975))
estdailyI2_lower_CI<-lower_CI[1:120,]
estdailyI2_upper_CI<-upper_CI[1:120,]

###wave3
seir<-function(time,y,params){

  q <- if (time <= 10) params$q1 else
    if (time <= 20) params$q2 else
      if (time <= 30) params$q3 else
        if (time <= 40) params$q4 else
          if (time <= 50) params$q5 else
            if (time <= 60) params$q6 else
              if (time <= 70) params$q7 else
                if (time <= 80) params$q8 else
                  if (time <= 90) params$q9 else
                    if (time <= 100) params$q10 else
                      if (time <= 110) params$q11 else
                        if (time <= 120) params$q12 else
                          params$q13
  beta_mat <- q * C

  epsilon=params[["epsilon"]]
  gamma=params[["gamma"]]
  alpha1<-params[["alpha1"]]
  alpha2<-params[["alpha2"]]
  alpha3<-params[["alpha3"]]
  alpha4<-params[["alpha4"]]
  alpha5<-params[["alpha5"]]
  alpha6<-params[["alpha6"]]
  alpha7<-params[["alpha7"]]
  alpha8<-params[["alpha8"]]
  alpha9<-params[["alpha9"]]
  alpha10<-params[["alpha10"]]
  alpha11<-params[["alpha11"]]
  alpha12<-params[["alpha12"]]
  alpha13<-params[["alpha13"]]
  alpha14<-params[["alpha14"]]

  alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14)

  S<-y[1:n_group]
  E<-y[(n_group+1):(2*n_group)]
  I<-y[(2*n_group+1):(3*n_group)]
  R<-y[(3*n_group+1):(4*n_group)]

  dS<-numeric(n_group)
  dE<-numeric(n_group)
  dI<-numeric(n_group)
  dR<-numeric(n_group)

  for(i in 1:14){
    infection_force <- sum(beta_mat[i,] * S[i] * alpha[i] * (3*I + E)/popN)
    dS[i] <- -infection_force
    dE[i] <- infection_force - epsilon * E[i]
    dI[i] <- epsilon * E[i] - gamma * I[i]
    dR[i] <- gamma * I[i]
  }

  solution=c(dS,dE,dI,dR)
  return(list(solution))
}
nonlinear_constraint <- function(par) {
  alpha <- par[14:27]
  q_vec <- par[1:13]
  sigma <- par[28]
  result <- calculate_R0(alpha, q_vec)
  R0_max <- result$R0_max
  constraint1 <- R0_max - 2.5
  constraint2 <- sigma - (sigma3+0.01)
  return(c(constraint1,constraint2))
}
init=initnew
times=seq(1,nrow(wave3df),by=1)
f<-function(par){

  paras<-list(q1=par[1],q2=par[2],q3=par[3],q4=par[4],q5=par[5],q6=par[6],q7=par[7],q8=par[8],q9=par[9],q10=par[10],
              q11=par[11],q12=par[12],q13=par[13],
              alpha1=par[14],alpha2=par[15],alpha3=par[16],alpha4=par[17],alpha5=par[18],alpha6=par[19],alpha7=par[20],
              alpha8=par[21],alpha9=par[22],alpha10=par[23],alpha11=par[24],alpha12=par[25],alpha13=par[26],alpha14=par[27],
              epsilon=1/4.5,gamma=1/1.5)
  sigma <- par[28]
  sigma_vec <- rep(sigma , 14)

  out <- tryCatch(ode(y=init, times=times, func=seir, parms=paras, method="lsoda", rtol=1e-8, atol=1e-8), error=function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) return(1e30)
  out=as.data.frame(out)
  out<-out[,-1]
  estE<-as.matrix(out[,(n_group+1):(2*n_group)])
  lambda=epsilon*estE

  Iobs <- as.matrix(wave3df)
  lnL <- numeric(14)
  for(i in 1:14){
    lambda_i <- pmax(lambda[, i], 1e-10)
    obs_i <- Iobs[, i]
    sigma_i <- sigma_vec[i]

    lnL[i] <- sum(
      lgamma(obs_i + 1/sigma_i) - lgamma(1/sigma_i) +
        obs_i * log(lambda_i) + obs_i * log(sigma_i) -
        (obs_i + 1/sigma_i) * log(1 + lambda_i * sigma_i)
    )
  }

  obj <- -sum(lnL)
  if (!is.finite(obj)) obj <- 1e30
  return(obj)
}
intipar<-c(initq3,0.01,initalpha3,sigma3)
model<-nloptr(
  x0 = intipar,
  eval_f = f,
  eval_g_ineq = nonlinear_constraint,
  opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-6, "maxeval" = 1000)
)
estpar<-model$solution
allpara<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
              q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],q12=estpar[12],q13=estpar[13],
              alpha1=estpar[14],alpha2=estpar[15],alpha3=estpar[16],alpha4=estpar[17],alpha5=estpar[18],alpha6=estpar[19],alpha7=estpar[20],
              alpha8=estpar[21],alpha9=estpar[22],alpha10=estpar[23],alpha11=estpar[24],alpha12=estpar[25],alpha13=estpar[26],alpha14=estpar[27],
              epsilon=1/4.5,gamma=1/1.5)
out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
out=as.data.frame(out)
out<-out[,-1]
estdailyI<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
alpha<-estpar[14:27]
q_vec<-estpar[1:12]
result <- calculate_R0(alpha, q_vec)
R0_values <- result$R0_values
#
q3<-q_vec
alpha3<-alpha
R03<-R0_values
estdailyI3<-estdailyI[1:120,]
estsigma3<-estpar[28]
initnew <- as.numeric(out[120,])
##CI
hess_mat<-hessian(func=f,x=estpar)
cov_mat<-get_inv(hess_mat)
boot_pars <- mvrnorm(1000,mu=estpar,Sigma=cov_mat)
q_boot <- boot_pars[, 1:13]
alpha_boot <- boot_pars[, 14:27]
sigma_boot <- boot_pars[, 28]
R0_values<-matrix(0,nrow=1000,ncol=length(R0_values))
for(i in 1:1000){
  result <- calculate_R0(alpha_boot[i,], q_boot[i,1:ncol(R0_values)])
  R0_values[i,]<-result$R0_values
}
estdailyI_boot<-list()
for(i in 1:1000){
  estpar<-c(q_boot[i,],alpha_boot[i,])
  allpara<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
                q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],q12=estpar[12],q13=estpar[13],
                alpha1=estpar[14],alpha2=estpar[15],alpha3=estpar[16],alpha4=estpar[17],alpha5=estpar[18],alpha6=estpar[19],alpha7=estpar[20],
                alpha8=estpar[21],alpha9=estpar[22],alpha10=estpar[23],alpha11=estpar[24],alpha12=estpar[25],alpha13=estpar[26],alpha14=estpar[27],
                epsilon=1/4.5,gamma=1/1.5)
  out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
  out=as.data.frame(out)
  out<-out[,-1]
  estdailyI_boot[[i]]<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
}
estdailyI_boot_nb <- vector("list", length(estdailyI_boot))
for(i in seq_along(estdailyI_boot)) {
  lambda_mat <- estdailyI_boot[[i]]
  sigma_vec <- rep(sigma_boot[i],14)
  obs_mat <- matrix(NA, nrow = nrow(lambda_mat), ncol = ncol(lambda_mat))
  for(g in 1:ncol(lambda_mat)) {
    mu_vec <- lambda_mat[, g]
    size_param <- 1 / sigma_vec[g]
    mu_vec <- pmax(mu_vec, 1e-6)
    size_param <- max(size_param, 1e-6)
    obs_mat[, g] <- rnbinom(length(mu_vec), size = size_param, mu = mu_vec)
  }
  estdailyI_boot_nb[[i]] <- obs_mat
}
array_boot <- array(unlist(estdailyI_boot_nb),
                    dim = c(nrow(estdailyI_boot_nb[[1]]),
                            ncol(estdailyI_boot_nb[[1]]),
                            length(estdailyI_boot_nb)))
lower_CI <- matrix(NA, nrow = dim(array_boot)[1], ncol = dim(array_boot)[2])
upper_CI <- matrix(NA, nrow = dim(array_boot)[1], ncol = dim(array_boot)[2])
for (t in 1:dim(array_boot)[1]) {
  for (g in 1:dim(array_boot)[2]) {
    vals <- array_boot[t, g, ]
    lower_CI[t, g] <- quantile(vals, 0.025)
    upper_CI[t, g] <- quantile(vals, 0.975)
  }
}
#
q_boot3<-q_boot[,1:12]
alpha_boot3<-alpha_boot
sigma_boot3<-sigma_boot
q_CI3 <- apply(q_boot, 2, quantile, probs = c(0.025, 0.975))[,1:12]
alpha_CI3 <- apply(alpha_boot, 2, quantile, probs = c(0.025, 0.975))
R0_CI3 <- apply(R0_values, 2, quantile, probs = c(0.025, 0.975))
estsigma_CI3 <- quantile(sigma_boot, probs = c(0.025, 0.975))
estdailyI3_lower_CI<-lower_CI[1:120,]
estdailyI3_upper_CI<-upper_CI[1:120,]

###wave4
seir<-function(time,y,params){

  q <- if (time <= 10) params$q1 else
    if (time <= 20) params$q2 else
      if (time <= 30) params$q3 else
        if (time <= 40) params$q4 else
          if (time <= 50) params$q5 else
            if (time <= 60) params$q6 else
              if (time <= 70) params$q7 else
                if (time <= 80) params$q8 else
                  if (time <= 90) params$q9 else
                    if (time <= 100) params$q10 else
                      if (time <= 110) params$q11 else
                        if (time <= 120) params$q12 else
                          if (time <= 130) params$q13 else
                            params$q14
  beta_mat <- q * C

  epsilon=params[["epsilon"]]
  gamma=params[["gamma"]]
  alpha1<-params[["alpha1"]]
  alpha2<-params[["alpha2"]]
  alpha3<-params[["alpha3"]]
  alpha4<-params[["alpha4"]]
  alpha5<-params[["alpha5"]]
  alpha6<-params[["alpha6"]]
  alpha7<-params[["alpha7"]]
  alpha8<-params[["alpha8"]]
  alpha9<-params[["alpha9"]]
  alpha10<-params[["alpha10"]]
  alpha11<-params[["alpha11"]]
  alpha12<-params[["alpha12"]]
  alpha13<-params[["alpha13"]]
  alpha14<-params[["alpha14"]]

  alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14)

  S<-y[1:n_group]
  E<-y[(n_group+1):(2*n_group)]
  I<-y[(2*n_group+1):(3*n_group)]
  R<-y[(3*n_group+1):(4*n_group)]

  dS<-numeric(n_group)
  dE<-numeric(n_group)
  dI<-numeric(n_group)
  dR<-numeric(n_group)

  for(i in 1:14){
    infection_force <- sum(beta_mat[i,] * S[i] * alpha[i] * (3*I + E)/popN)
    dS[i] <- -infection_force
    dE[i] <- infection_force - epsilon * E[i]
    dI[i] <- epsilon * E[i] - gamma * I[i]
    dR[i] <- gamma * I[i]
  }

  solution=c(dS,dE,dI,dR)
  return(list(solution))
}
nonlinear_constraint <- function(par) {
  alpha <- par[15:28]
  q_vec <- par[1:14]
  sigma <- par[29]
  result <- calculate_R0(alpha, q_vec)
  R0_max <- result$R0_max
  constraint1 <- R0_max - 2.5
  constraint2 <- sigma - (sigma4+0.01)
  return(c(constraint1,constraint2))
}
init=initnew
times=seq(1,nrow(wave4df),by=1)
f<-function(par){
  #
  q_vec  <- par[1:11]
  q_init <- initq4
  lambda_dev <- 10000
  pen_dev <- sum(((q_vec - q_init)/q_init)^2)
  #

  paras<-list(q1=par[1],q2=par[2],q3=par[3],q4=par[4],q5=par[5],q6=par[6],q7=par[7],q8=par[8],q9=par[9],q10=par[10],
              q11=par[11],q12=par[12],q13=par[13],q14=par[14],
              alpha1=par[15],alpha2=par[16],alpha3=par[17],alpha4=par[18],alpha5=par[19],alpha6=par[20],alpha7=par[21],
              alpha8=par[22],alpha9=par[23],alpha10=par[24],alpha11=par[25],alpha12=par[26],alpha13=par[27],alpha14=par[28],
              epsilon=1/4.5,gamma=1/1.5)
  sigma <- par[29]
  sigma_vec <- rep(sigma , 14)

  out <- tryCatch(ode(y=init, times=times, func=seir, parms=paras, method="lsoda", rtol=1e-8, atol=1e-8), error=function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) return(1e30)
  out=as.data.frame(out)
  out<-out[,-1]
  estE<-as.matrix(out[,(n_group+1):(2*n_group)])
  lambda=epsilon*estE

  Iobs <- as.matrix(wave4df)
  lnL <- numeric(14)
  for(i in 1:14){
    lambda_i <- pmax(lambda[, i], 1e-10)
    obs_i <- Iobs[, i]
    sigma_i <- sigma_vec[i]

    lnL[i] <- sum(
      lgamma(obs_i + 1/sigma_i) - lgamma(1/sigma_i) +
        obs_i * log(lambda_i) + obs_i * log(sigma_i) -
        (obs_i + 1/sigma_i) * log(1 + lambda_i * sigma_i)
    )
  }

  obj <- -sum(lnL) + lambda_dev * pen_dev
  if (!is.finite(obj)) obj <- 1e30
  return(obj)
}
intipar<-c(initq4,rep(0.01,3),initalpha4,sigma4)
model<-nloptr(
  x0 = intipar,
  eval_f = f,
  eval_g_ineq = nonlinear_constraint,
  opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-6, "maxeval" = 1000)
)
estpar<-model$solution
allpara<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
              q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],q12=estpar[12],q13=estpar[13],q14=estpar[14],
              alpha1=estpar[15],alpha2=estpar[16],alpha3=estpar[17],alpha4=estpar[18],alpha5=estpar[19],alpha6=estpar[20],alpha7=estpar[21],
              alpha8=estpar[22],alpha9=estpar[23],alpha10=estpar[24],alpha11=estpar[25],alpha12=estpar[26],alpha13=estpar[27],alpha14=estpar[28],
              epsilon=1/4.5,gamma=1/1.5)
out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
out=as.data.frame(out)
out<-out[,-1]
estdailyI<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
alpha<-estpar[15:28]
q_vec<-estpar[1:13]
result <- calculate_R0(alpha, q_vec)
R0_values <- result$R0_values
#
q4<-q_vec[1:11]
alpha4<-alpha
R04<-R0_values[1:11]
estdailyI4<-estdailyI[1:110,]
estsigma4<-estpar[29]
initnew <- as.numeric(out[110,])
##CI
hess_mat<-hessian(func=f,x=estpar)
cov_mat<-get_inv(hess_mat)
boot_pars <- mvrnorm(1000,mu=estpar,Sigma=cov_mat)
q_boot <- boot_pars[, 1:14]
alpha_boot <- boot_pars[, 15:28]
sigma_boot <- boot_pars[, 29]
R0_values<-matrix(0,nrow=1000,ncol=length(R0_values))
for(i in 1:1000){
  result <- calculate_R0(alpha_boot[i,], q_boot[i,1:ncol(R0_values)])
  R0_values[i,]<-result$R0_values
}
estdailyI_boot<-list()
for(i in 1:1000){
  estpar<-c(q_boot[i,],alpha_boot[i,])
  allpara<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
                q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],q12=estpar[12],q13=estpar[13],q14=estpar[14],
                alpha1=estpar[15],alpha2=estpar[16],alpha3=estpar[17],alpha4=estpar[18],alpha5=estpar[19],alpha6=estpar[20],alpha7=estpar[21],
                alpha8=estpar[22],alpha9=estpar[23],alpha10=estpar[24],alpha11=estpar[25],alpha12=estpar[26],alpha13=estpar[27],alpha14=estpar[28],
                epsilon=1/4.5,gamma=1/1.5)
  out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
  out=as.data.frame(out)
  out<-out[,-1]
  estdailyI_boot[[i]]<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
}
estdailyI_boot_nb <- vector("list", length(estdailyI_boot))
for(i in seq_along(estdailyI_boot)) {
  lambda_mat <- estdailyI_boot[[i]]
  sigma_vec <- rep(sigma_boot[i],14)
  obs_mat <- matrix(NA, nrow = nrow(lambda_mat), ncol = ncol(lambda_mat))
  for(g in 1:ncol(lambda_mat)) {
    mu_vec <- lambda_mat[, g]
    size_param <- 1 / sigma_vec[g]
    mu_vec <- pmax(mu_vec, 1e-6)
    size_param <- max(size_param, 1e-6)
    obs_mat[, g] <- rnbinom(length(mu_vec), size = size_param, mu = mu_vec)
  }
  estdailyI_boot_nb[[i]] <- obs_mat
}
array_boot <- array(unlist(estdailyI_boot_nb),
                    dim = c(nrow(estdailyI_boot_nb[[1]]),
                            ncol(estdailyI_boot_nb[[1]]),
                            length(estdailyI_boot_nb)))
lower_CI <- matrix(NA, nrow = dim(array_boot)[1], ncol = dim(array_boot)[2])
upper_CI <- matrix(NA, nrow = dim(array_boot)[1], ncol = dim(array_boot)[2])
for (t in 1:dim(array_boot)[1]) {
  for (g in 1:dim(array_boot)[2]) {
    vals <- array_boot[t, g, ]
    lower_CI[t, g] <- quantile(vals, 0.025)
    upper_CI[t, g] <- quantile(vals, 0.975)
  }
}
#
q_boot4<-q_boot[,1:11]
alpha_boot4<-alpha_boot
sigma_boot4<-sigma_boot
q_CI4 <- apply(q_boot, 2, quantile, probs = c(0.025, 0.975))[,1:11]
alpha_CI4 <- apply(alpha_boot, 2, quantile, probs = c(0.025, 0.975))
R0_CI4 <- apply(R0_values, 2, quantile, probs = c(0.025, 0.975))[,1:11]
estsigma_CI4 <- quantile(sigma_boot, probs = c(0.025, 0.975))
estdailyI4_lower_CI<-lower_CI[1:110,]
estdailyI4_upper_CI<-upper_CI[1:110,]

###wave5
seir<-function(time,y,params){

  q <- if (time <= 10) params$q1 else
    if (time <= 20) params$q2 else
      if (time <= 30) params$q3 else
        if (time <= 40) params$q4 else
          if (time <= 50) params$q5 else
            if (time <= 60) params$q6 else
              if (time <= 70) params$q7 else
                if (time <= 80) params$q8 else
                  if (time <= 90) params$q9 else
                    if (time <= 100) params$q10 else
                      if (time <= 110) params$q11 else
                        if (time <= 120) params$q12 else
                          if (time <= 130) params$q13 else
                            if (time <= 140) params$q14 else
                              if (time <= 150) params$q15 else
                                if (time <= 160) params$q16 else
                                  if (time <= 170) params$q17 else
                                    if (time <= 180) params$q18 else
                                      if (time <= 190) params$q19 else
                                        params$q20
  beta_mat <- q * C

  epsilon=params[["epsilon"]]
  gamma=params[["gamma"]]
  alpha1<-params[["alpha1"]]
  alpha2<-params[["alpha2"]]
  alpha3<-params[["alpha3"]]
  alpha4<-params[["alpha4"]]
  alpha5<-params[["alpha5"]]
  alpha6<-params[["alpha6"]]
  alpha7<-params[["alpha7"]]
  alpha8<-params[["alpha8"]]
  alpha9<-params[["alpha9"]]
  alpha10<-params[["alpha10"]]
  alpha11<-params[["alpha11"]]
  alpha12<-params[["alpha12"]]
  alpha13<-params[["alpha13"]]
  alpha14<-params[["alpha14"]]

  alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14)

  S<-y[1:n_group]
  E<-y[(n_group+1):(2*n_group)]
  I<-y[(2*n_group+1):(3*n_group)]
  R<-y[(3*n_group+1):(4*n_group)]

  dS<-numeric(n_group)
  dE<-numeric(n_group)
  dI<-numeric(n_group)
  dR<-numeric(n_group)

  for(i in 1:14){
    infection_force <- sum(beta_mat[i,] * S[i] * alpha[i] * (3*I + E)/popN)
    dS[i] <- -infection_force
    dE[i] <- infection_force - epsilon * E[i]
    dI[i] <- epsilon * E[i] - gamma * I[i]
    dR[i] <- gamma * I[i]
  }

  solution=c(dS,dE,dI,dR)
  return(list(solution))
}
nonlinear_constraint <- function(par) {
  alpha <- par[21:34]
  q_vec <- par[1:20]
  sigma <- par[35]
  result <- calculate_R0(alpha, q_vec)
  R0_max <- result$R0_max
  constraint1 <- R0_max - 2.5
  constraint2 <- sigma - (sigma5+0.01)
  return(c(constraint1,constraint2))
}
init=initnew
times=seq(1,nrow(wave5df),by=1)
f<-function(par){
  ##
  q_vec  <- par[1:20]
  q_init <- initq5
  pen_q <- sum(((q_vec - q_init)/q_init)^2)
  #
  alpha_vec  <- par[21:34]
  alpha_init <- initalpha5
  pen_alpha <- sum(((alpha_vec - alpha_init)/alpha_init)^2)
  #
  pen_dev <- pen_q + pen_alpha
  lambda_dev <- 10000
  ##

  paras<-list(q1=par[1],q2=par[2],q3=par[3],q4=par[4],q5=par[5],q6=par[6],q7=par[7],q8=par[8],q9=par[9],q10=par[10],
              q11=par[11],q12=par[12],q13=par[13],q14=par[14],q15=par[15],q16=par[16],q17=par[17],q18=par[18],q19=par[19],q20=par[20],
              alpha1=par[21],alpha2=par[22],alpha3=par[23],alpha4=par[24],alpha5=par[25],alpha6=par[26],alpha7=par[27],
              alpha8=par[28],alpha9=par[29],alpha10=par[30],alpha11=par[31],alpha12=par[32],alpha13=par[33],alpha14=par[34],
              epsilon=1/4.5,gamma=1/1.5)
  sigma <- par[35]
  sigma_vec <- rep(sigma , 14)

  out <- tryCatch(ode(y=init, times=times, func=seir, parms=paras, method="lsoda", rtol=1e-8, atol=1e-8), error=function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) return(1e30)
  out=as.data.frame(out)
  out<-out[,-1]
  estE<-as.matrix(out[,(n_group+1):(2*n_group)])
  lambda=epsilon*estE

  Iobs <- as.matrix(wave5df)
  lnL <- numeric(14)
  for(i in 1:14){
    lambda_i <- pmax(lambda[, i], 1e-10)
    obs_i <- Iobs[, i]
    sigma_i <- sigma_vec[i]

    lnL[i] <- sum(
      lgamma(obs_i + 1/sigma_i) - lgamma(1/sigma_i) +
        obs_i * log(lambda_i) + obs_i * log(sigma_i) -
        (obs_i + 1/sigma_i) * log(1 + lambda_i * sigma_i)
    )
  }

  obj <- -sum(lnL) + lambda_dev * pen_dev
  if (!is.finite(obj)) obj <- 1e30
  return(obj)
}
intipar<-c(initq5,initalpha5,sigma5)
model<-nloptr(
  x0 = intipar,
  eval_f = f,
  eval_g_ineq = nonlinear_constraint,
  opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-6, "maxeval" = 1000)
)
estpar<-model$solution
allpara<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],q8=estpar[8],q9=estpar[9],q10=estpar[10],
              q11=estpar[11],q12=estpar[12],q13=estpar[13],q14=estpar[14],q15=estpar[15],q16=estpar[16],q17=estpar[17],q18=estpar[18],q19=estpar[19],q20=estpar[20],
              alpha1=estpar[21],alpha2=estpar[22],alpha3=estpar[23],alpha4=estpar[24],alpha5=estpar[25],alpha6=estpar[26],alpha7=estpar[27],
              alpha8=estpar[28],alpha9=estpar[29],alpha10=estpar[30],alpha11=estpar[31],alpha12=estpar[32],alpha13=estpar[33],alpha14=estpar[34],
              epsilon=1/4.5,gamma=1/1.5)
out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
out=as.data.frame(out)
out<-out[,-1]
estdailyI<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
alpha<-estpar[21:34]
q_vec<-estpar[1:20]
result <- calculate_R0(alpha, q_vec)
R0_values <- result$R0_values
#
q5<-q_vec
alpha5<-alpha
R05<-R0_values
estdailyI5<-estdailyI
estsigma5<-estpar[35]
##CI
hess_mat<-hessian(func=f,x=estpar)
cov_mat<-get_inv(hess_mat)
boot_pars <- mvrnorm(1000,mu=estpar,Sigma=cov_mat)
q_boot <- boot_pars[, 1:20]
alpha_boot <- boot_pars[, 21:34]
sigma_boot <- boot_pars[, 35]
R0_values<-matrix(0,nrow=1000,ncol=length(R0_values))
for(i in 1:1000){
  result <- calculate_R0(alpha_boot[i,], q_boot[i,1:ncol(R0_values)])
  R0_values[i,]<-result$R0_values
}
estdailyI_boot<-list()
for(i in 1:1000){
  estpar<-c(q_boot[i,],alpha_boot[i,])
  allpara<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],q8=estpar[8],q9=estpar[9],q10=estpar[10],
                q11=estpar[11],q12=estpar[12],q13=estpar[13],q14=estpar[14],q15=estpar[15],q16=estpar[16],q17=estpar[17],q18=estpar[18],q19=estpar[19],q20=estpar[20],
                alpha1=estpar[21],alpha2=estpar[22],alpha3=estpar[23],alpha4=estpar[24],alpha5=estpar[25],alpha6=estpar[26],alpha7=estpar[27],
                alpha8=estpar[28],alpha9=estpar[29],alpha10=estpar[30],alpha11=estpar[31],alpha12=estpar[32],alpha13=estpar[33],alpha14=estpar[34],
                epsilon=1/4.5,gamma=1/1.5)
  out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
  out=as.data.frame(out)
  out<-out[,-1]
  estdailyI_boot[[i]]<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
}
estdailyI_boot_nb <- vector("list", length(estdailyI_boot))
for(i in seq_along(estdailyI_boot)) {
  lambda_mat <- estdailyI_boot[[i]]
  sigma_vec <- rep(sigma_boot[i],14)
  obs_mat <- matrix(NA, nrow = nrow(lambda_mat), ncol = ncol(lambda_mat))
  for(g in 1:ncol(lambda_mat)) {
    mu_vec <- lambda_mat[, g]
    size_param <- 1 / sigma_vec[g]
    mu_vec <- pmax(mu_vec, 1e-6)
    size_param <- max(size_param, 1e-6)
    obs_mat[, g] <- rnbinom(length(mu_vec), size = size_param, mu = mu_vec)
  }
  estdailyI_boot_nb[[i]] <- obs_mat
}
array_boot <- array(unlist(estdailyI_boot_nb),
                    dim = c(nrow(estdailyI_boot_nb[[1]]),
                            ncol(estdailyI_boot_nb[[1]]),
                            length(estdailyI_boot_nb)))
lower_CI <- matrix(NA, nrow = dim(array_boot)[1], ncol = dim(array_boot)[2])
upper_CI <- matrix(NA, nrow = dim(array_boot)[1], ncol = dim(array_boot)[2])
for (t in 1:dim(array_boot)[1]) {
  for (g in 1:dim(array_boot)[2]) {
    vals <- array_boot[t, g, ]
    lower_CI[t, g] <- quantile(vals, 0.025)
    upper_CI[t, g] <- quantile(vals, 0.975)
  }
}
#
q_boot5<-q_boot
alpha_boot5<-alpha_boot
sigma_boot5<-sigma_boot
q_CI5 <- apply(q_boot, 2, quantile, probs = c(0.025, 0.975))
alpha_CI5 <- apply(alpha_boot, 2, quantile, probs = c(0.025, 0.975))
R0_CI5 <- apply(R0_values, 2, quantile, probs = c(0.025, 0.975))
estsigma_CI5 <- quantile(sigma_boot, probs = c(0.025, 0.975))
estdailyI5_lower_CI<-lower_CI
estdailyI5_upper_CI<-upper_CI

###### Aggregation

### output
q_boot<-cbind(q_boot1,q_boot2,q_boot3,q_boot4,q_boot5)
alpha_boot<-cbind(alpha_boot1,alpha_boot2,alpha_boot3,alpha_boot4,alpha_boot5)
sigma_boot<-cbind(sigma_boot1,sigma_boot2,sigma_boot3,sigma_boot4,sigma_boot5)
q<-c(q1,q2,q3,q4,q5)
qCI<-cbind(q_CI1,q_CI2,q_CI3,q_CI4,q_CI5)
alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5)
alphaCI<-cbind(alpha_CI1,alpha_CI2,alpha_CI3,alpha_CI4,alpha_CI5)
estsigma<-c(estsigma1,estsigma2,estsigma3,estsigma4,estsigma5)
estsigma_CI<-cbind(estsigma_CI1,estsigma_CI2,estsigma_CI3,estsigma_CI4,estsigma_CI5)
R0<-c(R01,R02,R03,R04,R05)
R0CI<-cbind(R0_CI1,R0_CI2,R0_CI3,R0_CI4,R0_CI5)

estdailyIall<-data.frame(maleu20=c(estdailyI1[,1],estdailyI2[,1],estdailyI3[,1],estdailyI4[,1],estdailyI5[,1]),
                         male20=c(estdailyI1[,2],estdailyI2[,2],estdailyI3[,2],estdailyI4[,2],estdailyI5[,2]),
                         male30=c(estdailyI1[,3],estdailyI2[,3],estdailyI3[,3],estdailyI4[,3],estdailyI5[,3]),
                         male40=c(estdailyI1[,4],estdailyI2[,4],estdailyI3[,4],estdailyI4[,4],estdailyI5[,4]),
                         male50=c(estdailyI1[,5],estdailyI2[,5],estdailyI3[,5],estdailyI4[,5],estdailyI5[,5]),
                         male60=c(estdailyI1[,6],estdailyI2[,6],estdailyI3[,6],estdailyI4[,6],estdailyI5[,6]),
                         malea70=c(estdailyI1[,7],estdailyI2[,7],estdailyI3[,7],estdailyI4[,7],estdailyI5[,7]),
                         femaleu20=c(estdailyI1[,8],estdailyI2[,8],estdailyI3[,8],estdailyI4[,8],estdailyI5[,8]),
                         female20=c(estdailyI1[,9],estdailyI2[,9],estdailyI3[,9],estdailyI4[,9],estdailyI5[,9]),
                         female30=c(estdailyI1[,10],estdailyI2[,10],estdailyI3[,10],estdailyI4[,10],estdailyI5[,10]),
                         female40=c(estdailyI1[,11],estdailyI2[,11],estdailyI3[,11],estdailyI4[,11],estdailyI5[,11]),
                         female50=c(estdailyI1[,12],estdailyI2[,12],estdailyI3[,12],estdailyI4[,12],estdailyI5[,12]),
                         female60=c(estdailyI1[,13],estdailyI2[,13],estdailyI3[,13],estdailyI4[,13],estdailyI5[,13]),
                         femalea70=c(estdailyI1[,14],estdailyI2[,14],estdailyI3[,14],estdailyI4[,14],estdailyI5[,14]))
estdailyIall_lower<-data.frame(maleu20=c(estdailyI1_lower_CI[,1],estdailyI2_lower_CI[,1],estdailyI3_lower_CI[,1],estdailyI4_lower_CI[,1],estdailyI5_lower_CI[,1]),
                               male20=c(estdailyI1_lower_CI[,2],estdailyI2_lower_CI[,2],estdailyI3_lower_CI[,2],estdailyI4_lower_CI[,2],estdailyI5_lower_CI[,2]),
                               male30=c(estdailyI1_lower_CI[,3],estdailyI2_lower_CI[,3],estdailyI3_lower_CI[,3],estdailyI4_lower_CI[,3],estdailyI5_lower_CI[,3]),
                               male40=c(estdailyI1_lower_CI[,4],estdailyI2_lower_CI[,4],estdailyI3_lower_CI[,4],estdailyI4_lower_CI[,4],estdailyI5_lower_CI[,4]),
                               male50=c(estdailyI1_lower_CI[,5],estdailyI2_lower_CI[,5],estdailyI3_lower_CI[,5],estdailyI4_lower_CI[,5],estdailyI5_lower_CI[,5]),
                               male60=c(estdailyI1_lower_CI[,6],estdailyI2_lower_CI[,6],estdailyI3_lower_CI[,6],estdailyI4_lower_CI[,6],estdailyI5_lower_CI[,6]),
                               malea70=c(estdailyI1_lower_CI[,7],estdailyI2_lower_CI[,7],estdailyI3_lower_CI[,7],estdailyI4_lower_CI[,7],estdailyI5_lower_CI[,7]),
                               femaleu20=c(estdailyI1_lower_CI[,8],estdailyI2_lower_CI[,8],estdailyI3_lower_CI[,8],estdailyI4_lower_CI[,8],estdailyI5_lower_CI[,8]),
                               female20=c(estdailyI1_lower_CI[,9],estdailyI2_lower_CI[,9],estdailyI3_lower_CI[,9],estdailyI4_lower_CI[,9],estdailyI5_lower_CI[,9]),
                               female30=c(estdailyI1_lower_CI[,10],estdailyI2_lower_CI[,10],estdailyI3_lower_CI[,10],estdailyI4_lower_CI[,10],estdailyI5_lower_CI[,10]),
                               female40=c(estdailyI1_lower_CI[,11],estdailyI2_lower_CI[,11],estdailyI3_lower_CI[,11],estdailyI4_lower_CI[,11],estdailyI5_lower_CI[,11]),
                               female50=c(estdailyI1_lower_CI[,12],estdailyI2_lower_CI[,12],estdailyI3_lower_CI[,12],estdailyI4_lower_CI[,12],estdailyI5_lower_CI[,12]),
                               female60=c(estdailyI1_lower_CI[,13],estdailyI2_lower_CI[,13],estdailyI3_lower_CI[,13],estdailyI4_lower_CI[,13],estdailyI5_lower_CI[,13]),
                               femalea70=c(estdailyI1_lower_CI[,14],estdailyI2_lower_CI[,14],estdailyI3_lower_CI[,14],estdailyI4_lower_CI[,14],estdailyI5_lower_CI[,14]))
estdailyIall_upper<-data.frame(maleu20=c(estdailyI1_upper_CI[,1],estdailyI2_upper_CI[,1],estdailyI3_upper_CI[,1],estdailyI4_upper_CI[,1],estdailyI5_upper_CI[,1]),
                               male20=c(estdailyI1_upper_CI[,2],estdailyI2_upper_CI[,2],estdailyI3_upper_CI[,2],estdailyI4_upper_CI[,2],estdailyI5_upper_CI[,2]),
                               male30=c(estdailyI1_upper_CI[,3],estdailyI2_upper_CI[,3],estdailyI3_upper_CI[,3],estdailyI4_upper_CI[,3],estdailyI5_upper_CI[,3]),
                               male40=c(estdailyI1_upper_CI[,4],estdailyI2_upper_CI[,4],estdailyI3_upper_CI[,4],estdailyI4_upper_CI[,4],estdailyI5_upper_CI[,4]),
                               male50=c(estdailyI1_upper_CI[,5],estdailyI2_upper_CI[,5],estdailyI3_upper_CI[,5],estdailyI4_upper_CI[,5],estdailyI5_upper_CI[,5]),
                               male60=c(estdailyI1_upper_CI[,6],estdailyI2_upper_CI[,6],estdailyI3_upper_CI[,6],estdailyI4_upper_CI[,6],estdailyI5_upper_CI[,6]),
                               malea70=c(estdailyI1_upper_CI[,7],estdailyI2_upper_CI[,7],estdailyI3_upper_CI[,7],estdailyI4_upper_CI[,7],estdailyI5_upper_CI[,7]),
                               femaleu20=c(estdailyI1_upper_CI[,8],estdailyI2_upper_CI[,8],estdailyI3_upper_CI[,8],estdailyI4_upper_CI[,8],estdailyI5_upper_CI[,8]),
                               female20=c(estdailyI1_upper_CI[,9],estdailyI2_upper_CI[,9],estdailyI3_upper_CI[,9],estdailyI4_upper_CI[,9],estdailyI5_upper_CI[,9]),
                               female30=c(estdailyI1_upper_CI[,10],estdailyI2_upper_CI[,10],estdailyI3_upper_CI[,10],estdailyI4_upper_CI[,10],estdailyI5_upper_CI[,10]),
                               female40=c(estdailyI1_upper_CI[,11],estdailyI2_upper_CI[,11],estdailyI3_upper_CI[,11],estdailyI4_upper_CI[,11],estdailyI5_upper_CI[,11]),
                               female50=c(estdailyI1_upper_CI[,12],estdailyI2_upper_CI[,12],estdailyI3_upper_CI[,12],estdailyI4_upper_CI[,12],estdailyI5_upper_CI[,12]),
                               female60=c(estdailyI1_upper_CI[,13],estdailyI2_upper_CI[,13],estdailyI3_upper_CI[,13],estdailyI4_upper_CI[,13],estdailyI5_upper_CI[,13]),
                               femalea70=c(estdailyI1_upper_CI[,14],estdailyI2_upper_CI[,14],estdailyI3_upper_CI[,14],estdailyI4_upper_CI[,14],estdailyI5_upper_CI[,14]))

### plot
varnames <- c("Male <20", "Male 20-29", "Male 30-39", "Male 40-49", "Male 50-59", "Male 60-69", "Male ≥70",
              "Female <20", "Female 20-29", "Female 30-39", "Female 40-49", "Female 50-59", "Female 60-69", "Female ≥70")
colnames(estdailyIall)<-colnames(estdailyIall_lower)<-colnames(estdailyIall_upper)<-colnames(df)[-1]<-varnames
ymin_global <- min(c(unlist(estdailyIall), unlist(estdailyIall_lower)), na.rm = TRUE)
ymax_global <- max(c(unlist(estdailyIall), unlist(estdailyIall_upper)), na.rm = TRUE)
xmax_global <- 700
create_unified_theme <- function() {
  theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(size = 18),
      axis.title = element_blank(),
      plot.title = element_text(size = 24, hjust = 0.5, face = "bold", margin = margin(b = 5)),
      plot.margin = margin(5, 5, 5, 5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank()
    )
}
plot_single <- function(true_y, est_y, lower_y, upper_y, title) {
  df_plot <- data.frame(
    day = 1:xmax_global,
    true = true_y,
    est = est_y,
    lower = lower_y,
    upper = upper_y
  )
  
  ggplot(df_plot, aes(x = day)) +
    geom_point(aes(y = true), color = "black", size = 0.8, alpha = 0.6) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#CC0000", alpha = 0.3, color="red", linetype="dotted") +
    geom_line(aes(y = est), color = "red", linewidth = 0.6) +
    labs(title = title) +
    scale_y_continuous(limits = c(ymin_global, ymax_global)) +
    scale_x_continuous(limits = c(1, xmax_global),      
                       breaks = c(1,100,300,500,700),
                       labels = c("1","100","300","500","700")) +
    create_unified_theme()
}
plot_list <- lapply(varnames, function(var) {
  plot_single(df[[var]],
              estdailyIall[[var]],
              estdailyIall_lower[[var]],
              estdailyIall_upper[[var]],
              title = var)
})
plot_list <- lapply(plot_list, function(p) p + create_unified_theme())
big_plot <- wrap_plots(plot_list, nrow = 2, ncol = 7) &
  theme(plot.margin = margin(t = 10, r = 10, b = 15, l = 18))
final_plot_labeled <- ggdraw(big_plot, clip = "off") +
  draw_label("Days from February 1, 2020 to December 31, 2021", x = 0.5, y = -0.005, vjust = -0.5, size = 28) +
  draw_label("COVID-19 cases", x = 0.006, y = 0.75, angle = 90, hjust = 1.2, size = 28)
# ggsave("results/figures/Figure_SEIR_fitted_vs_observed.jpg",final_plot_labeled,width=28,height=8,dpi=300)

###### rate of change
dates  <- seq(as.Date("2020-02-01"), as.Date("2021-12-31"), by = "day")
months <- format(dates, "%Y-%m")
estdailyI_monthly <- aggregate(
  estdailyIall,
  by = list(month = months),
  FUN = sum
)
estdailyI_monthly <- estdailyI_monthly[, -1]
estdailyI_monthly <- rbind(
  rep(0, 14),
  rep(1, 14),
  estdailyI_monthly
)
estdailyI_monthly_all <- rowSums(estdailyI_monthly)
rate_of_change <- diff(estdailyI_monthly_all) /
  (head(estdailyI_monthly_all, -1) + 1e-1)
rate_of_change <- c(0, 0, rate_of_change)

###### control effort pt
alpha_mat<-rbind(matrix(rep(alpha[1:n_group],(nrow(wave1df)/10)-1),nrow=(nrow(wave1df)/10)-1,byrow=TRUE),
                 matrix(rep(alpha[(n_group+1):(2*n_group)],(nrow(wave2df)/10)-1),nrow=(nrow(wave2df)/10)-1,byrow=TRUE),
                 matrix(rep(alpha[(2*n_group+1):(3*n_group)],(nrow(wave3df)/10)-1),nrow=(nrow(wave3df)/10)-1,byrow=TRUE),
                 matrix(rep(alpha[(3*n_group+1):(4*n_group)],(nrow(wave4df)/10)-3),nrow=(nrow(wave4df)/10)-3,byrow=TRUE),
                 matrix(rep(alpha[(4*n_group+1):(5*n_group)],nrow(wave5df)/10),nrow=nrow(wave5df)/10,byrow=TRUE))
Rt<-matrix(0,nrow=length(q),ncol=n_group)
for(i in 1:ncol(Rt)){
  for(j in 1:nrow(Rt)){
    Rt[j,i]<-(4.5+3*1.5)*q[j]*sum(alpha_mat[j,]*C[,i]*popN_mat[,i])
  }
}
Rt_hat <- Rt 
Rt_max <- apply(Rt, 2, max)
pt<-1-sweep(Rt, 2, Rt_max, "/")
pt_hat <- pt
pt_expanded <- kronecker(pt, matrix(1, 10, 1))
pt_monthly_means <- aggregate(pt_expanded, by = list(month = months), FUN = mean)
pt_monthly_means<-rbind(rep(0,14),pt_monthly_means[,-1])

### Counterfactual analysis under scaled transmission intensity with bootstrap uncertainty
q_boot_1<-q_boot[,1:15]
q_boot_2<-q_boot[,(15+1):(15+12)]
q_boot_3<-q_boot[,(15+12+1):(15+12+12)]
q_boot_4<-q_boot[,(15+12+12+1):(15+12+12+11)]
q_boot_5<-q_boot[,(15+12+12+11+1):(15+12+12+11+20)]
alpha_boot_1<-alpha_boot[,1:14]
alpha_boot_2<-alpha_boot[,(14+1):(14+14)]
alpha_boot_3<-alpha_boot[,(14+14+1):(14+14+14)]
alpha_boot_4<-alpha_boot[,(14+14+14+1):(14+14+14+14)]
alpha_boot_5<-alpha_boot[,(14+14+14+14+1):(14+14+14+14+14)]
# additional effort
Rt_new<-list()
for(i in 1:length(40:100)){
  l<-c(rep((i+39)/100,15),rep((i+39)/100,12),rep((i+39)/100,12),rep((i+39)/100,11),rep((i+39)/100,20))
  Rt_new[[i]]<-sweep(Rt, 1, l, "*")
}
pt_new<-list()
for(i in 1:length(40:100)){
  Rtl<-Rt_new[[i]]
  pt_new[[i]]<-1-sweep(Rtl,2,Rt_max,"/")
}
pt_monthly_new<-list()
for(i in 1:length(40:100)){
  ptl<-pt_new[[i]]
  ptl_expanded <- kronecker(ptl, matrix(1, 10, 1))
  ptl_monthly_means <- aggregate(ptl_expanded, by = list(month = months), FUN = mean)
  pt_monthly_new[[i]]<-rbind(rep(0,14),ptl_monthly_means[,-1])
}
# bootstrap uncertainty
n_boot <- nrow(q_boot)
l_grid <- 40:100 / 100
Rt_boot_list <- vector("list", n_boot)
Rt_new_boot_list <- vector("list", n_boot)
pt_new_boot_list <- vector("list", n_boot)
pt_monthly_boot_list <- vector("list", n_boot)
for (b in 1:n_boot) {
  q_all <- c(as.numeric(q_boot_1[b,]), as.numeric(q_boot_2[b,]), as.numeric(q_boot_3[b,]), as.numeric(q_boot_4[b,]), as.numeric(q_boot_5[b,]))
  alpha_all <- c(as.numeric(alpha_boot_1[b,]), as.numeric(alpha_boot_2[b,]), as.numeric(alpha_boot_3[b,]), as.numeric(alpha_boot_4[b,]), as.numeric(alpha_boot_5[b,]))

  alpha_mat <- rbind(
    matrix(rep(alpha_all[1:14], nrow(wave1df) %/% 10 - 1), ncol=14, byrow=TRUE),
    matrix(rep(alpha_all[15:28], nrow(wave2df) %/% 10 - 1), ncol=14, byrow=TRUE),
    matrix(rep(alpha_all[29:42], nrow(wave3df) %/% 10 - 1), ncol=14, byrow=TRUE),
    matrix(rep(alpha_all[43:56], nrow(wave4df) %/% 10 - 3), ncol=14, byrow=TRUE),
    matrix(rep(alpha_all[57:70], nrow(wave5df) %/% 10), ncol=14, byrow=TRUE)
  )

  Rt <- matrix(0, nrow=length(q_all), ncol=n_group)
  for (i in 1:ncol(Rt)) {
    for (j in 1:nrow(Rt)) {
      Rt[j,i] <- (4.5+3*1.5) * q_all[j] * sum(alpha_mat[j,] * C[,i] * popN_mat[,i])
    }
  }

  Rt_boot_list[[b]] <- Rt

  Rt_new_list <- vector("list", length(l_grid))
  pt_new_list <- vector("list", length(l_grid))

  for (k in 1:length(l_grid)) {
    l <- c(rep(l_grid[k], 15), rep(l_grid[k], 12), rep(l_grid[k], 12), rep(l_grid[k], 11), rep(l_grid[k], 20))
    Rt_new <- sweep(Rt, 1, l, "*")
    pt <- 1 - sweep(Rt_new, 2, Rt_max, "/")

    Rt_new_list[[k]] <- Rt_new
    pt_new_list[[k]] <- pt
  }

  Rt_new_boot_list[[b]] <- Rt_new_list
  pt_new_boot_list[[b]] <- pt_new_list

  pt_monthly_list <- vector("list", length(l_grid))
  for (k in 1:length(l_grid)) {
    ptl <- pt_new_list[[k]]
    ptl_expanded <- ptl[rep(1:nrow(ptl), each=10), ]
    ptl_monthly <- aggregate(ptl_expanded, by = list(month=months), FUN=mean)
    ptl_monthly_matrix <- as.matrix(ptl_monthly[,-1])
    pt_monthly_list[[k]] <- ptl_monthly_matrix
  }
  pt_monthly_boot_list[[b]] <- pt_monthly_list
}
Rt_array <- array(NA, dim=c(nrow(Rt), ncol(Rt), n_boot))
for (b in 1:n_boot) Rt_array[,,b] <- Rt_boot_list[[b]]
Rt_lower <- apply(Rt_array, c(1,2), quantile, 0.025, na.rm = TRUE)
Rt_upper <- apply(Rt_array, c(1,2), quantile, 0.975, na.rm = TRUE)
pt_lower<-1-sweep(Rt_upper, 2, Rt_max, "/")
pt_upper<-1-sweep(Rt_lower, 2, Rt_max, "/")
Rt_new_CI <- vector("list", length(l_grid))
pt_new_CI <- vector("list", length(l_grid))
pt_monthly_CI <- vector("list", length(l_grid))
for (k in 1:length(l_grid)) {

  Rt_new_array <- array(NA, dim=c(nrow(Rt), ncol(Rt), n_boot))
  for (b in 1:n_boot) Rt_new_array[,,b] <- Rt_new_boot_list[[b]][[k]]
  Rt_new_lower <- apply(Rt_new_array, c(1,2), quantile, 0.025)
  Rt_new_upper <- apply(Rt_new_array, c(1,2), quantile, 0.975)
  Rt_new_CI[[k]] <- list(lower=Rt_new_lower, upper=Rt_new_upper)

  pt_new_array <- array(NA, dim=c(nrow(Rt), ncol(Rt), n_boot))
  for (b in 1:n_boot) pt_new_array[,,b] <- pt_new_boot_list[[b]][[k]]
  pt_new_lower <- apply(pt_new_array, c(1,2), quantile, 0.025)
  pt_new_upper <- apply(pt_new_array, c(1,2), quantile, 0.975)
  pt_new_CI[[k]] <- list(lower=pt_new_lower, upper=pt_new_upper)

  ptl_monthly_first <- pt_monthly_boot_list[[1]][[k]]
  monthly_rows <- nrow(ptl_monthly_first)
  monthly_cols <- ncol(ptl_monthly_first)
  pt_monthly_array <- array(NA, dim=c(monthly_rows, monthly_cols, n_boot))
  for (b in 1:n_boot) pt_monthly_array[,,b] <- pt_monthly_boot_list[[b]][[k]]
  pt_monthly_lower <- apply(pt_monthly_array, c(1,2), quantile, 0.025)
  pt_monthly_upper <- apply(pt_monthly_array, c(1,2), quantile, 0.975)
  pt_monthly_CI[[k]] <- list(lower=pt_monthly_lower, upper=pt_monthly_upper)
}
### plot
Rt <- Rt_hat
Rt_df <- as.data.frame(Rt)
colnames(Rt_df)[1:14] <- varnames
Rt_df$time <- seq_len(nrow(Rt_df))
Rt_lower <- as.data.frame(Rt_lower)
Rt_upper <- as.data.frame(Rt_upper)
colnames(Rt_lower)[1:14] <- varnames
colnames(Rt_upper)[1:14] <- varnames
Rt_lower$time <- seq_len(nrow(Rt_lower))
Rt_upper$time <- seq_len(nrow(Rt_upper))
Rt_long <- Rt_df %>%
  pivot_longer(cols = -time, names_to = "group", values_to = "Rt") %>%
  left_join(
    Rt_lower %>% pivot_longer(cols = -time, names_to = "group", values_to = "Rt_lower"),
    by = c("time", "group")
  ) %>%
  left_join(
    Rt_upper %>% pivot_longer(cols = -time, names_to = "group", values_to = "Rt_upper"),
    by = c("time", "group")
  )
pt <- pt_hat
pt_df <- as.data.frame(pt)
colnames(pt_df)[1:14] <- varnames
pt_df$time <- seq_len(nrow(pt_df))
pt_lower <- as.data.frame(pt_lower)
pt_upper <- as.data.frame(pt_upper)
colnames(pt_lower)[1:14] <- varnames
colnames(pt_upper)[1:14] <- varnames
pt_lower$time <- seq_len(nrow(pt_lower))
pt_upper$time <- seq_len(nrow(pt_upper))
pt_long <- pt_df %>%
  pivot_longer(cols = -time, names_to = "group", values_to = "pt") %>%
  left_join(
    pt_lower %>% pivot_longer(cols = -time, names_to = "group", values_to = "pt_lower"),
    by = c("time", "group")
  ) %>%
  left_join(
    pt_upper %>% pivot_longer(cols = -time, names_to = "group", values_to = "pt_upper"),
    by = c("time", "group")
  )

combined_long <- Rt_long %>%
  left_join(pt_long, by = c("time", "group")) %>%
  mutate(group = factor(group, levels = varnames)) %>%
  arrange(group, time)

create_unified_theme <- function() {
  theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text  = element_text(size = 15),
      axis.title = element_blank(),
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      plot.margin = margin(5, 5, 5, 5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

ymin_global_Rt <- 0
ymax_global_Rt <- max(combined_long$Rt_upper, na.rm = TRUE)
y_breaks <- seq(0, ceiling(ymax_global_Rt), by = 1)

plot_single_rtpt <- function(df_sub, title, xmax = 70) {
  
  df_sub <- df_sub %>%
    mutate(across(c(pt_lower, pt, pt_upper), ~ pmax(., 0)))
  
  y_max_rt <- max(df_sub$Rt_upper, na.rm = TRUE)
  y_breaks <- seq(0, ceiling(y_max_rt), by = 1)
  
  pt_min <- min(c(df_sub$pt_lower, df_sub$pt, df_sub$pt_upper), na.rm = TRUE)
  pt_max <- max(c(df_sub$pt_lower, df_sub$pt, df_sub$pt_upper), na.rm = TRUE)
  denom  <- ifelse(pt_max > pt_min, (pt_max - pt_min), 1)
  
  df_sub <- df_sub %>%
    mutate(across(c(pt_lower, pt, pt_upper), ~ pmin(pmax(., 0), 1))) %>%  
    mutate(
      pt_scaled       = pt       * ymax_global_Rt,
      pt_lower_scaled = pt_lower * ymax_global_Rt,
      pt_upper_scaled = pt_upper * ymax_global_Rt
    )
  
  ggplot(df_sub, aes(x = time)) +

    geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), fill = "#CC0000", alpha = 0.3) +
    geom_line(aes(y = Rt), color = "red", linewidth = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
    
    geom_ribbon(aes(ymin = pt_lower_scaled, ymax = pt_upper_scaled), fill = "grey10", alpha = 0.25) +
    geom_line(aes(y = pt_scaled), color = "black", linewidth = 0.5) +
    
    labs(title = title) +
    scale_x_continuous(
      limits = c(1, xmax),
      breaks = c(1, 10, 30, 50, 70),
      labels = c("1","100","300","500","700")
    ) +
    scale_y_continuous(
      limits = c(0, ymax_global_Rt),
      breaks = seq(0, ceiling(ymax_global_Rt), by = 1),
      sec.axis = sec_axis(
        transform = ~ . / ymax_global_Rt,
        breaks = seq(0, 1, 0.2),
        labels = scales::number_format(accuracy = 0.1)
      )
    ) +
    create_unified_theme()
}

plot_list_rtpt <- lapply(varnames, function(g) {
  df_sub <- filter(combined_long, group == g)
  plot_single_rtpt(df_sub, title = g, xmax = 70)
})

rtpt_plot_grid <- cowplot::plot_grid(
  plotlist = plot_list_rtpt,
  nrow = 2, ncol = 7,
  align = "hv", axis = "tblr"
)

rtpt_final_plot <- ggdraw(rtpt_plot_grid, xlim = c(-0.02, 1.02), ylim = c(-0.04, 1)) +
  draw_label("Days from February 1, 2020 to December 31, 2021",
             x = 0.5, y = -0.01, size = 24, hjust = 0.5) +
  draw_label("Effective reproduction number ℛ(t)",
             x = -0.01, y = 0.5, angle = 90, size = 24, hjust = 0.5) +
  draw_label("Required control effort 𝑝(𝑡)",
             x = 1.01, y = 0.5, angle = -90, size = 24, hjust = 0.5)
# ggsave("results/figures/Figure_SEIR_rtpt.jpg",rtpt_final_plot,width=28,height=8,dpi=300)

############################################################
## 10. Senario analysis
############################################################
q_1<-q[1:15]
q_2<-q[(15+1):(15+12)]
q_3<-q[(15+12+1):(15+12+12)]
q_4<-q[(15+12+12+1):(15+12+12+11)]
q_5<-q[(15+12+12+11+1):(15+12+12+11+20)]
q_boot_1<-q_boot[,1:15]
q_boot_2<-q_boot[,(15+1):(15+12)]
q_boot_3<-q_boot[,(15+12+1):(15+12+12)]
q_boot_4<-q_boot[,(15+12+12+1):(15+12+12+11)]
q_boot_5<-q_boot[,(15+12+12+11+1):(15+12+12+11+20)]
alpha_1<-alpha[1:14]
alpha_2<-alpha[(14+1):(14+14)]
alpha_3<-alpha[(14+14+1):(14+14+14)]
alpha_4<-alpha[(14+14+14+1):(14+14+14+14)]
alpha_5<-alpha[(14+14+14+14+1):(14+14+14+14+14)]
alpha_boot_1<-alpha_boot[,1:14]
alpha_boot_2<-alpha_boot[,(14+1):(14+14)]
alpha_boot_3<-alpha_boot[,(14+14+1):(14+14+14)]
alpha_boot_4<-alpha_boot[,(14+14+14+1):(14+14+14+14)]
alpha_boot_5<-alpha_boot[,(14+14+14+14+1):(14+14+14+14+14)]
sigma_boot_1<-sigma_boot[,1]
sigma_boot_2<-sigma_boot[,2]
sigma_boot_3<-sigma_boot[,3]
sigma_boot_4<-sigma_boot[,4]
sigma_boot_5<-sigma_boot[,5]

### wave1
seir<-function(time,y,params){

  q <- if (time <= 10) params$q1 else
    if (time <= 20) params$q2 else
      if (time <= 30) params$q3 else
        if (time <= 40) params$q4 else
          if (time <= 50) params$q5 else
            if (time <= 60) params$q6 else
              if (time <= 70) params$q7 else
                if (time <= 80) params$q8 else
                  if (time <= 90) params$q9 else
                    if (time <= 100) params$q10 else
                      if (time <= 110) params$q11 else
                        if (time <= 120) params$q12 else
                          if (time <= 130) params$q13 else
                            if (time <= 140) params$q14 else
                              params$q15
  beta_mat <- q * C

  epsilon=params[["epsilon"]]
  gamma=params[["gamma"]]
  alpha1<-params[["alpha1"]]
  alpha2<-params[["alpha2"]]
  alpha3<-params[["alpha3"]]
  alpha4<-params[["alpha4"]]
  alpha5<-params[["alpha5"]]
  alpha6<-params[["alpha6"]]
  alpha7<-params[["alpha7"]]
  alpha8<-params[["alpha8"]]
  alpha9<-params[["alpha9"]]
  alpha10<-params[["alpha10"]]
  alpha11<-params[["alpha11"]]
  alpha12<-params[["alpha12"]]
  alpha13<-params[["alpha13"]]
  alpha14<-params[["alpha14"]]
  l=params[["l"]]

  alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14)

  S<-y[1:n_group]
  E<-y[(n_group+1):(2*n_group)]
  I<-y[(2*n_group+1):(3*n_group)]
  R<-y[(3*n_group+1):(4*n_group)]

  dS<-numeric(n_group)
  dE<-numeric(n_group)
  dI<-numeric(n_group)
  dR<-numeric(n_group)

  for(i in 1:14){
    infection_force <- l*sum(beta_mat[i,] * S[i] * alpha[i] * (3*I + E)/popN)
    dS[i] <- -infection_force
    dE[i] <- infection_force - epsilon * E[i]
    dI[i] <- epsilon * E[i] - gamma * I[i]
    dR[i] <- gamma * I[i]
  }

  solution=c(dS,dE,dI,dR)
  return(list(solution))
}
init=c(population_group$male,population_group$female,rep(0.01,n_group),rep(0,n_group),rep(0,n_group))
times=seq(1,150,by=1)
estdailyI1_new<-vector("list", length(40:100))
estdailyI1_lower_CI<-vector("list", length(40:100))
estdailyI1_upper_CI<-vector("list", length(40:100))
initnew1<-matrix(0,nrow=length(40:100),ncol=56)
for(i in 40:100){
  allpara<-list(q1=q_1[1],q2=q_1[2],q3=q_1[3],q4=q_1[4],q5=q_1[5],q6=q_1[6],q7=q_1[7],
                q8=q_1[8],q9=q_1[9],q10=q_1[10],q11=q_1[11],q12=q_1[12],q13=q_1[13],q14=q_1[14],q15=q_1[15],
                epsilon=1/4.5,gamma=1/1.5,l=i/100,
                alpha1=alpha_1[1],alpha2=alpha_1[2],alpha3=alpha_1[3],alpha4=alpha_1[4],alpha5=alpha_1[5],alpha6=alpha_1[6],alpha7=alpha_1[7],
                alpha8=alpha_1[8],alpha9=alpha_1[9],alpha10=alpha_1[10],alpha11=alpha_1[11],alpha12=alpha_1[12],alpha13=alpha_1[13],alpha14=alpha_1[14])
  out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
  out=as.data.frame(out)
  out<-out[,-1]
  estdailyI1_new[[i-39]]<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
  initnew1[i-39,]<-as.numeric(out[150,])
  ##CI
  q_boot <- q_boot_1
  alpha_boot <- alpha_boot_1
  sigma_hat <- sigma_boot_1
  estdailyI_boot <- vector("list", 1000)
  for (j in 1:1000) {
    estpar <- c(as.numeric(q_boot[j,]), as.numeric(alpha_boot[j,]))
    allpara_boot <- list(
      q1=estpar[1], q2=estpar[2], q3=estpar[3], q4=estpar[4], q5=estpar[5],
      q6=estpar[6], q7=estpar[7], q8=estpar[8], q9=estpar[9], q10=estpar[10],
      q11=estpar[11], q12=estpar[12], q13=estpar[13], q14=estpar[14], q15=estpar[15],
      l=i/100,
      alpha1=estpar[16], alpha2=estpar[17], alpha3=estpar[18], alpha4=estpar[19],
      alpha5=estpar[20], alpha6=estpar[21], alpha7=estpar[22], alpha8=estpar[23],
      alpha9=estpar[24], alpha10=estpar[25], alpha11=estpar[26], alpha12=estpar[27],
      alpha13=estpar[28], alpha14=estpar[29],
      epsilon=1/4.5, gamma=1/1.5
    )
    out_boot <- ode(y=init, times=times, func=seir, parms=allpara_boot, method="lsoda",rtol = 1e-8,atol = 1e-8)
    out_boot <- as.data.frame(out_boot)[, -1]
    estdailyI_boot[[j]] <- epsilon * as.matrix(out_boot[,(n_group+1):(2*n_group)])
  }
  estdailyI_boot_nb <- vector("list", length(estdailyI_boot))
  for (j in seq_along(estdailyI_boot)) {
    lambda_mat <- estdailyI_boot[[j]]
    obs_mat <- matrix(NA, nrow = nrow(lambda_mat), ncol = ncol(lambda_mat))
    for (g in 1:ncol(lambda_mat)) {
      for (t in 1:nrow(lambda_mat)) {
        lambda_val <- lambda_mat[t, g]
        lambda_val <- pmax(lambda_val, 1e-6)
        size <- 1 / sigma_hat[j]
        obs_mat[t, g] <- rnbinom(1, size = size, mu = lambda_val)
      }
    }
    estdailyI_boot_nb[[j]] <- obs_mat
  }
  array_boot_nb <- array(unlist(estdailyI_boot_nb),
                         dim = c(nrow(estdailyI_boot_nb[[1]]),
                                 ncol(estdailyI_boot_nb[[1]]),
                                 length(estdailyI_boot_nb)))
  lower_CI <- apply(array_boot_nb, c(1,2), quantile, probs=0.025)
  upper_CI <- apply(array_boot_nb, c(1,2), quantile, probs=0.975)

  estdailyI1_lower_CI[[i - 39]] <- lower_CI
  estdailyI1_upper_CI[[i - 39]] <- upper_CI
}

### wave2
seir<-function(time,y,params){

  q <- if (time <= 10) params$q1 else
    if (time <= 20) params$q2 else
      if (time <= 30) params$q3 else
        if (time <= 40) params$q4 else
          if (time <= 50) params$q5 else
            if (time <= 60) params$q6 else
              if (time <= 70) params$q7 else
                if (time <= 80) params$q8 else
                  if (time <= 90) params$q9 else
                    if (time <= 100) params$q10 else
                      if (time <= 110) params$q11 else
                        params$q12
  beta_mat <- q * C

  epsilon=params[["epsilon"]]
  gamma=params[["gamma"]]
  alpha1<-params[["alpha1"]]
  alpha2<-params[["alpha2"]]
  alpha3<-params[["alpha3"]]
  alpha4<-params[["alpha4"]]
  alpha5<-params[["alpha5"]]
  alpha6<-params[["alpha6"]]
  alpha7<-params[["alpha7"]]
  alpha8<-params[["alpha8"]]
  alpha9<-params[["alpha9"]]
  alpha10<-params[["alpha10"]]
  alpha11<-params[["alpha11"]]
  alpha12<-params[["alpha12"]]
  alpha13<-params[["alpha13"]]
  alpha14<-params[["alpha14"]]
  l=params[["l"]]

  alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14)

  S<-y[1:n_group]
  E<-y[(n_group+1):(2*n_group)]
  I<-y[(2*n_group+1):(3*n_group)]
  R<-y[(3*n_group+1):(4*n_group)]

  dS<-numeric(n_group)
  dE<-numeric(n_group)
  dI<-numeric(n_group)
  dR<-numeric(n_group)

  for(i in 1:14){
    infection_force <- l*sum(beta_mat[i,] * S[i] * alpha[i] * (3*I + E)/popN)
    dS[i] <- -infection_force
    dE[i] <- infection_force - epsilon * E[i]
    dI[i] <- epsilon * E[i] - gamma * I[i]
    dR[i] <- gamma * I[i]
  }

  solution=c(dS,dE,dI,dR)
  return(list(solution))
}
times=seq(1,120,by=1)
estdailyI2_new<-vector("list", length(40:100))
estdailyI2_lower_CI<-vector("list", length(40:100))
estdailyI2_upper_CI<-vector("list", length(40:100))
initnew2<-matrix(0,nrow=length(40:100),ncol=56)
for(i in 40:100){
  allpara<-list(q1=q_2[1],q2=q_2[2],q3=q_2[3],q4=q_2[4],q5=q_2[5],q6=q_2[6],q7=q_2[7],
                q8=q_2[8],q9=q_2[9],q10=q_2[10],q11=q_2[11],q12=q_2[12],
                epsilon=1/4.5,gamma=1/1.5,l=i/100,
                alpha1=alpha_2[1],alpha2=alpha_2[2],alpha3=alpha_2[3],alpha4=alpha_2[4],alpha5=alpha_2[5],alpha6=alpha_2[6],alpha7=alpha_2[7],
                alpha8=alpha_2[8],alpha9=alpha_2[9],alpha10=alpha_2[10],alpha11=alpha_2[11],alpha12=alpha_2[12],alpha13=alpha_2[13],alpha14=alpha_2[14])
  init=initnew1[i-39,]
  out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
  out=as.data.frame(out)
  out<-out[,-1]
  estdailyI2_new[[i-39]]<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
  initnew2[i-39,]<-as.numeric(out[120,])
  ##CI
  q_boot <- q_boot_2
  alpha_boot <- alpha_boot_2
  sigma_hat <- sigma_boot_2
  estdailyI_boot <- vector("list", 1000)
  for (j in 1:1000) {
    estpar <- c(as.numeric(q_boot[j,]), as.numeric(alpha_boot[j,]))
    allpara_boot<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
                       q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],q12=estpar[12],l=i/100,
                       alpha1=estpar[13],alpha2=estpar[14],alpha3=estpar[15],alpha4=estpar[16],alpha5=estpar[17],alpha6=estpar[18],alpha7=estpar[19],
                       alpha8=estpar[20],alpha9=estpar[21],alpha10=estpar[22],alpha11=estpar[23],alpha12=estpar[24],alpha13=estpar[25],alpha14=estpar[26],
                       epsilon=1/4.5,gamma=1/1.5)
    out_boot <- ode(y=init, times=times, func=seir, parms=allpara_boot, method="lsoda",rtol = 1e-8,atol = 1e-8)
    out_boot <- as.data.frame(out_boot)[, -1]
    estdailyI_boot[[j]] <- epsilon * as.matrix(out_boot[,(n_group+1):(2*n_group)])
  }
  estdailyI_boot_nb <- vector("list", length(estdailyI_boot))
  for (j in seq_along(estdailyI_boot)) {
    lambda_mat <- estdailyI_boot[[j]]
    obs_mat <- matrix(NA, nrow = nrow(lambda_mat), ncol = ncol(lambda_mat))
    for (g in 1:ncol(lambda_mat)) {
      for (t in 1:nrow(lambda_mat)) {
        lambda_val <- lambda_mat[t, g]
        lambda_val <- pmax(lambda_val, 1e-6)
        size <- 1 / sigma_hat[j]
        obs_mat[t, g] <- rnbinom(1, size = size, mu = lambda_val)
      }
    }
    estdailyI_boot_nb[[j]] <- obs_mat
  }
  array_boot_nb <- array(unlist(estdailyI_boot_nb),
                         dim = c(nrow(estdailyI_boot_nb[[1]]),
                                 ncol(estdailyI_boot_nb[[1]]),
                                 length(estdailyI_boot_nb)))
  lower_CI <- apply(array_boot_nb, c(1,2), quantile, probs=0.025)
  upper_CI <- apply(array_boot_nb, c(1,2), quantile, probs=0.975)

  estdailyI2_lower_CI[[i - 39]] <- lower_CI
  estdailyI2_upper_CI[[i - 39]] <- upper_CI
}

### wave3
seir<-function(time,y,params){

  q <- if (time <= 10) params$q1 else
    if (time <= 20) params$q2 else
      if (time <= 30) params$q3 else
        if (time <= 40) params$q4 else
          if (time <= 50) params$q5 else
            if (time <= 60) params$q6 else
              if (time <= 70) params$q7 else
                if (time <= 80) params$q8 else
                  if (time <= 90) params$q9 else
                    if (time <= 100) params$q10 else
                      if (time <= 110) params$q11 else
                        params$q12
  beta_mat <- q * C

  epsilon=params[["epsilon"]]
  gamma=params[["gamma"]]
  alpha1<-params[["alpha1"]]
  alpha2<-params[["alpha2"]]
  alpha3<-params[["alpha3"]]
  alpha4<-params[["alpha4"]]
  alpha5<-params[["alpha5"]]
  alpha6<-params[["alpha6"]]
  alpha7<-params[["alpha7"]]
  alpha8<-params[["alpha8"]]
  alpha9<-params[["alpha9"]]
  alpha10<-params[["alpha10"]]
  alpha11<-params[["alpha11"]]
  alpha12<-params[["alpha12"]]
  alpha13<-params[["alpha13"]]
  alpha14<-params[["alpha14"]]
  l=params[["l"]]

  alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14)

  S<-y[1:n_group]
  E<-y[(n_group+1):(2*n_group)]
  I<-y[(2*n_group+1):(3*n_group)]
  R<-y[(3*n_group+1):(4*n_group)]

  dS<-numeric(n_group)
  dE<-numeric(n_group)
  dI<-numeric(n_group)
  dR<-numeric(n_group)

  for(i in 1:14){
    infection_force <- l*sum(beta_mat[i,] * S[i] * alpha[i] * (3*I + E)/popN)
    dS[i] <- -infection_force
    dE[i] <- infection_force - epsilon * E[i]
    dI[i] <- epsilon * E[i] - gamma * I[i]
    dR[i] <- gamma * I[i]
  }

  solution=c(dS,dE,dI,dR)
  return(list(solution))
}
times=seq(1,120,by=1)
estdailyI3_new<-vector("list", length(40:100))
estdailyI3_lower_CI<-vector("list", length(40:100))
estdailyI3_upper_CI<-vector("list", length(40:100))
initnew3<-matrix(0,nrow=length(40:100),ncol=56)
for(i in 40:100){
  allpara<-list(q1=q_3[1],q2=q_3[2],q3=q_3[3],q4=q_3[4],q5=q_3[5],q6=q_3[6],q7=q_3[7],
                q8=q_3[8],q9=q_3[9],q10=q_3[10],q11=q_3[11],q12=q_3[12],
                epsilon=1/4.5,gamma=1/1.5,l=i/100,
                alpha1=alpha_3[1],alpha2=alpha_3[2],alpha3=alpha_3[3],alpha4=alpha_3[4],alpha5=alpha_3[5],alpha6=alpha_3[6],alpha7=alpha_3[7],
                alpha8=alpha_3[8],alpha9=alpha_3[9],alpha10=alpha_3[10],alpha11=alpha_3[11],alpha12=alpha_3[12],alpha13=alpha_3[13],alpha14=alpha_3[14])
  init=initnew2[i-39,]
  out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
  out=as.data.frame(out)
  out<-out[,-1]
  estdailyI3_new[[i-39]]<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
  initnew3[i-39,]<-as.numeric(out[120,])
  ##CI
  q_boot <- q_boot_3
  alpha_boot <- alpha_boot_3
  sigma_hat <- sigma_boot_3
  estdailyI_boot <- vector("list", 1000)
  for (j in 1:1000) {
    estpar <- c(as.numeric(q_boot[j,]), as.numeric(alpha_boot[j,]))
    allpara_boot<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
                       q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],q12=estpar[12],l=i/100,
                       alpha1=estpar[13],alpha2=estpar[14],alpha3=estpar[15],alpha4=estpar[16],alpha5=estpar[17],alpha6=estpar[18],alpha7=estpar[19],
                       alpha8=estpar[20],alpha9=estpar[21],alpha10=estpar[22],alpha11=estpar[23],alpha12=estpar[24],alpha13=estpar[25],alpha14=estpar[26],
                       epsilon=1/4.5,gamma=1/1.5)
    out_boot <- ode(y=init, times=times, func=seir, parms=allpara_boot, method="lsoda",rtol = 1e-8,atol = 1e-8)
    out_boot <- as.data.frame(out_boot)[, -1]
    estdailyI_boot[[j]] <- epsilon * as.matrix(out_boot[,(n_group+1):(2*n_group)])
  }
  estdailyI_boot_nb <- vector("list", length(estdailyI_boot))
  for (j in seq_along(estdailyI_boot)) {
    lambda_mat <- estdailyI_boot[[j]]
    obs_mat <- matrix(NA, nrow = nrow(lambda_mat), ncol = ncol(lambda_mat))
    for (g in 1:ncol(lambda_mat)) {
      for (t in 1:nrow(lambda_mat)) {
        lambda_val <- lambda_mat[t, g]
        lambda_val <- pmax(lambda_val, 1e-6)
        size <- 1 / sigma_hat[j]
        obs_mat[t, g] <- rnbinom(1, size = size, mu = lambda_val)
      }
    }
    estdailyI_boot_nb[[j]] <- obs_mat
  }
  array_boot_nb <- array(unlist(estdailyI_boot_nb),
                         dim = c(nrow(estdailyI_boot_nb[[1]]),
                                 ncol(estdailyI_boot_nb[[1]]),
                                 length(estdailyI_boot_nb)))
  lower_CI <- apply(array_boot_nb, c(1,2), quantile, probs=0.025)
  upper_CI <- apply(array_boot_nb, c(1,2), quantile, probs=0.975)

  estdailyI3_lower_CI[[i - 39]] <- lower_CI
  estdailyI3_upper_CI[[i - 39]] <- upper_CI
}

### wave4
seir<-function(time,y,params){

  q <- if (time <= 10) params$q1 else
    if (time <= 20) params$q2 else
      if (time <= 30) params$q3 else
        if (time <= 40) params$q4 else
          if (time <= 50) params$q5 else
            if (time <= 60) params$q6 else
              if (time <= 70) params$q7 else
                if (time <= 80) params$q8 else
                  if (time <= 90) params$q9 else
                    if (time <= 100) params$q10 else
                      params$q11
  beta_mat <- q * C

  epsilon=params[["epsilon"]]
  gamma=params[["gamma"]]
  alpha1<-params[["alpha1"]]
  alpha2<-params[["alpha2"]]
  alpha3<-params[["alpha3"]]
  alpha4<-params[["alpha4"]]
  alpha5<-params[["alpha5"]]
  alpha6<-params[["alpha6"]]
  alpha7<-params[["alpha7"]]
  alpha8<-params[["alpha8"]]
  alpha9<-params[["alpha9"]]
  alpha10<-params[["alpha10"]]
  alpha11<-params[["alpha11"]]
  alpha12<-params[["alpha12"]]
  alpha13<-params[["alpha13"]]
  alpha14<-params[["alpha14"]]
  l=params[["l"]]

  alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14)

  S<-y[1:n_group]
  E<-y[(n_group+1):(2*n_group)]
  I<-y[(2*n_group+1):(3*n_group)]
  R<-y[(3*n_group+1):(4*n_group)]

  dS<-numeric(n_group)
  dE<-numeric(n_group)
  dI<-numeric(n_group)
  dR<-numeric(n_group)

  for(i in 1:14){
    infection_force <- l*sum(beta_mat[i,] * S[i] * alpha[i] * (3*I + E)/popN)
    dS[i] <- -infection_force
    dE[i] <- infection_force - epsilon * E[i]
    dI[i] <- epsilon * E[i] - gamma * I[i]
    dR[i] <- gamma * I[i]
  }

  solution=c(dS,dE,dI,dR)
  return(list(solution))
}
times=seq(1,110,by=1)
estdailyI4_new<-vector("list", length(40:100))
estdailyI4_lower_CI<-vector("list", length(40:100))
estdailyI4_upper_CI<-vector("list", length(40:100))
initnew4<-matrix(0,nrow=length(40:100),ncol=56)
for(i in 40:100){
  allpara<-list(q1=q_4[1],q2=q_4[2],q3=q_4[3],q4=q_4[4],q5=q_4[5],q6=q_4[6],q7=q_4[7],
                q8=q_4[8],q9=q_4[9],q10=q_4[10],q11=q_4[11],
                epsilon=1/4.5,gamma=1/1.5,l=i/100,
                alpha1=alpha_4[1],alpha2=alpha_4[2],alpha3=alpha_4[3],alpha4=alpha_4[4],alpha5=alpha_4[5],alpha6=alpha_4[6],alpha7=alpha_4[7],
                alpha8=alpha_4[8],alpha9=alpha_4[9],alpha10=alpha_4[10],alpha11=alpha_4[11],alpha12=alpha_4[12],alpha13=alpha_4[13],alpha14=alpha_4[14])
  init=initnew3[i-39,]
  out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
  out=as.data.frame(out)
  out<-out[,-1]
  estdailyI4_new[[i-39]]<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
  initnew4[i-39,]<-as.numeric(out[110,])
  ##CI
  q_boot <- q_boot_4
  alpha_boot <- alpha_boot_4
  sigma_hat <- sigma_boot_4
  estdailyI_boot <- vector("list", 1000)
  for (j in 1:1000) {
    estpar <- c(as.numeric(q_boot[j,]), as.numeric(alpha_boot[j,]))
    allpara_boot<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],
                       q8=estpar[8],q9=estpar[9],q10=estpar[10],q11=estpar[11],l=i/100,
                       alpha1=estpar[12],alpha2=estpar[13],alpha3=estpar[14],alpha4=estpar[15],alpha5=estpar[16],alpha6=estpar[17],alpha7=estpar[18],
                       alpha8=estpar[19],alpha9=estpar[20],alpha10=estpar[21],alpha11=estpar[22],alpha12=estpar[23],alpha13=estpar[24],alpha14=estpar[25],
                       epsilon=1/4.5,gamma=1/1.5)
    out_boot <- ode(y=init, times=times, func=seir, parms=allpara_boot, method="lsoda",rtol = 1e-8,atol = 1e-8)
    out_boot <- as.data.frame(out_boot)[, -1]
    estdailyI_boot[[j]] <- epsilon * as.matrix(out_boot[,(n_group+1):(2*n_group)])
  }
  estdailyI_boot_nb <- vector("list", length(estdailyI_boot))
  for (j in seq_along(estdailyI_boot)) {
    lambda_mat <- estdailyI_boot[[j]]
    obs_mat <- matrix(NA, nrow = nrow(lambda_mat), ncol = ncol(lambda_mat))
    for (g in 1:ncol(lambda_mat)) {
      for (t in 1:nrow(lambda_mat)) {
        lambda_val <- lambda_mat[t, g]
        lambda_val <- pmax(lambda_val, 1e-6)
        size <- 1 / sigma_hat[j]
        obs_mat[t, g] <- rnbinom(1, size = size, mu = lambda_val)
      }
    }
    estdailyI_boot_nb[[j]] <- obs_mat
  }
  array_boot_nb <- array(unlist(estdailyI_boot_nb),
                         dim = c(nrow(estdailyI_boot_nb[[1]]),
                                 ncol(estdailyI_boot_nb[[1]]),
                                 length(estdailyI_boot_nb)))
  lower_CI <- apply(array_boot_nb, c(1,2), quantile, probs=0.025)
  upper_CI <- apply(array_boot_nb, c(1,2), quantile, probs=0.975)

  estdailyI4_lower_CI[[i - 39]] <- lower_CI
  estdailyI4_upper_CI[[i - 39]] <- upper_CI
}

### wave5
seir<-function(time,y,params){

  q <- if (time <= 10) params$q1 else
    if (time <= 20) params$q2 else
      if (time <= 30) params$q3 else
        if (time <= 40) params$q4 else
          if (time <= 50) params$q5 else
            if (time <= 60) params$q6 else
              if (time <= 70) params$q7 else
                if (time <= 80) params$q8 else
                  if (time <= 90) params$q9 else
                    if (time <= 100) params$q10 else
                      if (time <= 110) params$q11 else
                        if (time <= 120) params$q12 else
                          if (time <= 130) params$q13 else
                            if (time <= 140) params$q14 else
                              if (time <= 150) params$q15 else
                                if (time <= 160) params$q16 else
                                  if (time <= 170) params$q17 else
                                    if (time <= 180) params$q18 else
                                      if (time <= 190) params$q19 else
                                        params$q20
  beta_mat <- q * C

  epsilon=params[["epsilon"]]
  gamma=params[["gamma"]]
  alpha1<-params[["alpha1"]]
  alpha2<-params[["alpha2"]]
  alpha3<-params[["alpha3"]]
  alpha4<-params[["alpha4"]]
  alpha5<-params[["alpha5"]]
  alpha6<-params[["alpha6"]]
  alpha7<-params[["alpha7"]]
  alpha8<-params[["alpha8"]]
  alpha9<-params[["alpha9"]]
  alpha10<-params[["alpha10"]]
  alpha11<-params[["alpha11"]]
  alpha12<-params[["alpha12"]]
  alpha13<-params[["alpha13"]]
  alpha14<-params[["alpha14"]]
  l=params[["l"]]

  alpha<-c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14)

  S<-y[1:n_group]
  E<-y[(n_group+1):(2*n_group)]
  I<-y[(2*n_group+1):(3*n_group)]
  R<-y[(3*n_group+1):(4*n_group)]

  dS<-numeric(n_group)
  dE<-numeric(n_group)
  dI<-numeric(n_group)
  dR<-numeric(n_group)

  for(i in 1:14){
    infection_force <- l*sum(beta_mat[i,] * S[i] * alpha[i] * (3*I + E)/popN)
    dS[i] <- -infection_force
    dE[i] <- infection_force - epsilon * E[i]
    dI[i] <- epsilon * E[i] - gamma * I[i]
    dR[i] <- gamma * I[i]
  }

  solution=c(dS,dE,dI,dR)
  return(list(solution))
}
times=seq(1,200,by=1)
estdailyI5_new<-vector("list", length(40:100))
estdailyI5_lower_CI<-vector("list", length(40:100))
estdailyI5_upper_CI<-vector("list", length(40:100))
for(i in 40:100){
  allpara<-list(q1=q_5[1],q2=q_5[2],q3=q_5[3],q4=q_5[4],q5=q_5[5],q6=q_5[6],q7=q_5[7],q8=q_5[8],q9=q_5[9],q10=q_5[10],
                q11=q_5[11],q12=q_5[12],q13=q_5[13],q14=q_5[14],q15=q_5[15],q16=q_5[16],q17=q_5[17],q18=q_5[18],q19=q_5[19],q20=q_5[20],
                epsilon=1/4.5,gamma=1/1.5,l=i/100,
                alpha1=alpha_5[1],alpha2=alpha_5[2],alpha3=alpha_5[3],alpha4=alpha_5[4],alpha5=alpha_5[5],alpha6=alpha_5[6],alpha7=alpha_5[7],
                alpha8=alpha_5[8],alpha9=alpha_5[9],alpha10=alpha_5[10],alpha11=alpha_5[11],alpha12=alpha_5[12],alpha13=alpha_5[13],alpha14=alpha_5[14])
  init=initnew4[i-39,]
  out=ode(y=init,times=times,func=seir,parms=allpara, method="lsoda",rtol = 1e-8,atol = 1e-8)
  out=as.data.frame(out)
  out<-out[,-1]
  estdailyI5_new[[i-39]]<-epsilon*as.matrix(out[,(n_group+1):(2*n_group)])
  ##CI
  q_boot <- q_boot_5
  alpha_boot <- alpha_boot_5
  sigma_hat <- sigma_boot_5
  estdailyI_boot <- vector("list", 1000)
  for (j in 1:1000) {
    estpar <- c(as.numeric(q_boot[j,]), as.numeric(alpha_boot[j,]))
    allpara_boot<-list(q1=estpar[1],q2=estpar[2],q3=estpar[3],q4=estpar[4],q5=estpar[5],q6=estpar[6],q7=estpar[7],q8=estpar[8],q9=estpar[9],q10=estpar[10],
                       q11=estpar[11],q12=estpar[12],q13=estpar[13],q14=estpar[14],q15=estpar[15],q16=estpar[16],q17=estpar[17],q18=estpar[18],q19=estpar[19],q20=estpar[20],l=i/100,
                       alpha1=estpar[21],alpha2=estpar[22],alpha3=estpar[23],alpha4=estpar[24],alpha5=estpar[25],alpha6=estpar[26],alpha7=estpar[27],
                       alpha8=estpar[28],alpha9=estpar[29],alpha10=estpar[30],alpha11=estpar[31],alpha12=estpar[32],alpha13=estpar[33],alpha14=estpar[34],
                       epsilon=1/4.5,gamma=1/1.5)
    out_boot <- ode(y=init, times=times, func=seir, parms=allpara_boot, method="lsoda",rtol = 1e-8,atol = 1e-8)
    out_boot <- as.data.frame(out_boot)[, -1]
    estdailyI_boot[[j]] <- epsilon * as.matrix(out_boot[,(n_group+1):(2*n_group)])
  }
  estdailyI_boot_nb <- vector("list", length(estdailyI_boot))
  for (j in seq_along(estdailyI_boot)) {
    lambda_mat <- estdailyI_boot[[j]]
    obs_mat <- matrix(NA, nrow = nrow(lambda_mat), ncol = ncol(lambda_mat))
    for (g in 1:ncol(lambda_mat)) {
      for (t in 1:nrow(lambda_mat)) {
        lambda_val <- lambda_mat[t, g]
        lambda_val <- pmax(lambda_val, 1e-6)
        size <- 1 / sigma_hat[j]
        obs_mat[t, g] <- rnbinom(1, size = size, mu = lambda_val)
      }
    }
    estdailyI_boot_nb[[j]] <- obs_mat
  }
  array_boot_nb <- array(unlist(estdailyI_boot_nb),
                         dim = c(nrow(estdailyI_boot_nb[[1]]),
                                 ncol(estdailyI_boot_nb[[1]]),
                                 length(estdailyI_boot_nb)))
  lower_CI <- apply(array_boot_nb, c(1,2), quantile, probs=0.025)
  upper_CI <- apply(array_boot_nb, c(1,2), quantile, probs=0.975)

  estdailyI5_lower_CI[[i - 39]] <- lower_CI
  estdailyI5_upper_CI[[i - 39]] <- upper_CI
}

### aggregate
estdailyIall_new<-list()
for(i in 1:length(40:100)){
  estdailyIall_new[[i]]<-rbind(estdailyI1_new[[i]],estdailyI2_new[[i]],estdailyI3_new[[i]],estdailyI4_new[[i]],estdailyI5_new[[i]])
}
estdailyIall_new_origin<-estdailyIall_new
# lower
estdailyIall_new<-list()
for(i in 1:length(40:100)){
  estdailyIall_new[[i]]<-rbind(estdailyI1_lower_CI[[i]],estdailyI2_lower_CI[[i]],estdailyI3_lower_CI[[i]],estdailyI4_lower_CI[[i]],estdailyI5_lower_CI[[i]])
}
estdailyIall_new_lower<-estdailyIall_new
# upper
estdailyIall_new<-list()
for(i in 1:length(40:100)){
  estdailyIall_new[[i]]<-rbind(estdailyI1_upper_CI[[i]],estdailyI2_upper_CI[[i]],estdailyI3_upper_CI[[i]],estdailyI4_upper_CI[[i]],estdailyI5_upper_CI[[i]])
}
estdailyIall_new_upper<-estdailyIall_new

### rate of change
estdailyIall_new<-estdailyIall_new_origin
estdailyI_monthly_new<-list()
for(i in 1:length(40:100)){
  estdailyI_monthly_n<-aggregate(estdailyIall_new[[i]], by = list(month = months), FUN = sum)
  estdailyI_monthly_new[[i]]<-rbind(rep(0,14),rep(1,14),estdailyI_monthly_n[,-1])
}
estdailyI_monthly_all_new<-matrix(0,nrow=25,ncol=length(40:100))
for(i in 1:length(40:100)){
  estdailyI_monthly_all_new[,i]<-rowSums(estdailyI_monthly_new[[i]])
}
rate_of_change_m<-matrix(0,nrow=24,ncol=length(40:100))
for(i in 1:length(40:100)){
  rate_of_change_m[,i]<-diff(estdailyI_monthly_all_new[,i])/(head(estdailyI_monthly_all_new[,i],-1)+1e-1)
}
rate_of_change_m<-rbind(rep(0,length(40:100)),rep(0,length(40:100)),rate_of_change_m)
# lower
estdailyIall_new<-estdailyIall_new_lower
estdailyI_monthly_new<-list()
for(i in 1:length(40:100)){
  estdailyI_monthly_n<-aggregate(estdailyIall_new[[i]], by = list(month = months), FUN = sum)
  estdailyI_monthly_new[[i]]<-rbind(rep(0,14),rep(1,14),estdailyI_monthly_n[,-1])
}
estdailyI_monthly_all_new<-matrix(0,nrow=25,ncol=length(40:100))
for(i in 1:length(40:100)){
  estdailyI_monthly_all_new[,i]<-rowSums(estdailyI_monthly_new[[i]])
}
rate_of_change_m<-matrix(0,nrow=24,ncol=length(40:100))
for(i in 1:length(40:100)){
  rate_of_change_m[,i]<-diff(estdailyI_monthly_all_new[,i])/(head(estdailyI_monthly_all_new[,i],-1)+1e+1)
}
rate_of_change_m_lower<-rbind(rep(0,length(40:100)),rep(0,length(40:100)),rate_of_change_m)
# upper
estdailyIall_new<-estdailyIall_new_upper
estdailyI_monthly_new<-list()
for(i in 1:length(40:100)){
  estdailyI_monthly_n<-aggregate(estdailyIall_new[[i]], by = list(month = months), FUN = sum)
  estdailyI_monthly_new[[i]]<-rbind(rep(0,14),rep(1,14),estdailyI_monthly_n[,-1])
}
estdailyI_monthly_all_new<-matrix(0,nrow=25,ncol=length(40:100))
for(i in 1:length(40:100)){
  estdailyI_monthly_all_new[,i]<-rowSums(estdailyI_monthly_new[[i]])
}
rate_of_change_m<-matrix(0,nrow=24,ncol=length(40:100))
for(i in 1:length(40:100)){
  rate_of_change_m[,i]<-diff(estdailyI_monthly_all_new[,i])/(head(estdailyI_monthly_all_new[,i],-1)+1e-1)
}
rate_of_change_m_upper<-rbind(rep(0,length(40:100)),rep(0,length(40:100)),rate_of_change_m)

############################################################
## 10. Case fatality risk (CFR)
############################################################
# death_Osaka
Deathdata <- read.csv(
  "data/deathdata_Osaka.csv",
  header = TRUE,
  row.names = 1
)
Deathdata[, "number"] <- rep(1,nrow(Deathdata))
rownames(Deathdata)<-1:nrow(Deathdata)
Deathdata_summary <- Deathdata %>%
  group_by(agedecade, sex, date_death) %>%
  summarise(number = sum(number, na.rm = TRUE), .groups = "drop") %>%
  arrange(date_death)
Deathdata_list <- split(Deathdata_summary, 
                        list(Deathdata_summary$sex, Deathdata_summary$agedecade), 
                        drop = TRUE)
Deathdata_list <- lapply(Deathdata_list, function(df) {
  df[order(df$date_death), ]
})
periods <- list(
  period1 = as.Date(c("2020-02-01", "2020-06-29")),
  period2 = as.Date(c("2020-06-30", "2020-10-27")),
  period3 = as.Date(c("2020-10-28", "2021-02-24")),
  period4 = as.Date(c("2021-02-25", "2021-06-14")),
  period5 = as.Date(c("2021-06-15", "2021-12-31"))
)
sum_by_period <- function(df, periods) {
  sapply(periods, function(dates) {
    sum(df$number[df$date_death >= dates[1] & df$date_death <= dates[2]], na.rm = TRUE)
  })
}
Deathdata_sums <- lapply(Deathdata_list, sum_by_period, periods = periods)
Deathdata_sums_df <- do.call(rbind, Deathdata_sums)
Deathdata_sums_df <- as.data.frame(Deathdata_sums_df)
Deathdata_sums_df$group <- rownames(Deathdata_sums_df)
rownames(Deathdata_sums_df) <- NULL
Deathdata_sums_df <- Deathdata_sums_df %>%
  separate(group, into = c("sex", "age_start"), sep = "\\.")
Deathdata_sums_df$sex <- ifelse(Deathdata_sums_df$sex == "M", "male", "female")
Deathdata_sums_df$agedecade <- paste0(Deathdata_sums_df$age_start, "-", as.numeric(Deathdata_sums_df$age_start) + 9)
Deathdata_sums_df$sex <- factor(Deathdata_sums_df$sex, levels = c("male", "female"))
Deathdata_sums_df$age_start <- as.numeric(Deathdata_sums_df$age_start)
Deathdata_sums_df <- Deathdata_sums_df %>%
  arrange(sex, age_start)
Deathdata_sums_df <- Deathdata_sums_df %>%
  dplyr::select(sex, agedecade, everything(), -age_start)
deathwave1<-c(Deathdata_sums_df$period1[1],0,Deathdata_sums_df$period1[2:5],sum(Deathdata_sums_df$period1[6:9]),0,Deathdata_sums_df$period1[10:14],sum(Deathdata_sums_df$period1[15:18]))
deathwave2<-c(Deathdata_sums_df$period2[1],0,Deathdata_sums_df$period2[2:5],sum(Deathdata_sums_df$period2[6:9]),0,Deathdata_sums_df$period2[10:14],sum(Deathdata_sums_df$period2[15:18]))
deathwave3<-c(Deathdata_sums_df$period3[1],0,Deathdata_sums_df$period3[2:5],sum(Deathdata_sums_df$period3[6:9]),0,Deathdata_sums_df$period3[10:14],sum(Deathdata_sums_df$period3[15:18]))
deathwave4<-c(Deathdata_sums_df$period4[1],0,Deathdata_sums_df$period4[2:5],sum(Deathdata_sums_df$period4[6:9]),0,Deathdata_sums_df$period4[10:14],sum(Deathdata_sums_df$period4[15:18]))
deathwave5<-c(Deathdata_sums_df$period5[1],0,Deathdata_sums_df$period5[2:5],sum(Deathdata_sums_df$period5[6:9]),0,Deathdata_sums_df$period5[10:14],sum(Deathdata_sums_df$period5[15:18]))

### case_Osaka
Casedata <- read.csv(
  "data/casedata_Osaka.csv",
  header = TRUE,
  row.names = 1
)
df_summary <- Casedata %>%
  mutate(age = as.numeric(age),  
         age_group = case_when(
           age < 20 ~ "0-19",
           age >= 20 & age < 30 ~ "20-29",
           age >= 30 & age < 40 ~ "30-39",
           age >= 40 & age < 50 ~ "40-49",
           age >= 50 & age < 60 ~ "50-59",
           age >= 60 & age < 70 ~ "60-69",
           age >= 70 ~ "70+"
         )) %>%
  group_by(Date, age_group) %>%
  summarise(case = sum(case, na.rm = TRUE), .groups = "drop") %>%
  arrange(Date, age_group)
get_period <- function(date, periods) {
  for (nm in names(periods)) {
    range <- periods[[nm]]
    if (date >= range[1] && date <= range[2]) {
      return(nm)
    }
  }
  return(NA_character_)
}
df_summary$period <- sapply(df_summary$Date, get_period, periods = periods)
df_period_summary <- df_summary %>%
  filter(!is.na(period)) %>%
  group_by(age_group, period) %>%
  summarise(case = sum(case, na.rm = TRUE), .groups = "drop")
df_wide <- df_period_summary %>%
  pivot_wider(names_from = period, values_from = case, values_fill = 0) %>%
  arrange(factor(age_group, levels = c("0-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70+")))
wave1df<-wave1df[1:150,]
wave2df<-wave2df[1:120,]
wave3df<-wave3df[1:120,]
wave4df<-wave4df[1:110,]
wave5df<-wave5df[1:200,]
age_groups <- c("u20" = "0-19", "20" = "20-29", "30" = "30-39", "40" = "40-49", 
                "50" = "50-59", "60" = "60-69", "a70" = "70+")
get_gender_ratio <- function(wavedf, period_name) {
  male_cols <- grep("^male", names(wavedf), value = TRUE)
  female_cols <- grep("^female", names(wavedf), value = TRUE)
  male_sum <- colSums(wavedf[male_cols])
  female_sum <- colSums(wavedf[female_cols])
  
  age_keys <- gsub("male|female", "", names(male_sum))
  
  data.frame(
    age_group = age_groups[age_keys],
    male = as.numeric(male_sum),
    female = as.numeric(female_sum),
    period = period_name,
    stringsAsFactors = FALSE
  )
}
gender_ratios <- bind_rows(
  get_gender_ratio(wave1df, "period1"),
  get_gender_ratio(wave2df, "period2"),
  get_gender_ratio(wave3df, "period3"),
  get_gender_ratio(wave4df, "period4"),
  get_gender_ratio(wave5df, "period5")
) %>%
  mutate(
    total = male + female,
    male_ratio = ifelse(total > 0, male / total, 0.5),
    female_ratio = 1 - male_ratio
  )
df_long <- df_wide %>%
  pivot_longer(cols = starts_with("period"), names_to = "period", values_to = "total_cases") %>%
  left_join(gender_ratios, by = c("age_group", "period")) %>%
  mutate(
    male_cases = total_cases * male_ratio,
    female_cases = total_cases * female_ratio
  )
df_split_by_sex <- df_long %>%
  dplyr::select(age_group, period, male_cases, female_cases) %>%
  pivot_longer(cols = c(male_cases, female_cases),
               names_to = "sex", values_to = "cases") %>%
  mutate(
    sex = ifelse(sex == "male_cases", "male", "female")
  ) %>%
  arrange(sex, factor(age_group, levels = c("0-19", "20-29", "30-39", "40-49", 
                                            "50-59", "60-69", "70+")), period)
infwave1 <- round(c(
  df_split_by_sex$cases[df_split_by_sex$period=="period1" & df_split_by_sex$sex=="male"],
  df_split_by_sex$cases[df_split_by_sex$period=="period1" & df_split_by_sex$sex=="female"]
), 0)
infwave2 <- round(c(
  df_split_by_sex$cases[df_split_by_sex$period=="period2" & df_split_by_sex$sex=="male"],
  df_split_by_sex$cases[df_split_by_sex$period=="period2" & df_split_by_sex$sex=="female"]
), 0)
infwave3 <- round(c(
  df_split_by_sex$cases[df_split_by_sex$period=="period3" & df_split_by_sex$sex=="male"],
  df_split_by_sex$cases[df_split_by_sex$period=="period3" & df_split_by_sex$sex=="female"]
), 0)
infwave4 <- round(c(
  df_split_by_sex$cases[df_split_by_sex$period=="period4" & df_split_by_sex$sex=="male"],
  df_split_by_sex$cases[df_split_by_sex$period=="period4" & df_split_by_sex$sex=="female"]
), 0)
infwave5 <- round(c(
  df_split_by_sex$cases[df_split_by_sex$period=="period5" & df_split_by_sex$sex=="male"],
  df_split_by_sex$cases[df_split_by_sex$period=="period5" & df_split_by_sex$sex=="female"]
), 0)

### cfr
calc_cfr <- function(deaths, infections) {
  deaths / infections
}

calc_cfr_ci <- function(deaths, infections) {
  ci <- matrix(NA_real_, nrow = length(deaths), ncol = 2)
  for (k in seq_along(deaths)) {
    ci[k, ] <- as.numeric(
      prop.test(deaths[k], infections[k], correct = FALSE)$conf.int
    )
  }
  ci
}

cfrwave1   <- calc_cfr(deathwave1, infwave1)
cfrwave2   <- calc_cfr(deathwave2, infwave2)
cfrwave3   <- calc_cfr(deathwave3, infwave3)
cfrwave4   <- calc_cfr(deathwave4, infwave4)
cfrwave5   <- calc_cfr(deathwave5, infwave5)

cfrwave1CI <- calc_cfr_ci(deathwave1, infwave1)
cfrwave2CI <- calc_cfr_ci(deathwave2, infwave2)
cfrwave3CI <- calc_cfr_ci(deathwave3, infwave3)
cfrwave4CI <- calc_cfr_ci(deathwave4, infwave4)
cfrwave5CI <- calc_cfr_ci(deathwave5, infwave5)

### estimate COVID-19 deaths
n_rho <- 61 
est   <- numeric(n_rho)
lower <- numeric(n_rho)
upper <- numeric(n_rho)

# inf
for (i in seq_len(n_rho)) {
  est[i]   <- sum(as.matrix(estdailyIall_new_origin[[i]]))
  lower[i] <- sum(as.matrix(estdailyIall_new_lower[[i]]))
  upper[i] <- sum(as.matrix(estdailyIall_new_upper[[i]]))
}

# death
wave_slices <- list(
  wave1 = 1:150,
  wave2 = 151:270,
  wave3 = 271:390,
  wave4 = 391:500,
  wave5 = 501:700
)

estimate_deaths_by_rho <- function(estdailyI_list,
                                   cfr1, cfr2, cfr3, cfr4, cfr5) {
  out <- matrix(0, nrow = length(estdailyI_list), ncol = 14)
  
  for (i in seq_along(estdailyI_list)) {
    mat <- estdailyI_list[[i]]
    
    I1 <- mat[wave_slices$wave1, , drop = FALSE]
    I2 <- mat[wave_slices$wave2, , drop = FALSE]
    I3 <- mat[wave_slices$wave3, , drop = FALSE]
    I4 <- mat[wave_slices$wave4, , drop = FALSE]
    I5 <- mat[wave_slices$wave5, , drop = FALSE]
    
    out[i, ] <-
      as.numeric(colSums(I1)) * cfr1 +
      as.numeric(colSums(I2)) * cfr2 +
      as.numeric(colSums(I3)) * cfr3 +
      as.numeric(colSums(I4)) * cfr4 +
      as.numeric(colSums(I5)) * cfr5
  }
  out
}

estdeath_df <- estimate_deaths_by_rho(
  estdailyI_list = estdailyIall_new_origin,
  cfr1 = cfrwave1, cfr2 = cfrwave2, cfr3 = cfrwave3, cfr4 = cfrwave4, cfr5 = cfrwave5
)
estdeath_df_lower <- estimate_deaths_by_rho(
  estdailyI_list = estdailyIall_new_lower,
  cfr1 = cfrwave1CI[, 1], cfr2 = cfrwave2CI[, 1], cfr3 = cfrwave3CI[, 1],
  cfr4 = cfrwave4CI[, 1], cfr5 = cfrwave5CI[, 1]
)
estdeath_df_upper <- estimate_deaths_by_rho(
  estdailyI_list = estdailyIall_new_upper,
  cfr1 = cfrwave1CI[, 2], cfr2 = cfrwave2CI[, 2], cfr3 = cfrwave3CI[, 2],
  cfr4 = cfrwave4CI[, 2], cfr5 = cfrwave5CI[, 2]
)

estdeathall       <- rowSums(estdeath_df)
estdeathall_lower <- rowSums(estdeath_df_lower)
estdeathall_upper <- rowSums(estdeath_df_upper)

### plot
plot_df <- data.frame(
  rho_value = 0:60,
  est   = rev(estdeathall),
  lower = rev(estdeathall_lower),
  upper = rev(estdeathall_upper)
)
y_min <- min(plot_df$lower, na.rm = TRUE)
y_max <- max(plot_df$upper, na.rm = TRUE)
ggplot(plot_df, aes(x = rho_value, y = est)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "gray30",
    alpha = 0.3
  ) +
  geom_line(color = "black", linewidth = 1) +
  geom_line(aes(y = lower), color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = upper), color = "black", linetype = "dashed", linewidth = 0.8) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    labels = label_comma()
  ) +
  scale_x_continuous(
    breaks = c(0, 10, 20, 30, 40, 50, 60)
  ) +
  labs(
    x = "Additional relative control effort (ρ)",
    y = "Estimated COVID-19 deaths"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 28),
    axis.text  = element_text(size = 24),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  )
# ggsave("results/figures/est_COVIDdeath.jpg",width=10,height=8,dpi=300)

############################################################
# End of script
############################################################