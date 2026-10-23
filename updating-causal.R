#clean file with model updating and validation
#also interval creating function

###############################################################################

#libraries needed:
library(survival)
library(tidyverse)
#needed for validation functions:
library(riskRegression)
library(MatrixModels)
library(geepack)

# For plots:
library(reshape2) # For melt
library(ggpubr)   # For ggarrange

#############################################################################
#load data 
load("data.rda")

#load functions for validation of model
source("validation-functions.R") 

################################################################################
# -------- First wave model ----------------------------------------------------
################################################################################

#covariates:
#all covariates to include in cox model
var <- c("age", "sexM", "numb_comorb", "rr_min", "sp_o2", "crp_mg_l", "creatinine_ser_mmol_l")

#cox model fitted on first wave data
fmla <- formula(paste("Surv(survtime, status) ~ ", paste(var, collapse = "+")))
model_w1 <- coxph(fmla, data=datasub_w1_cc, method="breslow", x=T, y=T)


################################################################################
# ------- Second wave ----------------------------------------------------------
################################################################################

#linear predictor for new data (wave 2) based on original model
datasub_w2_cc$lp_origmod <- predict(model_w1, newdata = datasub_w2_cc, type = "lp", reference = "zero") 
# datasub_w2_cc$lp_origmod <- predict(model_w1, newdata = datasub_w2_cc, type = "lp", reference = "sample") ### NEW

################################################################################
# --------- Intervals/periods --------------------------------------------------
################################################################################
startdate <- as.Date("2020-10-15") # Start validation dataset (and date of first update)
enddate <- as.Date("2020-12-31")   # Predict for patients who are admitted until this date
InDate <- as.Date("2020-08-01")    # Introduction date: Here dexamethason is introduced
OutDate <- as.Date("2020-11-01")   # From here we have the "new standard care"
N_period <- 18                     # Overall number of updates
newevents <- 25 
oldevents <- 115

# -----------
# creating intervals with 180 events:
# Intervals for model fitting: 10 new events, 180 old events
df_int_fit0 <- lapply(1:N_period, function(n_int) {
  make_intervals(data = datasub_w2_cc, startdate, n_events = newevents, n_interval = n_int, n_oldevents=oldevents)
  }
)

# Intervals for prediction: 0 old events
df_int_predd <- lapply(1:(N_period+1), function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = newevents, n_interval = n_int, n_nextevents=100 - newevents)
}
)

# Check the starting date of each period
startdates_fit <- lapply(1:N_period, function(a) as.Date(min(df_int_fit0[[a]]$admitted_date)))
startdates_pred <- lapply(1:(N_period+1), function(a) as.Date(min(df_int_predd[[a]]$admitted_date))) 
# we go to N_100+1 to see when the next update would have been 
N <- sum(do.call(c, startdates_pred) <=  enddate)
df_int_fit0 <- df_int_fit0[1:N]
df_int_predd <- df_int_predd[1:N]

# Match each period to an "interventional" updating phase: Before new treatment 
# introduction ("Before"), Transitional phase ("During"), New standard care ("After")
update_periods <- lapply(1:N, function(a) {
  case_when(startdates_pred[[a]] < InDate ~ "Before", 
            startdates_fit[[a]] >= OutDate ~ "After",
            TRUE ~ "During")})

interventional_N <- sum(do.call(c, update_periods) == "During")

period_starts <- do.call(c,startdates_pred)
datasub_w2_cc$update_interval <- findInterval(as.Date(datasub_w2_cc$admitted_date), period_starts)
df_pred <- filter(datasub_w2_cc, update_interval != c(0,(N+1)))

df_int_pred <- lapply(1:N, function(a) filter(df_pred, update_interval == a & admitted_date <= enddate))
df_int_pred_all <- do.call(rbind, df_int_pred)

periods <- data.frame(
  period = factor(paste0("Period ", 1:N), levels = paste0("Period ", N:1)),
  type = c(rep("Interventional update", interventional_N), rep("Factual update", N - interventional_N)),
  start  = do.call(c, startdates_fit[1:N]),
  end    = do.call(c, startdates_pred[1:N])
)


################################################################################
# --------------- Applying standard methods ------------------------------------
################################################################################

#no updating
no_update <- getValidation(model_w1, df_int_pred_all)
# no_update_int <- getValidation(model_w1, df_int_pred_int)

#refitting: 100 events
refit_cfit <- lapply(df_int_fit0, function(a) my_cox(fmla, a, method = "breslow", x=T, y=T))
refit_cfit <- c(list(model_w1), refit_cfit[-N])
refit_pred <- lapply(1:N, function(a) predCox(refit_cfit[[a]], df_int_pred[[a]], time = 28))
refit <- getValidation(do.call(c, refit_pred), df_int_pred_all)
# refit_int <- getValidation(do.call(c, refit_pred[1:interventional_N]), df_int_pred_int)

#recalibration: 100 events
recal_fmlaO <- formula("Surv(survtime, status) ~ offset(lp_origmod)")
recal_cfit <- lapply(df_int_fit0, function(a) my_cox(recal_fmlaO, a, method = "breslow", x=T, y=T))
recal_cfit <- c(list(model_w1), recal_cfit[-N])
recal_pred <- lapply(1:N, function(a) predCox(recal_cfit[[a]], df_int_pred[[a]], time = 28))
recal <- getValidation(do.call(c, recal_pred), df_int_pred_all)
# recal_int <- getValidation(do.call(c, recal_pred[1:interventional_N]), df_int_pred_int)


################################################################################
# --------------- Applying interventional methods ------------------------------
################################################################################

# Treatment effect HR on death
# By groups: Invasive Mechanical Ventilation 0.64 (N=1007), 
# Oxygen Only 0.82 (N=3883), 
# No Oxygen Received 1.19 (N=1535)

# Multiply treatment by treatment effect and create lp_origmod_dex (sum of lp_origmodel 
# and dexamethason effect) in both development data and prediction data
# For the small number of patients who received dexamethasone prior to August A is set to 0

mutate_dex <- function(data) {
  data <- data |> 
  mutate(dexamethason_base = if_else(admitted_date < InDate, 0, dexamethason_base)) |>
  mutate(dex = case_when(
    MV == "No oxygen" ~ dexamethason_base*log(1.19), 
    MV == "Oxygen" ~ dexamethason_base*log(0.82), 
    MV == "Invasive" ~ dexamethason_base*log(0.64) 
    )) |>
  mutate(dex_marg = dexamethason_base * log(0.83)) |>
  mutate(lp_origmod_dex = lp_origmod + dex,
         lp_origmod_dex_marg = lp_origmod + dex_marg) |>
  select(-c(dexamethason))
  
  return(data)
}

refit_fit0_dex <- lapply(df_int_fit0, mutate_dex)
df_int_pred_dex <- lapply(df_int_pred, mutate_dex)

model_w1_i <- add_coeff(model_w1, "dex", 1)

# Interventional refitting
refit_fmlai <- lapply(update_periods, function(a) refit_fmla(model_w1, a, "dex"))
refit_cfiti <- lapply(1:N, function(a) my_cox(refit_fmlai[[a]], refit_fit0_dex[[a]], method = "breslow", x=T, y=T))
refit_cfiti <- c(list(model_w1_i),refit_cfiti[-N])
refiti_pred <- lapply(1:N, function(a) predCox(refit_cfiti[[a]], df_int_pred_dex[[a]], time = 28))
refit_i <- getValidation(do.call(c, refiti_pred), df_int_pred_all)

# Interventional recalibration
recal_fmlai <- lapply(update_periods, function(a) recal_fmla(model_w1, a, "lp_origmod", "lp_origmod_dex"))
recal_cfiti <- lapply(1:N, function(a) my_cox(recal_fmlai[[a]], refit_fit0_dex[[a]], method = "breslow", x=T, y=T))
recal_cfiti <- c(list(model_w1_i),recal_cfiti[-N])
recali_pred <- lapply(1:N, function(a) predCox(recal_cfiti[[a]], df_int_pred_dex[[a]], time = 28))
recal_i <- getValidation(do.call(c, recali_pred), df_int_pred_all)

# Inverventional refitting by subgroups
refitgr_fmlai <- lapply(update_periods, function(a) refit_fmla(model_w1, a, "dex", "MV_Invasive + MV_Oxygen"))
refitgr_cfiti <- lapply(1:N, function(a) my_cox(refitgr_fmlai[[a]], refit_fit0_dex[[a]], method = "breslow", x=T, y=T))
refitgr_cfiti <- c(list(model_w1_i), refitgr_cfiti[-N])
refitgri_pred <- lapply(1:N, function(a) predCox(refitgr_cfiti[[a]], df_int_pred_dex[[a]], time = 28))
refitgr_i <- getValidation(do.call(c, refitgri_pred), df_int_pred_all)

# Interevntional recalibration by subgroups 
recalgr_fmlai <- lapply(update_periods, function(a) recal_fmla(model_w1, a, "lp_origmod", "lp_origmod_dex", "MV_Invasive + MV_Oxygen"))
recalgr_cfiti <- lapply(1:N, function(a) my_cox(recalgr_fmlai[[a]], refit_fit0_dex[[a]], method = "breslow", x=T, y=T))
recalgr_cfiti <- c(list(model_w1_i), recalgr_cfiti[-N])
recalgri_pred <- lapply(1:N, function(a) predCox(recalgr_cfiti[[a]], df_int_pred_dex[[a]], time = 28))
recalgr_i <- getValidation(do.call(c, recalgri_pred), df_int_pred_all)


################################################################################
# --------------- Numerical summaries ------------------------------------------
################################################################################

ll <- list(no_update, refit, refit_i, refitgr_i, recal, recal_i, recalgr_i)
names(ll)<- c("Never update", "Refit", "Interventional Refit", "Interventional Refit (groups)",
              "Recalibration", "Interventional Recalibration",  "Interventional Recalibration (groups)")

# Overall
summm1 <- do.call(rbind, lapply(ll, function(a) {
  x <- paste(a[,1], " [", a[,2], "-", a[,3], "]", sep = "")
  names(x) <- rownames(a)
  x
  }))
summm1[, 2:7]

filename <- paste("Results/summ", newevents, "n", oldevents, "o.rda", sep = "")
save(summm1, periods, startdate, enddate, InDate, OutDate, file = filename)