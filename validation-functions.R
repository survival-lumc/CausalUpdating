################################################################################
# -------------------- CREATE INTERVALS ----------------------------------------
################################################################################

make_intervals <- function(data, startdate, n_events, n_interval, n_oldevents = 0,
                           n_nextevents = 0, sliding_win = 0) {
  
  # n_oldevents == include these many old events
  # sliding window == include a fixed amount of old data, eg 1 month (1/12) or 2 weeks (1/24)
  
  if (n_oldevents != 0 & sliding_win != 0) {
    return("Only set either n_oldevents != 0 or sliding_win != 0")
  }
  
  # Arrange data by event date
  data <- arrange(data, date_event)
  enddate <- NULL
  
  # Find startdate and enddate for interval number n
  # Here I assume the new interval starts the day after the previous one ends
  for (i in 1:n_interval) { 
    data$event_count <- cumsum((data$admitted_date >= startdate) & (data$status == 1))
    enddate <- as.Date(first(data$date_event[data$event_count == n_events])) 
    if (i != n_interval) startdate <- enddate + 1 
  }
  
  if (n_nextevents != 0) enddate <- as.Date(data$date_event[data$event_count == n_events + n_nextevents][1])
  
  if (n_oldevents != 0) {
    # From startdate look in the past until you see the correct number of old events
    data$oldevent_count <- sum((data$admitted_date < startdate) & (data$status == 1)) -
      cumsum((data$admitted_date < startdate) & (data$status == 1))
    startdate <- as.Date(last(data$date_event[data$oldevent_count == n_oldevents]))
    # The last line will return NA if the condition is never met and there
    # simply are not that many old events in the dataset
    
    # Then further look back to find the correct starting date for those events
    for (i in 1:50) { # I must have seen the events within 28 days, which is time horizon
      df <- subset(data, admitted_date >= startdate & admitted_date < enddate+1)
      df$status <- if_else(df$date_event >= enddate + 1, 0, df$status)
      evnts <- (sum(df$status) < n_events + n_oldevents + n_nextevents)
      if (evnts) startdate <- startdate - 1
    }
    if (is.na(startdate)) startdate <- as.Date(min(data$admitted_date))
    data <- select(data, -c(oldevent_count))
  }
  
  # Select data between starting date (minus sliding window) and end date 
  data <- subset(data, admitted_date >= startdate - sliding_win & admitted_date < enddate + 1)
  
  # censor event after end date
  data$status <- if_else(data$date_event >= enddate + 1, 0, data$status)
  # max survival time put to date of enddate
  data$survtime <- if_else(data$date_event >= enddate + 1,
                           pmin(28.001,round(as.numeric(difftime(enddate + 1, data$admitted_date, units = "days")),2)), 
                                data$survtime)
  
  data <- select(data, -c(event_count))
  
  return(data)
}


################################################################################
# ---------------- MODEL FITTING -----------------------------------------------
################################################################################

remove_offset <- function(model, offsetvar, mean_offsetvar) {
  model$linear.predictors <- model$linear.predictors - mean_offsetvar   # model$linear.predictor is the uncentered offset + centered lp
  lp <- model$linear.predictors # Linear predictors is equal to a centered lp prediction, already including the offset term
  newnames <- c(names(model$coefficients), offsetvar)
  model$coefficients <- c(model$coefficients, 1)
  names(model$coefficients) <- newnames
  model$means <- c(model$means, mean_offsetvar) 
  names(model$means) <- newnames
  L <- length(newnames)
  fmlatxt <- paste(
    paste(model$formula[2]), "~",
    gsub("\\)", "",(gsub("\\(", "", gsub("offset", "", paste(model$formula)[3])))))
  terms <- formula(fmlatxt)
  attributes(terms) <- attributes(model$terms)
  attr(terms, "variables") <- parse(text = paste("list(", paste(model$formula[2]), ",",
                                                 paste(newnames, collapse = ", "), ")", sep = ""))[[1]]
  newmat <- matrix(0L, nrow = L + 1, ncol = L + 1, 
                   dimnames = list(c("Surv(survtime, status)", newnames),
                                   c("Surv(survtime, status)", newnames)))
  diag(newmat) <- 1L
  attr(terms, "factors") <- newmat[,-1, drop = FALSE]
  attr(terms, "term.labels") <- newnames
  attr(terms, "order") <- rep(1L, L)
  attr(terms, "offset") <- NULL
  attr(terms, "predvars") <- parse(text = paste("list(", paste(model$formula[2]), ",",
                                                paste(newnames, collapse = ", "), ")", sep = ""))[[1]]
  attr(attr(terms, "dataClasses"), "names") <- c("Surv(survtime, status)", newnames)
  model$terms <- terms
  model$assign[[L]] <- as.integer(L)
  names(model$assign) <- newnames
  if (is.null(model$x)) {
    x <- matrix(c(lp), dimnames = list(attr(model$x, "dimnames")[[1]], offsetvar))
    attr(x, "assign") = c(1L)
  }
  if (!is.null(model$x)) {
    newvar <- round(model$linear.predictors -
                      rowSums(t(t(model$x)*model$coefficients[-L])) + sum(model$means*model$coefficients), 3)
    x <- cbind(model$x, newvar)
    colnames(x) <- newnames
    attr(x, "assign") = 1:L
  }
  model$x <- x
  model$formula <- formula(fmlatxt, env = globalenv())
  model$offset <- NULL
  model$call <- parse(text = paste('coxph(', fmlatxt, ', data = ', paste(model$call)[3], 
                                   ', x = T, y = T, method = "breslow")', sep = ""))[[1]]
  model$loglik <- coxph(model$y ~ 1,  method = "breslow")$loglik
  if (!is.null(model$var)) model$var <- rbind(cbind(model$var, 0), 0)
  if (is.null(model$var)) model$var <- 0
  return(model)
}

add_coeff <- function(model, namevar, coeffvar) {
  newnames <- c(names(model$coefficients), namevar)
  model$coefficients <- c(model$coefficients, coeffvar)
  names(model$coefficients) <- newnames
  model$means <- c(model$means, 0) 
  names(model$means) <- newnames
  L <- length(newnames)
  fmlatxt <- paste(model$formula[2], "~", model$formula[3], "+", namevar)
  terms <- formula(fmlatxt)
  attributes(terms) <- attributes(model$terms)
  attr(terms, "variables") <- parse(text = paste("list(", paste(model$formula[2]), ",",
                                                 paste(newnames, collapse = ", "), ")", sep = ""))[[1]]
  newmat <- matrix(0L, nrow = L + 1, ncol = L + 1, 
                   dimnames = list(c("Surv(survtime, status)", newnames),
                                   c("Surv(survtime, status)", newnames)))
  diag(newmat) <- 1L
  attr(terms, "factors") <- newmat[,-1, drop = FALSE]
  attr(terms, "term.labels") <- newnames
  attr(terms, "order") <- rep(1L, L)
  attr(terms, "offset") <- NULL
  attr(terms, "predvars") <- parse(text = paste("list(", paste(model$formula[2]), ",",
                                                paste(newnames, collapse = ", "), ")", sep = ""))[[1]]
  attr(terms, "dataClasses") <- c(attr(terms, "dataClasses"), newvar = "numeric")
  attr(attr(terms, "dataClasses"), "names") <- c("Surv(survtime, status)", newnames)
  model$terms <- terms
  model$assign[[L]] <- as.integer(L)
  names(model$assign) <- newnames
  if (is.null(model$x)) {
    x <- matrix(0, dimnames = list(attr(model$x, "dimnames")[[1]], namevar))
    attr(x, "assign") = c(1L)
  }
  if (!is.null(model$x)) {
    newvar <- rep(0, nrow(model$x))
    x <- cbind(model$x, newvar)
    colnames(x) <- newnames
    attr(x, "assign") = 1:L
  }
  model$x <- x
  model$formula <- formula(fmlatxt, env = globalenv())
  model$offset <- NULL
  model$call <- parse(text = paste('coxph(', fmlatxt, ', data = ', paste(model$call)[3], 
                                   ', x = T, y = T, method = "breslow")', sep = ""))[[1]]
  model$loglik <- coxph(model$y ~ 1,  method = "breslow")$loglik
  if (!is.null(model$var)) model$var <- rbind(cbind(model$var, 0), 0)
  if (is.null(model$var)) model$var <- 0
  return(model)
}

my_cox<- function(fmla, data, x = TRUE, y = TRUE, method = "breslow") { 
  
  cfit <- coxph(fmla, data, x = x, y = y, method = method)
  
  if (!is.null(cfit$offset)) {
    offvarname <- gsub("\\)", "", gsub(".*\\(", "", paste(fmla)[3]))
    offvar <- data[, offvarname]
    cfit <- remove_offset(cfit, offvarname, mean(offvar))
  }
  
  return(cfit)
}


################################################################################
# ---- INTERVENTIONAL UPDATES: FORMULAS PER INTERVAL TYPE ----------------------
################################################################################


refit_fmla <- function (origmodel, period_type, newvar, adjvars = NULL) {
  
  if (!period_type %in% c("Before", "During", "After")) {
    return("Specify period type : Before, During of After")
  }
  
  fmla <- origmodel$formula
  
  if (period_type == c("During")) {
    lhs <- paste(fmla)[2]                               # Left-handside of the formula
    rhs <- paste(paste(fmla)[3], adjvars, sep = " + ")  # Right-handside of the formula
    fmla <- formula(paste(lhs, " ~ ", rhs, " + offset(", newvar, ")", sep = "")) # Add treatment as offset
  }
  
  return(fmla)
}

recal_fmla <- function (origmodel, period_type, var1, var2, adjvars = NULL) {
  
  if (!period_type %in% c("Before", "During", "After")) {
    return("Specify period type : Before, During of After")
  }
  
  if (period_type %in% c("Before", "After")) {
    fmla <- formula(paste(paste(origmodel$formula)[2], " ~ offset(", var1, ")", sep = ""))
    }
  
  if (period_type == c("During")) {
    fmla <- formula(paste(paste(origmodel$formula)[2], "~" , adjvars, if(!is.null(adjvars)){"+"}, " offset(", var2, ")", sep = ""))
  }
  return(fmla)
  
}


refit_newvar_fmla <- function (origmodel, period_type, newvar) {
  
  if (!period_type %in% c("Before", "During", "After")) {
    return("Specify period type : Before, During of After")
  }
  
  fmla <- origmodel$formula
  
  if (period_type == c("During")) {
    lhs <- paste(fmla)[2]  # Left-handside of the formula
    rhs <- paste(fmla)[3]  # Right-handside of the formula
    fmla <- formula(paste(lhs, " ~ ", rhs, " + ", newvar, sep = "")) # Add treatment as covariate
  }
  
  return(fmla)
}

################################################################################
# -------------------- PREDICT -------------------------------------------------
################################################################################

predCox <- function (coxModel, newdata, time) {
  
  newdata <- data.frame(newdata)
  if(length(time)!=1) "Function implemented for one time-point only"
  
  bh <- basehaz(coxModel, centered = FALSE) |>
    merge(data.frame(time = time), all  = TRUE) |>
    arrange(time) |>
    fill(hazard)
  
  coefs <- coxModel$coefficients
  lp <- rowSums(t(t(newdata[, names(coefs)])*coefs))
  predRisk <- 1- exp(-bh$hazard[bh$time==time]*exp(lp))
  
  return(predRisk)
  
}


################################################################################
# ------------------ VALIDATION FUNCTIONS --------------------------------------
################################################################################

################################################################################
# In this file:
# Functions for model validation are provided.
#
# - Overall performance:
#   - Brier score
# - Discrimination:
#   - Uno's c-index
#   - AUCt
# - Calibration
#   - observed-expected ratio
#   - calibration plot
#   - calibration intercept and slope

# Code based on
# - https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal.R
# - https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_CSC_minimal.R


# ------- Overall performance -------------------------------------------------
# --------- Brier score -------------------------------------------------------

#function that computes Brier score given the model and data
getBrier <- function(model, data){
  
  if(any(is.numeric(model))) data$risk <- model # always provide the risk
  
  if(any(grepl("coxph", is(model)))) data$risk <- predCox(model, data, time = 28)
  
  bs <- riskRegression::Score(list("numeric" = data$risk), #model
                                formula = Surv(survtime, status) ~ 1,
                                data = data, #new data, only complete cases
                                conf.int = TRUE,
                                times = 28, #prediction horizon
                                cens.model = "km", #using KM to estimate observed
                                cens.method = c("ipcw", "pseudo"),
                                metrics = "brier",
                                summary = "ipa")
    
  brier <- c(as.numeric(bs$Brier$score[2,"Brier"]), #estimate Brier score
             as.numeric(bs$Brier$score[2,"lower"]), #lower
             as.numeric(bs$Brier$score[2,"upper"])) #upper
  
  nullmod <- survfit(Surv(survtime, status) ~ 1, data = data)
  nullpred <- rep(1 - summary(nullmod, times = 28)$surv, nrow(data))
  
  nullbs <- riskRegression::Score(list("numeric" = nullpred), #model
                              formula = Surv(survtime, status) ~ 1,
                              data = data, #new data, only complete cases
                              conf.int = TRUE,
                              times = 28, #prediction horizon
                              cens.model = "km", #using KM to estimate observed
                              cens.method = c("ipcw", "pseudo"), ### NEW
                              metrics = "brier",
                              summary = "ipa")
  
  
  #estimated scale Brier score (IPA) (using formula)
  ipa <- c((1-(as.numeric(bs$Brier$score[2,"Brier"])/as.numeric(nullbs$Brier$score[2,"Brier"]))), #estimate
           (1-(as.numeric(bs$Brier$score[2,"upper"])/as.numeric(nullbs$Brier$score[2,"upper"]))), #lower
           (1-(as.numeric(bs$Brier$score[2,"lower"])/as.numeric(nullbs$Brier$score[2,"lower"])))) #upper
  return(rbind(brier, ipa))
}



# --------- Discrimination ----------------------------------------------------
# ------ Uno's C index --------------------------------------------------------

#function that computes Uno's c-index with 95% CI
getCindex <- function(model, data){
  if(any(is.numeric(model))) data$lp <- log(-log(1 - model)) # always provide the risk
  
  if(any(grepl("coxph", is(model))))  data$lp <- log(-log(1 - predCox(model, data, time = 28)))
  
  Uno_C <- concordance(Surv(survtime, status) ~ lp, 
                       data = data, 
                       reverse = TRUE,
                       timewt = "n/G2") #weights?
  alpha <- 0.05
  res_C <- matrix(
    c(Uno_C$concordance,
      Uno_C$concordance - 
        qnorm(1 - alpha/2) * sqrt(Uno_C$var),
      Uno_C$concordance + 
        qnorm(1 - alpha/2) * sqrt(Uno_C$var)), 
    nrow = 1, ncol = 3, byrow = T,
    dimnames = list(c("Uno C"),
                    c("Estimate", "2.5 %", "97.5 %")))
  return(res_C)
}

# ------- AUCt ----------------------------------------------------------------

#function that computes time-dependent ROC curve estimation with 95% CI
getAUCt <- function(model, data){
  #add linear predictor (bX)
  if(any(is.numeric(model))) data$lp <- log(-log(1 - model)) # always provide the risk
  
  if(any(grepl("coxph", is(model)))) if(any(grepl("coxph", is(model))))  data$lp <- log(-log(1 - predCox(model, data, time = 28)))

  #IPCW estimation of time-dependent ROC curve
  Uno <- timeROC::timeROC(
    T = data$survtime,  # vector of event-times
    delta = data$status, # vector of event indicator
    marker = data$lp,  # marker values for which we want to compute the time-dependent ROC curves
    cause = 1, # value that indicates non-censored values
    weighting = "marginal", # -> KM estimator for censoring distribution
    times = 28, # time horizon
    ROC = TRUE,
    iid = TRUE # want to compute the iid-representation of the area under 
    #time-dependent ROC curve estimator
  )
  alpha <- 0.05
  #95% CI table of estimate
  Uno_AUC_res <- c(
    "time-dependent Uno AUC" = unname(Uno$AUC[2]),
    "2.5 %" = unname(Uno$AUC["t=28"] -
                       qnorm(1 - alpha / 2) * Uno$inference$vect_sd_1["t=28"]),
    "97. 5 %" = unname(Uno$AUC["t=28"] +
                         qnorm(1 - alpha / 2) * Uno$inference$vect_sd_1["t=28"])
  )
  return(Uno_AUC_res)
}

# ---- Calibration ------------------------------------------------------------
# ----- Observed-Expected ratio -----------------------------------------------

#function that computes observed-expected ratio with 95% CI
getOEratio <- function(model, data){
  
  time_horizon <- 28
  #observed risk:
  observed <- summary(survfit(Surv(survtime, status)~1, data = data), 
                      times = time_horizon) #KM survival curve
  observed_t <- 1 - observed$surv #observed risk (=1-surival)
  
  #predicted risk:
  if(any(is.numeric(model))) data$predrisk <- model

  if(any(grepl("coxph", is(model)))) {
    data$predrisk <-  predCox(model, data, time = time_horizon) 
  }

  expected_t <- mean(data$predrisk) #expected
  
  #ratio:
  OE_t <- observed_t / expected_t
  alpha <- .05
  OE_summary <- c(
    "OE" = OE_t,
    "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / observed$n.event)),
    "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / observed$n.event))
  )
  return(OE_summary)
}

# --- Calibration slope and intercept -------------------------------------------------------

#function that computes the calibration slope and intercept with 95% CI
#based on: https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_CSC_minimal.R
getCalslopeint <- function(model, data){
  
  if(any(is.numeric(model))) data$risk <- model # always provide the risk
  
  if(any(grepl("coxph", is(model)))) data$risk <- predCox(model, data, time = 28)
    
  cal <-riskRegression::Score(list("numeric" = data$risk), #model
                              formula = Surv(survtime, status) ~ 1,
                              data = data, #new data
                              conf.int = TRUE,
                              times = 28, #prediction horizon
                              cens.model = "km", #using KM to estimate observed
                              cens.method = c("ipcw", "pseudo"),
                              metrics = c("auc", "brier"),
                              summary = c("ipa"),
                              cause = 1,
                              plots="calibration")
  
  pseudos <- data.frame(cal$Calibration$plotframe) |> filter(model == "numeric")
  pseudos$risk[pseudos$risk == 1] <- 1 - 10^-6   #sometimes the pseudovalue is 1 and that creates issues
  pseudos$cll_pred <- log(-log(1-pseudos$risk)) #change scale

  
  #Fit model for calibration intercept
  fit_cal_int <- geese(
    pseudovalue ~ offset(cll_pred),
    data = pseudos,
    id = riskRegression_ID,
    scale.fix = TRUE,
    family = gaussian,
    mean.link = "cloglog", #back to original scale
    corstr = "independence", #correlation structure
    jack = TRUE)
  
  # Fit model for calibration slope
  fit_cal_slope <- geese(
    pseudovalue ~ offset(cll_pred) + cll_pred,
    data = pseudos,
    id = riskRegression_ID,
    scale.fix = TRUE,
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence",
    jack = TRUE)
  
  #Value and CI
  sum1 <- with(
    summary(fit_cal_slope)$mean["cll_pred", ],
    c("slope" = 1 + estimate,
      "2.5 %" = 1 + (estimate - qnorm(0.975) * san.se),
      "97.5 %" = 1 + (estimate + qnorm(0.975) * san.se)))
  #Value and CI
  sum2 <- with(
    summary(fit_cal_int)$mean,
    c("intercept" = estimate,
      "2.5 %" = estimate - qnorm(0.975) * san.se,
      "97.5 %" = estimate + qnorm(0.975) * san.se))
  output <- rbind(sum2, sum1)
  colnames(output)[1] <- "estimate"
  rownames(output) <- c("intercept", "slope")
  return(output)
}


#------------------------------------------------------------------------------
# -------- Function that returns information of all validation functions -------

getValidation <- function(model, data){
  b1 <- getBrier(model = model, data = data) #Brier+IPA score
  b1_t <- t(b1)
  cindex1 <- getCindex(model = model, data = data) #C-index Uno
  auct1 <- getAUCt(model = model, data = data) #AUCt
  oe1 <- getOEratio(model = model, data = data) #OE ratio
  #calibration slope and intercept
  slopeintercept1 <- getCalslopeint(model = model, data = data)
  slopeintercept1_t <- t(slopeintercept1)

  #summary of numerical estimates (with CI)
  summary_model <- data.frame("Brier score" = round(b1_t[,1], digits = 3),
                              "Scaled Brier" = round(b1_t[,2], 3),
                              "Unos c-index" = round(cindex1[1:3], 3),
                              "AUCt" = round(auct1, 3),
                              "OE ratio" = round(oe1, 3),
                              "Calibration intercept" = round(slopeintercept1_t[,1], 3),
                              "Calibration slope" = round(slopeintercept1_t[,2], 3))
  rownames(summary_model) <- c("estimate", "2.5%", "97.5%")
  summary_model_t <- t(summary_model)
  return(summary_model_t)
}

