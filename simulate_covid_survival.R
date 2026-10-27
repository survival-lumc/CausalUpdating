## ------------------------------------------------------------------
## Fake COVID-19 survival data (admissions Mar-Dec 2020)
## 28-day mortality, administrative censoring at 28.001 days
## ------------------------------------------------------------------
library(survival)

set.seed(2020)
n <- 5000

## Independent baseline covariates -----------------------------------
age      <- pmin(pmax(round(rnorm(n, mean = 65, sd = 13)), 18), 95)
sexM     <- rbinom(n, 1, 0.60)
n_comorb <- pmin(rpois(n, lambda = 1.5), 6)

admitted_date <- as.Date("2020-03-01") +
  sample(0:as.numeric(as.Date("2020-12-31") - as.Date("2020-03-01")), n, replace = TRUE)

## Dexamethasone at baseline ----------------------------------------
## 0 before 1 Aug 2020; then a logistic ramp up to ~70% by Nov/Dec,
## plus a small effect of age, sex and comorbidities
days_aug <- as.numeric(admitted_date - as.Date("2020-08-01"))
p_time   <- 0.70 * plogis((days_aug - 50) / 15)
p_dexa   <- plogis(qlogis(pmax(p_time, 1e-9)) +
                     0.01 * (age - 65) + 0.10 * sexM + 0.10 * (n_comorb - 1.5))
p_dexa[days_aug < 0] <- 0
dexamethason_base <- rbinom(n, 1, p_dexa)

## 6. Respiratory support at baseline (multinomial logit) -------------
eta_oxy <-  0.5 + 0.02 * (age - 65) + 0.2 * sexM + 0.15 * n_comorb
eta_inv <- -1.2 + 0.02 * (age - 65) + 0.4 * sexM + 0.20 * n_comorb
probs   <- cbind(1, exp(eta_oxy), exp(eta_inv))
probs   <- probs / rowSums(probs)
MV <- apply(probs, 1, function(p) sample(c("No oxygen", "Oxygen", "Invasive"), 1, prob = p))
MV <- factor(MV, levels = c("No oxygen", "Oxygen", "Invasive"))
MV_Oxygen   <- as.integer(MV == "Oxygen")
MV_Invasive <- as.integer(MV == "Invasive")

## 7. Event time: Weibull proportional hazards ------------------------
months_since_mar <- as.numeric(admitted_date - as.Date("2020-03-01")) / 30.44
day  <- as.numeric(admitted_date - as.Date("2020-03-01"))   # 0 = 1 March
wave <- 1.0 * exp(-0.5 *((day -  35)/20)^2) +  # 1st wave, peak ~5 Apr  
        1.0*exp(-0.5*((day - 245)/35)^2)       # 2nd wave, peak ~1 Nov
  
lp <- 0.06 * (age - 65) +
     -0.30 * sexM +
      0.15 * n_comorb +
      0.05 * n_comorb * months_since_mar +
      0.50 * MV_Oxygen +
      1.20 * MV_Invasive +
      0.2 * dexamethason_base * (MV == "oxygen") +
     -0.2 * dexamethason_base * MV_Oxygen +   
     -0.5 * dexamethason_base * MV_Invasive +
      0.50 * wave     

lambda <- 0.005   
k      <- 0.9
T_event <- (-log(runif(n)) / (lambda * exp(lp)))^(1 / k)

## 8. Status, survtime, date_event ------------------------------------
status     <- as.integer(T_event <= 28)
survtime   <- round(ifelse(status == 1, T_event, 28.001),3)
date_event <- ifelse(status == 1, as.Date(admitted_date + floor(T_event)), as.Date("2021-01-01"))

dat <- data.frame(id = 1:n, age, sexM, n_comorb, admitted_date, 
                  dexamethason_base, MV, MV_Oxygen, MV_Invasive,
                  date_event, status, survtime)

## 9. Save -------------------------------------------------------------
save(dat, file = "simulated_data.rda")

## Quick checks --------------------------------------------------------
str(dat)
mean(dat$status)                                                                 # 28-day mortality
round(tapply(dat$dexamethason_base, format(dat$admitted_date, "%m"), mean), 2)  # dexa by month
prop.table(table(dat$MV))


