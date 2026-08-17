# Forecasting Multivariado - ARIMAX
library(vars)
library(rminer)
library(forecast)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R")
source("../multi-utils.R")

# Configuração inicial do ficheiro
H <- 7
K <- 7

# Carregar os CSVs das quatro lojas
bt <- load_store("baltimore_clean.csv",    "Baltimore")
lc <- load_store("lancaster_clean.csv",    "Lancaster")
ph <- load_store("philadelphia_clean.csv", "Philadelphia")
rc <- load_store("richmond_clean.csv",     "Richmond")

YR_bt <- bt$YR; YR_lc <- lc$YR
YR_ph <- ph$YR; YR_rc <- rc$YR

L   <- length(bt$S)
LTR <- L - H

Y_bt <- bt$S[(LTR + 1):L]
Y_lc <- lc$S[(LTR + 1):L]
Y_ph <- ph$S[(LTR + 1):L]
Y_rc <- rc$S[(LTR + 1):L]

# Preparar variáveis exógenas para treino e teste
prep_exog <- function(df, LTR, H) {
  tourist <- as.numeric(df$TouristEvent == "Yes")
  promo   <- df$Pct_On_Sale
  empl    <- df$Num_Employees
  exotr   <- cbind(TouristEvent = tourist[1:LTR], Pct_On_Sale = promo[1:LTR], Num_Employees = empl[1:LTR])
  exots   <- cbind(TouristEvent = tourist[(LTR+1):(LTR+H)], Pct_On_Sale = promo[(LTR+1):(LTR+H)], Num_Employees = empl[(LTR+1):(LTR+H)])
  list(exotr = exotr, exots = exots)
}

# Fase 1 - Cenário A
cat("FASE 1 — CENÁRIO A: ARIMAX por loja (Customers + Sales)\n")

fit_arimax_A <- function(store, store_name, YR, Y_real, LTR, H, K) {
  cust  <- store$df$Num_Customers; sales <- store$df$Sales
  cdata <- cbind(Customers = cust[1:LTR], Sales = sales[1:LTR])
  mtr   <- ts(cdata, frequency = K)
  exog  <- prep_exog(store$df, LTR, H)
  arimax <- autoARIMAX(mtr, frequency = K, exogen = exog$exotr)
  FA     <- forecastARIMAX(arimax, h = H, exogen = exog$exots)
  Pred_cust <- FA[[1]]
  nmae <- mmetric(Y_real, Pred_cust, metric = "NMAE", val = YR)
  rmse <- mmetric(Y_real, Pred_cust, metric = "RMSE")
  r2   <- mmetric(Y_real, Pred_cust, metric = "R22")
  cat(sprintf("  %-15s NMAE: %6.2f%%  RMSE: %8.2f  R2: %.4f\n", store_name, nmae, rmse, r2))
  list(Pred_cust = Pred_cust, model = arimax, nmae = nmae, rmse = rmse, r2 = r2)
}

aa_bt <- fit_arimax_A(bt, "Baltimore",    YR_bt, Y_bt, LTR, H, K)
aa_lc <- fit_arimax_A(lc, "Lancaster",    YR_lc, Y_lc, LTR, H, K)
aa_ph <- fit_arimax_A(ph, "Philadelphia", YR_ph, Y_ph, LTR, H, K)
aa_rc <- fit_arimax_A(rc, "Richmond",     YR_rc, Y_rc, LTR, H, K)

plot_with_legend("ARIMAX-A", Y_bt, aa_bt$Pred_cust, "Baltimore",    YR_bt)
plot_with_legend("ARIMAX-A", Y_lc, aa_lc$Pred_cust, "Lancaster",    YR_lc)
plot_with_legend("ARIMAX-A", Y_ph, aa_ph$Pred_cust, "Philadelphia", YR_ph)
plot_with_legend("ARIMAX-A", Y_rc, aa_rc$Pred_cust, "Richmond",     YR_rc)

# Fase 1 - Cenário D
cat("FASE 1 — CENÁRIO D: ARIMAX univariado (Customers + exógenas)\n")

fit_arimax_D <- function(store, store_name, YR, Y_real, LTR, H, K) {
  cust <- store$df$Num_Customers
  ytr  <- ts(cust[1:LTR], frequency = K)
  exog <- prep_exog(store$df, LTR, H)
  model <- auto.arima(ytr, xreg = exog$exotr)
  fc <- forecast(model, h = H, xreg = exog$exots)
  Pred <- as.numeric(fc$mean)
  nmae <- mmetric(Y_real, Pred, metric = "NMAE", val = YR)
  rmse <- mmetric(Y_real, Pred, metric = "RMSE")
  r2   <- mmetric(Y_real, Pred, metric = "R22")
  cat(sprintf("  %-15s NMAE: %6.2f%%  RMSE: %8.2f  R2: %.4f\n", store_name, nmae, rmse, r2))
  list(Pred = Pred, model = model, nmae = nmae, rmse = rmse, r2 = r2)
}

ad_bt <- fit_arimax_D(bt, "Baltimore",    YR_bt, Y_bt, LTR, H, K)
ad_lc <- fit_arimax_D(lc, "Lancaster",    YR_lc, Y_lc, LTR, H, K)
ad_ph <- fit_arimax_D(ph, "Philadelphia", YR_ph, Y_ph, LTR, H, K)
ad_rc <- fit_arimax_D(rc, "Richmond",     YR_rc, Y_rc, LTR, H, K)

plot_with_legend("ARIMAX-D", Y_bt, ad_bt$Pred, "Baltimore",    YR_bt)
plot_with_legend("ARIMAX-D", Y_lc, ad_lc$Pred, "Lancaster",    YR_lc)
plot_with_legend("ARIMAX-D", Y_ph, ad_ph$Pred, "Philadelphia", YR_ph)
plot_with_legend("ARIMAX-D", Y_rc, ad_rc$Pred, "Richmond",     YR_rc)


# Fase 2 - com growing window
cat("FASE 2 — GROWING WINDOW (ARIMAX)\n")

RUNS   <- 20
stores <- list(bt, lc, ph, rc)
names_s <- c("Baltimore", "Lancaster", "Philadelphia", "Richmond")
YR_s   <- c(YR_bt, YR_lc, YR_ph, YR_rc)

# Growing window para o Cenário A
cat("GW Cenário A: ARIMAX por loja (Customers + Sales)\n")

gw_nmae_A <- matrix(NA, nrow = RUNS, ncol = 4, dimnames = list(NULL, names_s))
gw_r2_A   <- matrix(NA, nrow = RUNS, ncol = 4, dimnames = list(NULL, names_s))
gw_allY_A <- list(c(), c(), c(), c())
gw_allP_A <- list(c(), c(), c(), c())

for (run in 1:RUNS) {
  test_end   <- L - (RUNS - run) * H
  test_start <- test_end - H + 1
  train_end  <- test_start - 1
  cat(sprintf("  Iteração %2d/%d\n", run, RUNS))
  
  for (s in 1:4) {
    store <- stores[[s]]
    cust  <- store$df$Num_Customers; sales <- store$df$Sales
    Y_real <- cust[test_start:test_end]
    
    cdata <- cbind(Customers = cust[1:train_end], Sales = sales[1:train_end])
    mtr   <- ts(cdata, frequency = K)
    exog  <- prep_exog(store$df, train_end, H)
    
    pred <- tryCatch({
      arimax <- autoARIMAX(mtr, frequency = K, exogen = exog$exotr)
      FA     <- forecastARIMAX(arimax, h = H, exogen = exog$exots)
      FA[[1]]
    }, error = function(e) { rep(NA, H) })
    
    gw_nmae_A[run, s] <- mmetric(Y_real, pred, metric = "NMAE", val = YR_s[s])
    gw_r2_A[run, s]   <- mmetric(Y_real, pred, metric = "R22")
    gw_allY_A[[s]] <- c(gw_allY_A[[s]], Y_real)
    gw_allP_A[[s]] <- c(gw_allP_A[[s]], pred)
  }
}

cat("\nMedianas — GW Cenário A\n")
for (s in 1:4) cat(sprintf("  %-15s Mediana NMAE: %5.2f%%  Mediana R2: %.4f\n",
                           names_s[s], median(gw_nmae_A[,s], na.rm=T), median(gw_r2_A[,s], na.rm=T)))

for (s in 1:4) plot_gw_legend(
  gw_allY_A[[s]], list(gw_allP_A[[s]]),
  "ARIMAX-A", "blue",
  names_s[s], "ARIMAX-A")

# Growing window para o Cenário D
cat("GW Cenário D: ARIMAX univariado (Customers + exógenas)\n")

gw_nmae_D <- matrix(NA, nrow = RUNS, ncol = 4, dimnames = list(NULL, names_s))
gw_r2_D   <- matrix(NA, nrow = RUNS, ncol = 4, dimnames = list(NULL, names_s))
gw_allY_D <- list(c(), c(), c(), c())
gw_allP_D <- list(c(), c(), c(), c())

for (run in 1:RUNS) {
  test_end   <- L - (RUNS - run) * H
  test_start <- test_end - H + 1
  train_end  <- test_start - 1
  cat(sprintf("  Iteração %2d/%d\n", run, RUNS))
  
  for (s in 1:4) {
    store <- stores[[s]]
    cust  <- store$df$Num_Customers
    Y_real <- cust[test_start:test_end]
    
    ytr  <- ts(cust[1:train_end], frequency = K)
    exog <- prep_exog(store$df, train_end, H)
    
    pred <- tryCatch({
      model <- auto.arima(ytr, xreg = exog$exotr)
      fc    <- forecast(model, h = H, xreg = exog$exots)
      as.numeric(fc$mean)
    }, error = function(e) { rep(NA, H) })
    
    gw_nmae_D[run, s] <- mmetric(Y_real, pred, metric = "NMAE", val = YR_s[s])
    gw_r2_D[run, s]   <- mmetric(Y_real, pred, metric = "R22")
    gw_allY_D[[s]] <- c(gw_allY_D[[s]], Y_real)
    gw_allP_D[[s]] <- c(gw_allP_D[[s]], pred)
  }
}

cat("\nMedianas — GW Cenário D\n")
for (s in 1:4) cat(sprintf("  %-15s Mediana NMAE: %5.2f%%  Mediana R2: %.4f\n",
                           names_s[s], median(gw_nmae_D[,s], na.rm=T), median(gw_r2_D[,s], na.rm=T)))

for (s in 1:4) plot_gw_legend(
  gw_allY_D[[s]], list(gw_allP_D[[s]]),
  "ARIMAX-D", "darkgreen",
  names_s[s], "ARIMAX-D")

# Tabela final de comparação
cat("COMPARAÇÃO FINAL — Mediana NMAE e R2 (GW 20 iter.)\n")

cat(sprintf("  %-15s  %10s  %10s  %8s  %8s\n", "Loja", "ARIMAX-A", "ARIMAX-D", "R2-A", "R2-D"))
for (s in 1:4) cat(sprintf("  %-15s  %9.2f%%  %9.2f%%  %8.4f  %8.4f\n", names_s[s],
                           median(gw_nmae_A[,s], na.rm=T), median(gw_nmae_D[,s], na.rm=T),
                           median(gw_r2_A[,s], na.rm=T), median(gw_r2_D[,s], na.rm=T)))