# Forecasting Multivariado - VAR
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

# Valores reais de Customers da semana de teste 
Y_bt <- bt$S[(LTR + 1):L]
Y_lc <- lc$S[(LTR + 1):L]
Y_ph <- ph$S[(LTR + 1):L]
Y_rc <- rc$S[(LTR + 1):L]

# Preparar variáveis exógenas
prep_exog <- function(df, LTR, H) {
  tourist <- as.numeric(df$TouristEvent == "Yes")
  promo   <- df$Pct_On_Sale
  empl    <- df$Num_Employees
  exotr   <- cbind(TouristEvent  = tourist[1:LTR],
                   Pct_On_Sale   = promo[1:LTR],
                   Num_Employees = empl[1:LTR])
  exots   <- cbind(TouristEvent  = tourist[(LTR + 1):(LTR + H)],
                   Pct_On_Sale   = promo[(LTR + 1):(LTR + H)],
                   Num_Employees = empl[(LTR + 1):(LTR + H)])
  list(exotr = exotr, exots = exots)
}

# Fase 1 - Cenário A
cat("FASE 1 — CENÁRIO A: VAR por loja (Customers + Sales endógenas)\n")

fit_var_A <- function(store, store_name, YR, Y_real, LTR, H, K) {
  cust  <- store$df$Num_Customers
  sales <- store$df$Sales
  cdata <- cbind(Customers = cust[1:LTR], Sales = sales[1:LTR])
  mtr   <- ts(cdata, frequency = K)
  exog <- prep_exog(store$df, LTR, H)
  exogen <<- exog$exotr
  mvar <- autoVAR(mtr, season = K, exogen = exogen)
  FV   <- forecastVAR(mvar, h = H, exogen = exog$exots)
  Pred_cust  <- FV[[1]]
  nmae <- mmetric(Y_real, Pred_cust, metric = "NMAE", val = YR)
  rmse <- mmetric(Y_real, Pred_cust, metric = "RMSE")
  r2   <- mmetric(Y_real, Pred_cust, metric = "R22")
  cat(sprintf("  %-15s NMAE: %6.2f%%  RMSE: %8.2f  R2: %.4f\n", store_name, nmae, rmse, r2))
  list(Pred_cust = Pred_cust, model = mvar, nmae = nmae, rmse = rmse, r2 = r2)
}

va_bt <- fit_var_A(bt, "Baltimore",    YR_bt, Y_bt, LTR, H, K)
va_lc <- fit_var_A(lc, "Lancaster",    YR_lc, Y_lc, LTR, H, K)
va_ph <- fit_var_A(ph, "Philadelphia", YR_ph, Y_ph, LTR, H, K)
va_rc <- fit_var_A(rc, "Richmond",     YR_rc, Y_rc, LTR, H, K)

plot_with_legend("VAR-A Cust+Sales", Y_bt, va_bt$Pred_cust, "Baltimore",    YR_bt)
plot_with_legend("VAR-A Cust+Sales", Y_lc, va_lc$Pred_cust, "Lancaster",    YR_lc)
plot_with_legend("VAR-A Cust+Sales", Y_ph, va_ph$Pred_cust, "Philadelphia", YR_ph)
plot_with_legend("VAR-A Cust+Sales", Y_rc, va_rc$Pred_cust, "Richmond",     YR_rc)

# Fase 1 - Cenário C
cat("FASE 1 — CENÁRIO C: VAR global (4 lojas Customers + exógenas)\n")

# Agregar variáveis exógenas das quatro lojas
prep_exog_agg <- function(train_end, H) {
  ex_bt <- prep_exog(bt$df, train_end, H)
  ex_lc <- prep_exog(lc$df, train_end, H)
  ex_ph <- prep_exog(ph$df, train_end, H)
  ex_rc <- prep_exog(rc$df, train_end, H)
  agg <- function(part) {
    cbind(
      Tourist = rowMeans(cbind(ex_bt[[part]][,1], ex_lc[[part]][,1], ex_ph[[part]][,1], ex_rc[[part]][,1])),
      Promo   = apply(cbind(ex_bt[[part]][,2], ex_lc[[part]][,2], ex_ph[[part]][,2], ex_rc[[part]][,2]), 1, median),
      Empl    = apply(cbind(ex_bt[[part]][,3], ex_lc[[part]][,3], ex_ph[[part]][,3], ex_rc[[part]][,3]), 1, median)
    )
  }
  list(exotr = agg("exotr"), exots = agg("exots"))
}

cdata_4 <- cbind(Baltimore = bt$S[1:LTR], Lancaster = lc$S[1:LTR],
                 Philadelphia = ph$S[1:LTR], Richmond = rc$S[1:LTR])
mtr_4 <- ts(cdata_4, frequency = K)

exog_agg <- prep_exog_agg(LTR, H)
exogen <- exog_agg$exotr

mvar_4 <- autoVAR(mtr_4, season = K, exogen = exogen)
FV_4 <- forecastVAR(mvar_4, h = H, exogen = exog_agg$exots)

Pred_C <- list(FV_4[[1]], FV_4[[2]], FV_4[[3]], FV_4[[4]])
Y_s    <- list(Y_bt, Y_lc, Y_ph, Y_rc)
names_s <- c("Baltimore", "Lancaster", "Philadelphia", "Richmond")
YR_s   <- c(YR_bt, YR_lc, YR_ph, YR_rc)

cat("\nRESULTADOS CENÁRIO C - VAR Global\n")
for (s in 1:4) {
  cat(sprintf("  %-15s NMAE: %6.2f%%  RMSE: %8.2f  R2: %.4f\n", names_s[s],
              mmetric(Y_s[[s]], Pred_C[[s]], metric = "NMAE", val = YR_s[s]),
              mmetric(Y_s[[s]], Pred_C[[s]], metric = "RMSE"),
              mmetric(Y_s[[s]], Pred_C[[s]], metric = "R22")))
}

for (s in 1:4) plot_with_legend("VAR-C Global", Y_s[[s]], Pred_C[[s]], names_s[s], YR_s[s])

# Fase 2 - com growing window
cat("FASE 2 — GROWING WINDOW (VAR)\n")

RUNS   <- 20
stores <- list(bt, lc, ph, rc)

# Growing window para o Cenário A
cat("GW Cenário A: VAR por loja (Customers + Sales)\n")

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
    cust  <- store$df$Num_Customers
    sales <- store$df$Sales
    Y_real <- cust[test_start:test_end]
    
    cdata <- cbind(Customers = cust[1:train_end], Sales = sales[1:train_end])
    mtr   <- ts(cdata, frequency = K)
    exog  <- prep_exog(store$df, train_end, H)
    exogen <<- exog$exotr
    
    pred <- tryCatch({
      mvar <- autoVAR(mtr, season = K, exogen = exogen)
      FV   <- forecastVAR(mvar, h = H, exogen = exog$exots)
      FV[[1]]
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
  "VAR-A", "blue",
  names_s[s], "VAR-A")

# Growing window para o Cenário C
cat("GW Cenário C: VAR global (4 lojas)\n")

gw_nmae_C <- matrix(NA, nrow = RUNS, ncol = 4, dimnames = list(NULL, names_s))
gw_r2_C   <- matrix(NA, nrow = RUNS, ncol = 4, dimnames = list(NULL, names_s))
gw_allY_C <- list(c(), c(), c(), c())
gw_allP_C <- list(c(), c(), c(), c())

for (run in 1:RUNS) {
  test_end   <- L - (RUNS - run) * H
  test_start <- test_end - H + 1
  train_end  <- test_start - 1
  
  cat(sprintf("  Iteração %2d/%d\n", run, RUNS))
  
  cdata_gw <- cbind(Baltimore = bt$S[1:train_end], Lancaster = lc$S[1:train_end],
                    Philadelphia = ph$S[1:train_end], Richmond = rc$S[1:train_end])
  mtr_gw <- ts(cdata_gw, frequency = K)
  
  exog_agg <- prep_exog_agg(train_end, H)
  exogen <<- exog_agg$exotr
  
  res_all <- tryCatch({
    mvar_gw <- autoVAR(mtr_gw, season = K, exogen = exogen)
    forecastVAR(mvar_gw, h = H, exogen = exog_agg$exots)
  }, error = function(e) { list(rep(NA,H), rep(NA,H), rep(NA,H), rep(NA,H)) })
  
  for (s in 1:4) {
    Y_real <- stores[[s]]$S[test_start:test_end]
    gw_nmae_C[run, s] <- mmetric(Y_real, res_all[[s]], metric = "NMAE", val = YR_s[s])
    gw_r2_C[run, s]   <- mmetric(Y_real, res_all[[s]], metric = "R22")
    gw_allY_C[[s]] <- c(gw_allY_C[[s]], Y_real)
    gw_allP_C[[s]] <- c(gw_allP_C[[s]], res_all[[s]])
  }
}

cat("\nMedianas — GW Cenário C\n")
for (s in 1:4) cat(sprintf("  %-15s Mediana NMAE: %5.2f%%  Mediana R2: %.4f\n",
  names_s[s], median(gw_nmae_C[,s], na.rm=T), median(gw_r2_C[,s], na.rm=T)))

for (s in 1:4) plot_gw_legend(
  gw_allY_C[[s]], list(gw_allP_C[[s]]),
  "VAR-C", "red",
  names_s[s], "VAR-C")

# Tabela final de comparação
cat("COMPARAÇÃO FINAL — Mediana NMAE e R2 (GW 20 iter.)\n")

cat(sprintf("  %-15s  %8s  %8s  %8s  %8s\n", "Loja", "VAR-A", "VAR-C", "R2-A", "R2-C"))
for (s in 1:4) cat(sprintf("  %-15s  %7.2f%%  %7.2f%%  %8.4f  %8.4f\n", names_s[s],
                           median(gw_nmae_A[,s], na.rm=T), median(gw_nmae_C[,s], na.rm=T),
                           median(gw_r2_A[,s], na.rm=T), median(gw_r2_C[,s], na.rm=T)))