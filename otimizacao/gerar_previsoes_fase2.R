# SCRIPT DE GERAÇÃO DE PREVISÕES — GROWING WINDOW (Para a Otimização)
# Modelo: Random Forest - Cenário C (Global - 4 Lojas)
# Lags: 1, 3, 7, 14, 21, 28
# Fase 2: guarda as 20 iterações do Growing Window

library(rminer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R")
source("../multi-utils.R")


# Config
H    <- 7
RUNS <- 20


# Carregar dados
cat("A carregar os datasets...\n")
bt <- load_store("baltimore_clean.csv",    "Baltimore")
lc <- load_store("lancaster_clean.csv",    "Lancaster")
ph <- load_store("philadelphia_clean.csv", "Philadelphia")
rc <- load_store("richmond_clean.csv",     "Richmond")

L   <- length(bt$S)
LTR <- L - H

stores  <- list(bt, lc, ph, rc)
names_s <- c("Baltimore", "Lancaster", "Philadelphia", "Richmond")


# Lags
lag_atual <- c(1, 3, 7, 14, 21, 28)

VINP_C_dyn <- list(
  list(lag_atual, 1, 1, 1),
  list(1, lag_atual, 1, 1),
  list(1, 1, lag_atual, 1),
  list(1, 1, 1, lag_atual)
)

# Exógenas
prep_exog_full <- function(df) {
  cbind(
    TouristEvent  = as.numeric(df$TouristEvent == "Yes"),
    Pct_On_Sale   = df$Pct_On_Sale / 100,
    Num_Employees = df$Num_Employees / max(df$Num_Employees[1:LTR])
  )
}

prep_exog_all <- function(train_end, H) {
  exf_bt <- prep_exog_full(bt$df)
  exf_lc <- prep_exog_full(lc$df)
  exf_ph <- prep_exog_full(ph$df)
  exf_rc <- prep_exog_full(rc$df)
  exotr <- cbind(
    BT_Tourist=exf_bt[1:train_end,1], BT_Promo=exf_bt[1:train_end,2], BT_Empl=exf_bt[1:train_end,3],
    LC_Tourist=exf_lc[1:train_end,1], LC_Promo=exf_lc[1:train_end,2], LC_Empl=exf_lc[1:train_end,3],
    PH_Tourist=exf_ph[1:train_end,1], PH_Promo=exf_ph[1:train_end,2], PH_Empl=exf_ph[1:train_end,3],
    RC_Tourist=exf_rc[1:train_end,1], RC_Promo=exf_rc[1:train_end,2], RC_Empl=exf_rc[1:train_end,3]
  )
  exots <- cbind(
    BT_Tourist=exf_bt[(train_end+1):(train_end+H),1], BT_Promo=exf_bt[(train_end+1):(train_end+H),2], BT_Empl=exf_bt[(train_end+1):(train_end+H),3],
    LC_Tourist=exf_lc[(train_end+1):(train_end+H),1], LC_Promo=exf_lc[(train_end+1):(train_end+H),2], LC_Empl=exf_lc[(train_end+1):(train_end+H),3],
    PH_Tourist=exf_ph[(train_end+1):(train_end+H),1], PH_Promo=exf_ph[(train_end+1):(train_end+H),2], PH_Empl=exf_ph[(train_end+1):(train_end+H),3],
    RC_Tourist=exf_rc[(train_end+1):(train_end+H),1], RC_Promo=exf_rc[(train_end+1):(train_end+H),2], RC_Empl=exf_rc[(train_end+1):(train_end+H),3]
  )
  list(exotr = exotr, exots = exots)
}


# Função SAVE ajustada para incluir valores reais
save_forecasts_gw <- function(Pred_final, Real_final, store_name, test_dates, run) {
  
  Pred_final <- round(pmax(Pred_final, 0))  # garante não-negativos e arredonda
  
  out <- data.frame(
    Run                 = run,
    Store               = store_name,
    Date                = test_dates,
    Weekday             = weekdays(test_dates),
    Predicted_Customers = Pred_final,
    Actual_Customers    = Real_final        
  )
  
  invisible(out)
}


# Growing Window
cat("\nA correr Growing Window RF-C (20 iterações)...\n")

# Uma lista por loja para acumular as linhas
acum <- list(
  Baltimore    = data.frame(),
  Lancaster    = data.frame(),
  Philadelphia = data.frame(),
  Richmond     = data.frame()
)

for (run in 1:RUNS) {
  test_end   <- L - (RUNS - run) * H
  test_start <- test_end - H + 1
  train_end  <- test_start - 1
  
  cat(sprintf("  Iteração %2d/%d | treino: 1:%d | teste: %d:%d\n",
              run, RUNS, train_end, test_start, test_end))
  
  mtr_gw <- cbind(
    Baltimore    = bt$S[1:train_end],
    Lancaster    = lc$S[1:train_end],
    Philadelphia = ph$S[1:train_end],
    Richmond     = rc$S[1:train_end]
  )
  
  exog_gw <- prep_exog_all(train_end, H)
  
  res <- tryCatch({
    set.seed(42)
    MNN <- mfit(mtr_gw, "randomForest", VINP_C_dyn, exogen = exog_gw$exotr)
    lforecastm(MNN, h = H, exogen = exog_gw$exots)
  }, error = function(e) {
    cat(sprintf("    ERRO na iteração %d: %s\n", run, e$message))
    list(rep(NA, H), rep(NA, H), rep(NA, H), rep(NA, H))
  })
  
  # datas reais do período de teste (já existem no df, não são "futuras")
  test_dates <- stores[[1]]$df$Date[test_start:test_end]
  
  for (s in 1:4) {
    # extrair valores reais
    real_vals <- stores[[s]]$S[test_start:test_end]
    # Passamos os valores reais (real_vals) como o segundo argumento
    out <- save_forecasts_gw(res[[s]], real_vals, names_s[s], test_dates, run)
    acum[[names_s[s]]] <- rbind(acum[[names_s[s]]], out)
  }
}


# Exportar um CSV por loja
cat("\nA exportar CSVs...\n")
for (nome in names_s) {
  output_csv <- paste0("previsoes_gw_", tolower(nome), ".csv")
  write.csv(acum[[nome]], output_csv, row.names = FALSE)
  cat("Previsões guardadas em:", output_csv, "\n")
}

cat("\nSucesso! CSVs gerados para as 4 lojas (20 iterações cada).\n")