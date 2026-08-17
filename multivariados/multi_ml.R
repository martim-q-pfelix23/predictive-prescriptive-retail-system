# Forecasting Multivariado - Machine Learning
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

# Preparar variáveis exógenas
prep_exog_full <- function(df) {
  cbind(TouristEvent  = as.numeric(df$TouristEvent == "Yes"),
        Pct_On_Sale   = df$Pct_On_Sale / 100,
        Num_Employees = df$Num_Employees / max(df$Num_Employees[1:LTR]))
}

# Configuração dos lags para cada cenário
VINP_A <- vector("list", length = 2)
VINP_A[[1]] <- list(1:7, 1:7)
VINP_A[[2]] <- list(1:7, 1:7)

VINP_B <- vector("list", length = 1)
VINP_B[[1]] <- list(1:7)

VINP_C <- vector("list", length = 4)
VINP_C[[1]] <- list(1:7, 1, 1, 1)
VINP_C[[2]] <- list(1, 1:7, 1, 1)
VINP_C[[3]] <- list(1, 1, 1:7, 1)
VINP_C[[4]] <- list(1, 1, 1, 1:7)

stores  <- list(bt, lc, ph, rc)
names_s <- c("Baltimore", "Lancaster", "Philadelphia", "Richmond")
YR_s    <- c(YR_bt, YR_lc, YR_ph, YR_rc)
Y_s     <- list(Y_bt, Y_lc, Y_ph, Y_rc)

# Fase 1 - Cenário A
cat("FASE 1 — CENÁRIO A: ML por loja (Customers + Sales)\n")
fit_ml_A <- function(store, store_name, YR, Y_real, LTR, H, VINP, ml_model) {
  cust  <- store$df$Num_Customers; sales <- store$df$Sales
  mtr <- cbind(Customers = cust[1:LTR], Sales = sales[1:LTR])
  exog_all <- prep_exog_full(store$df)
  exotr <- exog_all[1:LTR, ]; exots <- exog_all[(LTR+1):(LTR+H), ]
  set.seed(42)
  MNN  <- mfit(mtr, ml_model, VINP, exogen = exotr)
  Pred <- lforecastm(MNN, h = H, exogen = exots)
  Pred_cust <- Pred[[1]]
  nmae <- mmetric(Y_real, Pred_cust, metric = "NMAE", val = YR)
  rmse <- mmetric(Y_real, Pred_cust, metric = "RMSE")
  r2   <- mmetric(Y_real, Pred_cust, metric = "R22")
  cat(sprintf("  %-15s NMAE: %6.2f%%  RMSE: %8.2f  R2: %.4f\n", store_name, nmae, rmse, r2))
  list(Pred_cust = Pred_cust, nmae = nmae, rmse = rmse, r2 = r2)
}

cat("MLPE (Cenário A):\n")
for (i in 1:4) fit_ml_A(stores[[i]], names_s[i], YR_s[i], Y_s[[i]], LTR, H, VINP_A, "mlpe")
cat("\nRF (Cenário A):\n")
for (i in 1:4) fit_ml_A(stores[[i]], names_s[i], YR_s[i], Y_s[[i]], LTR, H, VINP_A, "randomForest")


# Fase 1 - Cenário B
cat("FASE 1 — CENÁRIO B: ML por loja (só Customers)\n")
fit_ml_B <- function(store, store_name, YR, Y_real, LTR, H, VINP, ml_model) {
  cust <- store$df$Num_Customers
  mtr <- cbind(Customers = cust[1:LTR])
  exog_all <- prep_exog_full(store$df)
  exotr <- exog_all[1:LTR, ]; exots <- exog_all[(LTR+1):(LTR+H), ]
  set.seed(42)
  MNN  <- mfit(mtr, ml_model, VINP, exogen = exotr)
  Pred <- lforecastm(MNN, h = H, exogen = exots)
  Pred_cust <- Pred[[1]]
  nmae <- mmetric(Y_real, Pred_cust, metric = "NMAE", val = YR)
  rmse <- mmetric(Y_real, Pred_cust, metric = "RMSE")
  r2   <- mmetric(Y_real, Pred_cust, metric = "R22")
  cat(sprintf("  %-15s NMAE: %6.2f%%  RMSE: %8.2f  R2: %.4f\n", store_name, nmae, rmse, r2))
  list(Pred_cust = Pred_cust, nmae = nmae, rmse = rmse, r2 = r2)
}

cat("MLPE (Cenário B):\n")
for (i in 1:4) fit_ml_B(stores[[i]], names_s[i], YR_s[i], Y_s[[i]], LTR, H, VINP_B, "mlpe")
cat("\nRF (Cenário B):\n")
for (i in 1:4) fit_ml_B(stores[[i]], names_s[i], YR_s[i], Y_s[[i]], LTR, H, VINP_B, "randomForest")


# Fase 1 - Cenário C
cat("FASE 1 — CENÁRIO C: ML global (4 lojas)\n")
mtr_4 <- cbind(Baltimore = bt$S[1:LTR], Lancaster = lc$S[1:LTR],
               Philadelphia = ph$S[1:LTR], Richmond = rc$S[1:LTR])

# Preparar variáveis exógenas das quatro lojas
prep_exog_all <- function(train_end, H) {
  exf_bt <- prep_exog_full(bt$df); exf_lc <- prep_exog_full(lc$df)
  exf_ph <- prep_exog_full(ph$df); exf_rc <- prep_exog_full(rc$df)
  exotr <- cbind(BT_Tourist=exf_bt[1:train_end,1], BT_Promo=exf_bt[1:train_end,2], BT_Empl=exf_bt[1:train_end,3],
                 LC_Tourist=exf_lc[1:train_end,1], LC_Promo=exf_lc[1:train_end,2], LC_Empl=exf_lc[1:train_end,3],
                 PH_Tourist=exf_ph[1:train_end,1], PH_Promo=exf_ph[1:train_end,2], PH_Empl=exf_ph[1:train_end,3],
                 RC_Tourist=exf_rc[1:train_end,1], RC_Promo=exf_rc[1:train_end,2], RC_Empl=exf_rc[1:train_end,3])
  exots <- cbind(BT_Tourist=exf_bt[(train_end+1):(train_end+H),1], BT_Promo=exf_bt[(train_end+1):(train_end+H),2], BT_Empl=exf_bt[(train_end+1):(train_end+H),3],
                 LC_Tourist=exf_lc[(train_end+1):(train_end+H),1], LC_Promo=exf_lc[(train_end+1):(train_end+H),2], LC_Empl=exf_lc[(train_end+1):(train_end+H),3],
                 PH_Tourist=exf_ph[(train_end+1):(train_end+H),1], PH_Promo=exf_ph[(train_end+1):(train_end+H),2], PH_Empl=exf_ph[(train_end+1):(train_end+H),3],
                 RC_Tourist=exf_rc[(train_end+1):(train_end+H),1], RC_Promo=exf_rc[(train_end+1):(train_end+H),2], RC_Empl=exf_rc[(train_end+1):(train_end+H),3])
  list(exotr = exotr, exots = exots)
}

exog_C <- prep_exog_all(LTR, H)

for (ml_model in c("mlpe", "randomForest")) {
  cat(sprintf("%s (Cenário C):\n", toupper(ml_model)))
  set.seed(42)
  MNN <- mfit(mtr_4, ml_model, VINP_C, exogen = exog_C$exotr)
  Pred <- lforecastm(MNN, h = H, exogen = exog_C$exots)
  for (s in 1:4) {
    nmae <- mmetric(Y_s[[s]], Pred[[s]], metric = "NMAE", val = YR_s[s])
    r2   <- mmetric(Y_s[[s]], Pred[[s]], metric = "R22")
    cat(sprintf("  %-15s NMAE: %6.2f%%  R2: %.4f\n", names_s[s], nmae, r2))
  }
  cat("\n")
}


# Gráficos da Fase 1
cat("GRÁFICOS — FASE 1\n")

# Gráficos do Cenário A
cat("MLPE (Gráficos Cenário A):\n")
mlpe_a_bt <- fit_ml_A(bt, "Baltimore",    YR_bt, Y_bt, LTR, H, VINP_A, "mlpe")
mlpe_a_lc <- fit_ml_A(lc, "Lancaster",    YR_lc, Y_lc, LTR, H, VINP_A, "mlpe")
mlpe_a_ph <- fit_ml_A(ph, "Philadelphia", YR_ph, Y_ph, LTR, H, VINP_A, "mlpe")
mlpe_a_rc <- fit_ml_A(rc, "Richmond",     YR_rc, Y_rc, LTR, H, VINP_A, "mlpe")

plot_with_legend("MLPE-A Cust+Sales", Y_bt, mlpe_a_bt$Pred_cust, "Baltimore",    YR_bt)
plot_with_legend("MLPE-A Cust+Sales", Y_lc, mlpe_a_lc$Pred_cust, "Lancaster",    YR_lc)
plot_with_legend("MLPE-A Cust+Sales", Y_ph, mlpe_a_ph$Pred_cust, "Philadelphia", YR_ph)
plot_with_legend("MLPE-A Cust+Sales", Y_rc, mlpe_a_rc$Pred_cust, "Richmond",     YR_rc)


cat("\nRF (Gráficos Cenário A):\n")
rf_a_bt <- fit_ml_A(bt, "Baltimore",    YR_bt, Y_bt, LTR, H, VINP_A, "randomForest")
rf_a_lc <- fit_ml_A(lc, "Lancaster",    YR_lc, Y_lc, LTR, H, VINP_A, "randomForest")
rf_a_ph <- fit_ml_A(ph, "Philadelphia", YR_ph, Y_ph, LTR, H, VINP_A, "randomForest")
rf_a_rc <- fit_ml_A(rc, "Richmond",     YR_rc, Y_rc, LTR, H, VINP_A, "randomForest")

plot_with_legend("RF-A Cust+Sales", Y_bt, rf_a_bt$Pred_cust, "Baltimore",    YR_bt)
plot_with_legend("RF-A Cust+Sales", Y_lc, rf_a_lc$Pred_cust, "Lancaster",    YR_lc)
plot_with_legend("RF-A Cust+Sales", Y_ph, rf_a_ph$Pred_cust, "Philadelphia", YR_ph)
plot_with_legend("RF-A Cust+Sales", Y_rc, rf_a_rc$Pred_cust, "Richmond",     YR_rc)


# Gráficos do Cenário B
cat("\nMLPE (Gráficos Cenário B):\n")
mlpe_b_bt <- fit_ml_B(bt, "Baltimore",    YR_bt, Y_bt, LTR, H, VINP_B, "mlpe")
mlpe_b_lc <- fit_ml_B(lc, "Lancaster",    YR_lc, Y_lc, LTR, H, VINP_B, "mlpe")
mlpe_b_ph <- fit_ml_B(ph, "Philadelphia", YR_ph, Y_ph, LTR, H, VINP_B, "mlpe")
mlpe_b_rc <- fit_ml_B(rc, "Richmond",     YR_rc, Y_rc, LTR, H, VINP_B, "mlpe")

plot_with_legend("MLPE-B", Y_bt, mlpe_b_bt$Pred_cust, "Baltimore",    YR_bt)
plot_with_legend("MLPE-B", Y_lc, mlpe_b_lc$Pred_cust, "Lancaster",    YR_lc)
plot_with_legend("MLPE-B", Y_ph, mlpe_b_ph$Pred_cust, "Philadelphia", YR_ph)
plot_with_legend("MLPE-B", Y_rc, mlpe_b_rc$Pred_cust, "Richmond",     YR_rc)


cat("\nRF (Gráficos Cenário B):\n")
rf_b_bt <- fit_ml_B(bt, "Baltimore",    YR_bt, Y_bt, LTR, H, VINP_B, "randomForest")
rf_b_lc <- fit_ml_B(lc, "Lancaster",    YR_lc, Y_lc, LTR, H, VINP_B, "randomForest")
rf_b_ph <- fit_ml_B(ph, "Philadelphia", YR_ph, Y_ph, LTR, H, VINP_B, "randomForest")
rf_b_rc <- fit_ml_B(rc, "Richmond",     YR_rc, Y_rc, LTR, H, VINP_B, "randomForest")

plot_with_legend("RF-B", Y_bt, rf_b_bt$Pred_cust, "Baltimore",    YR_bt)
plot_with_legend("RF-B", Y_lc, rf_b_lc$Pred_cust, "Lancaster",    YR_lc)
plot_with_legend("RF-B", Y_ph, rf_b_ph$Pred_cust, "Philadelphia", YR_ph)
plot_with_legend("RF-B", Y_rc, rf_b_rc$Pred_cust, "Richmond",     YR_rc)



# Gráficos do Cenário C
cat("\nMLPE (Gráficos Cenário C):\n")
set.seed(42)
MNN_mlpe_C <- mfit(mtr_4, "mlpe", VINP_C, exogen = exog_C$exotr)
Pred_mlpe_C <- lforecastm(MNN_mlpe_C, h = H, exogen = exog_C$exots)

for (s in 1:4) plot_with_legend("MLPE-C", Y_s[[s]], Pred_mlpe_C[[s]], names_s[s], YR_s[s])


cat("\nRF (Gráficos Cenário C):\n")
set.seed(42)
MNN_rf_C <- mfit(mtr_4, "randomForest", VINP_C, exogen = exog_C$exotr)
Pred_rf_C <- lforecastm(MNN_rf_C, h = H, exogen = exog_C$exots)

for (s in 1:4) plot_with_legend("RF-C", Y_s[[s]], Pred_rf_C[[s]], names_s[s], YR_s[s])

# Fase 2 - com growing window
cat("FASE 2 — GROWING WINDOW (ML)\n")
RUNS <- 20

# Função para growing window dos cenários por loja
gw_ml_loja <- function(cenario, VINP, endo_cols, ml_model) {
  label <- sprintf("%s-%s", ml_model, cenario)
  cat(sprintf("  A correr GW %s...\n", label))
  
  nmae_mat <- matrix(NA, nrow = RUNS, ncol = 4, dimnames = list(NULL, names_s))
  r2_mat   <- matrix(NA, nrow = RUNS, ncol = 4, dimnames = list(NULL, names_s))
  allY <- list(c(), c(), c(), c())
  allP <- list(c(), c(), c(), c())
  
  for (run in 1:RUNS) {
    test_end   <- L - (RUNS - run) * H
    test_start <- test_end - H + 1
    train_end  <- test_start - 1
    
    for (s in 1:4) {
      store <- stores[[s]]
      cust  <- store$df$Num_Customers
      Y_real <- cust[test_start:test_end]
      
      mtr <- NULL
      if ("Customers" %in% endo_cols) mtr <- cbind(mtr, Customers = cust[1:train_end])
      if ("Sales" %in% endo_cols) mtr <- cbind(mtr, Sales = store$df$Sales[1:train_end])
      
      exog_all <- prep_exog_full(store$df)
      exotr <- exog_all[1:train_end, ]
      exots <- exog_all[(train_end+1):(train_end+H), ]
      
      target_col <- which(colnames(mtr) == "Customers")
      
      pred <- tryCatch({
        set.seed(42)
        MNN  <- mfit(mtr, ml_model, VINP, exogen = exotr)
        Pred <- lforecastm(MNN, h = H, exogen = exots)
        Pred[[target_col]]
      }, error = function(e) { rep(NA, H) })
      
      nmae_mat[run, s] <- mmetric(Y_real, pred, metric = "NMAE", val = YR_s[s])
      r2_mat[run, s]   <- mmetric(Y_real, pred, metric = "R22")
      allY[[s]] <- c(allY[[s]], Y_real)
      allP[[s]] <- c(allP[[s]], pred)
    }
  }
  
  list(nmae = nmae_mat, r2 = r2_mat, allY = allY, allP = allP, label = label)
}

# Função para growing window do cenário global
gw_ml_global <- function(ml_model, VINP) {
  label <- sprintf("%s-C", ml_model)
  cat(sprintf("  A correr GW %s...\n", label))
  
  nmae_mat <- matrix(NA, nrow = RUNS, ncol = 4, dimnames = list(NULL, names_s))
  r2_mat   <- matrix(NA, nrow = RUNS, ncol = 4, dimnames = list(NULL, names_s))
  allY <- list(c(), c(), c(), c())
  allP <- list(c(), c(), c(), c())
  
  for (run in 1:RUNS) {
    test_end   <- L - (RUNS - run) * H
    test_start <- test_end - H + 1
    train_end  <- test_start - 1
    
    mtr_gw <- cbind(Baltimore = bt$S[1:train_end], Lancaster = lc$S[1:train_end],
                    Philadelphia = ph$S[1:train_end], Richmond = rc$S[1:train_end])
    
    exog_gw <- prep_exog_all(train_end, H)
    
    res <- tryCatch({
      set.seed(42)
      MNN  <- mfit(mtr_gw, ml_model, VINP, exogen = exog_gw$exotr)
      lforecastm(MNN, h = H, exogen = exog_gw$exots)
    }, error = function(e) { list(rep(NA,H), rep(NA,H), rep(NA,H), rep(NA,H)) })
    
    for (s in 1:4) {
      Y_real <- stores[[s]]$S[test_start:test_end]
      nmae_mat[run, s] <- mmetric(Y_real, res[[s]], metric = "NMAE", val = YR_s[s])
      r2_mat[run, s]   <- mmetric(Y_real, res[[s]], metric = "R22")
      allY[[s]] <- c(allY[[s]], Y_real)
      allP[[s]] <- c(allP[[s]], res[[s]])
    }
  }
  
  list(nmae = nmae_mat, r2 = r2_mat, allY = allY, allP = allP, label = label)
}


cat("A correr Growing Window para todos os lags em análise...\n")

# Configurações de lags testadas
configuracoes_lags <- list(
  "Lag_1_7"          = c(1, 7),
  "Lag_1_7_28"       = c(1, 7, 28),
  "Lag_1_2_3_7"      = c(1, 2, 3, 7),
  "Lag_1_3_7_14_21_28" = c(1, 3, 7, 14, 21, 28)
)

# Guardar os resultados de todas as combinações
all_gw <- list()

for (nome_lag in names(configuracoes_lags)) {
  lag_atual <- configuracoes_lags[[nome_lag]]
  
  cat(sprintf("\n\n>>> A TESTAR CONFIGURAÇÃO: %s <<<\n", nome_lag))
  
  # VINP dinâmico para o Cenário A
  VINP_A_dyn <- list(lag_atual, lag_atual)
  
  # VINP dinâmico para o Cenário B
  VINP_B_dyn <- list(lag_atual)
  
  # VINP dinâmico para o Cenário C
  VINP_C_dyn <- list(
    list(lag_atual, 1, 1, 1),
    list(1, lag_atual, 1, 1),
    list(1, 1, lag_atual, 1),
    list(1, 1, 1, lag_atual)
  )
  
  # Correr Cenário A
  gw_mlpe_A <- gw_ml_loja(sprintf("A-%s", nome_lag), VINP_A_dyn, c("Customers","Sales"), "mlpe")
  gw_rf_A   <- gw_ml_loja(sprintf("A-%s", nome_lag), VINP_A_dyn, c("Customers","Sales"), "randomForest")
  
  # Correr Cenário B
  gw_mlpe_B <- gw_ml_loja(sprintf("B-%s", nome_lag), VINP_B_dyn, "Customers", "mlpe")
  gw_rf_B   <- gw_ml_loja(sprintf("B-%s", nome_lag), VINP_B_dyn, "Customers", "randomForest")
  
  # Correr Cenário C
  gw_mlpe_C <- gw_ml_global("mlpe", VINP_C_dyn)
  gw_mlpe_C$label <- sprintf("mlpe-C-%s", nome_lag)
  
  gw_rf_C   <- gw_ml_global("randomForest", VINP_C_dyn)
  gw_rf_C$label <- sprintf("RF-C-%s", nome_lag)
  
  # Guardar resultados da configuração atual
  all_gw <- c(all_gw, list(gw_mlpe_A, gw_rf_A, gw_mlpe_B, gw_rf_B, gw_mlpe_C, gw_rf_C))
}


# Resultados globais
cat("RESULTADOS - Mediana NMAE e R2 (GW 20 iter.) - TODOS OS LAGS\n")


res_df <- data.frame()
for (g in all_gw) {
  for (s in 1:4) {
    nmae_med = median(g$nmae[,s], na.rm = TRUE)
    r2_med   = median(g$r2[,s], na.rm = TRUE)
    res_df <- rbind(res_df, data.frame(
      Configuracao = g$label,
      Loja = names_s[s],
      NMAE = round(nmae_med, 2),
      R2 = round(r2_med, 4)
    ))
  }
}

print(res_df)

# Exportar gráficos da growing window para PDF
cat("A gerar PDF com todos os gráficos comparativos...\n")

pdf("Graficos_GW_Multivariados.pdf", width = 12, height = 9)

cenarios <- c("A", "B", "C")

for (cenario in cenarios) {
  for (nome_lag in names(configuracoes_lags)) {
    
    lbl_mlpe <- sprintf("mlpe-%s-%s", cenario, nome_lag) 
    lbl_rf   <- sprintf("RF-%s-%s", cenario, nome_lag)   
    
    gw_mlpe_atual <- NULL
    gw_rf_atual   <- NULL
    
    for (g in all_gw) {
      if (grepl(sprintf("(?i)mlpe.*%s.*%s", cenario, nome_lag), g$label, perl=TRUE)) gw_mlpe_atual <- g
      if (grepl(sprintf("(?i)rf|randomforest.*%s.*%s", cenario, nome_lag), g$label, perl=TRUE)) gw_rf_atual <- g
    }
    
    if (!is.null(gw_mlpe_atual) && !is.null(gw_rf_atual)) {
      
      par(mfrow = c(2, 1), mar = c(4, 4, 6, 2)) 
      
      for (s in 1:4) {
        titulo <- sprintf("GW: MLPE + RF | Lags %s | Cenário %s", nome_lag, cenario)
        
        plot_window(gw_mlpe_atual$allY[[s]],
                    list(gw_mlpe_atual$allP[[s]], gw_rf_atual$allP[[s]]),
                    c("MLPE", "RandomForest"), 
                    c("blue", "red"),
                    names_s[s], 
                    titulo)
      }
    }
  }
}

dev.off() 
cat("Ficheiro 'Graficos_GW_Multivariados.pdf' guardado com sucesso.\n")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)