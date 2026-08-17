# SCRIPT DE GERAÇÃO DE PREVISÕES FINAIS (Para a Otimização)
# Modelo: Random Forest - Cenário C (Global - 4 Lojas)
# Lags: 1, 3, 7, 14, 21, 28


library(rminer)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R")
source("../multi-utils.R")


# Config

H <- 7


# Carregar dados

cat("A carregar os datasets...\n")

bt <- load_store("baltimore_clean.csv",    "Baltimore")
lc <- load_store("lancaster_clean.csv",    "Lancaster")
ph <- load_store("philadelphia_clean.csv", "Philadelphia")
rc <- load_store("richmond_clean.csv",     "Richmond")

L <- length(bt$S)
LTR <- L - H   # último ponto de treino



# Matriz global
mtr_4 <- cbind(
  Baltimore    = bt$S[1:LTR],
  Lancaster    = lc$S[1:LTR],
  Philadelphia = ph$S[1:LTR],
  Richmond     = rc$S[1:LTR]
)

# Lags (VINP) 
lags_principais <- c(1,3,7,14,21,28)

VINP_C <- vector("list", length = 4)

VINP_C[[1]] <- list(lags_principais, 1, 1, 1)
VINP_C[[2]] <- list(1, lags_principais, 1, 1)
VINP_C[[3]] <- list(1, 1, lags_principais, 1)
VINP_C[[4]] <- list(1, 1, 1, lags_principais)


# Exógenas
prep_exog_full <- function(df) {
  cbind(
    TouristEvent  = as.numeric(df$TouristEvent == "Yes"),
    Pct_On_Sale   = df$Pct_On_Sale / 100,
    Num_Employees = df$Num_Employees / max(df$Num_Employees[1:LTR])
  )
}

exf_bt <- prep_exog_full(bt$df)
exf_lc <- prep_exog_full(lc$df)
exf_ph <- prep_exog_full(ph$df)
exf_rc <- prep_exog_full(rc$df)

# treino
exotr <- cbind(
  BT_Tourist=exf_bt[1:LTR,1], BT_Promo=exf_bt[1:LTR,2], BT_Empl=exf_bt[1:LTR,3],
  LC_Tourist=exf_lc[1:LTR,1], LC_Promo=exf_lc[1:LTR,2], LC_Empl=exf_lc[1:LTR,3],
  PH_Tourist=exf_ph[1:LTR,1], PH_Promo=exf_ph[1:LTR,2], PH_Empl=exf_ph[1:LTR,3],
  RC_Tourist=exf_rc[1:LTR,1], RC_Promo=exf_rc[1:LTR,2], RC_Empl=exf_rc[1:LTR,3]
)

# "futuro" (proxy: últimos dados disponíveis)
exots <- cbind(
  BT_Tourist=exf_bt[(LTR+1):L,1], BT_Promo=exf_bt[(LTR+1):L,2], BT_Empl=exf_bt[(LTR+1):L,3],
  LC_Tourist=exf_lc[(LTR+1):L,1], LC_Promo=exf_lc[(LTR+1):L,2], LC_Empl=exf_lc[(LTR+1):L,3],
  PH_Tourist=exf_ph[(LTR+1):L,1], PH_Promo=exf_ph[(LTR+1):L,2], PH_Empl=exf_ph[(LTR+1):L,3],
  RC_Tourist=exf_rc[(LTR+1):L,1], RC_Promo=exf_rc[(LTR+1):L,2], RC_Empl=exf_rc[(LTR+1):L,3]
)


# Treinar modelo final
cat("\nA treinar Random Forest (cenário C global)...\n")

set.seed(42)
MNN_rf_C <- mfit(mtr_4, "randomForest", VINP_C, exogen = exotr)


# Previsões
cat("A gerar previsões para os próximos 7 dias...\n")

Pred_rf_C <- lforecastm(MNN_rf_C, h = H, exogen = exots)


# Exportar CSVs
cat("A exportar previsões...\n")

save_forecasts(round(Pred_rf_C[[1]]), bt$df[1:LTR, ], "Baltimore", H)
save_forecasts(round(Pred_rf_C[[2]]), lc$df[1:LTR, ], "Lancaster", H)
save_forecasts(round(Pred_rf_C[[3]]), ph$df[1:LTR, ], "Philadelphia", H)
save_forecasts(round(Pred_rf_C[[4]]), rc$df[1:LTR, ], "Richmond", H)

cat("\nSucesso! CSVs gerados para as 4 lojas.\n")