# Forecasting Univariado - NNETAR + ETS

library(forecast)
library(rminer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R")

# Configuração inicial do ficheiro
store_file <- "richmond_clean.csv"
store_name <- "Richmond"
H <- 7
K <- 7

# Carregar o CSV e preparar
store <- load_store(store_file, store_name)
S  <- store$S
L  <- store$L
YR <- store$YR
df <- store$df

# Divisão em treino e teste
split <- train_test_split(S, L, H, K)
TR_ts <- split$TR_ts
Y     <- split$Y

# Fase 1 - sem growing window
cat(" FASE I -", store_name, "\n")

# NNETAR
set.seed(42)
NN1      <- nnetar(TR_ts, P = 1, repeats = 5)
print(NN1)
Pred_nn1 <- forecast(NN1, h = H)$mean[1:H]

r_nn1 <- plot_with_legend(
  "NNETAR",
  Y,
  Pred_nn1,
  store_name,
  YR
)

# ETS
ETS      <- ets(TR_ts)
print(ETS)
Pred_ets <- forecast(ETS, h = H)$mean[1:H]

r_ets <- plot_with_legend(
  "ETS",
  Y,
  Pred_ets,
  store_name,
  YR
)

cat("\n RESUMO FASE I -", store_name, "\n")
print(summary_table(
  list(r_nn1, r_ets),
  c("NNETAR", "ETS")
))


# Fase 2 - com growing window
RUNS <- 20
S_step <- 7


cat("\n FASE II - GROWING WINDOW -", store_name, "\n")
cat("RUNS:", RUNS, "| H:", H, "| Step:", S_step, "\n\n")

W <- (L - H) - (RUNS - 1) * S_step

ev_nn1  <- numeric(RUNS)
ev_ets  <- numeric(RUNS)
all_Y   <- numeric(0)
all_nn1 <- numeric(0)
all_ets <- numeric(0)

for (b in 1:RUNS) {
  Hgw <- holdout(S, ratio = H, mode = "incremental",
                 iter = b, window = W, increment = S_step)
  
  dtr <- ts(S[Hgw$tr], frequency = K)
  Y_b <- S[Hgw$ts]
  
  # NNETAR
  set.seed(42)
  Pred_nn1_b <- tryCatch(
    forecast(nnetar(dtr, P = 1, repeats = 5), h = H)$mean[1:H],
    error = function(e) rep(mean(S[Hgw$tr]), H)
  )
  
  # ETS
  Pred_ets_b <- tryCatch(
    forecast(ets(dtr), h = H)$mean[1:H],
    error = function(e) rep(mean(S[Hgw$tr]), H)
  )
  
  ev_nn1[b] <- mmetric(Y_b, Pred_nn1_b, metric = "NMAE", val = YR)
  ev_ets[b] <- mmetric(Y_b, Pred_ets_b, metric = "NMAE", val = YR)
  
  all_Y   <- c(all_Y,   Y_b)
  all_nn1 <- c(all_nn1, Pred_nn1_b)
  all_ets <- c(all_ets, Pred_ets_b)
  
  cat(sprintf("  iter %2d | TR: %d-%d (%d dias) | TS: %d-%d | NNETAR=%.4f | ETS=%.4f\n",
              b,
              Hgw$tr[1], Hgw$tr[length(Hgw$tr)], length(Hgw$tr),
              Hgw$ts[1], Hgw$ts[length(Hgw$ts)],
              ev_nn1[b], ev_ets[b]))
}

cat("\n RESULTADOS FASE II -", store_name, "\n")
cat(sprintf("  %-15s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "NNETAR",
            median(ev_nn1), mean(ev_nn1), sd(ev_nn1),
            mmetric(all_Y, all_nn1, metric = "R22")))
cat(sprintf("  %-15s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "ETS",
            median(ev_ets), mean(ev_ets), sd(ev_ets),
            mmetric(all_Y, all_ets, metric = "R22")))


# Boxplot das métricas da FASE II
dados_box <- data.frame(
  NMAE   = c(ev_nn1, ev_ets),
  Modelo = rep(c("NNETAR", "ETS"), each = RUNS)
)
boxplot(NMAE ~ Modelo, data = dados_box,
        col  = c("steelblue", "tomato"),
        main = paste("Distribuição do NMAE por Iteração |", store_name),
        ylab = "NMAE (%)", xlab = "Modelo",
        outline = TRUE)

plot_gw_legend(
  all_Y,
  list(all_nn1, all_ets),
  c("NNETAR", "ETS"),
  c("blue", "red"),
  store_name,
  "NNETAR + ETS"
)