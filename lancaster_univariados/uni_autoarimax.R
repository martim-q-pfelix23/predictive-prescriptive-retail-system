# Forecasting Univariado - AUTO.ARIMA (SARIMA)
library(forecast)
library(rminer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R")

# Configuração inicial do ficheiro
store_file <- "lancaster_clean.csv"
store_name <- "Lancaster"
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
cat("FASE I -", store_name, "\n")

# AUTO.ARIMA
AR <- suppressWarnings(auto.arima(TR_ts, seasonal = TRUE))
print(AR)  # mostrar o modelo selecionado
Pred_ar <- forecast(AR, h = H)$mean[1:H]

r_ar <- plot_with_legend(
  "AutoARIMA",
  Y,
  Pred_ar,
  store_name,
  YR
)


cat("\nRESUMO FASE I -", store_name, "\n")
print(summary_table(
  list(r_ar),
  c("AutoARIMA")
))


# Fase 2 - com growing window
RUNS <- 20
S_step <- 7


cat("\nFASE II - GROWING WINDOW -", store_name, "\n")
cat("RUNS:", RUNS, "| H:", H, "| Step:", S_step, "\n\n")

W <- (L - H) - (RUNS - 1) * S_step

ev_ar  <- numeric(RUNS)
all_Y  <- numeric(0)
all_ar <- numeric(0)

for (b in 1:RUNS) {
  Hgw <- holdout(S, ratio = H, mode = "incremental",
                 iter = b, window = W, increment = S_step)
  
  dtr <- ts(S[Hgw$tr], frequency = K)
  Y_b <- S[Hgw$ts]
  
  Pred_ar_b <- tryCatch(
    forecast(suppressWarnings(auto.arima(dtr, seasonal = TRUE)), h = H)$mean[1:H],
    error = function(e) rep(mean(S[Hgw$tr]), H)
  )
  
  ev_ar[b] <- mmetric(Y_b, Pred_ar_b, metric = "NMAE", val = YR)
  
  all_Y  <- c(all_Y,  Y_b)
  all_ar <- c(all_ar, Pred_ar_b)
  
  cat(sprintf("  iter %2d | TR: %d-%d (%d dias) | TS: %d-%d | AutoARIMA=%.4f\n",
              b,
              Hgw$tr[1], Hgw$tr[length(Hgw$tr)], length(Hgw$tr),
              Hgw$ts[1], Hgw$ts[length(Hgw$ts)],
              ev_ar[b]))
}

cat("\nRESULTADOS FASE II -", store_name, "\n")
cat(sprintf("  %-15s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "AutoARIMA",
            median(ev_ar), mean(ev_ar), sd(ev_ar),
            mmetric(all_Y, all_ar, metric = "R22")))

# Boxplot das métricas da Fase 2
boxplot(ev_ar,
        col  = "steelblue",
        main = paste("Distribuição do NMAE - AutoARIMA |", store_name),
        ylab = "NMAE (%)", xlab = "AutoARIMA",
        outline = TRUE)

# Gráfico com valores reais e previsões acumuladas da growing window
plot_gw_legend(
  all_Y,
  list(all_ar),
  c("AutoARIMA"),
  c("blue"),
  store_name,
  "AutoARIMA"
)