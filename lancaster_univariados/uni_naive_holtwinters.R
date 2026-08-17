# Forecasting Univariado - Seasonal Naive + Holt-Winters
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
store  <- load_store(store_file, store_name)
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

# Seasonal Naive
NAIVE      <- snaive(TR_ts, h = H)
Pred_naive <- as.numeric(NAIVE$mean)


r_naive <- plot_with_legend(
  "SeasonalNaive",
  Y,
  Pred_naive,
  store_name,
  YR
)



# Holt-Winters
HW      <- suppressWarnings(HoltWinters(TR_ts))
Pred_hw <- forecast(HW, h = H)$mean[1:H]



r_hw <- plot_with_legend(
  "HoltWinters",
  Y,
  Pred_hw,
  store_name,
  YR
)



cat("\nRESUMO FASE I -", store_name, "\n")
print(summary_table(
  list(r_naive, r_hw),
  c("SeasonalNaive", "HoltWinters")
))

# Fase 2 - com growing window
RUNS <- 20
S_step <- 7


cat("\n FASE II - GROWING WINDOW -", store_name, "\n")
cat("RUNS:", RUNS, "| H:", H, "| Step:", S_step, "\n\n")

W <- (L - H) - (RUNS - 1) * S_step

ev_naive  <- numeric(RUNS)
ev_hw     <- numeric(RUNS)
all_Y     <- numeric(0)
all_naive <- numeric(0)
all_hw    <- numeric(0)

for (b in 1:RUNS) {
  Hgw <- holdout(S, ratio = H, mode = "incremental",
                 iter = b, window = W, increment = S_step)
  
  dtr <- ts(S[Hgw$tr], frequency = K)
  Y_b <- S[Hgw$ts]
  
  Pred_naive_b <- as.numeric(snaive(dtr, h = H)$mean)
  
  Pred_hw_b <- tryCatch(
    forecast(suppressWarnings(HoltWinters(dtr)), h = H)$mean[1:H],
    error = function(e) rep(mean(S[Hgw$tr]), H)
  )
  
  ev_naive[b] <- mmetric(Y_b, Pred_naive_b, metric = "NMAE", val = YR)
  ev_hw[b]    <- mmetric(Y_b, Pred_hw_b,    metric = "NMAE", val = YR)
  
  all_Y     <- c(all_Y,     Y_b)
  all_naive <- c(all_naive, Pred_naive_b)
  all_hw    <- c(all_hw,    Pred_hw_b)
  
  cat(sprintf("  iter %2d | TR: %d-%d (%d dias) | TS: %d-%d | Naive=%.4f | HW=%.4f\n",
              b,
              Hgw$tr[1], Hgw$tr[length(Hgw$tr)], length(Hgw$tr),
              Hgw$ts[1], Hgw$ts[length(Hgw$ts)],
              ev_naive[b], ev_hw[b]))
}

cat("\n RESULTADOS FASE II -", store_name, "\n")
cat(sprintf("  %-15s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "SeasonalNaive",
            median(ev_naive), mean(ev_naive), sd(ev_naive),
            mmetric(all_Y, all_naive, metric = "R22")))
cat(sprintf("  %-15s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "HoltWinters",
            median(ev_hw), mean(ev_hw), sd(ev_hw),
            mmetric(all_Y, all_hw, metric = "R22")))

# Boxplot das métricas da Fase 2
dados_box <- data.frame(
  NMAE   = c(ev_naive, ev_hw),
  Modelo = rep(c("SeasonalNaive", "HoltWinters"), each = RUNS)
)
boxplot(NMAE ~ Modelo, data = dados_box,
        col  = c("steelblue", "tomato"),
        main = paste("Distribuição do NMAE por Iteração |", store_name),
        ylab = "NMAE (%)", xlab = "Modelo",
        outline = TRUE)

# Gráfico com valores reais e previsões acumuladas da growing window
plot_gw_legend(
  all_Y,
  list(all_naive, all_hw),
  c("SeasonalNaive", "HoltWinters"),
  c("blue", "red"),
  store_name,
  "SeasonalNaive + HoltWinters"
)
