# Forecasting Univariado - MLPE + Random Forest


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
Y     <- split$Y

# Preparar dados com lags
D  <- CasesSeries(S, c(1, 7))
hd <- holdout(D$y, ratio = H, mode = "order")

# Fase 1 - sem growing window
cat(" FASE I -", store_name, "\n")

# MLPE
set.seed(42)
M_mlpe    <- fit(y~., D[hd$tr,], model = "mlpe")
Pred_mlpe <- lforecast(M_mlpe, D, start = hd$ts[1], horizon = H)

r_mlpe <- plot_with_legend(
  "MLPE",
  Y,
  Pred_mlpe,
  store_name,
  YR
)

# Random Forest
set.seed(42)
M_rf    <- fit(y~., D[hd$tr,], model = "randomForest")
Pred_rf <- lforecast(M_rf, D, start = hd$ts[1], horizon = H)


r_rf <- plot_with_legend(
  "RandomForest",
  Y,
  Pred_rf,
  store_name,
  YR
)


cat("\n RESUMO FASE I -", store_name, "\n")
print(summary_table(
  list(r_mlpe, r_rf),
  c("MLPE", "RandomForest")
))



# Fase 2 - com growing window
RUNS   <- 20
S_step <- 7

# Conjuntos de lags testados
timelags_17   <- c(1, 7)
timelags_1728 <- c(1, 7, 28)
timelags_1237      <- c(1, 2, 3, 7)
timelags_137142128 <- c(1, 3, 7, 14, 21, 28)

# Criar séries com os diferentes conjuntos de lags
D2 <- CasesSeries(S, timelags_1728)
D3 <- CasesSeries(S, timelags_1237)
D4 <- CasesSeries(S, timelags_137142128)

cat("\nFASE II - GROWING WINDOW -", store_name, "\n")
cat("RUNS:", RUNS, "| H:", H, "| Step:", S_step, "\n\n")

W   <- (L - H) - (RUNS - 1) * S_step
W2  <- W - max(timelags_17)
W2b <- W - max(timelags_1728)
W3 <- W - max(timelags_1237)
W4 <- W - max(timelags_137142128)

ev_mlpe    <- numeric(RUNS); ev_rf    <- numeric(RUNS)
ev_mlpe_28 <- numeric(RUNS); ev_rf_28 <- numeric(RUNS)
all_Y <- numeric(0)
all_mlpe <- numeric(0); all_rf    <- numeric(0)
all_mlpe_28 <- numeric(0); all_rf_28 <- numeric(0)

ev_mlpe_1237 <- numeric(RUNS)
ev_rf_1237   <- numeric(RUNS)

ev_mlpe_137 <- numeric(RUNS)
ev_rf_137   <- numeric(RUNS)

all_mlpe_1237 <- numeric(0)
all_rf_1237   <- numeric(0)

all_mlpe_137 <- numeric(0)
all_rf_137   <- numeric(0)

for (b in 1:RUNS) {
  Hgw   <- holdout(S,    ratio = H, mode = "incremental",
                   iter = b, window = W,   increment = S_step)
  H2gw  <- holdout(D$y,  ratio = H, mode = "incremental",
                   iter = b, window = W2,  increment = S_step)
  H2gw2 <- holdout(D2$y, ratio = H, mode = "incremental",
                   iter = b, window = W2b, increment = S_step)
  
  H3gw <- holdout(D3$y, ratio = H, mode = "incremental",
                  iter = b, window = W3, increment = S_step)
  
  H4gw <- holdout(D4$y, ratio = H, mode = "incremental",
                  iter = b, window = W4, increment = S_step)
  
  Y_b <- S[Hgw$ts]
  
  # Modelos com lags c(1,7)
  set.seed(42)
  M_mlpe_b    <- fit(y~., D[H2gw$tr,], model = "mlpe")
  Pred_mlpe_b <- lforecast(M_mlpe_b, D, start = (length(H2gw$tr) + 1), horizon = H)
  
  set.seed(42)
  M_rf_b    <- fit(y~., D[H2gw$tr,], model = "randomForest")
  Pred_rf_b <- lforecast(M_rf_b, D, start = (length(H2gw$tr) + 1), horizon = H)
  
  # Modelos com lags c(1,7,28)
  set.seed(42)
  M_mlpe_b2    <- fit(y~., D2[H2gw2$tr,], model = "mlpe")
  Pred_mlpe_b2 <- lforecast(M_mlpe_b2, D2, start = (length(H2gw2$tr) + 1), horizon = H)
  
  set.seed(42)
  M_rf_b2    <- fit(y~., D2[H2gw2$tr,], model = "randomForest")
  Pred_rf_b2 <- lforecast(M_rf_b2, D2, start = (length(H2gw2$tr) + 1), horizon = H)
  
  # Modelos com lags c(1,2,3,7)
  set.seed(42)
  M_mlpe_1237 <- fit(y~., D3[H3gw$tr,], model = "mlpe")
  Pred_mlpe_1237 <- lforecast(M_mlpe_1237, D3,
                              start = (length(H3gw$tr) + 1), horizon = H)
  
  set.seed(42)
  M_rf_1237 <- fit(y~., D3[H3gw$tr,], model = "randomForest")
  Pred_rf_1237 <- lforecast(M_rf_1237, D3,
                            start = (length(H3gw$tr) + 1), horizon = H)
  
  # Modelos com lags c(1,3,7,14,21,28)
  set.seed(42)
  M_mlpe_137 <- fit(y~., D4[H4gw$tr,], model = "mlpe")
  Pred_mlpe_137 <- lforecast(M_mlpe_137, D4,
                             start = (length(H4gw$tr) + 1), horizon = H)
  
  set.seed(42)
  M_rf_137 <- fit(y~., D4[H4gw$tr,], model = "randomForest")
  Pred_rf_137 <- lforecast(M_rf_137, D4,
                           start = (length(H4gw$tr) + 1), horizon = H)
  
  
  ev_mlpe[b]    <- mmetric(Y_b, Pred_mlpe_b,  metric = "NMAE", val = YR)
  ev_rf[b]      <- mmetric(Y_b, Pred_rf_b,    metric = "NMAE", val = YR)
  ev_mlpe_28[b] <- mmetric(Y_b, Pred_mlpe_b2, metric = "NMAE", val = YR)
  ev_rf_28[b]   <- mmetric(Y_b, Pred_rf_b2,   metric = "NMAE", val = YR)
  ev_mlpe_1237[b] <- mmetric(Y_b, Pred_mlpe_1237, metric = "NMAE", val = YR)
  ev_rf_1237[b]   <- mmetric(Y_b, Pred_rf_1237,   metric = "NMAE", val = YR)
  
  ev_mlpe_137[b] <- mmetric(Y_b, Pred_mlpe_137, metric = "NMAE", val = YR)
  ev_rf_137[b]   <- mmetric(Y_b, Pred_rf_137,   metric = "NMAE", val = YR)
  
  
  all_Y       <- c(all_Y,       Y_b)
  all_mlpe    <- c(all_mlpe,    Pred_mlpe_b)
  all_rf      <- c(all_rf,      Pred_rf_b)
  all_mlpe_28 <- c(all_mlpe_28, Pred_mlpe_b2)
  all_rf_28   <- c(all_rf_28,   Pred_rf_b2)
  all_mlpe_1237 <- c(all_mlpe_1237, Pred_mlpe_1237)
  all_rf_1237   <- c(all_rf_1237,   Pred_rf_1237)
  
  all_mlpe_137 <- c(all_mlpe_137, Pred_mlpe_137)
  all_rf_137   <- c(all_rf_137,   Pred_rf_137)
  
  
  cat(sprintf("  iter %2d | MLPE(1,7)=%.4f | RF(1,7)=%.4f | MLPE(1,7,28)=%.4f | RF(1,7,28)=%.4f | MLPE(1,2,3,7)=%.4f | RF(1,2,3,7)=%.4f | MLPE(1,3,7,14,21,28)=%.4f | RF(1,3,7,14,21,28)=%.4f\n",
              b,
              ev_mlpe[b], ev_rf[b],
              ev_mlpe_28[b], ev_rf_28[b],
              ev_mlpe_1237[b], ev_rf_1237[b],
              ev_mlpe_137[b], ev_rf_137[b]))
}

cat("\nRESULTADOS FASE II -", store_name, "\n")
cat(sprintf("  %-20s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "MLPE c(1,7)",
            median(ev_mlpe), mean(ev_mlpe), sd(ev_mlpe),
            mmetric(all_Y, all_mlpe, metric = "R22")))
cat(sprintf("  %-20s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "RF c(1,7)",
            median(ev_rf), mean(ev_rf), sd(ev_rf),
            mmetric(all_Y, all_rf, metric = "R22")))
cat(sprintf("  %-20s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "MLPE c(1,7,28)",
            median(ev_mlpe_28), mean(ev_mlpe_28), sd(ev_mlpe_28),
            mmetric(all_Y, all_mlpe_28, metric = "R22")))
cat(sprintf("  %-20s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "RF c(1,7,28)",
            median(ev_rf_28), mean(ev_rf_28), sd(ev_rf_28),
            mmetric(all_Y, all_rf_28, metric = "R22")))

cat(sprintf("  %-20s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "MLPE c(1,2,3,7)",
            median(ev_mlpe_1237), mean(ev_mlpe_1237), sd(ev_mlpe_1237),
            mmetric(all_Y, all_mlpe_1237, metric = "R22")))
cat(sprintf("  %-20s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "RF c(1,2,3,7)",
            median(ev_rf_1237), mean(ev_rf_1237), sd(ev_rf_1237),
            mmetric(all_Y, all_rf_1237, metric = "R22")))
cat(sprintf("  %-20s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "MLPE c(1,3,7,14,21,28)",
            median(ev_mlpe_137), mean(ev_mlpe_137), sd(ev_mlpe_137),
            mmetric(all_Y, all_mlpe_137, metric = "R22")))
cat(sprintf("  %-20s Mediana NMAE: %.4f%% | Media: %.4f%% | SD: %.4f%% | R2: %.4f\n",
            "RF c(1,3,7,14,21,28)",
            median(ev_rf_137), mean(ev_rf_137), sd(ev_rf_137),
            mmetric(all_Y, all_rf_137, metric = "R22")))

# BOXPLOT DAS MÉTRICAS - FASE II
dados_box <- data.frame(
  NMAE = c(ev_mlpe, ev_rf,
           ev_mlpe_28, ev_rf_28,
           ev_mlpe_1237, ev_rf_1237,
           ev_mlpe_137, ev_rf_137),
  
  Modelo = rep(c(
    "MLPE c(1,7)", "RF c(1,7)",
    "MLPE c(1,7,28)", "RF c(1,7,28)",
    "MLPE c(1,2,3,7)", "RF c(1,2,3,7)",
    "MLPE c(1,3,7,14,21,28)", "RF c(1,3,7,14,21,28)"
  ), each = RUNS)
)
par(mar = c(10, 4, 4, 2))
boxplot(NMAE ~ Modelo, data = dados_box,
        col  = c("steelblue", "tomato", "lightblue", "salmon"),
        main = paste("Distribuição do NMAE por Iteração |", store_name),
        ylab = "NMAE (%)", xlab = "Modelo",
        outline = TRUE,
        las = 2,
        cex.axis = 0.8)

# Gráfico com valores reais e previsões para lags c(1,7)
plot_gw_legend(
  all_Y,
  list(all_mlpe, all_rf),
  c("MLPE", "RandomForest"),
  c("blue", "red"),
  store_name,
  "MLPE + RF | lags c(1,7)"
)

# Gráfico com valores reais e previsões para lags c(1,7,28)
plot_gw_legend(
  all_Y,
  list(all_mlpe_28, all_rf_28),
  c("MLPE", "RandomForest"),
  c("blue", "red"),
  store_name,
  "MLPE + RF | lags c(1,7,28)"
)

# Gráfico com valores reais e previsões para lags c(1,2,3,7)
plot_gw_legend(
  all_Y,
  list(all_mlpe_1237, all_rf_1237),
  c("MLPE", "RandomForest"),
  c("blue", "red"),
  store_name,
  "MLPE + RF | lags c(1,2,3,7)"
)

# Gráfico com valores reais e previsões para lags c(1,3,7,14,21,28)
plot_gw_legend(
  all_Y,
  list(all_mlpe_137, all_rf_137),
  c("MLPE", "RandomForest"),
  c("blue", "red"),
  store_name,
  "MLPE + RF | lags c(1,3,7,14,21,28)"
)
     
# Gráfico comparativo dos modelos MLPE
old_par <- par(no.readonly = TRUE)
par(oma = c(5, 0, 0, 0), mar = c(4, 4, 4, 2))

n <- length(all_Y)

ylim_all <- range(c(
  all_Y,
  all_mlpe,
  all_mlpe_28,
  all_mlpe_1237,
  all_mlpe_137
))

plot(1:n, all_Y,
     type = "l",
     lwd = 2,
     col = "black",
     ylim = ylim_all,
     main = paste("Comparação MLPE - Todos os Lags |", store_name),
     xlab = "Observações (períodos de teste acumulados)",
     ylab = "Num_Customers")

lines(all_mlpe,      col = "red",    lty = 2)
lines(all_mlpe_28,   col = "green",  lty = 2)
lines(all_mlpe_1237, col = "blue",   lty = 2)
lines(all_mlpe_137,  col = "purple", lty = 2)

par(fig = c(0,1,0,1), oma = c(0,0,0,0),
    mar = c(0,0,0,0), new = TRUE)
plot.new()

legend("bottom",
       legend = c(
         "Real",
         "c(1,7)",
         "c(1,7,28)",
         "c(1,2,3,7)",
         "c(1,3,7,14,21,28)"
       ),
       col  = c("black", "red", "green", "blue", "purple"),
       lty  = c(1, 2, 2, 2, 2),
       lwd  = 2,
       horiz = TRUE,
       bty = "n",
       cex = 0.75)

par(old_par) 

# Gráfico comparativo RF - Todos os Lags
old_par <- par(no.readonly = TRUE)
par(oma = c(5, 0, 0, 0), mar = c(4, 4, 4, 2))

ylim_all <- range(c(
  all_Y,
  all_rf,
  all_rf_28,
  all_rf_1237,
  all_rf_137
))

plot(1:n, all_Y,
     type = "l",
     lwd = 2,
     col = "black",
     ylim = ylim_all,
     main = paste("Comparação RF - Todos os Lags |", store_name),
     xlab = "Observações (períodos de teste acumulados)",
     ylab = "Num_Customers")

lines(all_rf,      col = "red",    lty = 2)
lines(all_rf_28,   col = "green",  lty = 2)
lines(all_rf_1237, col = "blue",   lty = 2)
lines(all_rf_137,  col = "purple", lty = 2)

par(fig = c(0,1,0,1), oma = c(0,0,0,0),
    mar = c(0,0,0,0), new = TRUE)
plot.new()

legend("bottom",
       legend = c(
         "Real",
         "c(1,7)",
         "c(1,7,28)",
         "c(1,2,3,7)",
         "c(1,3,7,14,21,28)"
       ),
       col  = c("black", "red", "green", "blue", "purple"),
       lty  = c(1, 2, 2, 2, 2),
       lwd  = 2,
       horiz = TRUE,
       bty = "n",
       cex = 0.75)

par(old_par)