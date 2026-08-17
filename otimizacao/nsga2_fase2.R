# NSGA-II — Otimização Multi-Objetivo (O3)
# Fase 2: 20 semanas do Growing Window


library(mco)
source("avaliacao.R")


# 1. Carregar previsões GW
cat("A carregar previsões GW do modelo RF Cenário C...\n")

bt_pred <- read.csv("previsoes_gw_baltimore.csv",    header = TRUE)
lc_pred <- read.csv("previsoes_gw_lancaster.csv",    header = TRUE)
ph_pred <- read.csv("previsoes_gw_philadelphia.csv", header = TRUE)
rv_pred <- read.csv("previsoes_gw_richmond.csv",     header = TRUE)

RUNS <- max(bt_pred$Run)


# 2. Funções auxiliares
is_weekend <- function(weekday_vec) {
  weekday_vec %in% c("sábado", "domingo", "Saturday", "Sunday")
}

calc_bounds <- function(C) {
  lower_store <- c()
  upper_store <- c()
  for (d in 1:7) {
    lower_store <- c(lower_store, 0, 0, 0.00)
    upper_store <- c(upper_store,
                     ceiling(C[d] / 6),
                     ceiling(C[d] / 7),
                     0.30)
  }
  list(lower = lower_store, upper = upper_store)
}

REF_PROFIT <- 11
REF_HR     <- 55

calc_hv <- function(profit_vec, hr_vec, ref_profit, ref_HR) {
  points <- cbind(-profit_vec, hr_vec)
  ref    <- c(-ref_profit, ref_HR)
  
  ord    <- order(points[, 1])
  points <- points[ord, , drop = FALSE]
  
  n  <- nrow(points)
  hv <- 0
  
  current_ref2 <- ref[2]
  
  for (i in 1:n) {
    width  <- ref[1] - points[i, 1]
    height <- current_ref2 - points[i, 2]
    if (width > 0 && height > 0) {
      hv <- hv + width * height
    }
    current_ref2 <- points[i, 2]
  }
  return(hv)
}


# 3. NSGA-II — 20 SEMANAS
popsize <- 200
gens    <- 600

media_profit_vec   <- numeric(RUNS)
mediana_profit_vec <- numeric(RUNS)
media_HR_vec       <- numeric(RUNS)
mediana_HR_vec     <- numeric(RUNS)
max_profit_vec     <- numeric(RUNS)
n_pareto_vec       <- numeric(RUNS)
hv_vec             <- numeric(RUNS)
best_units_vec     <- numeric(RUNS)

all_pareto <- list()

all_profits_global <- c()
all_HR_global      <- c()

cat(sprintf("\nNSGA-II O3 - %d semanas (GW)\n", RUNS))
cat(sprintf("População: %d | Gerações: %d\n\n", popsize, gens))


# 4. Loop principal — 20 Semanas
for (run in 1:RUNS) {
  
  bt_run <- bt_pred[bt_pred$Run == run, ]
  lc_run <- lc_pred[lc_pred$Run == run, ]
  ph_run <- ph_pred[ph_pred$Run == run, ]
  rv_run <- rv_pred[rv_pred$Run == run, ]
  
  C_baltimore    <<- bt_run$Predicted_Customers
  C_lancaster    <<- lc_run$Predicted_Customers
  C_philadelphia <<- ph_run$Predicted_Customers
  C_richmond     <<- rv_run$Predicted_Customers
  
  weekend_baltimore    <<- is_weekend(bt_run$Weekday)
  weekend_lancaster    <<- is_weekend(lc_run$Weekday)
  weekend_philadelphia <<- is_weekend(ph_run$Weekday)
  weekend_richmond     <<- is_weekend(rv_run$Weekday)
  
  b_bt <- calc_bounds(C_baltimore)
  b_lc <- calc_bounds(C_lancaster)
  b_ph <- calc_bounds(C_philadelphia)
  b_rv <- calc_bounds(C_richmond)
  
  lower <- c(b_bt$lower, b_lc$lower, b_ph$lower, b_rv$lower)
  upper <- c(b_bt$upper, b_lc$upper, b_ph$upper, b_rv$upper)
  
  cat(sprintf("--- Run %2d | %s a %s ---\n",
              run, bt_run$Date[1], bt_run$Date[7]))
  
  set.seed(42)
  NSGA <- nsga2(
    fn           = eval_O3,
    idim         = length(lower),
    odim         = 2,
    lower.bounds = lower,
    upper.bounds = upper,
    popsize      = popsize,
    generations  = gens
  )
  
  pareto_idx <- which(NSGA$pareto.optimal)
  raw_profit <- -NSGA$value[pareto_idx, 1]
  raw_HR     <-  NSGA$value[pareto_idx, 2]
  raw_sols   <-  NSGA$par[pareto_idx, , drop = FALSE]
  
  valid_idx <- which(raw_profit > 0)
  
  if (length(valid_idx) == 0) {
    cat("  Sem soluções Pareto válidas (Lucro >= 0) nesta semana.\n\n")
    all_pareto[[run]] <- NULL
    next
  }
  
  pareto_profit <- raw_profit[valid_idx]
  pareto_HR     <- raw_HR[valid_idx]
  pareto_sols   <- raw_sols[valid_idx, , drop = FALSE]
  
  best_idx_run   <- which.max(pareto_profit)
  best_sol_run   <- pareto_sols[best_idx_run, ]
  best_units_run <- eval_O3_restricao(best_sol_run) + 10000
  
  all_profits_global <- c(all_profits_global, pareto_profit)
  all_HR_global      <<- c(all_HR_global,      pareto_HR)
  
  hv_vec[run] <- calc_hv(pareto_profit, pareto_HR, REF_PROFIT, REF_HR)
  
  n_pareto_vec[run]       <- length(pareto_profit)
  max_profit_vec[run]     <- max(pareto_profit)
  media_profit_vec[run]   <- mean(pareto_profit)
  mediana_profit_vec[run] <- median(pareto_profit)
  media_HR_vec[run]       <- mean(pareto_HR)
  mediana_HR_vec[run]     <- median(pareto_HR)
  best_units_vec[run]     <- best_units_run
  
  all_pareto[[run]] <- list(
    profit      = pareto_profit,
    HR          = pareto_HR,
    data_inicio = bt_run$Date[1],
    data_fim    = bt_run$Date[7]
  )
  
  cat(sprintf("  Pareto: %d soluções limpas\n", n_pareto_vec[run]))
  cat(sprintf("  Lucro  — Média: %d  Mediana: %d  Máx: %d USD (com %d Unidades)\n",
              round(media_profit_vec[run]),
              round(mediana_profit_vec[run]),
              round(max_profit_vec[run]),
              round(best_units_run)))
  cat(sprintf("  HR     — Média: %d  Mediana: %d\n",
              round(media_HR_vec[run]),
              round(mediana_HR_vec[run])))
  cat(sprintf("  Hipervolume (ref global): %.2f\n\n", hv_vec[run]))
}


# 5. Hipervolume global
hv_global_vec <- numeric(RUNS)
for (run in 1:RUNS) {
  if (!is.null(all_pareto[[run]])) {
    hv_global_vec[run] <- calc_hv(
      all_pareto[[run]]$profit,
      all_pareto[[run]]$HR,
      REF_PROFIT,
      REF_HR
    )
  }
}


# 6. Agregação final
valid_runs  <- which(!sapply(all_pareto, is.null))
units_valid <- best_units_vec[valid_runs]

cat("RESUMO NSGA-II O3 — FASE 2 (20 semanas)\n")
cat(sprintf("%-30s %10s %10s\n", "Métrica", "Média", "Mediana"))
cat(rep("-", 52), "\n", sep = "")
cat(sprintf("%-30s %10d %10d\n", "Nº soluções Pareto",
            round(mean(n_pareto_vec[n_pareto_vec>0])),
            round(median(n_pareto_vec[n_pareto_vec>0]))))
cat(sprintf("%-30s %10d %10d\n", "Máx lucro da curva (USD)",
            round(mean(max_profit_vec[max_profit_vec>0])),
            round(median(max_profit_vec[max_profit_vec>0]))))
cat(sprintf("%-30s %10d %10d\n", "Média lucro da curva (USD)",
            round(mean(media_profit_vec[media_profit_vec>0])),
            round(median(media_profit_vec[media_profit_vec>0]))))
cat(sprintf("%-30s %10d %10d\n", "Mediana lucro da curva (USD)",
            round(mean(mediana_profit_vec[mediana_profit_vec>0])),
            round(median(mediana_profit_vec[mediana_profit_vec>0]))))
cat(sprintf("%-30s %10d %10d\n", "Média HR da curva",
            round(mean(media_HR_vec[media_HR_vec>0])),
            round(median(media_HR_vec[media_HR_vec>0]))))
cat(sprintf("%-30s %10d %10d\n", "Mediana HR da curva",
            round(mean(mediana_HR_vec[mediana_HR_vec>0])),
            round(median(mediana_HR_vec[mediana_HR_vec>0]))))
cat(sprintf("%-30s %10d %10d\n", "Unidades vendidas (melhor sol.)",
            round(mean(units_valid)),
            round(median(units_valid))))
cat(rep("-", 52), "\n", sep = "")
cat(sprintf("%-30s %10d\n", ">> Mediana da mediana lucro",
            round(median(mediana_profit_vec[mediana_profit_vec>0]))))
cat(sprintf("%-30s %10d\n", ">> Mediana da mediana HR",
            round(median(mediana_HR_vec[mediana_HR_vec>0]))))
cat(sprintf("%-30s %10d\n", ">> Média das unidades vendidas",
            round(mean(units_valid))))
cat(sprintf("%-30s %10d\n", ">> Mediana das unidades vendidas",
            round(median(units_valid))))
cat(rep("-", 52), "\n", sep = "")
cat(sprintf("%-30s %10.2f %10.2f\n", "Hipervolume (ref global)",
            mean(hv_global_vec[hv_global_vec>0]),
            median(hv_global_vec[hv_global_vec>0])))
cat(rep("=", 52), "\n", sep = "")

# 7. Tabela final
resultados_df <- data.frame(
  Run           = 1:RUNS,
  Semana        = sapply(1:RUNS, function(r) ifelse(is.null(all_pareto[[r]]), "N/A", all_pareto[[r]]$data_inicio)),
  N_Pareto      = n_pareto_vec,
  Max_Lucro     = round(max_profit_vec),
  Media_Lucro   = round(media_profit_vec),
  Mediana_Lucro = round(mediana_profit_vec),
  Media_HR      = round(media_HR_vec),
  Mediana_HR    = round(mediana_HR_vec),
  Unidades      = round(best_units_vec),
  HV_global     = round(hv_global_vec, 2)
)
cat("\nTabela completa por semana:\n")
print(resultados_df, row.names = FALSE)

# 8. Gráficos — 20 Curvas de Pareto
pdf("NSGA2_Pareto_GW_20semanas.pdf", width = 20, height = 16)
par(mfrow = c(4, 5), mar = c(4, 3, 4, 1), oma = c(2, 1, 4, 1))

for (run in 1:RUNS) {
  
  if (is.null(all_pareto[[run]])) {
    plot(1, type = "n", main = sprintf("Run %d\n(sem soluções válidas)", run),
         xlab = "", ylab = "", cex.main = 0.85)
    next
  }
  
  p   <- all_pareto[[run]]
  ord <- order(p$HR)
  
  plot(p$HR[ord], p$profit[ord],
       type = "b",
       pch  = 16, cex = 0.7,
       col  = "steelblue",
       xlab = "HR", ylab = "Lucro (USD)",
       main = sprintf("Run %d | %s", run, p$data_inicio),
       cex.main = 0.85)
  
  abline(h = media_profit_vec[run],   col = "orange", lty = 2, lwd = 1.5)
  abline(h = mediana_profit_vec[run], col = "green",  lty = 3, lwd = 1.5)
  
  legend("top", bty = "n", cex = 0.7, horiz = TRUE,
         legend = c(sprintf("HV=%.0f",    hv_global_vec[run]),
                    sprintf("Média=%d",   round(media_profit_vec[run])),
                    sprintf("Mediana=%d", round(mediana_profit_vec[run]))),
         col    = c("steelblue", "orange", "green"),
         lty    = c(1, 2, 3),
         lwd    = c(1, 1.5, 1.5))
  
  grid()
}

mtext("Curvas de Pareto NSGA-II — O3 (20 semanas GW)",
      outer = TRUE, cex = 1.3, font = 2)

dev.off()
cat("\nFicheiro 'NSGA2_Pareto_GW_20semanas.pdf' guardado com sucesso.\n")