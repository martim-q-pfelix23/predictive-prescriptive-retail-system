# SIMULATED ANNEALING O3 — Soma Ponderada
# Fase 2: 20 semanas do Growing Window


source("hill.R")
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

get_metrics_O3 <- function(s) {
  stores <- make_stores()
  total_profit <- 0
  total_hr     <- 0
  for (i in 1:4) {
    idx <- ((i-1)*21 + 1):(i*21)
    res <- eval(s[idx], stores[[i]])
    total_profit <- total_profit + res$profit
    J_vals <- round(s[idx][seq(1, 21, by = 3)])
    X_vals <- round(s[idx][seq(2, 21, by = 3)])
    total_hr <- total_hr + sum(J_vals) + sum(X_vals)
  }
  list(profit = total_profit, hr = total_hr)
}

get_units_total <- function(s) {
  stores <- make_stores()
  total_units <- 0
  for (i in 1:4) {
    idx <- ((i-1)*21 + 1):(i*21)
    res <- eval(s[idx], stores[[i]])
    total_units <- total_units + res$units
  }
  return(total_units)
}

rchange_sann <- function(par, ...) {
  lim_inf <- get("lower", envir = .GlobalEnv)
  lim_sup <- get("upper", envir = .GlobalEnv)
  hchange(par, lower = lim_inf, upper = lim_sup,
          operator = "*", dist = rnorm, mean = 1, sd = 0.05,
          round = FALSE)
}

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
    if (width > 0 && height > 0) hv <- hv + width * height
    current_ref2 <- points[i, 2]
  }
  return(hv)
}

is_dominated <- function(df) {
  n   <- nrow(df)
  dom <- logical(n)
  if (n < 2) return(dom)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      if (df$Lucro[j] >= df$Lucro[i] &&
          df$HR[j]    <= df$HR[i]    &&
          (df$Lucro[j] > df$Lucro[i] || df$HR[j] < df$HR[i])) {
        dom[i] <- TRUE
        break
      }
    }
  }
  return(dom)
}

# 3. Configuração
REF_PROFIT <- 11
REF_HR     <- 55

pesos_lucro_testar <- seq(0.1, 0.9, by = 0.1)
N_O3 <- 10000
temp <- 0.05

hv_vec         <- numeric(RUNS)
med_lucro_run  <- numeric(RUNS)
mean_lucro_run <- numeric(RUNS)
med_hr_run     <- numeric(RUNS)
mean_hr_run    <- numeric(RUNS)
max_lucro_run  <- numeric(RUNS)
n_pareto_run   <- numeric(RUNS)
best_units_vec <- numeric(RUNS)

todas_pareto <- data.frame(
  Run        = integer(),
  Semana     = character(),
  Peso_Lucro = numeric(),
  Peso_HR    = numeric(),
  Lucro      = numeric(),
  HR         = numeric()
)

# guardar também os pares (par, lucro) para extrair unidades depois
all_pareto      <- list()
all_pareto_sols <- list()

cat(sprintf("\nSANN O3 - %d semanas (GW)\n\n", RUNS))

# 4. Loop Principal — 20 Semanas
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
  
  lower <<- c(b_bt$lower, b_lc$lower, b_ph$lower, b_rv$lower)
  upper <<- c(b_bt$upper, b_lc$upper, b_ph$upper, b_rv$upper)
  
  s_mid <- lower + (upper - lower) / 2
  
  cat(sprintf("Run %2d | %s a %s\n",
              run, bt_run$Date[1], bt_run$Date[7]))
  
  pareto_run      <- data.frame(
    Peso_Lucro = numeric(),
    Peso_HR    = numeric(),
    Lucro      = numeric(),
    HR         = numeric()
  )
  pareto_run_sols <- list()
  
  for (w_lucro in pesos_lucro_testar) {
    w_hr <- 1 - w_lucro
    
    set.seed(42)
    SA_O3 <- optim(
      par     = s_mid,
      fn      = neg_eval_O3_pesos,
      gr      = rchange_sann,
      method  = "SANN",
      control = list(maxit = N_O3, temp = temp, trace = FALSE),
      w_lucro = w_lucro,
      w_hr    = w_hr
    )
    
    metrics <- get_metrics_O3(SA_O3$par)
    pareto_run <- rbind(pareto_run, data.frame(
      Peso_Lucro = w_lucro,
      Peso_HR    = w_hr,
      Lucro      = metrics$profit,
      HR         = metrics$hr
    ))
    pareto_run_sols[[length(pareto_run_sols) + 1]] <- SA_O3$par
    
    cat(sprintf("  w=%.1f | Lucro: %d USD | HR: %d\n",
                w_lucro, round(metrics$profit), round(metrics$hr)))
  }
  
  # Filtrar negativos e dominadas
  p_valido     <- pareto_run[pareto_run$Lucro > 0, ]
  valid_idx_pv <- which(pareto_run$Lucro > 0)
  
  if (nrow(p_valido) > 1) {
    dom_mask     <- !is_dominated(p_valido)
    p_valido     <- p_valido[dom_mask, ]
    valid_idx_pv <- valid_idx_pv[dom_mask]
  }
  
  n_pareto_run[run]  <- nrow(p_valido)
  max_lucro_run[run] <- if (nrow(p_valido) > 0) max(p_valido$Lucro) else 0
  hv_vec[run]        <- if (nrow(p_valido) > 0)
    calc_hv(p_valido$Lucro, p_valido$HR, REF_PROFIT, REF_HR) else 0
  
  # Unidades da solução com maior lucro no Pareto limpo
  if (nrow(p_valido) > 0) {
    best_local_idx     <- valid_idx_pv[which.max(p_valido$Lucro)]
    best_sol_run       <- pareto_run_sols[[best_local_idx]]
    best_units_vec[run] <- get_units_total(best_sol_run)
  } else {
    best_units_vec[run] <- NA
  }
  
  all_pareto[[run]] <- if (nrow(p_valido) > 0) list(
    profit      = p_valido$Lucro,
    HR          = p_valido$HR,
    data_inicio = bt_run$Date[1],
    data_fim    = bt_run$Date[7]
  ) else NULL
  
  med_lucro_run[run]  <- if (nrow(p_valido) > 0) median(p_valido$Lucro) else NA
  mean_lucro_run[run] <- if (nrow(p_valido) > 0) mean(p_valido$Lucro)   else NA
  med_hr_run[run]     <- if (nrow(p_valido) > 0) median(p_valido$HR)    else NA
  mean_hr_run[run]    <- if (nrow(p_valido) > 0) mean(p_valido$HR)      else NA
  
  cat(sprintf("  Pareto limpo: %d soluções\n", n_pareto_run[run]))
  cat(sprintf("  Lucro  — Média: %d  Mediana: %d  Máx: %d USD (com %d Unidades)\n",
              round(mean_lucro_run[run]), round(med_lucro_run[run]),
              round(max_lucro_run[run]),  round(best_units_vec[run])))
  cat(sprintf("  HR     — Média: %d  Mediana: %d\n",
              round(mean_hr_run[run]), round(med_hr_run[run])))
  cat(sprintf("  Hipervolume (ref global): %.2f\n\n", hv_vec[run]))
  
  todas_pareto <- rbind(todas_pareto, data.frame(
    Run        = run,
    Semana     = bt_run$Date[1],
    pareto_run
  ))
}


# 5. Agregação final
valid_runs  <- which(!sapply(all_pareto, is.null))
units_valid <- best_units_vec[valid_runs]

cat("RESUMO SANN O3 — FASE 2 (20 semanas)\n")
cat(sprintf("%-30s %10s %10s\n", "Métrica", "Média", "Mediana"))
cat(rep("-", 52), "\n", sep = "")
cat(sprintf("%-30s %10d %10d\n", "Nº soluções Pareto",
            round(mean(n_pareto_run[n_pareto_run > 0])),
            round(median(n_pareto_run[n_pareto_run > 0]))))
cat(sprintf("%-30s %10d %10d\n", "Máx lucro da curva (USD)",
            round(mean(max_lucro_run[max_lucro_run > 0])),
            round(median(max_lucro_run[max_lucro_run > 0]))))
cat(sprintf("%-30s %10d %10d\n", "Média lucro da curva (USD)",
            round(mean(mean_lucro_run)), round(median(mean_lucro_run))))
cat(sprintf("%-30s %10d %10d\n", "Mediana lucro da curva (USD)",
            round(mean(med_lucro_run)),  round(median(med_lucro_run))))
cat(sprintf("%-30s %10d %10d\n", "Média HR da curva",
            round(mean(mean_hr_run)),    round(median(mean_hr_run))))
cat(sprintf("%-30s %10d %10d\n", "Mediana HR da curva",
            round(mean(med_hr_run)),     round(median(med_hr_run))))
cat(sprintf("%-30s %10d %10d\n", "Unidades vendidas (melhor sol.)",
            round(mean(units_valid)),
            round(median(units_valid))))
cat(rep("-", 52), "\n", sep = "")
cat(sprintf("%-30s %10d\n", ">> Mediana da mediana lucro",
            round(median(med_lucro_run))))
cat(sprintf("%-30s %10d\n", ">> Mediana da mediana HR",
            round(median(med_hr_run))))
cat(sprintf("%-30s %10d\n", ">> Média das unidades vendidas",
            round(mean(units_valid))))
cat(sprintf("%-30s %10d\n", ">> Mediana das unidades vendidas",
            round(median(units_valid))))
cat(rep("-", 52), "\n", sep = "")
cat(sprintf("%-30s %10.2f %10.2f\n", "Hipervolume (ref global)",
            mean(hv_vec[hv_vec > 0]), median(hv_vec[hv_vec > 0])))
cat(rep("=", 52), "\n", sep = "")



# 6. Tabela final
resultados_df <- data.frame(
  Run           = 1:RUNS,
  Semana        = sapply(1:RUNS, function(r) bt_pred$Date[bt_pred$Run == r][1]),
  N_Pareto      = n_pareto_run,
  Max_Lucro     = round(max_lucro_run),
  Media_Lucro   = round(mean_lucro_run),
  Mediana_Lucro = round(med_lucro_run),
  Media_HR      = round(mean_hr_run),
  Mediana_HR    = round(med_hr_run),
  Unidades      = ifelse(is.na(best_units_vec), NA, round(best_units_vec)),
  HV_global     = round(hv_vec, 2)
)
cat("\nTabela completa por semana:\n")
print(resultados_df, row.names = FALSE)

# 7. Gráficos — 20 Curvas de Pareto
pdf("SANN_O3_Pareto_GW.pdf", width = 20, height = 16)
par(mfrow = c(4, 5), mar = c(4, 3, 4, 1), oma = c(2, 1, 4, 1))

for (run in 1:RUNS) {
  p   <- todas_pareto[todas_pareto$Run == run, ]
  sem <- bt_pred$Date[bt_pred$Run == run][1]
  
  p <- p[p$Lucro >= 0, ]
  if (nrow(p) > 1) p <- p[!is_dominated(p), ]
  p <- p[order(p$HR), ]
  
  if (nrow(p) == 0) {
    plot(1, type = "n",
         main = sprintf("Run %d | %s\n(sem soluções válidas)", run, sem),
         xlab = "", ylab = "", cex.main = 0.85)
    next
  }
  
  plot(p$HR, p$Lucro,
       type = "b", pch = 16, cex = 0.8,
       col  = "darkorange", lwd = 1.5,
       xlab = "HR", ylab = "Lucro (USD)",
       main = sprintf("Run %d | %s", run, sem),
       cex.main = 0.85)
  
  abline(h = median(p$Lucro), col = "green", lty = 3, lwd = 1.5)
  abline(h = mean(p$Lucro),   col = "blue",  lty = 2, lwd = 1.5)
  
  legend("top", bty = "n", cex = 0.7, horiz = TRUE,
         legend = c(sprintf("HV=%.0f",    hv_vec[run]),
                    sprintf("Média=%d",   round(mean(p$Lucro))),
                    sprintf("Mediana=%d", round(median(p$Lucro)))),
         col    = c("darkorange", "blue", "green"),
         lty    = c(1, 2, 3), lwd = c(1.5, 1.5, 1.5))
  
  grid()
}

mtext("Fronteiras de Pareto SANN O3 — 20 semanas GW",
      outer = TRUE, cex = 1.3, font = 2)

dev.off()
cat("Ficheiro 'SANN_O3_Pareto_GW.pdf' guardado com sucesso.\n")