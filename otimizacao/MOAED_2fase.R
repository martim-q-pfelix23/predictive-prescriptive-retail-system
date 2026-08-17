# MOEA/D — Otimização Multi-Objetivo (O3)
# Fase 2: 20 semanas do Growing Window

library(MOEADr)
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

moead_eval <- function(X, ...) {
  Y <- matrix(0, nrow = nrow(X), ncol = 2)
  for (i in 1:nrow(X)) {
    s <- X[i, ]
    stores <- make_stores()
    total_profit <- 0
    total_units  <- 0
    total_hr     <- 0
    
    for (j in 1:4) {
      idx <- ((j-1)*21 + 1):(j*21)
      res <- eval(s[idx], stores[[j]])
      total_profit <- total_profit + res$profit
      total_units  <- total_units  + res$units
      
      J_vals <- round(s[idx][seq(1, 21, by = 3)])
      X_vals <- round(s[idx][seq(2, 21, by = 3)])
      total_hr <- total_hr + sum(J_vals) + sum(X_vals)
    }
    
    if (total_units > 10000) {
      total_profit <- -1e9
    }
    
    Y[i, 1] <- -total_profit
    Y[i, 2] <- total_hr
  }
  return(Y)
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

is_dominated <- function(df) {
  n   <- nrow(df)
  dom <- logical(n)
  if(n < 2) return(dom)
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


# 3. Configuração MOEA/D e variáveis
gens <- 200

REF_PROFIT <- 11
REF_HR     <- 55

media_profit_vec   <- numeric(RUNS)
mediana_profit_vec <- numeric(RUNS)
media_HR_vec       <- numeric(RUNS)
mediana_HR_vec     <- numeric(RUNS)
max_profit_vec     <- numeric(RUNS)
n_pareto_vec       <- numeric(RUNS)
hv_vec             <- numeric(RUNS)
best_units_vec     <- numeric(RUNS)

all_pareto         <- list()
all_profits_global <- c()
all_HR_global      <- c()

cat(sprintf("\nMOEA/D O3 - %d semanas (GW)\n", RUNS))
cat(sprintf("Gerações: %d | Decomposição: 100 subproblemas\n", gens))
flush.console()


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
  
  lower <<- c(b_bt$lower, b_lc$lower, b_ph$lower, b_rv$lower)
  upper <<- c(b_bt$upper, b_lc$upper, b_ph$upper, b_rv$upper)
  
  cat(sprintf("Run %2d | %s a %s\n", run, bt_run$Date[1], bt_run$Date[7]))
  flush.console()
  
  problem <- list(name = "moead_eval", xmin = lower, xmax = upper, m = 2)
  decomp  <- list(name = "sld", H = 99)
  stopcrt <- list(list(name = "maxiter", maxiter = gens))
  
  set.seed(42)
  
  MOEAD_out <- moead(
    problem  = problem,
    decomp   = decomp,
    stopcrit = stopcrt,
    preset   = preset_moead("original"),
    showpars = list(show.iters = "dots", showevery = 10)
  )
  
  df_moead <- data.frame(
    ID    = 1:nrow(MOEAD_out$Y),
    Lucro = -MOEAD_out$Y[, 1],
    HR    = MOEAD_out$Y[, 2]
  )
  
  df_moead <- df_moead[df_moead$Lucro >= 0, ]
  
  if (nrow(df_moead) > 1) {
    df_moead <- df_moead[!is_dominated(df_moead), ]
  }
  
  # Remover soluções duplicadas no espaço objetivo
  df_moead <- df_moead[!duplicated(df_moead[, c("Lucro", "HR")]), ]
  
  df_moead <- df_moead[order(df_moead$HR), ]
  
  if (nrow(df_moead) == 0) {
    cat("\n  Sem soluções Pareto válidas (Lucro >= 0) nesta semana.\n\n")
    all_pareto[[run]]    <- NULL
    best_units_vec[run]  <- NA
    flush.console()
    next
  }
  
  best_idx_run        <- df_moead$ID[which.max(df_moead$Lucro)]
  best_sol_run        <- MOEAD_out$X[best_idx_run, ]
  best_units_run      <- get_units_total(best_sol_run)
  best_units_vec[run] <- best_units_run
  
  all_profits_global <- c(all_profits_global, df_moead$Lucro)
  all_HR_global      <- c(all_HR_global,      df_moead$HR)
  
  # HV com referência global fixa (igual ao NSGA-II e SANN)
  hv_vec[run] <- calc_hv(df_moead$Lucro, df_moead$HR, REF_PROFIT, REF_HR)
  
  n_pareto_vec[run]       <- nrow(df_moead)
  max_profit_vec[run]     <- max(df_moead$Lucro)
  media_profit_vec[run]   <- mean(df_moead$Lucro)
  mediana_profit_vec[run] <- median(df_moead$Lucro)
  media_HR_vec[run]       <- mean(df_moead$HR)
  mediana_HR_vec[run]     <- median(df_moead$HR)
  
  all_pareto[[run]] <- list(
    profit      = df_moead$Lucro,
    HR          = df_moead$HR,
    data_inicio = bt_run$Date[1],
    data_fim    = bt_run$Date[7]
  )
  
  cat(sprintf("\n  Pareto Limpo: %d soluções\n", n_pareto_vec[run]))
  cat(sprintf("  Lucro  — Média: %d  Mediana: %d  Máx: %d USD (com %d Unidades)\n",
              round(media_profit_vec[run]), round(mediana_profit_vec[run]),
              round(max_profit_vec[run]),   round(best_units_run)))
  cat(sprintf("  HR     — Média: %d  Mediana: %d\n",
              round(media_HR_vec[run]), round(mediana_HR_vec[run])))
  cat(sprintf("  Hipervolume (ref global): %.2f\n\n", hv_vec[run]))
  flush.console()
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

cat("RESUMO MOEA/D O3 — FASE 2 (20 semanas)\n")
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
  Unidades      = ifelse(is.na(best_units_vec), NA, round(best_units_vec)),
  HV_global     = round(hv_global_vec, 2)
)
cat("\nTabela completa por semana:\n")
print(resultados_df, row.names = FALSE)

# 8. Gráficos — 20 Curvas de pareto
pdf("MOEAD_O3_Pareto_GW.pdf", width = 20, height = 16)
par(mfrow = c(4, 5), mar = c(4, 3, 4, 1), oma = c(2, 1, 4, 1))

for (run in 1:RUNS) {
  if (is.null(all_pareto[[run]])) {
    plot(1, type = "n", main = sprintf("Run %d\n(sem soluções válidas)", run),
         xlab = "", ylab = "", cex.main = 0.85)
    next
  }
  
  p <- all_pareto[[run]]
  plot(p$HR, p$profit,
       type = "b", pch = 16, cex = 0.7, col = "purple",
       xlab = "HR", ylab = "Lucro (USD)",
       main = sprintf("Run %d | %s", run, p$data_inicio),
       cex.main = 0.85)
  
  abline(h = media_profit_vec[run],   col = "blue",  lty = 2, lwd = 1.5)
  abline(h = mediana_profit_vec[run], col = "green", lty = 3, lwd = 1.5)
  
  legend("top", bty = "n", cex = 0.7, horiz = TRUE,
         legend = c(sprintf("HV=%.0f",    hv_global_vec[run]),
                    sprintf("Média=%d",   round(media_profit_vec[run])),
                    sprintf("Mediana=%d", round(mediana_profit_vec[run]))),
         col    = c("purple", "blue", "green"),
         lty    = c(1, 2, 3), lwd = c(1.5, 1.5, 1.5))
  
  grid()
}

mtext("Curvas de Pareto MOEA/D — O3 (20 semanas GW)", outer = TRUE, cex = 1.3, font = 2)
dev.off()

cat("\nFicheiro 'MOEAD_O3_Pareto_GW.pdf' guardado com sucesso.\n")