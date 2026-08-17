# GENETIC ALGORITHM — Otimização Fase 2 (O1 e O2)
# 20 semanas do Growing Window

library(GA)
source("avaliacao.R")


# 1. Cqrregar previsões GW
cat("A carregar previsões GW do modelo RF Cenário C...\n")

bt_pred <- read.csv("previsoes_gw_baltimore.csv",    header = TRUE)
lc_pred <- read.csv("previsoes_gw_lancaster.csv",    header = TRUE)
ph_pred <- read.csv("previsoes_gw_philadelphia.csv", header = TRUE)
rv_pred <- read.csv("previsoes_gw_richmond.csv",     header = TRUE)

RUNS <- max(bt_pred$Run)


# 2. Funções originais
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

get_units_total <- function(s) {
  stores      <- make_stores()
  total_units <- 0
  for (i in 1:4) {
    idx <- ((i-1)*21 + 1):(i*21)
    res <- eval(s[idx], stores[[i]])
    total_units <- total_units + res$units
  }
  return(total_units)
}

# 3. Configuração
Runs_GA     <- 5
popSize     <- 200
maxiter     <- 600
pmutation   <- 0.15
pcrossover  <- 0.65
elitism_val <- round(popSize * 0.10)

lucros_O1_gw   <- numeric(RUNS)
lucros_O2_gw   <- numeric(RUNS)
unidades_O2_gw <- numeric(RUNS)

curvas_O1_gw <- list()
curvas_O2_gw <- list()

cat(sprintf("\nGenetic Algorithm O1 + O2 - %d semanas (GW)\n", RUNS))
cat(sprintf("Corridas por semana: %d | Gerações máx.: %d | População: %d\n\n",
            Runs_GA, maxiter, popSize))


# 4. Loop principal - 20 semanas
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
  

  # O1

  cat("O1\n")
  
  best_global_O1     <- -Inf
  best_sol_global_O1 <- NULL
  curvas_O1_runs     <- list()
  
  for (r in 1:Runs_GA) {
    
    set.seed(r * 100)
    GA_O1 <- ga(
      type       = "real-valued",
      fitness    = eval_O1,
      lower      = lower,
      upper      = upper,
      popSize    = popSize,
      maxiter    = maxiter,
      pmutation  = pmutation,
      pcrossover = pcrossover,
      elitism    = elitism_val,
      run        = 100,
      monitor    = FALSE
    )
    
    lucro <- GA_O1@fitnessValue
    cat(sprintf("    Run %d: %d USD\n", r, round(lucro)))
    
    curvas_O1_runs[[r]] <- GA_O1@summary[, 1]
    
    if (lucro > best_global_O1) {
      best_global_O1     <- lucro
      best_sol_global_O1 <- GA_O1@solution[1, ]
    }
  }
  
  lucros_O1_gw[run]   <- best_global_O1
  curvas_O1_gw[[run]] <- curvas_O1_runs
  
  cat(sprintf("  >> Melhor O1: %d USD\n\n", round(best_global_O1)))
  
  # O2 — ponto inicial via reparação do O1
  cat("O2\n")
  
  s_smart         <- best_sol_global_O1
  unidades_atuais <- get_units_total(s_smart)
  iteracoes_corte <- 0
  
  while (unidades_atuais > 10000) {
    for (i in 1:4) {
      idx    <- ((i-1)*21 + 1):(i*21)
      idx_J  <- seq(1, 21, by = 3)
      idx_X  <- seq(2, 21, by = 3)
      idx_PR <- seq(3, 21, by = 3)
      s_smart[idx[idx_J]]  <- floor(s_smart[idx[idx_J]]  * 0.95)
      s_smart[idx[idx_X]]  <- floor(s_smart[idx[idx_X]]  * 0.95)
      s_smart[idx[idx_PR]] <- s_smart[idx[idx_PR]] * 0.95
    }
    unidades_atuais <- get_units_total(s_smart)
    iteracoes_corte <- iteracoes_corte + 1
  }
  
  cat(sprintf("  Reparação: %d passos | Unidades: %d (VÁLIDO)\n",
              iteracoes_corte, round(unidades_atuais)))
  
  best_global_O2     <- -Inf
  best_sol_global_O2 <- NULL
  curvas_O2_runs     <- list()
  
  for (r in 1:Runs_GA) {
    
    suggestions <- rbind(s_smart, runif(length(lower), lower, upper), best_sol_global_O1)
    
    set.seed(r * 200)
    GA_O2 <- ga(
      type        = "real-valued",
      fitness     = eval_O2_ga,
      lower       = lower,
      upper       = upper,
      popSize     = popSize,
      maxiter     = maxiter,
      pmutation   = pmutation,
      pcrossover  = pcrossover,
      elitism     = elitism_val,
      suggestions = suggestions,
      run         = 100,
      monitor     = FALSE
    )
    
    lucro <- GA_O2@fitnessValue
    cat(sprintf("    Run %d: %d USD\n", r, round(lucro)))
    
    curvas_O2_runs[[r]] <- GA_O2@summary[, 1]
    
    if (lucro > best_global_O2) {
      best_global_O2     <- lucro
      best_sol_global_O2 <- GA_O2@solution[1, ]
    }
  }
  
  lucros_O2_gw[run]   <- best_global_O2
  curvas_O2_gw[[run]] <- curvas_O2_runs
  
  # unidades da melhor solução O2
  if (!is.null(best_sol_global_O2) && is.finite(best_global_O2)) {
    unidades_O2_gw[run] <- get_units_total(best_sol_global_O2)
  } else {
    unidades_O2_gw[run] <- NA
  }
  
  cat(sprintf("  >> Melhor O2: %d USD | Unidades: %s\n\n",
              round(best_global_O2),
              ifelse(is.na(unidades_O2_gw[run]), "N/A", as.character(round(unidades_O2_gw[run])))))
}


# 5. Agregação final
cat("RESUMO GENETIC ALGORITHM — FASE 2 (20 semanas)\n")

lucros_O1_fin <- lucros_O1_gw[is.finite(lucros_O1_gw)]
lucros_O2_fin <- lucros_O2_gw[is.finite(lucros_O2_gw)]
unidades_fin  <- unidades_O2_gw[!is.na(unidades_O2_gw)]

cat(sprintf("O1 - Média   dos lucros: %d USD\n", round(mean(lucros_O1_fin))))
cat(sprintf("O1 - Mediana dos lucros: %d USD\n", round(median(lucros_O1_fin))))
cat(sprintf("O1 - Mínimo:             %d USD\n", round(min(lucros_O1_fin))))
cat(sprintf("O1 - Máximo:             %d USD\n", round(max(lucros_O1_fin))))

cat("\n")

cat(sprintf("O2 - Média   dos lucros: %d USD\n", round(mean(lucros_O2_fin))))
cat(sprintf("O2 - Mediana dos lucros: %d USD\n", round(median(lucros_O2_fin))))
cat(sprintf("O2 - Mínimo:             %d USD\n", round(min(lucros_O2_fin))))
cat(sprintf("O2 - Máximo:             %d USD\n", round(max(lucros_O2_fin))))

cat("\n")

cat(sprintf("O2 - Média   das unidades vendidas: %d\n", round(mean(unidades_fin))))
cat(sprintf("O2 - Mediana das unidades vendidas: %d\n", round(median(unidades_fin))))


# 6. Tabela final
resultados_df <- data.frame(
  Run         = 1:RUNS,
  Semana      = sapply(1:RUNS, function(r) bt_pred$Date[bt_pred$Run == r][1]),
  GA_O1       = round(lucros_O1_gw),
  GA_O2       = ifelse(is.finite(lucros_O2_gw), round(lucros_O2_gw), NA),
  Unidades_O2 = ifelse(is.na(unidades_O2_gw), NA, round(unidades_O2_gw))
)
cat("\nTabela completa:\n")
print(resultados_df, row.names = FALSE)

# 7. Gráficos - 1 por objetivo
cores_grafico <- c("red", "blue", "darkgreen", "orange", "purple")

# PDF O1
pdf("GA_convergencia_O1_GW.pdf", width = 20, height = 16)
par(mfrow = c(4, 5), mar = c(3, 3, 3, 1), oma = c(2, 1, 4, 1))

for (run in 1:RUNS) {
  curvas      <- curvas_O1_gw[[run]]
  vals        <- unlist(curvas)
  y_min       <- min(vals, na.rm = TRUE)
  y_max       <- max(vals, na.rm = TRUE)
  max_gens    <- max(sapply(curvas, length))
  data_semana <- bt_pred$Date[bt_pred$Run == run][1]
  
  plot(1, type = "n", xlim = c(1, max_gens), ylim = c(y_min, y_max),
       main = sprintf("Run %d | %s", run, data_semana),
       xlab = "Geração", ylab = "Lucro (USD)", cex.main = 0.85)
  
  for (r in 1:Runs_GA) {
    lines(curvas[[r]], col = cores_grafico[r], lwd = 1.5)
  }
  
  abline(h = lucros_O1_gw[run], col = "black", lty = 3, lwd = 1)
  grid()
}

mtext("GA Convergência O1 — 20 semanas GW", outer = TRUE, cex = 1.3, font = 2)
par(fig = c(0, 1, 0, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
legend("bottom", horiz = TRUE, bty = "n", cex = 0.85,
       legend = paste("Run", 1:Runs_GA),
       col    = cores_grafico, lwd = 2)
dev.off()
cat("Ficheiro 'GA_convergencia_O1_GW.pdf' guardado.\n")

# PDF O2
pdf("GA_convergencia_O2_GW.pdf", width = 20, height = 16)
par(mfrow = c(4, 5), mar = c(3, 3, 3, 1), oma = c(2, 1, 4, 1))

for (run in 1:RUNS) {
  curvas      <- curvas_O2_gw[[run]]
  vals        <- unlist(curvas)
  data_semana <- bt_pred$Date[bt_pred$Run == run][1]
  max_gens    <- max(sapply(curvas, length))
  
  if (length(vals) > 0) {
    y_min <- min(vals, na.rm = TRUE)
    y_max <- max(vals, na.rm = TRUE)
    
    plot(1, type = "n", xlim = c(1, max_gens), ylim = c(y_min, y_max),
         main = sprintf("Run %d | %s", run, data_semana),
         xlab = "Geração", ylab = "Lucro (USD)", cex.main = 0.85)
    
    for (r in 1:Runs_GA) {
      lines(curvas[[r]], col = cores_grafico[r], lwd = 1.5)
    }
    
    if (is.finite(lucros_O2_gw[run])) {
      abline(h = lucros_O2_gw[run], col = "black", lty = 3, lwd = 1)
    }
    grid()
    
  } else {
    plot(1, type = "n",
         main = sprintf("Run %d | %s\n(sem solução válida)", run, data_semana),
         xlab = "", ylab = "", cex.main = 0.85)
  }
}

mtext("GA Convergência O2 — 20 semanas GW", outer = TRUE, cex = 1.3, font = 2)
par(fig = c(0, 1, 0, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
legend("bottom", horiz = TRUE, bty = "n", cex = 0.85,
       legend = paste("Run", 1:Runs_GA),
       col    = cores_grafico, lwd = 2)
dev.off()
cat("Ficheiro 'GA_convergencia_O2_GW.pdf' guardado.\n")