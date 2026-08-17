# SIMULATED ANNEALING — Otimização Fase 2 (O1 e O2)
# 20 semanas do Growing Window

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

rchange_sann <- function(par) {
  hchange(par, lower = lower, upper = upper,
          operator = "*", dist = rnorm, mean = 1, sd = 0.05,
          round = FALSE)
}


# 3. Configuração
Runs_SA            <- 5
N                  <- 10000
temperaturas_teste <- c(10, 250, 500, 750, 1000)
cores_grafico      <- c("red", "blue", "darkgreen", "orange", "purple")

lucros_O1_gw    <- numeric(RUNS)
lucros_O2_gw    <- numeric(RUNS)
unidades_O2_gw  <- numeric(RUNS)
best_temp_O1_gw <- numeric(RUNS)
best_temp_O2_gw <- numeric(RUNS)

curvas_O1_gw <- list()
curvas_O2_gw <- list()

cat(sprintf("\nSimulated Annealing O1 + O2 - %d semanas (GW)\n", RUNS))
cat(sprintf("Temperaturas: %s | Corridas por temp: %d | Iterações: %d\n\n",
            paste(temperaturas_teste, collapse=", "), Runs_SA, N))


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
  
  cat(sprintf("Run %2d | %s a %s\n", run, bt_run$Date[1], bt_run$Date[7]))

  # Funções de Monitorização
  EV      <- 0
  BEST    <- -Inf
  F_curva <- numeric()
  
  m_neg_eval_O1 <- function(s) {
    v <- eval_O1(s)
    EV    <<- EV + 1
    if (v > BEST) BEST <<- v
    if (EV <= N)  F_curva[EV] <<- BEST
    return(-v)
  }
  
  m_neg_eval_O2 <- function(s) {
    v <- eval_O2(s)
    EV <<- EV + 1
    if (!is.infinite(v) && v > BEST) BEST <<- v
    if (EV <= N) F_curva[EV] <<- BEST
    if (is.infinite(v) && v == -Inf) return(Inf)
    return(-v)
  }
  
  # O1 — estudo de temperaturas
  cat("O1\n")
  
  resultados_O1       <- list()
  curvas_O1_temp      <- list()
  best_global_O1      <- -Inf
  best_sol_global_O1  <- NULL
  best_temp_global_O1 <- NA
  
  for (t_idx in 1:length(temperaturas_teste)) {
    temp_atual      <- temperaturas_teste[t_idx]
    best_temp_O1    <- -Inf
    best_curva_temp <- NULL
    
    for (i in 1:Runs_SA) {
      set.seed(i * 10)
      EV      <<- 0
      BEST    <<- -Inf
      F_curva <<- rep(NA, N)
      
      s_init <- runif(length(lower), lower, upper)
      
      SA <- optim(par     = s_init,
                  fn      = m_neg_eval_O1,
                  gr      = rchange_sann,
                  method  = "SANN",
                  control = list(maxit = N, temp = temp_atual, trace = FALSE))
      
      lucro <- -SA$value
      if (lucro > best_temp_O1) {
        best_temp_O1    <- lucro
        best_curva_temp <- F_curva
      }
      if (lucro > best_global_O1) {
        best_global_O1      <- lucro
        best_sol_global_O1  <- SA$par
        best_temp_global_O1 <- temp_atual
      }
    }
    
    resultados_O1[[t_idx]]  <- best_temp_O1
    curvas_O1_temp[[t_idx]] <- best_curva_temp
    cat(sprintf("    Temp %4d: %d USD\n", temp_atual, round(best_temp_O1)))
  }
  
  lucros_O1_gw[run]    <- best_global_O1
  best_temp_O1_gw[run] <- best_temp_global_O1
  curvas_O1_gw[[run]]  <- curvas_O1_temp
  
  cat(sprintf("  >> Melhor O1: %d USD (Temp %d)\n\n",
              round(best_global_O1), best_temp_global_O1))
  
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
  
  resultados_O2       <- list()
  curvas_O2_temp      <- list()
  best_global_O2      <- -Inf
  best_sol_global_O2  <- NULL
  best_temp_global_O2 <- NA
  
  for (t_idx in 1:length(temperaturas_teste)) {
    temp_atual      <- temperaturas_teste[t_idx]
    best_temp_O2    <- -Inf
    best_curva_temp <- NULL
    
    for (i in 1:Runs_SA) {
      set.seed(i * 10)
      EV      <<- 0
      BEST    <<- -Inf
      F_curva <<- rep(NA, N)
      
      SA <- optim(par     = s_smart,
                  fn      = m_neg_eval_O2,
                  gr      = rchange_sann,
                  method  = "SANN",
                  control = list(maxit = N, temp = temp_atual, trace = FALSE))
      
      lucro <- -SA$value
      if (is.infinite(lucro)) lucro <- -Inf
      if (lucro > best_temp_O2) {
        best_temp_O2    <- lucro
        best_curva_temp <- F_curva
      }
      if (lucro > best_global_O2) {
        best_global_O2      <- lucro
        best_sol_global_O2  <- SA$par
        best_temp_global_O2 <- temp_atual
      }
    }
    
    resultados_O2[[t_idx]]  <- best_temp_O2
    curvas_O2_temp[[t_idx]] <- best_curva_temp
    cat(sprintf("    Temp %4d: %s USD\n", temp_atual,
                ifelse(is.infinite(best_temp_O2), "N/A", as.character(round(best_temp_O2)))))
  }
  
  lucros_O2_gw[run]    <- best_global_O2
  best_temp_O2_gw[run] <- best_temp_global_O2
  curvas_O2_gw[[run]]  <- curvas_O2_temp
  
  # unidades da melhor solução O2 encontrada
  if (!is.null(best_sol_global_O2) && is.finite(best_global_O2)) {
    unidades_O2_gw[run] <- get_units_total(best_sol_global_O2)
  } else {
    unidades_O2_gw[run] <- NA
  }
  
  cat(sprintf("  >> Melhor O2: %s USD (Temp %s) | Unidades: %s\n\n",
              ifelse(is.infinite(best_global_O2), "N/A", as.character(round(best_global_O2))),
              ifelse(is.na(best_temp_global_O2),  "N/A", as.character(best_temp_global_O2)),
              ifelse(is.na(unidades_O2_gw[run]),  "N/A", as.character(round(unidades_O2_gw[run])))))
}

# 5. Agregação final
cat("RESUMO SIMULATED ANNEALING — FASE 2 (20 semanas)\n")

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
  SA_O1       = round(lucros_O1_gw),
  Temp_O1     = best_temp_O1_gw,
  SA_O2       = ifelse(is.finite(lucros_O2_gw), round(lucros_O2_gw), NA),
  Temp_O2     = best_temp_O2_gw,
  Unidades_O2 = ifelse(is.na(unidades_O2_gw), NA, round(unidades_O2_gw))
)
cat("\nTabela completa:\n")
print(resultados_df, row.names = FALSE)

# 7. Gráficos - Um PDF por objetivo

# PDF O1
pdf("SANN_convergencia_O1_GW.pdf", width = 20, height = 16)
par(mfrow = c(4, 5), mar = c(3, 3, 3, 1), oma = c(2, 1, 4, 1))

for (run in 1:RUNS) {
  curvas <- curvas_O1_gw[[run]]
  vals   <- unlist(curvas)
  y_min  <- min(vals, na.rm = TRUE)
  y_max  <- max(vals, na.rm = TRUE)
  
  data_semana <- bt_pred$Date[bt_pred$Run == run][1]
  
  plot(1, type = "n", xlim = c(1, N), ylim = c(y_min, y_max),
       main = sprintf("Run %d | %s", run, data_semana),
       xlab = "Avaliações", ylab = "Lucro (USD)", cex.main = 0.85)
  
  for (t_idx in 1:length(temperaturas_teste)) {
    lines(curvas[[t_idx]], col = cores_grafico[t_idx], lwd = 1.5)
  }
  
  abline(h = lucros_O1_gw[run], col = "black", lty = 3, lwd = 1)
  grid()
}

mtext("SA Convergência O1 — 20 semanas GW", outer = TRUE, cex = 1.3, font = 2)
par(fig = c(0, 1, 0, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
legend("bottom", horiz = TRUE, bty = "n", cex = 0.85,
       legend = paste("Temp =", temperaturas_teste),
       col    = cores_grafico, lwd = 2)
dev.off()
cat("Ficheiro 'SANN_convergencia_O1_GW.pdf' guardado.\n")

# PDF O2
pdf("SANN_convergencia_O2_GW.pdf", width = 20, height = 16)
par(mfrow = c(4, 5), mar = c(3, 3, 3, 1), oma = c(2, 1, 4, 1))

for (run in 1:RUNS) {
  curvas <- curvas_O2_gw[[run]]
  vals   <- unlist(curvas)
  vals   <- vals[is.finite(vals)]
  
  data_semana <- bt_pred$Date[bt_pred$Run == run][1]
  
  if (length(vals) > 0) {
    y_min <- min(vals, na.rm = TRUE)
    y_max <- max(vals, na.rm = TRUE)
    
    plot(1, type = "n", xlim = c(1, N), ylim = c(y_min, y_max),
         main = sprintf("Run %d | %s", run, data_semana),
         xlab = "Avaliações", ylab = "Lucro (USD)", cex.main = 0.85)
    
    for (t_idx in 1:length(temperaturas_teste)) {
      curva <- curvas[[t_idx]]
      if (!is.null(curva) && any(is.finite(curva))) {
        lines(curva, col = cores_grafico[t_idx], lwd = 1.5)
      }
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

mtext("SA Convergência O2 — 20 semanas GW", outer = TRUE, cex = 1.3, font = 2)
par(fig = c(0, 1, 0, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
legend("bottom", horiz = TRUE, bty = "n", cex = 0.85,
       legend = paste("Temp =", temperaturas_teste),
       col    = cores_grafico, lwd = 2)
dev.off()
cat("Ficheiro 'SANN_convergencia_O2_GW.pdf' guardado.\n")