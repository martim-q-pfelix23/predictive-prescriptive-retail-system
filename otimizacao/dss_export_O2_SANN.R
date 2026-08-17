# dss_export_O2_SANN.R
# Corre SANN O1 + SANN O2 e grava apenas os resultados do O2 para o DSS


source("hill.R")
source("avaliacao.R")

cat("A carregar previsões GW...\n")
bt_pred <- read.csv("previsoes_gw_baltimore.csv",    header = TRUE)
lc_pred <- read.csv("previsoes_gw_lancaster.csv",    header = TRUE)
ph_pred <- read.csv("previsoes_gw_philadelphia.csv", header = TRUE)
rv_pred <- read.csv("previsoes_gw_richmond.csv",     header = TRUE)

RUNS <- max(bt_pred$Run)


# Funções auxiliares
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

extrair_detalhes <- function(s) {
  lojas_nomes <- c("Baltimore", "Lancaster", "Philadelphia", "Richmond")
  stores      <- make_stores()
  rows        <- list()
  for (i in 1:4) {
    idx   <- ((i-1)*21 + 1):(i*21)
    store <- stores[[i]]
    C     <- store$C
    FJ    <- store$FJ
    FX    <- store$FX
    for (d in 1:7) {
      J  <- round(s[idx[(d-1)*3 + 1]])
      X  <- round(s[idx[(d-1)*3 + 2]])
      PR <- s[idx[(d-1)*3 + 3]]
      AX  <- min(7*X, C[d])
      AJ  <- min(6*J, max(0, C[d] - AX))
      UXc <- round(FX * 10 / log(2 - PR))
      UJc <- round(FJ * 10 / log(2 - PR))
      UX  <- AX * UXc
      UJ  <- AJ * UJc
      SX  <- round(UX * (1 - PR) * 1.07)
      SJ  <- round(UJ * (1 - PR) * 1.07)
      revenue <- SX + SJ
      cost_hr <- if (store$weekend[d]) J*70 + X*95 else J*60 + X*80
      rows[[length(rows) + 1]] <- data.frame(
        Loja       = lojas_nomes[i],
        Dia        = d,
        J          = J,
        X          = X,
        PR         = round(PR, 2),
        Clientes   = C[d],
        Unidades   = UX + UJ,
        Vendas_USD = revenue,
        Custo_HR   = cost_hr,
        Lucro_Dia  = revenue - cost_hr
      )
    }
  }
  do.call(rbind, rows)
}

rchange_sann <- function(par) {
  hchange(par, lower = lower, upper = upper,
          operator = "*", dist = rnorm, mean = 1, sd = 0.05,
          round = FALSE)
}


# Configuração
Runs_SA            <- 5
N                  <- 10000
temperaturas_teste <- c(10, 250, 500, 750, 1000)

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

resultados_O2 <- list()

cat(sprintf("\nDSS Export - SANN O2 | %d semanas\n\n", RUNS))


# Loop Principal — 20 Semanas
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
  
  # SANN O1
  cat("O1 (SANN interno para reparação)\n")
  
  best_global_O1     <- -Inf
  best_sol_global_O1 <- NULL
  
  for (t_idx in 1:length(temperaturas_teste)) {
    temp_atual   <- temperaturas_teste[t_idx]
    best_temp_O1 <- -Inf
    
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
      if (lucro > best_temp_O1) best_temp_O1 <- lucro
      if (lucro > best_global_O1) {
        best_global_O1     <- lucro
        best_sol_global_O1 <- SA$par
      }
    }
    cat(sprintf("    Temp %4d: %d USD\n", temp_atual, round(best_temp_O1)))
  }
  
  cat(sprintf("  >> Melhor O1 interno: %d USD\n\n", round(best_global_O1)))
  

  # Reparação — igual ao SANN original
  cat("Reparação\n")
  
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
  cat(sprintf("  Lucro inicial O2: %d USD\n\n", round(eval_O2(s_smart))))
  

  # SANN O2 — estudo de temperaturas
  cat("O2 (SANN)\n")
  
  best_global_O2      <- -Inf
  best_sol_global_O2  <- NULL
  best_temp_global_O2 <- NA
  curvas_O2_temps     <- list()
  
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
    
    curvas_O2_temps[[t_idx]] <- best_curva_temp
    cat(sprintf("    Temp %4d: %s USD\n", temp_atual,
                ifelse(is.infinite(best_temp_O2), "N/A",
                       as.character(round(best_temp_O2)))))
  }
  
  cat(sprintf("  >> Melhor O2: %d USD (Temp %d)\n\n",
              round(best_global_O2), best_temp_global_O2))
  
  resultados_O2[[run]] <- list(
    semana     = bt_run$Date[1],
    previsoes  = data.frame(
      Data         = bt_run$Date,
      Weekday      = bt_run$Weekday,
      Baltimore    = C_baltimore,
      Lancaster    = C_lancaster,
      Philadelphia = C_philadelphia,
      Richmond     = C_richmond
    ),
    lucro      = best_global_O2,
    unidades   = get_units_total(best_sol_global_O2),
    temp_usada = best_temp_global_O2,
    plano      = best_sol_global_O2,
    curvas     = curvas_O2_temps,
    detalhes   = extrair_detalhes(best_sol_global_O2)
  )
  
  # Gravar após cada semana
  saveRDS(resultados_O2, "resultados_O2_SANN.rds")
  cat(sprintf("  RDS gravado após run %d.\n\n", run))
}

cat("Export concluído.\n")
cat("  - resultados_O2_SANN.rds\n")