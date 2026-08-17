# dss_export_O3_SANN.R
# Corre SANN O3 com soma ponderada e grava resultados para o DSS


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

get_metrics_O3 <- function(s) {
  stores       <- make_stores()
  total_profit <- 0
  total_hr     <- 0
  for (i in 1:4) {
    idx <- ((i-1)*21 + 1):(i*21)
    res <- eval(s[idx], stores[[i]])
    total_profit <- total_profit + res$profit
    J_vals   <- round(s[idx][seq(1, 21, by = 3)])
    X_vals   <- round(s[idx][seq(2, 21, by = 3)])
    total_hr <- total_hr + sum(J_vals) + sum(X_vals)
  }
  list(profit = total_profit, hr = total_hr)
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

is_dominated <- function(df) {
  n   <- nrow(df)
  dom <- logical(n)
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

# Ponto de referência fixo, o pior ponto de ambos os algortimos para o O3
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
    width  <- ref[1] - points[i, 1]          # Equivalente ao Lucro do ponto
    height <- current_ref2 - points[i, 2]   # Diferença de HR para a fatia anterior
    if (width > 0 && height > 0) {
      hv <- hv + width * height
    }
    current_ref2 <- points[i, 2]
  }
  return(hv)
}

MAX_PROFIT_EXPECTED <- 10000
MAX_HR_EXPECTED     <- 500

neg_eval_O3_pesos <- function(s, w_lucro, w_hr) {
  stores       <- make_stores()
  total_profit <- 0
  total_units  <- 0
  total_hr     <- 0
  for (i in 1:4) {
    idx <- ((i-1)*21 + 1):(i*21)
    res <- eval(s[idx], stores[[i]])
    total_profit <- total_profit + res$profit
    total_units  <- total_units  + res$units
    J_vals   <- round(s[idx][seq(1, 21, by = 3)])
    X_vals   <- round(s[idx][seq(2, 21, by = 3)])
    total_hr <- total_hr + sum(J_vals) + sum(X_vals)
  }
  if (total_units > 10000) {
    total_profit <- total_profit - (total_units - 10000) * 100
  }
  norm_profit <- total_profit / MAX_PROFIT_EXPECTED
  norm_hr     <- total_hr     / MAX_HR_EXPECTED
  fitness     <- (w_lucro * norm_profit) - (w_hr * norm_hr)
  return(-fitness)
}

rchange_sann <- function(par, ...) {
  lim_inf <- get("lower", envir = .GlobalEnv)
  lim_sup <- get("upper", envir = .GlobalEnv)
  hchange(par, lower = lim_inf, upper = lim_sup,
          operator = "*", dist = rnorm, mean = 1, sd = 0.05,
          round = FALSE)
}


# Configuração
pesos_lucro_testar <- seq(0.1, 0.9, by = 0.1)
N_O3 <- 10000
temp <- 0.05

resultados_O3_sann <- list()

cat(sprintf("\nDSS Export — SANN O3 | %d semanas\n\n", RUNS))

# Loop principal — 20 semanas
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
  
  cat(sprintf("Run %2d | %s a %s\n", run, bt_run$Date[1], bt_run$Date[7]))
  
  pontos_pareto <- list()
  
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
    
    cat(sprintf("  w=%.1f | Lucro: %d USD | HR: %d\n",
                w_lucro, round(metrics$profit), round(metrics$hr)))
    
    pontos_pareto[[length(pontos_pareto) + 1]] <- list(
      peso_lucro = w_lucro,
      peso_hr    = w_hr,
      lucro      = metrics$profit,
      hr         = metrics$hr,
      plano      = SA_O3$par,
      detalhes   = extrair_detalhes(SA_O3$par)
    )
  }
  
  lucros <- sapply(pontos_pareto, function(p) p$lucro)
  hrs    <- sapply(pontos_pareto, function(p) p$hr)
  
  # filtrar negativos
  validos       <- lucros >= 0
  pontos_pareto <- pontos_pareto[validos]
  lucros        <- lucros[validos]
  hrs           <- hrs[validos]
  
  # filtrar dominados
  if (length(pontos_pareto) > 1) {
    df_check  <- data.frame(Lucro = lucros, HR = hrs)
    nao_dom   <- !is_dominated(df_check)
    pontos_pareto <- pontos_pareto[nao_dom]
    lucros        <- lucros[nao_dom]
    hrs           <- hrs[nao_dom]
  }
  
  n_pontos <- length(pontos_pareto)
  cat(sprintf("  %d pontos Pareto (após filtros)\n", n_pontos))
  
  if (n_pontos > 0) {
    media_lucro   <- mean(lucros)
    mediana_lucro <- median(lucros)
    
    if (max(hrs) >= REF_HR) {
      cat(sprintf("  HR máximo (%d) >= REF_HR (%d)\n",
                  round(max(hrs)), REF_HR))
    }
    
    hv <- calc_hv(lucros, hrs, REF_PROFIT, REF_HR)
    
    cat(sprintf("  Lucro máx: %d USD | HR mín: %d | HV: %.2f\n\n",
                round(max(lucros)), round(min(hrs)), hv))
  } else {
    cat("  Sem soluções válidas nesta semana.\n\n")
    media_lucro   <- NA
    mediana_lucro <- NA
    hv            <- NA
  }
  
  resultados_O3_sann[[run]] <- list(
    semana        = bt_run$Date[1],
    previsoes     = data.frame(
      Data         = bt_run$Date,
      Weekday      = bt_run$Weekday,
      Baltimore    = C_baltimore,
      Lancaster    = C_lancaster,
      Philadelphia = C_philadelphia,
      Richmond     = C_richmond
    ),
    pareto_profit = lucros,
    pareto_HR     = hrs,
    media_lucro   = media_lucro,
    mediana_lucro = mediana_lucro,
    hv            = hv,
    pontos_pareto = pontos_pareto
  )
  
  saveRDS(resultados_O3_sann, "resultados_O3_SANN.rds")
  cat(sprintf("  RDS gravado após run %d.\n\n", run))
}

cat("Export concluído.\n")
cat("  - resultados_O3_SANN.rds\n")