# dss_export_O3_NSGA2.R
# Corre NSGA-II O3 e grava resultados para o DSS

library(mco)
source("avaliacao.R")

cat("A carregar previsões GW...\n")
bt_pred <- read.csv("previsoes_gw_baltimore.csv",    header = TRUE)
lc_pred <- read.csv("previsoes_gw_lancaster.csv",    header = TRUE)
ph_pred <- read.csv("previsoes_gw_philadelphia.csv", header = TRUE)
rv_pred <- read.csv("previsoes_gw_richmond.csv",     header = TRUE)

RUNS <- max(bt_pred$Run)

# Ponto de referência fixo, o pior ponto de ambos os algortimos para o O3
REF_PROFIT <- 11
REF_HR     <- 55

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

# Configuração
popsize <- 200
gens    <- 600

resultados_O3 <- list()

cat(sprintf("\n=== DSS Export — NSGA-II O3 | %d semanas ===\n\n", RUNS))
cat(sprintf("Ponto de referência fixo: REF_PROFIT = %d | REF_HR = %d\n\n",
            REF_PROFIT, REF_HR))


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
  
  lower <- c(b_bt$lower, b_lc$lower, b_ph$lower, b_rv$lower)
  upper <- c(b_bt$upper, b_lc$upper, b_ph$upper, b_rv$upper)
  
  cat(sprintf("Run %2d | %s a %s\n", run, bt_run$Date[1], bt_run$Date[7]))
  
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
  
  # Filtrar soluções com lucro negativo
  valid_idx <- which(raw_profit > 0)
  
  if (length(valid_idx) == 0) {
    cat("  Sem soluções Pareto válidas (Lucro > 0) nesta semana.\n\n")
    resultados_O3[[run]] <- list(
      semana        = bt_run$Date[1],
      previsoes     = data.frame(
        Data         = bt_run$Date,
        Weekday      = bt_run$Weekday,
        Baltimore    = C_baltimore,
        Lancaster    = C_lancaster,
        Philadelphia = C_philadelphia,
        Richmond     = C_richmond
      ),
      pareto_profit = numeric(0),
      pareto_HR     = numeric(0),
      media_lucro   = NA,
      mediana_lucro = NA,
      hv            = NA,
      pontos_pareto = list()
    )
    saveRDS(resultados_O3, "resultados_O3_NSGA2.rds")
    next
  }
  
  pareto_profit <- raw_profit[valid_idx]
  pareto_HR     <- raw_HR[valid_idx]
  pareto_sols   <- raw_sols[valid_idx, , drop = FALSE]
  
  media_lucro   <- mean(pareto_profit)
  mediana_lucro <- median(pareto_profit)
  
  # Verificação: alertar se HR excede o ponto de referência
  if (max(pareto_HR) >= REF_HR) {
    cat(sprintf("  HR máximo (%d) >= REF_HR (%d)\n",
                round(max(pareto_HR)), REF_HR))
  }
  
  hv <- calc_hv(pareto_profit, pareto_HR, REF_PROFIT, REF_HR)
  
  n_pareto <- length(pareto_profit)
  cat(sprintf("  %d soluções Pareto | Lucro máx: %d | Média: %d | Mediana: %d | HR mín: %d | HV: %.2f\n\n",
              n_pareto, round(max(pareto_profit)), round(media_lucro),
              round(mediana_lucro), round(min(pareto_HR)), hv))
  
  pontos_pareto <- vector("list", n_pareto)
  for (k in 1:n_pareto) {
    s_k <- as.numeric(pareto_sols[k, ])
    pontos_pareto[[k]] <- list(
      lucro    = pareto_profit[k],
      hr       = pareto_HR[k],
      plano    = s_k,
      detalhes = extrair_detalhes(s_k)
    )
  }
  
  resultados_O3[[run]] <- list(
    semana        = bt_run$Date[1],
    previsoes     = data.frame(
      Data         = bt_run$Date,
      Weekday      = bt_run$Weekday,
      Baltimore    = C_baltimore,
      Lancaster    = C_lancaster,
      Philadelphia = C_philadelphia,
      Richmond     = C_richmond
    ),
    pareto_profit = pareto_profit,
    pareto_HR     = pareto_HR,
    media_lucro   = media_lucro,
    mediana_lucro = mediana_lucro,
    hv            = hv,
    pontos_pareto = pontos_pareto
  )
  
  saveRDS(resultados_O3, "resultados_O3_NSGA2.rds")
  cat(sprintf("  RDS gravado após run %d.\n\n", run))
}


cat("Export concluído.\n")
cat("  - resultados_O3_NSGA2.rds\n")