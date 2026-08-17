# dss_export_O1_GA.R
# Corre GA Fase 2 e grava resultados para o DSS

library(GA)
source("hill.R")
source("avaliacao.R")

cat("A carregar previsões GW...\n")
bt_pred <- read.csv("previsoes_gw_baltimore.csv",    header = TRUE)
lc_pred <- read.csv("previsoes_gw_lancaster.csv",    header = TRUE)
ph_pred <- read.csv("previsoes_gw_philadelphia.csv", header = TRUE)
rv_pred <- read.csv("previsoes_gw_richmond.csv",     header = TRUE)

RUNS <- max(bt_pred$Run)

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

# Função que extrai o data.frame de detalhes diários de um plano
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
        Loja        = lojas_nomes[i],
        Dia         = d,
        J           = J,
        X           = X,
        PR          = round(PR, 2),
        Clientes    = C[d],
        Unidades    = UX + UJ,
        Vendas_USD  = revenue,
        Custo_HR    = cost_hr,
        Lucro_Dia   = revenue - cost_hr
      )
    }
  }
  do.call(rbind, rows)
}


# Configuração do GA
Runs_GA     <- 5
popSize     <- 200
maxiter     <- 600
pmutation   <- 0.15
pcrossover  <- 0.65
elitism_val <- round(popSize * 0.10)

resultados_O1 <- list()

cat(sprintf("\nDSS Export — GA O1 | %d semanas\n\n", RUNS))

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
    lucro               <- GA_O1@fitnessValue
    curvas_O1_runs[[r]] <- GA_O1@summary[, 1]
    cat(sprintf("  Run %d: %d USD\n", r, round(lucro)))
    
    if (lucro > best_global_O1) {
      best_global_O1     <- lucro
      best_sol_global_O1 <- GA_O1@solution[1, ]
    }
  }
  
  cat(sprintf("  >> Melhor: %d USD\n\n", round(best_global_O1)))
  
  # Guardar tudo o que o DSS precisa para esta semana
  resultados_O1[[run]] <- list(
    semana    = bt_run$Date[1],
    previsoes = data.frame(
      Data      = bt_run$Date,
      Weekday   = bt_run$Weekday,
      Baltimore = C_baltimore,
      Lancaster = C_lancaster,
      Philadelphia = C_philadelphia,
      Richmond  = C_richmond
    ),
    lucro     = best_global_O1,
    plano     = best_sol_global_O1,
    curvas    = curvas_O1_runs,
    detalhes  = extrair_detalhes(best_sol_global_O1)
  )
}

saveRDS(resultados_O1, "resultados_O1_GA.rds")
cat("Ficheiro 'resultados_O1_GA.rds' gravado.\n")