# MONTE CARLO - Otimização Fase 2 (O1 + O2)
# 20 semanas do Growing Window

source("blind.R")
source("montecarlo.R")
source("avaliacao.R")

# 1. Leitura das previsões GW
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
    
    upper_store <- c(
      upper_store,
      ceiling(C[d] / 6),
      ceiling(C[d] / 7),
      0.30
    )
  }
  
  list(lower = lower_store, upper = upper_store)
}


# total de unidades produzidas
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


# 3. MONTE CARLO — 20 SEMANAS
N <- 10000

lucros_O1   <- numeric(RUNS)
lucros_O2   <- numeric(RUNS)
unidades_O2 <- numeric(RUNS)

cat(sprintf("\nMonte Carlo O1 + O2 - %d semanas (GW)\n", RUNS))
cat("Número de amostras por semana:", N, "\n\n")

for (run in 1:RUNS) {
  
  # filtrar semana desta iteração
  bt_run <- bt_pred[bt_pred$Run == run, ]
  lc_run <- lc_pred[lc_pred$Run == run, ]
  ph_run <- ph_pred[ph_pred$Run == run, ]
  rv_run <- rv_pred[rv_pred$Run == run, ]
  

  # atualizar clientes globais
  C_baltimore    <<- bt_run$Predicted_Customers
  C_lancaster    <<- lc_run$Predicted_Customers
  C_philadelphia <<- ph_run$Predicted_Customers
  C_richmond     <<- rv_run$Predicted_Customers
  
  weekend_baltimore    <<- is_weekend(bt_run$Weekday)
  weekend_lancaster    <<- is_weekend(lc_run$Weekday)
  weekend_philadelphia <<- is_weekend(ph_run$Weekday)
  weekend_richmond     <<- is_weekend(rv_run$Weekday)
  

  # bounds
  b_bt <- calc_bounds(C_baltimore)
  b_lc <- calc_bounds(C_lancaster)
  b_ph <- calc_bounds(C_philadelphia)
  b_rv <- calc_bounds(C_richmond)
  
  lower <- c(b_bt$lower, b_lc$lower, b_ph$lower, b_rv$lower)
  upper <- c(b_bt$upper, b_lc$upper, b_ph$upper, b_rv$upper)
  

  # O1 — Monte Carlo original
  set.seed(123)
  
  MC_O1 <- mcsearch(
    fn    = eval_O1,
    lower = lower,
    upper = upper,
    N     = N,
    type  = "max"
  )
  
  lucros_O1[run] <- MC_O1$eval
  
  cat(sprintf(
    "Run %2d | Semana: %s a %s | O1: %d USD\n",
    run,
    bt_run$Date[1],
    bt_run$Date[7],
    round(MC_O1$eval)
  ))
  
  # O2 — repair da melhor solução do O1
  s_repair <- MC_O1$sol
  
  unidades_atuais <- get_units_total(s_repair)
  
  iteracoes_corte <- 0
  
  while (unidades_atuais > 10000) {
    
    for (i in 1:4) {
      
      idx <- ((i-1)*21 + 1):(i*21)
      
      idx_J  <- seq(1, 21, by = 3)
      idx_X  <- seq(2, 21, by = 3)
      idx_PR <- seq(3, 21, by = 3)
      
      s_repair[idx[idx_J]]  <- floor(s_repair[idx[idx_J]] * 0.95)
      s_repair[idx[idx_X]]  <- floor(s_repair[idx[idx_X]] * 0.95)
      s_repair[idx[idx_PR]] <- s_repair[idx[idx_PR]] * 0.95
    }
    
    unidades_atuais <- get_units_total(s_repair)
    
    iteracoes_corte <- iteracoes_corte + 1
  }
  
  lucro_O2        <- eval_O2(s_repair)
  lucros_O2[run]   <- lucro_O2
  unidades_O2[run] <- unidades_atuais
  
  cat(sprintf(
    "  Repair: %d passos | Unidades: %d | O2: %d USD\n\n",
    iteracoes_corte,
    round(unidades_atuais),
    round(lucro_O2)
  ))
}


# 4. Agregação dos resultados
cat("\nResultados agregados (Monte Carlo O1 + O2)\n")

cat(sprintf("O1 - Média   dos lucros: %d USD\n", round(mean(lucros_O1))))
cat(sprintf("O1 - Mediana dos lucros: %d USD\n", round(median(lucros_O1))))
cat(sprintf("O1 - Mínimo:             %d USD\n", round(min(lucros_O1))))
cat(sprintf("O1 - Máximo:             %d USD\n", round(max(lucros_O1))))

cat("\n")

cat(sprintf("O2 - Média   dos lucros: %d USD\n", round(mean(lucros_O2[is.finite(lucros_O2)]))))
cat(sprintf("O2 - Mediana dos lucros: %d USD\n", round(median(lucros_O2[is.finite(lucros_O2)]))))
cat(sprintf("O2 - Mínimo:             %d USD\n", round(min(lucros_O2[is.finite(lucros_O2)]))))
cat(sprintf("O2 - Máximo:             %d USD\n", round(max(lucros_O2[is.finite(lucros_O2)]))))

cat("\n")

cat(sprintf("O2 - Média   das unidades vendidas: %d\n", round(mean(unidades_O2))))
cat(sprintf("O2 - Mediana das unidades vendidas: %d\n", round(median(unidades_O2))))

# 5. Tabela final
resultados_df <- data.frame(
  Run         = 1:RUNS,
  Semana      = sapply(1:RUNS, function(r) bt_pred$Date[bt_pred$Run == r][1]),
  MC_O1       = round(lucros_O1),
  MC_O2       = round(lucros_O2),
  Unidades_O2 = round(unidades_O2)
)

cat("\nTabela completa:\n")
print(resultados_df, row.names = FALSE)