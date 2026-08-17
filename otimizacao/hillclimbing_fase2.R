# HILL CLIMBING — Otimização Fase 2 (O1 e O2)

source("blind.R")
source("montecarlo.R")
source("hill.R")
source("avaliacao.R")


# 1. Carregar previsões GW
cat("A carregar previsões GW do modelo RF Cenário C...\n")

bt_pred <- read.csv("previsoes_gw_baltimore.csv",    header = TRUE)
lc_pred <- read.csv("previsoes_gw_lancaster.csv",    header = TRUE)
ph_pred <- read.csv("previsoes_gw_philadelphia.csv", header = TRUE)
rv_pred <- read.csv("previsoes_gw_richmond.csv",     header = TRUE)

RUNS <- max(bt_pred$Run)


# 2. FUNÇÕES AUXILIARES 
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

rchange_mult <- function(par, lower, upper) {
  hchange(par, lower = lower, upper = upper,
          operator = "*", dist = rnorm, mean = 1, sd = 0.05,
          round = FALSE)
}

# 3. HILL CLIMBING — 20 Semanas
N   <- 10000
REP <- N / 10

lucros_O1   <- numeric(RUNS)
lucros_O2   <- numeric(RUNS)
unidades_O2 <- numeric(RUNS)

cat(sprintf("\nHill Climbing O1 + O2 - %d semanas (GW)\n", RUNS))
cat("Iterações por semana:", N, "\n\n")

for (run in 1:RUNS) {
  
  # filtrar a semana desta iteração
  bt_run <- bt_pred[bt_pred$Run == run, ]
  lc_run <- lc_pred[lc_pred$Run == run, ]
  ph_run <- ph_pred[ph_pred$Run == run, ]
  rv_run <- rv_pred[rv_pred$Run == run, ]
  
  # atualizar ambiente global (make_stores() lê daqui)
  C_baltimore    <<- bt_run$Predicted_Customers
  C_lancaster    <<- lc_run$Predicted_Customers
  C_philadelphia <<- ph_run$Predicted_Customers
  C_richmond     <<- rv_run$Predicted_Customers
  
  weekend_baltimore    <<- is_weekend(bt_run$Weekday)
  weekend_lancaster    <<- is_weekend(lc_run$Weekday)
  weekend_philadelphia <<- is_weekend(ph_run$Weekday)
  weekend_richmond     <<- is_weekend(rv_run$Weekday)
  
  # bounds para esta semana
  b_bt <- calc_bounds(C_baltimore)
  b_lc <- calc_bounds(C_lancaster)
  b_ph <- calc_bounds(C_philadelphia)
  b_rv <- calc_bounds(C_richmond)
  
  lower <- c(b_bt$lower, b_lc$lower, b_ph$lower, b_rv$lower)
  upper <- c(b_bt$upper, b_lc$upper, b_ph$upper, b_rv$upper)
  
  cat(sprintf("Run %2d | %s a %s\n",
              run, bt_run$Date[1], bt_run$Date[7]))
  
  # O1
  set.seed(123)
  MC_init <- mcsearch(fn = eval_O1, lower = lower, upper = upper,
                      N = 5000, type = "max")
  s0 <- MC_init$sol
  
  set.seed(42)
  HC_O1 <- hclimbing(
    par     = s0,
    fn      = eval_O1,
    change  = rchange_mult,
    lower   = lower,
    upper   = upper,
    type    = "max",
    control = list(maxit = N, REPORT = 0, digits = 0)
  )
  lucros_O1[run] <- HC_O1$eval
  cat(sprintf("  O1 Lucro: %d USD\n", round(HC_O1$eval)))
  
  # O2
  s_smart <- HC_O1$sol
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
  s0_O2 <- s_smart
  
  set.seed(42)
  HC_O2 <- hclimbing(
    par     = s0_O2,
    fn      = eval_O2,
    change  = rchange_mult,
    lower   = lower,
    upper   = upper,
    type    = "max",
    control = list(maxit = N, REPORT = 0, digits = 0)
  )
  
  lucros_O2[run]   <- HC_O2$eval
  unidades_O2[run] <- get_units_total(HC_O2$sol)
  
  cat(sprintf("  O2 Lucro: %d USD | Unidades: %d | Reparação: %d passos\n\n",
              round(HC_O2$eval), round(unidades_O2[run]), iteracoes_corte))
}

# 4. Agregação dos resultados
cat("RESUMO HILL CLIMBING — FASE 2 (20 semanas)\n")
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
  HC_O1       = round(lucros_O1),
  HC_O2       = round(lucros_O2),
  Unidades_O2 = round(unidades_O2)
)
cat("\nTabela completa:\n")
print(resultados_df, row.names = FALSE)

cat(sprintf("\n%-6s %-12s %8d %8d %8d\n", "Média",   "",
            round(mean(lucros_O1)),   round(mean(lucros_O2)),   round(mean(unidades_O2))))
cat(sprintf("%-6s %-12s %8d %8d %8d\n",  "Mediana", "",
            round(median(lucros_O1)), round(median(lucros_O2)), round(median(unidades_O2))))


# 6. Gráfico de desempenho ao longo das 20 semanas
cat("\nA gerar gráfico de lucros ao longo das 20 semanas...\n")

limite_y <- range(c(resultados_df$HC_O1, resultados_df$HC_O2))

plot(resultados_df$Run, resultados_df$HC_O1,
     type = "b",
     pch  = 16,
     col  = "blue",
     ylim = limite_y,
     main = "Lucro Semanal - Hill Climbing (Growing Window)",
     xlab = "Semanas Avaliadas (Runs)",
     ylab = "Lucro Final (USD)",
     lwd  = 2)

lines(resultados_df$Run, resultados_df$HC_O2,
      type = "b",
      pch  = 15,
      col  = "red",
      lwd  = 2)

grid()

legend("topleft", bty = "n",
       legend = c("O1 (Sem restrições)", "O2 (<= 10.000 unidades)"),
       col    = c("blue", "red"),
       pch    = c(16, 15),
       lty    = 1, lwd = 2)