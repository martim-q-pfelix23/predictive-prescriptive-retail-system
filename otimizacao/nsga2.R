# NSGA-II — Otimização Multi-Objetivo (O3)

library(mco)
source("avaliacao.R")

# 1. Carregar previsões
cat("A carregar previsões...\n")

bt_pred <- read.csv("previsoes_baltimore.csv",    header = TRUE)
lc_pred <- read.csv("previsoes_lancaster.csv",    header = TRUE)
ph_pred <- read.csv("previsoes_philadelphia.csv", header = TRUE)
rv_pred <- read.csv("previsoes_richmond.csv",     header = TRUE)

C_baltimore    <- bt_pred$Predicted_Customers
C_lancaster    <- lc_pred$Predicted_Customers
C_philadelphia <- ph_pred$Predicted_Customers
C_richmond     <- rv_pred$Predicted_Customers

is_weekend <- function(weekday_vec) {
  weekday_vec %in% c("sábado", "domingo", "Saturday", "Sunday")
}

weekend_baltimore    <- is_weekend(bt_pred$Weekday)
weekend_lancaster    <- is_weekend(lc_pred$Weekday)
weekend_philadelphia <- is_weekend(ph_pred$Weekday)
weekend_richmond     <- is_weekend(rv_pred$Weekday)


# 2. Bounds
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

b_bt <- calc_bounds(C_baltimore)
b_lc <- calc_bounds(C_lancaster)
b_ph <- calc_bounds(C_philadelphia)
b_rv <- calc_bounds(C_richmond)

lower <- c(b_bt$lower, b_lc$lower, b_ph$lower, b_rv$lower)
upper <- c(b_bt$upper, b_lc$upper, b_ph$upper, b_rv$upper)



# 3. NSGA-II — O3
cat("NSGA-II O3 - Maximizar lucro e Minimizar HR\n\n")

popsize <- 200   # tamanho da população
gens    <- 600   # número de gerações

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

# 4. ANALISAR CURVA DE PARETO
# 4.1. Extrair soluções Pareto-ótimas brutas
pareto_idx    <- which(NSGA$pareto.optimal)
raw_profit    <- -NSGA$value[pareto_idx, 1] 
raw_HR        <- NSGA$value[pareto_idx, 2]
raw_sols      <- NSGA$par[pareto_idx, ]

# 4.2. Filtro: Manter apenas lucros > 0
valid_idx     <- which(raw_profit > 0)

# Se não houver lucros positivos, avisar (ajuda a detetar erros de convergência)
if(length(valid_idx) == 0) {
  stop("ERRO: Nenhuma solução com lucro positivo encontrada. Aumenta as gerações ou a população.")
}

pareto_profit <- raw_profit[valid_idx]
pareto_HR     <- raw_HR[valid_idx]
pareto_sols   <- raw_sols[valid_idx, ]

# Ordenar por lucro decrescente
ord <- order(pareto_profit, decreasing = TRUE)
cat(sprintf("%-6s %-12s %-10s\n", "Rank", "Lucro (USD)", "Total HR"))
cat(rep("-", 32), "\n", sep="")
for (k in 1:min(10, length(ord))) {
  i <- ord[k]
  cat(sprintf("%-6d %-12d %-10d\n", k, round(pareto_profit[i]), round(pareto_HR[i])))
}


# 5. Solução recomendada
best_idx    <- ord[1]
best_profit <- pareto_profit[best_idx]
best_HR     <- pareto_HR[best_idx]
best_sol_O3 <- pareto_sols[best_idx, ]

# Calcular as unidades reais vendidas para esta solução
# Como a restrição devolve (unidades - 10000), somamos 10000 de volta
best_units <- eval_O3_restricao(best_sol_O3) + 10000

best_units_reais <- min(10000, best_units)

cat("\nSolução recomendada (maior lucro Pareto)\n")
cat(sprintf("Lucro:     %d USD\n", round(best_profit)))
cat(sprintf("Total HR:  %d\n",    round(best_HR)))
cat(sprintf("Unidades:  %d / 10000\n\n", round(best_units)))

print_plan(best_sol_O3)


# 6. Gráfico da curva de Pareto
opar <- par(mar = c(5, 4, 7, 2), xpd = TRUE) 

plot(pareto_HR, pareto_profit,
     xlab = "Total HR (recursos humanos)",
     ylab = "Lucro total (USD)",
     main = "Curva de Pareto — O3 (NSGA-II)",
     pch  = 16, col = "steelblue")

points(best_HR, best_profit, pch = 8, col = "red", cex = 1.5)
grid()

legend("topleft", inset = c(0, -0.3), 
       legend = c("Soluções Pareto", "Solução recomendada"),
       pch    = c(16, 8),
       col    = c("steelblue", "red"),
       bty    = "n") # Sem caixa à volta para não cortar linhas

par(opar)


# 7. Resumo
cat("RESUMO NSGA-II (O3)\n")
cat(sprintf("Soluções Pareto encontradas: %d\n", length(pareto_profit)))
cat(sprintf("Melhor lucro Pareto:  %8d USD\n",   round(best_profit)))
cat(sprintf("HR correspondente:    %8d\n",        round(best_HR)))
