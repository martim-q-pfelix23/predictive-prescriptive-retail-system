# MONTE CARLO - Otimização (O1 + O2)
# Usando previsões do modelo RF Cenário C

source("blind.R")
source("montecarlo.R")
source("avaliacao.R")


# 1. Leitura das previsões
cat("A carregar previsões do modelo RF Cenário C...\n")

bt_pred <- read.csv("previsoes_baltimore.csv",    header = TRUE)
lc_pred <- read.csv("previsoes_lancaster.csv",    header = TRUE)
ph_pred <- read.csv("previsoes_philadelphia.csv", header = TRUE)
rv_pred <- read.csv("previsoes_richmond.csv",     header = TRUE)


# 2. Mostrar previsões

cat("Previsões carregadas\n")

cat("\nBaltimore:\n")
print(bt_pred[, c("Date", "Weekday", "Predicted_Customers")])

cat("\nLancaster:\n")
print(lc_pred[, c("Date", "Weekday", "Predicted_Customers")])

cat("\nPhiladelphia:\n")
print(ph_pred[, c("Date", "Weekday", "Predicted_Customers")])

cat("\nRichmond:\n")
print(rv_pred[, c("Date", "Weekday", "Predicted_Customers")])


# 3. Clientes previstos
C_baltimore    <- bt_pred$Predicted_Customers
C_lancaster    <- lc_pred$Predicted_Customers
C_philadelphia <- ph_pred$Predicted_Customers
C_richmond     <- rv_pred$Predicted_Customers


# 4. Identificar fins de semana
is_weekend <- function(weekday_vec) {
  weekday_vec %in% c("sábado", "domingo", "Saturday", "Sunday")
}

weekend_baltimore    <- is_weekend(bt_pred$Weekday)
weekend_lancaster    <- is_weekend(lc_pred$Weekday)
weekend_philadelphia <- is_weekend(ph_pred$Weekday)
weekend_richmond     <- is_weekend(rv_pred$Weekday)

cat("\nFins de semana identificados\n")

cat("Baltimore:   ", weekend_baltimore,    "\n")
cat("Lancaster:   ", weekend_lancaster,    "\n")
cat("Philadelphia:", weekend_philadelphia, "\n")
cat("Richmond:    ", weekend_richmond,     "\n")


# 5. Bounds
calc_bounds <- function(C) {
  
  lower_store <- c()
  upper_store <- c()
  
  for (d in 1:7) {
    
    lower_store <- c(lower_store, 0, 0, 0.00)
    
    upper_store <- c(
      upper_store,
      ceiling(C[d] / 6),  # J_max
      ceiling(C[d] / 7),  # X_max
      0.30                # PR_max
    )
  }
  
  list(lower = lower_store, upper = upper_store)
}

b_bt <- calc_bounds(C_baltimore)
b_lc <- calc_bounds(C_lancaster)
b_ph <- calc_bounds(C_philadelphia)
b_rv <- calc_bounds(C_richmond)

lower <- c(b_bt$lower, b_lc$lower, b_ph$lower, b_rv$lower)
upper <- c(b_bt$upper, b_lc$upper, b_ph$upper, b_rv$upper)


# 6. Função auxiliar - total de unidades
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


# 7. MONTE CARLO — O1
N <- 10000

cat("\nMonte Carlo O1 - Maximizar lucro total\n")
cat("Número de amostras:", N, "\n\n")

set.seed(123)

MC <- mcsearch(
  fn    = eval_O1,
  lower = lower,
  upper = upper,
  N     = N,
  type  = "max"
)

cat("Melhor lucro encontrado:", round(MC$eval), "USD\n")
cat("Encontrado na iteração:", MC$index, "\n\n")


# 8. O2 — Repair da solução do O1
cat("Monte Carlo O2 - Repair da solução O1\n")

# melhor solução do O1
s_repair <- MC$sol

# calcular unidades atuais
unidades_atuais <- get_units_total(s_repair)

iteracoes_corte <- 0

# repair até cumprir <= 10000 unidades
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

# avaliar solução reparada no O2
lucro_O2 <- eval_O2(s_repair)

cat("Repair concluído em", iteracoes_corte, "iterações\n")
cat("Total de unidades:", round(unidades_atuais), "\n")
cat("Melhor lucro O2:", round(lucro_O2), "USD\n\n")


# 9. Mostrar plano final do O1
cat("PLANO FINAL — O1\n")

s_best <- MC$sol

lojas <- c(
  "Baltimore",
  "Lancaster",
  "Philadelphia",
  "Richmond"
)

for (i in 1:4) {
  
  idx <- ((i - 1) * 21 + 1):(i * 21)
  
  s_loja <- s_best[idx]
  
  J_vals  <- round(s_loja[seq(1, 21, by = 3)])
  X_vals  <- round(s_loja[seq(2, 21, by = 3)])
  PR_vals <- round(s_loja[seq(3, 21, by = 3)], 2)
  
  cat("---", lojas[i], "---\n")
  
  cat(sprintf(
    "%-4s %5s %5s %5s %5s %5s %5s %5s\n",
    "", "d1", "d2", "d3", "d4", "d5", "d6", "d7"
  ))
  
  cat(sprintf(
    "%-4s %5d %5d %5d %5d %5d %5d %5d\n",
    "J",
    J_vals[1], J_vals[2], J_vals[3], J_vals[4],
    J_vals[5], J_vals[6], J_vals[7]
  ))
  
  cat(sprintf(
    "%-4s %5d %5d %5d %5d %5d %5d %5d\n",
    "X",
    X_vals[1], X_vals[2], X_vals[3], X_vals[4],
    X_vals[5], X_vals[6], X_vals[7]
  ))
  
  cat(sprintf(
    "%-4s %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n\n",
    "PR",
    PR_vals[1], PR_vals[2], PR_vals[3], PR_vals[4],
    PR_vals[5], PR_vals[6], PR_vals[7]
  ))
}


# 10. Mostrar plano final do O2
cat("PLANO FINAL — O2 (REPARADO)\n")

for (i in 1:4) {
  
  idx <- ((i - 1) * 21 + 1):(i * 21)
  
  s_loja <- s_repair[idx]
  
  J_vals  <- round(s_loja[seq(1, 21, by = 3)])
  X_vals  <- round(s_loja[seq(2, 21, by = 3)])
  PR_vals <- round(s_loja[seq(3, 21, by = 3)], 2)
  
  cat("---", lojas[i], "---\n")
  
  cat(sprintf(
    "%-4s %5s %5s %5s %5s %5s %5s %5s\n",
    "", "d1", "d2", "d3", "d4", "d5", "d6", "d7"
  ))
  
  cat(sprintf(
    "%-4s %5d %5d %5d %5d %5d %5d %5d\n",
    "J",
    J_vals[1], J_vals[2], J_vals[3], J_vals[4],
    J_vals[5], J_vals[6], J_vals[7]
  ))
  
  cat(sprintf(
    "%-4s %5d %5d %5d %5d %5d %5d %5d\n",
    "X",
    X_vals[1], X_vals[2], X_vals[3], X_vals[4],
    X_vals[5], X_vals[6], X_vals[7]
  ))
  
  cat(sprintf(
    "%-4s %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n\n",
    "PR",
    PR_vals[1], PR_vals[2], PR_vals[3], PR_vals[4],
    PR_vals[5], PR_vals[6], PR_vals[7]
  ))
}