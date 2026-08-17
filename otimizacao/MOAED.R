# MOEA/D — Otimização Multi-Objetivo (O3)

#install.packages("MOEADr")
library(MOEADr)
source("hill.R")
source("avaliacao.R")


# 1. Carregar previsões

cat("A carregar previsões...\n")

bt_pred <- read.csv("previsoes_baltimore.csv",    header = TRUE)
lc_pred <- read.csv("previsoes_lancaster.csv",    header = TRUE)
ph_pred <- read.csv("previsoes_philadelphia.csv", header = TRUE)
rv_pred <- read.csv("previsoes_richmond.csv",     header = TRUE)

C_baltimore    <<- bt_pred$Predicted_Customers
C_lancaster    <<- lc_pred$Predicted_Customers
C_philadelphia <<- ph_pred$Predicted_Customers
C_richmond     <<- rv_pred$Predicted_Customers

is_weekend <- function(weekday_vec) {
  weekday_vec %in% c("sábado", "domingo", "Saturday", "Sunday")
}

weekend_baltimore    <<- is_weekend(bt_pred$Weekday)
weekend_lancaster    <<- is_weekend(lc_pred$Weekday)
weekend_philadelphia <<- is_weekend(ph_pred$Weekday)
weekend_richmond     <<- is_weekend(rv_pred$Weekday)


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

lower <<- c(b_bt$lower, b_lc$lower, b_ph$lower, b_rv$lower)
upper <<- c(b_bt$upper, b_lc$upper, b_ph$upper, b_rv$upper)

# 3. Função de avaliação para o MOEA/D

moead_eval <- function(X, ...) {
  # Preparar a matriz de saída
  Y <- matrix(0, nrow = nrow(X), ncol = 2)
  
  for (i in 1:nrow(X)) {
    s <- X[i, ]
    stores <- make_stores()
    total_profit <- 0
    total_units  <- 0
    total_hr     <- 0
    
    for (j in 1:4) {
      idx <- ((j-1)*21 + 1):(j*21)
      res <- eval(s[idx], stores[[j]])
      total_profit <- total_profit + res$profit
      total_units  <- total_units  + res$units
      
      # Contar HR
      J_vals <- round(s[idx][seq(1, 21, by = 3)])
      X_vals <- round(s[idx][seq(2, 21, by = 3)])
      total_hr <- total_hr + sum(J_vals) + sum(X_vals)
    }
    
    # Aplicar a restrição das 10.000 unidades do O2 via penalização
    if (total_units > 10000) {
      multa <- (total_units - 10000) * 100
      total_profit <- total_profit - multa
    }
    
    # Guardar na matriz Y 
    Y[i, 1] <- -total_profit  
    Y[i, 2] <- total_hr   
  }
  
  return(Y)
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


# 4. Execução do MOEA/D

cat("\nA executar MOEA/D O3 - Maximizar Lucro e Minimizar HR\n")
cat("A gerar subproblemas e avaliar populações...\n")

# Configuração do problema para o pacote MOEADr
problem <- list(
  name = "moead_eval",
  xmin = lower,
  xmax = upper,
  m    = 2 # Número de objetivos
)

# Configuração da decomposição
decomp <- list(name = "sld", H = 99)

# Condição de paragem
stopcrt <- list(list(name = "maxiter", maxiter = 200))

set.seed(42)

# Executar o algoritmo
MOEAD_out <- moead(
  problem  = problem,
  decomp   = decomp,
  stopcrit = stopcrt,
  preset   = preset_moead("original"),
  showpars = list(show.iters = "dots", showevery = 10)
)

# 5. Analisar e extrair curva de pareto
lucros_finais <- -MOEAD_out$Y[, 1]
hr_finais     <- MOEAD_out$Y[, 2]
sols_finais   <- MOEAD_out$X

df_moead <- data.frame(
  ID    = 1:length(lucros_finais),
  Lucro = lucros_finais,
  HR    = hr_finais
)

df_moead <- df_moead[df_moead$Lucro > 0, ]
df_moead <- df_moead[order(df_moead$HR, -df_moead$Lucro), ]

pareto_front <- data.frame()
max_lucro_visto <- -Inf

if (nrow(df_moead) > 0) {
  for (i in 1:nrow(df_moead)) {
    if (df_moead$Lucro[i] > max_lucro_visto) {
      pareto_front <- rbind(pareto_front, df_moead[i, ])
      max_lucro_visto <- df_moead$Lucro[i]
    }
  }
}

cat("\n\nResultados do MOEA/D (Fronteira de Pareto - Lucro Positivo):\n")
cat(sprintf("%-6s %-12s %-10s %-12s\n", "Rank", "Lucro (USD)", "Total HR", "Unidades"))
cat(rep("-", 44), "\n", sep="")

if (nrow(pareto_front) > 0) {
  pareto_print <- pareto_front[order(-pareto_front$Lucro), ]
  for (k in 1:min(10, nrow(pareto_print))) {
    sol_idx <- pareto_print$ID[k]
    sol_k   <- sols_finais[sol_idx, ]
    unidades_k <- get_units_total(sol_k)
    
    cat(sprintf("%-6d %-12d %-10d %-12d\n", k, round(pareto_print$Lucro[k]), round(pareto_print$HR[k]), round(unidades_k)))
  }
  
  best_idx    <- pareto_print$ID[1]
  best_profit <- pareto_print$Lucro[1]
  best_hr     <- pareto_print$HR[1]
  best_sol_O3 <- sols_finais[best_idx, ]
  best_unidades <- get_units_total(best_sol_O3)
  
  cat("\nSolução Recomendada (Maior Lucro no MOEA/D)\n")
  cat(sprintf("Lucro:    %d USD\n", round(best_profit)))
  cat(sprintf("Total HR: %d\n",   round(best_hr)))
  cat(sprintf("Unidades: %d\n\n", round(best_unidades)))
  
} else {
  cat("Não foram encontradas soluções com lucro positivo e não dominadas.\n")
}

# 6. Gráfico de fronteira de Pareto
cat("A gerar gráfico da Fronteira de Pareto limpa...\n")

if (nrow(pareto_front) > 0) {
  # Garantir que está ordenado por HR para a linha não cruzar
  pareto_front <- pareto_front[order(pareto_front$HR), ]
  
  opar <- par(mar = c(5, 4, 7, 2), xpd = TRUE) 
  
  plot(pareto_front$HR, pareto_front$Lucro,
       type = "b", pch = 16, col = "purple", lwd = 2,
       xlab = "Total HR (Recursos Humanos)",
       ylab = "Lucro Total (USD)",
       main = "Curva de Pareto — O3 (MOEA/D)")
  
  # Destacar a solução recomendada
  points(best_hr, best_profit, pch = 8, col = "red", cex = 1.5)
  grid()
  
  legend("topleft", inset = c(0, -0.3), 
         legend = c("Soluções Não Dominadas", "Solução Recomendada"),
         pch    = c(16, 8),
         col    = c("purple", "red"),
         bty    = "n")
  
  par(opar)
} else {
  cat("Não há dados válidos para desenhar o gráfico da Fronteira de Pareto.\n")
}