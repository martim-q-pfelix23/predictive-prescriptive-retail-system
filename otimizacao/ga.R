# GENETIC ALGORITHM — VERSÃO MELHORADA

library(GA)
source("avaliacao.R")


# 1. Carregar previsões
cat("A carregar previsões...\n")

bt_pred <- read.csv("previsoes_baltimore.csv",    header = TRUE)
lc_pred <- read.csv("previsoes_lancaster.csv",    header = TRUE)
ph_pred <- read.csv("previsoes_philadelphia.csv", header = TRUE)
rv_pred <- read.csv("previsoes_richmond.csv",     header = TRUE)

is_weekend <- function(weekday_vec) {
  weekday_vec %in% c("sábado", "domingo", "Saturday", "Sunday")
}

C_baltimore    <<- bt_pred$Predicted_Customers
C_lancaster    <<- lc_pred$Predicted_Customers
C_philadelphia <<- ph_pred$Predicted_Customers
C_richmond     <<- rv_pred$Predicted_Customers

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
    
    upper_store <- c(
      upper_store,
      ceiling(C[d] / 6),
      ceiling(C[d] / 7),
      0.30
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


# 3. Função auxiliar
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

# 4. Configuração do GA
Runs        <- 5

popSize     <- 200
maxiter     <- 600

pmutation   <- 0.15
pcrossover  <- 0.65

elitism_val <- round(popSize * 0.10)


# 5. GA — O1
cat("GA — O1 (Maximizar Lucro)\n")

best_global_O1 <- -Inf
best_solution_O1 <- NULL
curvas_O1 <- list()

for(r in 1:Runs){
  
  cat(sprintf("\n>>> RUN %d <<<\n", r))
  
  set.seed(r * 100)
  GA_O1 <- ga(
    type        = "real-valued",
    fitness     = eval_O1,
    lower       = lower,
    upper       = upper,
    popSize     = popSize,
    maxiter     = maxiter,
    pmutation   = pmutation,
    pcrossover  = pcrossover,
    elitism     = elitism_val,
    run         = 100,
    monitor     = FALSE
  )
  
  lucro <- GA_O1@fitnessValue
  
  cat(sprintf("Lucro obtido: %d USD\n", round(lucro)))
  
  curvas_O1[[r]] <- GA_O1@summary[,1]
  
  if(lucro > best_global_O1){
    
    best_global_O1 <- lucro
    best_solution_O1 <- GA_O1@solution[1,]
  }
}

cat(sprintf("MELHOR O1 GLOBAL: %d USD\n", round(best_global_O1)))

print_plan(best_solution_O1)

# 6. Gerar solução inteligente para O2
cat("\nA gerar solução inicial inteligente para O2...\n")

s_smart <- best_solution_O1
units_now <- get_units_total(s_smart)
repair_steps <- 0

while(units_now > 10000){
  
  for(i in 1:4){
    
    idx <- ((i-1)*21 + 1):(i*21)
    
    idx_J  <- seq(1, 21, by = 3)
    idx_X  <- seq(2, 21, by = 3)
    idx_PR <- seq(3, 21, by = 3)
    
    s_smart[idx[idx_J]]  <- floor(s_smart[idx[idx_J]] * 0.95)
    s_smart[idx[idx_X]]  <- floor(s_smart[idx[idx_X]] * 0.95)
    s_smart[idx[idx_PR]] <- s_smart[idx[idx_PR]] * 0.95
  }
  
  units_now <- get_units_total(s_smart)
  repair_steps <- repair_steps + 1
}

cat(sprintf("Reparação concluída em %d passos\n", repair_steps))
cat(sprintf("Unidades finais: %d\n", round(units_now)))
cat(sprintf("Lucro inicial O2: %d USD\n\n", round(eval_O2(s_smart))))

# 7. GA — O2
cat("GA — O2 (Restrição <=10000 unidades)\n")

best_global_O2 <- -Inf
best_solution_O2 <- NULL
curvas_O2 <- list()

for(r in 1:Runs){
  
  cat(sprintf("\n>>> RUN %d <<<\n", r))
  
  # diversidade inicial
  suggestions <- rbind(s_smart, runif(length(lower), lower, upper), best_solution_O1)
  
  set.seed(r * 200)
  GA_O2 <- ga(
    type        = "real-valued",
    fitness     = eval_O2_ga,
    lower       = lower,
    upper       = upper,
    popSize     = popSize,
    maxiter     = maxiter,
    pmutation   = pmutation,
    pcrossover  = pcrossover,
    elitism     = elitism_val,
    suggestions = suggestions,
    run         = 100,
    monitor     = FALSE
  )
  
  lucro <- GA_O2@fitnessValue
  
  cat(sprintf("Lucro obtido: %d USD\n", round(lucro)))
  
  curvas_O2[[r]] <- GA_O2@summary[,1]
  
  if(lucro > best_global_O2){
    
    best_global_O2 <- lucro
    best_solution_O2 <- GA_O2@solution[1,]
  }
}

cat(sprintf("MELHOR O2 GLOBAL: %d USD\n", round(best_global_O2)))
cat(sprintf("UNIDADES TOTAIS: %d\n",round(get_units_total(best_solution_O2))))

print_plan(best_solution_O2)

# 8. Gráficos de convergência
cat("\nA gerar gráficos...\n")

par(mfrow = c(1,2))


# O1
valores_O1 <- unlist(curvas_O1)
max_gens <- max(sapply(curvas_O1, length))

plot(
  1,
  type = "n",
  #xlim = c(1, maxiter),
  xlim = c(1, max_gens),
  ylim = range(valores_O1, na.rm = TRUE),
  main = "GA O1 — Convergência",
  xlab = "Geração",
  ylab = "Melhor Fitness"
)

cores <- c("red", "blue", "darkgreen", "orange", "purple")

for(i in 1:Runs){
  
  lines(curvas_O1[[i]],
        col = cores[i],
        lwd = 2)
}

legend(
  "bottomright",
  legend = paste("Run", 1:Runs),
  col = cores,
  lwd = 2,
  cex = 0.8
)

grid()


# O2
valores_O2 <- unlist(curvas_O2)
max_gens_o2 <- max(sapply(curvas_O2, length))

plot(
  1,
  type = "n",
  #xlim = c(1, maxiter),
  xlim = c(1, max_gens_o2),
  ylim = range(valores_O2, na.rm = TRUE),
  main = "GA O2 — Convergência",
  xlab = "Geração",
  ylab = "Melhor Fitness"
)

for(i in 1:Runs){
  
  lines(curvas_O2[[i]],
        col = cores[i],
        lwd = 2)
}

legend(
  "bottomright",
  legend = paste("Run", 1:Runs),
  col = cores,
  lwd = 2,
  cex = 0.8
)

grid()

par(mfrow = c(1,1))


# 9. Resumo final
cat("RESUMO FINAL — GENETIC ALGORITHM\n")

cat(sprintf("O1 -> Melhor lucro: %d USD\n",
            round(best_global_O1)))

cat(sprintf("O2 -> Melhor lucro: %d USD\n",
            round(best_global_O2)))

cat(sprintf("O2 -> Total unidades: %d\n",
            round(get_units_total(best_solution_O2))))