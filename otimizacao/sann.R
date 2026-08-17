# SIMULATED ANNEALING — Otimização (O1 e O2)

source("hill.R")       # precisa de hchange para o rchange
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



# 3. Funções de avaliação para SANN 
EV <- 0
BEST <- -Inf
F_curva <- numeric()

m_neg_eval_O1 <- function(s) {
  v <- eval_O1(s)
  
  EV <<- EV + 1 
  if (v > BEST) { BEST <<- v } 
  if (EV <= N) { F_curva[EV] <<- BEST } 
  
  return(-v)
}

m_neg_eval_O2 <- function(s) {
  v <- eval_O2(s)
  
  EV <<- EV + 1
  if (!is.infinite(v) && v > BEST) { BEST <<- v } 
  if (EV <= N) { F_curva[EV] <<- BEST }

  if (is.infinite(v) && v == -Inf) return(Inf)
  return(-v)
}

# 4. Função de perturbação
rchange_sann <- function(par) {
  hchange(par, lower = lower, upper = upper,
          operator = "*", dist = rnorm, mean = 1, sd = 0.05,
          round = FALSE)
}



# 5. SImulated Annealing - O1 (Teste de temperaturas)
Runs  <- 5
N     <- 10000

temperaturas_teste <- c(10, 250, 500, 750, 1000)
cores_grafico      <- c("red", "blue", "darkgreen", "orange", "purple")

cat("Simulated Annealing O1 - Estudo de Temperaturas \n")
cat("Corridas:", Runs, "| Iterações por corrida:", N, "\n\n")

resultados_O1 <- list()
curvas_O1     <- list()

best_global_O1      <- -Inf
best_sol_global_O1  <- NULL
best_temp_global_O1 <- NA 

for (t_idx in 1:length(temperaturas_teste)) {
  temp_atual <- temperaturas_teste[t_idx]
  cat(sprintf("A testar Temperatura: %d\n", temp_atual))
  
  best_temp_O1    <- -Inf
  best_curva_temp <- NULL
  
  for (i in 1:Runs) {
    set.seed(i * 10)
    
    EV <<- 0
    BEST <<- -Inf
    F_curva <<- rep(NA, N) 
    
    s_init <- runif(length(lower), lower, upper)
    
    SA <- optim(
      par     = s_init,
      fn      = m_neg_eval_O1, 
      gr      = rchange_sann,
      method  = "SANN",
      control = list(maxit = N, temp = temp_atual, trace = FALSE)
    )
    
    lucro <- -SA$value
    cat(sprintf("  Corrida %d: lucro = %d USD\n", i, round(lucro)))
    
    if (lucro > best_temp_O1) {
      best_temp_O1    <- lucro
      best_curva_temp <- F_curva 
    }
    
    # Guarda o melhor absoluto e qual foi a temperatura que o conseguiu!
    if (lucro > best_global_O1) {
      best_global_O1      <- lucro
      best_sol_global_O1  <- SA$par
      best_temp_global_O1 <- temp_atual
    }
  }
  
  resultados_O1[[t_idx]] <- best_temp_O1
  curvas_O1[[t_idx]]     <- best_curva_temp
  cat(sprintf("Melhor lucro com Temp %d: %d USD\n\n", temp_atual, round(best_temp_O1)))
}

cat(sprintf("Melhor lucro O1 Global (SA): %d USD (alcançado com Temperatura %d)\n\n", 
            round(best_global_O1), best_temp_global_O1))

cat("PLANO SEMANAL O1\n\n")
print_plan(as.vector(best_sol_global_O1))


# 6. Simulated Annealing - O2 (Teste de temperaturas)
cat("Simulated Annealing O2 - Lucro com restrição <=10000 unidades\n")
cat("A iniciar o estudo de temperaturas também para o O2...\n\n")

# 6.1 Função auxiliar: Calcular total de unidades de um plano
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

# 6.2 Estratégia: Gerar Ponto Inicial Válido
cat("A gerar ponto inicial inteligente (reparação do melhor O1 global)...\n")

s_smart <- best_sol_global_O1
unidades_atuais <- get_units_total(s_smart)

iteracoes_corte <- 0
while(unidades_atuais > 10000) {
  for (i in 1:4) {
    idx <- ((i - 1) * 21 + 1):(i * 21)
    idx_J  <- seq(1, 21, by = 3)
    idx_X  <- seq(2, 21, by = 3)
    idx_PR <- seq(3, 21, by = 3)
    
    s_smart[idx[idx_J]] <- floor(s_smart[idx[idx_J]] * 0.95)
    s_smart[idx[idx_X]] <- floor(s_smart[idx[idx_X]] * 0.95)
    s_smart[idx[idx_PR]] <- s_smart[idx[idx_PR]] * 0.95
  }
  unidades_atuais <- get_units_total(s_smart)
  iteracoes_corte <- iteracoes_corte + 1
}

cat("Reparação concluída em", iteracoes_corte, "passos.\n")
cat("Unidades do novo ponto inicial:", unidades_atuais, "(VÁLIDO)\n")
cat("Lucro inicial no O2:", round(eval_O2(s_smart)), "USD\n\n")


# 6.3 Execução do Simulated Annealing (O2)
resultados_O2 <- list()
curvas_O2     <- list()

best_global_O2      <- -Inf
best_sol_global_O2  <- NULL
best_temp_global_O2 <- NA

for (t_idx in 1:length(temperaturas_teste)) {
  temp_atual <- temperaturas_teste[t_idx]
  cat(sprintf("A testar Temperatura (O2): %d\n", temp_atual))
  
  best_temp_O2    <- -Inf
  best_curva_temp <- NULL
  
  for (i in 1:Runs) {
    set.seed(i * 10)
    
    EV <<- 0
    BEST <<- -Inf
    F_curva <<- rep(NA, N)
    
    SA <- optim(
      par     = s_smart,  
      fn      = m_neg_eval_O2,
      gr      = rchange_sann,
      method  = "SANN",
      control = list(maxit = N, temp = temp_atual, trace = FALSE)
    )
    
    lucro <- -SA$value
    if (is.infinite(lucro)) lucro <- -Inf
    cat(sprintf("  Corrida %d: lucro = %s USD\n", i,
                ifelse(is.infinite(lucro), "inválido", as.character(round(lucro)))))
    
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
  
  resultados_O2[[t_idx]] <- best_temp_O2
  curvas_O2[[t_idx]]     <- best_curva_temp
  cat(sprintf("Melhor lucro (O2) com Temp %d: %s USD\n\n", 
              temp_atual, ifelse(is.infinite(best_temp_O2), "N/A", as.character(round(best_temp_O2)))))
}

if (!is.null(best_sol_global_O2)) {
  cat(sprintf("Melhor lucro O2 Global (SA): %d USD (alcançado com Temperatura %d)\n", 
              round(best_global_O2), best_temp_global_O2))
  cat("Total unidades vendidas do campeão:", round(get_units_total(best_sol_global_O2)), "\n\n")
  
  cat("PLANO SEMANAL O2\n\n")
  print_plan(as.vector(best_sol_global_O2))
} else {
  cat("Não foi encontrada nenhuma solução válida para o cenário O2 em nenhuma temperatura.\n")
}

# 7. Resumo
cat("RESUMO SIMULATED ANNEALING\n")
cat(sprintf("O1 - Lucro máximo:            %8d USD (Temp %d)\n", 
            round(best_global_O1), best_temp_global_O1))
cat(sprintf("O2 - Lucro c/ restrição:      %8s USD (Temp %s)\n",
            ifelse(is.infinite(best_global_O2), "N/A", as.character(round(best_global_O2))), 
            ifelse(is.infinite(best_global_O2), "N/A", as.character(best_temp_global_O2))))


# 8. Gráficos de aprendizagem (Convergência dupla)
cat("A gerar gráficos de convergência...\n")

par(mfrow = c(1, 2)) 

# Gráfico O1
y_min_O1 <- min(unlist(curvas_O1), na.rm = TRUE)
y_max_O1 <- max(unlist(curvas_O1), na.rm = TRUE)

plot(1, type = "n", xlim = c(1, N), ylim = c(y_min_O1, y_max_O1),
     main = "Curva de Aprendizagem O1\n(Comparação de Temperaturas)",
     xlab = "Avaliações", ylab = "Melhor Lucro (USD)")

for (t_idx in 1:length(temperaturas_teste)) {
  lines(curvas_O1[[t_idx]], col = cores_grafico[t_idx], lwd = 2)
}
legend("bottomright", legend = paste("Temp =", temperaturas_teste),
       col = cores_grafico, lwd = 2, cex = 0.8)

# Gráfico O2
valores_O2 <- unlist(curvas_O2)
valores_O2 <- valores_O2[is.finite(valores_O2)]

if(length(valores_O2) > 0) {
  y_min_O2 <- min(valores_O2, na.rm = TRUE)
  y_max_O2 <- max(valores_O2, na.rm = TRUE)
  
  plot(1, type = "n", xlim = c(1, N), ylim = c(y_min_O2, y_max_O2),
       main = "Curva de Aprendizagem O2\n(Comparação de Temperaturas)",
       xlab = "Avaliações", ylab = "Melhor Lucro (USD)")
  
  for (t_idx in 1:length(temperaturas_teste)) {
    curva_limpa <- curvas_O2[[t_idx]]
    if(!is.null(curva_limpa) && any(is.finite(curva_limpa))) {
      lines(curva_limpa, col = cores_grafico[t_idx], lwd = 2)
    }
  }
  legend("bottomright", legend = paste("Temp =", temperaturas_teste),
         col = cores_grafico, lwd = 2, cex = 0.8)
  
} else {
  plot(1, type="n", main="Sem convergência válida no O2", xlab="", ylab="")
}

par(mfrow = c(1, 1))