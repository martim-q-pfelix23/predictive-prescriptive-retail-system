
# Script de Avaliação (avaliacao.R)
eval = function(s, store){
  
  total_profit = 0
  total_units = 0
  
  C  = store$C
  FJ = store$FJ
  FX = store$FX
  cost_fixed = store$cost_fixed
  weekend = store$weekend
  
  for(d in 1:length(C)){
    
    J = round(s[(d-1)*3 + 1])
    X = round(s[(d-1)*3 + 2])
    PR = s[(d-1)*3 + 3]
    
    # clientes atendidos
    AX = min(7*X, C[d])
    remaining = C[d] - AX
    AJ = min(6*J, max(0, remaining))
    
    # unidades por cliente
    UXc = round(FX * 10 / log(2 - PR))
    UJc = round(FJ * 10 / log(2 - PR))
    
    UX = AX * UXc
    UJ = AJ * UJc
    
    # vendas
    SX = round(UX * (1 - PR) * 1.07)
    SJ = round(UJ * (1 - PR) * 1.07)
    
    revenue = SX + SJ
    
    # custos HR
    if(weekend[d]){
      cost = J*70 + X*95
    } else {
      cost = J*60 + X*80
    }
    
    total_profit = total_profit + (revenue - cost)
    total_units  = total_units  + UX + UJ
  }
  
  total_profit = total_profit - cost_fixed
  
  return(list(profit = total_profit, units = total_units))
}


# Função auxiliar: constrói as 4 store lists
make_stores = function() {
  list(
    baltimore = list(
      C = C_baltimore, weekend = weekend_baltimore,
      FJ = 1.00, FX = 1.15, cost_fixed = 700
    ),
    lancaster = list(
      C = C_lancaster, weekend = weekend_lancaster,
      FJ = 1.05, FX = 1.20, cost_fixed = 730
    ),
    philadelphia = list(
      C = C_philadelphia, weekend = weekend_philadelphia,
      FJ = 1.10, FX = 1.15, cost_fixed = 760
    ),
    richmond = list(
      C = C_richmond, weekend = weekend_richmond,
      FJ = 1.15, FX = 1.25, cost_fixed = 800
    )
  )
}


# O1 - Maximizar lucro total (sem restrições)
eval_O1 = function(s) {
  stores = make_stores()
  total  = 0
  for (i in 1:4) {
    idx   = ((i-1)*21 + 1):(i*21)
    res   = eval(s[idx], stores[[i]])
    total = total + res$profit
  }
  return(total)
}


# O2 - Maximizar lucro com restrição:
#       total de unidades vendidas <= 10000
#       Death penalty se restrição violada
eval_O2 = function(s) {
  stores      = make_stores()
  total_profit = 0
  total_units  = 0
  
  for (i in 1:4) {
    idx = ((i-1)*21 + 1):(i*21)
    res = eval(s[idx], stores[[i]])
    total_profit = total_profit + res$profit
    total_units  = total_units  + res$units
  }
  
  if (total_units > 10000) return(-Inf)
  return(total_profit)

}

eval_O2_ga = function(s) {
  
  stores       = make_stores()
  total_profit = 0
  total_units  = 0
  
  for (i in 1:4) {
    idx = ((i-1)*21 + 1):(i*21)
    res = eval(s[idx], stores[[i]])
    
    total_profit = total_profit + res$profit
    total_units  = total_units  + res$units
  }
  
  # Penalização suave
  if(total_units > 10000){
    
    excess <- total_units - 10000
    
    # penalização quadrática
    #penalty <- excess^2 * 0.15
    penalty <- excess * 8 + excess^2 * 0.01
    
    return(total_profit - penalty)
  }
  
  return(total_profit)
}


# O3 - Multi-objetivo para NSGA-II:
#       Minimizar (-profit) e minimizar HR total
#       Devolve vetor c(-profit, total_HR)
eval_O3 = function(s) {
  stores       = make_stores()
  total_profit = 0
  total_HR     = 0
  total_units  = 0
  
  for (i in 1:4) {
    idx = ((i-1)*21 + 1):(i*21)
    res = eval(s[idx], stores[[i]])
    total_profit = total_profit + res$profit
    total_units  = total_units  + res$units
    
    J_vals   = round(s[idx][seq(1, 21, by = 3)])
    X_vals   = round(s[idx][seq(2, 21, by = 3)])
    total_HR = total_HR + sum(J_vals) + sum(X_vals)
  }
  
  # HARD CONSTRAINT: Gradient Penalty
  if (total_units > 10000) {
    excess <- total_units - 10000
    penalty_profit <- 100000 + (excess * 1000)
    penalty_HR     <- 10000 + excess
    return(c(penalty_profit, penalty_HR))
  }
  
  return(c(-total_profit, total_HR))
}


# Avaliar a Restrição do O3
eval_O3_restricao = function(s) {
  stores      = make_stores()
  total_units = 0
  
  for (i in 1:4) {
    idx = ((i-1)*21 + 1):(i*21)
    res = eval(s[idx], stores[[i]])
    total_units = total_units + res$units
  }
  
  # g(x) <= 0 (Se vender menos de 10000, dá negativo -> Válido)
  return(total_units - 10000) # se for positivo o algortimo sabe que viola a restrição
}


# Função de Avaliação SANN O3
neg_eval_O3_pesos <- function(s, w_lucro, w_hr, max_profit = 10000, max_hr = 500) {
  stores       <- make_stores()
  total_profit <- 0
  total_units  <- 0
  total_hr     <- 0
  
  for (i in 1:4) {
    idx <- ((i-1)*21 + 1):(i*21)
    res <- eval(s[idx], stores[[i]])
    total_profit <- total_profit + res$profit
    total_units  = total_units  + res$units
    J_vals <- round(s[idx][seq(1, 21, by = 3)])
    X_vals <- round(s[idx][seq(2, 21, by = 3)])
    total_hr <- total_hr + sum(J_vals) + sum(X_vals)
  }
  
  if (total_units > 10000) {
    total_profit <- total_profit - (total_units - 10000) * 100
  }
  
  norm_profit <- total_profit / max_profit
  norm_hr     <- total_hr     / max_hr
  fitness     <- (w_lucro * norm_profit) - (w_hr * norm_hr)
  
  return(-fitness)
}


# Função auxiliar: mostra o plano final
print_plan = function(s) {
  lojas = c("Baltimore", "Lancaster", "Philadelphia", "Richmond")
  for (i in 1:4) {
    idx    = ((i-1)*21 + 1):(i*21)
    s_loja = s[idx]
    J_vals  = round(s_loja[seq(1, 21, by = 3)])
    X_vals  = round(s_loja[seq(2, 21, by = 3)])
    PR_vals = round(s_loja[seq(3, 21, by = 3)], 2)
    cat("---", lojas[i], "---\n")
    cat(sprintf("%-4s %5s %5s %5s %5s %5s %5s %5s\n",
                "", "d1","d2","d3","d4","d5","d6","d7"))
    cat(sprintf("%-4s %5d %5d %5d %5d %5d %5d %5d\n",
                "J", J_vals[1],J_vals[2],J_vals[3],J_vals[4],
                J_vals[5],J_vals[6],J_vals[7]))
    cat(sprintf("%-4s %5d %5d %5d %5d %5d %5d %5d\n",
                "X", X_vals[1],X_vals[2],X_vals[3],X_vals[4],
                X_vals[5],X_vals[6],X_vals[7]))
    cat(sprintf("%-4s %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n\n",
                "PR", PR_vals[1],PR_vals[2],PR_vals[3],PR_vals[4],
                PR_vals[5],PR_vals[6],PR_vals[7]))
  }
}



# Teste dos slides
s1 = c(
  0,4,0.00,
  10,0,0.05,
  4,8,0.10,
  0,20,0.15,
  5,0,0.20,
  5,4,0.25,
  4,3,0.30
)

store_data = list(
  C = c(97, 61, 65, 71, 65, 89, 125),
  weekend = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  FJ = 1.00,
  FX = 1.15,
  cost_fixed = 700
)

store_philadelphia = list(
  C = c(230, 144, 154, 168, 154, 211, 298),
  weekend = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  FJ = 1.10,
  FX = 1.15,
  cost_fixed = 760
)


#resultado_baltimor <- eval(s1, store_data)
#resultado_philadelphia <- eval(s1, store_philadelphia)

#print_res_bal <- resultado_baltimor$profit
#vendas_num <- resultado_baltimor$units
#print_res_phil <- resultado_philadelphia$profit

#cat("Lucro calculado para exemplo Baltimor: ", print_res_bal, "\n")
#cat("num de vendas para exemplo Baltimor: ", vendas_num, "\n")
#cat("Lucro calculado para exemplo Philadelphia: ", print_res_phil, "\n")

