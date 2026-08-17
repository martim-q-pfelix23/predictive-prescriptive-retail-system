# HILL CLIMBING — Otimização (O1 e O2)
source("blind.R")
source("montecarlo.R")
source("hill.R")
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


# 3. Ponto Inicial
cat("A gerar ponto inicial com Monte Carlo...\n")
set.seed(123)
MC_init <- mcsearch(fn = eval_O1, lower = lower, upper = upper,
                    N = 5000, type = "max")
s0 <- MC_init$sol
cat("Lucro ponto inicial (MC):", round(MC_init$eval), "USD\n\n")



# 4. Função de perturbação (multiplicativa)
rchange_mult <- function(par, lower, upper) {
  hchange(par, lower = lower, upper = upper,
          operator = "*", dist = rnorm, mean = 1, sd = 0.05,
          round = FALSE)
}


# 5. HILL CLIMBING — O1
N    <- 10000
REP  <- N / 10

cat("Hill Climbing O1 - Maximizar lucro total\n")
cat("Iterações:", N, "\n\n")

set.seed(42)
HC_O1 <- hclimbing(
  par     = s0,
  fn      = eval_O1,
  change  = rchange_mult,
  lower   = lower,
  upper   = upper,
  type    = "max",
  control = list(maxit = N, REPORT = REP, digits = 0)
)

cat("\nMelhor lucro O1 (Hill Climbing):", round(HC_O1$eval), "USD\n\n")
print_plan(HC_O1$sol)


# 6. HILL CLIMBING — O2
cat("Hill Climbing O2 - Lucro com restrição <=10000 unidades\n")
cat("Iterações:", N, "\n\n")

# Função auxiliar: Calcular total de unidades de um plano
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

# Estratégia de reparação heurística
cat("A gerar ponto inicial inteligente a partir do Hill Climbing O1...\n")

# Pegamos no plano do cenário O1 (que maximiza lucro sem limites)
s_smart <- HC_O1$sol
unidades_atuais <- get_units_total(s_smart)

iteracoes_corte <- 0
while(unidades_atuais > 10000) {
  for (i in 1:4) {
    idx <- ((i - 1) * 21 + 1):(i * 21)
    
    # Posições de J, X e PR no vetor da loja
    idx_J  <- seq(1, 21, by = 3)
    idx_X  <- seq(2, 21, by = 3)
    idx_PR <- seq(3, 21, by = 3)
    
    # Cortes de 5% (pessoal arredondado para baixo, promoção como decimal)
    s_smart[idx[idx_J]] <- floor(s_smart[idx[idx_J]] * 0.95)
    s_smart[idx[idx_X]] <- floor(s_smart[idx[idx_X]] * 0.95)
    s_smart[idx[idx_PR]] <- s_smart[idx[idx_PR]] * 0.95
  }
  unidades_atuais <- get_units_total(s_smart)
  iteracoes_corte <- iteracoes_corte + 1
}

# O nosso ponto inicial é agora o plano reparado e validado
s0_O2 <- s_smart

cat("Reparação concluída em", iteracoes_corte, "passos.\n")
cat("Unidades do novo ponto inicial:", round(unidades_atuais), "(VÁLIDO)\n")
cat("Valor ponto inicial O2:", round(eval_O2(s0_O2)), "USD\n\n")


set.seed(42)
HC_O2 <- hclimbing(
  par     = s0_O2,
  fn      = eval_O2,
  change  = rchange_mult,
  lower   = lower,
  upper   = upper,
  type    = "max",
  control = list(maxit = N, REPORT = REP, digits = 0)
)

cat("\nMelhor lucro O2 (Hill Climbing):", round(HC_O2$eval), "USD\n")
cat("Total de unidades vendidas (O2):", round(get_units_total(HC_O2$sol)), "\n\n")
print_plan(HC_O2$sol)


# 7. Comparação final
cat("RESUMO HILL CLIMBING\n")
cat(sprintf("O1 - Lucro máximo:            %8d USD\n", round(HC_O1$eval)))
cat(sprintf("O2 - Lucro c/ restrição:      %8d USD\n", round(HC_O2$eval)))
