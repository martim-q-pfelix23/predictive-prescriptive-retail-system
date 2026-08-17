# SIMULATED ANNEALING O3 — Soma Ponderada
# Fase 1: última semana

source("hill.R")
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


# 3. Funções auxiliares
get_metrics_O3 <- function(s) {
  stores <- make_stores()
  total_profit <- 0
  total_hr     <- 0
  for (i in 1:4) {
    idx <- ((i-1)*21 + 1):(i*21)
    res <- eval(s[idx], stores[[i]])
    total_profit <- total_profit + res$profit
    J_vals <- round(s[idx][seq(1, 21, by = 3)])
    X_vals <- round(s[idx][seq(2, 21, by = 3)])
    total_hr <- total_hr + sum(J_vals) + sum(X_vals)
  }
  list(profit = total_profit, hr = total_hr)
}



rchange_sann <- function(par, ...) {
  lim_inf <- get("lower", envir = .GlobalEnv)
  lim_sup <- get("upper", envir = .GlobalEnv)
  hchange(par, lower = lim_inf, upper = lim_sup,
          operator = "*", dist = rnorm, mean = 1, sd = 0.05,
          round = FALSE)
}


cat("\nSANN O3 - Fase 1 (última semana)\n\n")

pesos_lucro_testar <- seq(0.1, 0.9, by = 0.1)
N_O3   <- 10000
temp   <- 0.05
s_mid  <- lower + (upper - lower) / 2

pareto_sann <- data.frame(
  Peso_Lucro = numeric(),
  Peso_HR    = numeric(),
  Lucro      = numeric(),
  HR         = numeric()
)

for (w_lucro in pesos_lucro_testar) {
  w_hr <- 1 - w_lucro
  cat(sprintf("Peso Lucro: %2.0f%% | Peso HR: %2.0f%% ...\n",
              w_lucro * 100, w_hr * 100))
  
  set.seed(42)
  SA_O3 <- optim(
    par     = s_mid,
    fn      = neg_eval_O3_pesos,
    gr      = rchange_sann,
    method  = "SANN",
    control = list(maxit = N_O3, temp = temp, trace = FALSE),
    w_lucro = w_lucro,
    w_hr    = w_hr
  )
  
  metrics <- get_metrics_O3(SA_O3$par)
  pareto_sann <- rbind(pareto_sann, data.frame(
    Peso_Lucro = w_lucro,
    Peso_HR    = w_hr,
    Lucro      = metrics$profit,
    HR         = metrics$hr
  ))
}

# 4. Resultados e gráfico
cat("RESULTADOS DA FRONTEIRA DE PARETO (Fase 1)\n")
print(pareto_sann)

# Filtrar soluções inválidas (lucro negativo)
pareto_sann <- pareto_sann[pareto_sann$Lucro >= 0, ]

# Filtrar soluções dominadas 
is_dominated <- function(df) {
  n   <- nrow(df)
  dom <- logical(n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      if (df$Lucro[j] >= df$Lucro[i] &&
          df$HR[j]    <= df$HR[i]    &&
          (df$Lucro[j] > df$Lucro[i] || df$HR[j] < df$HR[i])) {
        dom[i] <- TRUE
        break
      }
    }
  }
  return(dom)
}

pareto_sann <- pareto_sann[!is_dominated(pareto_sann), ]

cat("\nSoluções na fronteira de Pareto (após filtros):\n")
print(pareto_sann)

# Gráfico
pareto_sann <- pareto_sann[order(pareto_sann$HR), ]

plot(pareto_sann$HR, pareto_sann$Lucro,
     type = "b", pch = 16, col = "darkorange", lwd = 2,
     xlab = "Total HR (Número de Colaboradores)",
     ylab = "Lucro Total (USD)",
     main = "Fronteira de Pareto — O3 (SANN) | Fase 1")
grid()