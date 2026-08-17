
# MONTE CARLO - Otimização (O1)

source("blind.R")
source("montecarlo.R")
source("avaliacao.R")


# 1. Leitura dos CSVs
bt <- read.table("baltimore_clean.csv",    header = TRUE, sep = ",")
lc <- read.table("lancaster_clean.csv",    header = TRUE, sep = ",")
ph <- read.table("philadelphia_clean.csv", header = TRUE, sep = ",")
rv <- read.table("richmond_clean.csv",     header = TRUE, sep = ",")

# Extrair os últimos 7 dias (última semana real do dataset)
last7 <- function(df) tail(df, 7)

bt_week <- last7(bt)
lc_week <- last7(lc)
ph_week <- last7(ph)
rv_week <- last7(rv)

# Mostrar as últimas semanas lidas para confirmar
cat("Últimos 7 dias lidos\n")
cat("\nBaltimore:\n")
print(bt_week[, c("Date", "DayOfWeek", "Weekend", "Num_Customers")])
cat("\nLancaster:\n")
print(lc_week[, c("Date", "DayOfWeek", "Weekend", "Num_Customers")])
cat("\nPhiladelphia:\n")
print(ph_week[, c("Date", "DayOfWeek", "Weekend", "Num_Customers")])
cat("\nRichmond:\n")
print(rv_week[, c("Date", "DayOfWeek", "Weekend", "Num_Customers")])

# Vetores de clientes para a função eval
C_baltimore    <- bt_week$Num_Customers
C_lancaster    <- lc_week$Num_Customers
C_philadelphia <- ph_week$Num_Customers
C_richmond     <- rv_week$Num_Customers

# Identificar quais os dias são fim de semana (para custo de HR)
weekend_baltimore    <- bt_week$Weekend
weekend_lancaster    <- lc_week$Weekend
weekend_philadelphia <- ph_week$Weekend
weekend_richmond     <- rv_week$Weekend

# 2. Bounds
calc_bounds <- function(C) {
  lower_store <- c()
  upper_store <- c()
  for (d in 1:7) {
    lower_store <- c(lower_store, 0, 0, 0.00)
    upper_store <- c(upper_store,
                     ceiling(C[d] / 6),   # J_max
                     ceiling(C[d] / 7),   # X_max
                     0.30)                # PR_max
  }
  return(list(lower = lower_store, upper = upper_store))
}

b_bt <- calc_bounds(C_baltimore)
b_lc <- calc_bounds(C_lancaster)
b_ph <- calc_bounds(C_philadelphia)
b_rv <- calc_bounds(C_richmond)

lower <- c(b_bt$lower, b_lc$lower, b_ph$lower, b_rv$lower)
upper <- c(b_bt$upper, b_lc$upper, b_ph$upper, b_rv$upper)



# 3. MONTE CARLO
N <- 10000  # número de soluções aleatórias testadas

cat("\nMonte Carlo O1 - Maximizar lucro total\n")
cat("Número de amostras:", N, "\n\n")

set.seed(123)
MC <- mcsearch(fn    = eval_O1,
               lower = lower,
               upper = upper,
               N     = N,
               type  = "max")

cat("Melhor lucro encontrado:", round(MC$eval), "USD\n")
cat("Encontrado na iteração: ", MC$index, "\n\n")


# 5. Mostrar plano final
s_best <- MC$sol
lojas  <- c("Baltimore", "Lancaster", "Philadelphia", "Richmond")

for (i in 1:4) {
  idx    <- ((i - 1) * 21 + 1):(i * 21)
  s_loja <- s_best[idx]

  J_vals  <- round(s_loja[seq(1, 21, by = 3)])
  X_vals  <- round(s_loja[seq(2, 21, by = 3)])
  PR_vals <- round(s_loja[seq(3, 21, by = 3)], 2)

  cat("---", lojas[i], "---\n")
  cat(sprintf("%-4s %5s %5s %5s %5s %5s %5s %5s\n",
              "", "d1", "d2", "d3", "d4", "d5", "d6", "d7"))
  cat(sprintf("%-4s %5d %5d %5d %5d %5d %5d %5d\n",
              "J",  J_vals[1],  J_vals[2],  J_vals[3],  J_vals[4],
              J_vals[5],  J_vals[6],  J_vals[7]))
  cat(sprintf("%-4s %5d %5d %5d %5d %5d %5d %5d\n",
              "X",  X_vals[1],  X_vals[2],  X_vals[3],  X_vals[4],
              X_vals[5],  X_vals[6],  X_vals[7]))
  cat(sprintf("%-4s %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n\n",
              "PR", PR_vals[1], PR_vals[2], PR_vals[3], PR_vals[4],
              PR_vals[5], PR_vals[6], PR_vals[7]))
}
