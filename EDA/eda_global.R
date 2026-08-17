library(tidyverse)
library(corrplot)
library(forecast)

# Carregar CSVs
baltimore <- read.csv("baltimore.csv")
lancaster <- read.csv("lancaster.csv")
philadelphia <- read.csv("philadelphia.csv")
richmond <- read.csv("richmond.csv")

# Converter datas
baltimore$Date <- as.Date(baltimore$Date)
lancaster$Date <- as.Date(lancaster$Date)
philadelphia$Date <- as.Date(philadelphia$Date)
richmond$Date <- as.Date(richmond$Date)

# Ordenar por data
baltimore <- baltimore %>% arrange(Date)
lancaster  <- lancaster  %>% arrange(Date)
philadelphia <- philadelphia %>% arrange(Date)
richmond   <- richmond   %>% arrange(Date)

# Tratar missing values (imputação por zero)
baltimore$Pct_On_Sale[is.na(baltimore$Pct_On_Sale)]     <- 0
lancaster$Pct_On_Sale[is.na(lancaster$Pct_On_Sale)]     <- 0
philadelphia$Pct_On_Sale[is.na(philadelphia$Pct_On_Sale)] <- 0
richmond$Pct_On_Sale[is.na(richmond$Pct_On_Sale)]       <- 0

# Identificar loja
baltimore$Store    <- "Baltimore"
lancaster$Store    <- "Lancaster"
philadelphia$Store <- "Philadelphia"
richmond$Store     <- "Richmond"

# Juntar datasets das lojas todas
df_all <- bind_rows(baltimore, lancaster, philadelphia, richmond)

# Ordenar por data
df_all <- df_all %>% arrange(Date)

# Ordem correta dos dias da semana
day_levels <- c("segunda-feira", "terça-feira", "quarta-feira",
                "quinta-feira", "sexta-feira", "sábado", "domingo")

# Confirmar junção
head(df_all)


# Série temporal de clientes por loja
p1 <- ggplot(df_all, aes(x = Date, y = Num_Customers, color = Store)) +
  geom_line(alpha = 0.75, linewidth = 0.45) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  facet_wrap(~Store, ncol = 1, scales = "free_y") +
  labs(
    title    = "Série Temporal de Clientes por Loja",
    subtitle = "Variável alvo do forecasting: Num_Customers",
    x = NULL, y = "Nº de Clientes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )
print(p1)

# Média de clientes por dia da semana e loja 
df_all %>%
  mutate(DayOfWeek = factor(weekdays(Date), levels = day_levels)) %>%
  group_by(Store, DayOfWeek) %>%
  summarise(avg = mean(Num_Customers), .groups = "drop") %>%
  ggplot(aes(DayOfWeek, avg, fill = Store)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Average Customers by Day of Week and Store",
       x = "Day of Week", y = "Average Customers")


# Distribuição entre lojas
ggplot(df_all, aes(Store, Num_Customers, fill = Store)) +
  geom_boxplot() +
  labs(title = "Distribution of Customers by Store",
       x = "Store", y = "Num_Customers")


# Correlação entre lojas
df_wide <- df_all %>%
  select(Date, Store, Num_Customers) %>%
  pivot_wider(names_from = Store, values_from = Num_Customers)

cor_matrix <- cor(df_wide[,-1], use = "complete.obs")

corrplot(cor_matrix,
         method = "color",
         addCoef.col = "black",
         tl.col = "black",
         title = "Correlation Between Stores - Num_Customers")


# ACF e PACF por loja
# ACF
par(mfrow = c(2, 2))
acf(baltimore$Num_Customers,    lag.max = 35, main = "ACF - Baltimore") # lag.max = 35 para ver picos em 7, 14, 21 e 28
acf(lancaster$Num_Customers,    lag.max = 35, main = "ACF - Lancaster")
acf(philadelphia$Num_Customers, lag.max = 35, main = "ACF - Philadelphia")
acf(richmond$Num_Customers,     lag.max = 35, main = "ACF - Richmond")
par(mfrow = c(1, 1))

# PACF
par(mfrow = c(2, 2))
pacf(baltimore$Num_Customers,    lag.max = 35, main = "PACF - Baltimore")
pacf(lancaster$Num_Customers,    lag.max = 35, main = "PACF - Lancaster")
pacf(philadelphia$Num_Customers, lag.max = 35, main = "PACF - Philadelphia")
pacf(richmond$Num_Customers,     lag.max = 35, main = "PACF - Richmond")
par(mfrow = c(1, 1))



# CCF entre pares de lojas (apenas lags positivos)
plot_ccf_pos <- function(x, y, lag.max = 14, main = "") {
  cc <- ccf(x, y, lag.max = lag.max, plot = FALSE)
  # filtrar apenas lags >= 0
  idx <- cc$lag[,,1] >= 0
  lags_pos <- cc$lag[,,1][idx]
  acf_pos  <- cc$acf[,,1][idx]
  # intervalo de confiança
  ci <- qnorm(0.975) / sqrt(length(x))
  plot(lags_pos, acf_pos,
       type = "h", lwd = 2,
       xlab = "Lag (dias)", ylab = "CCF",
       main = main,
       ylim = c(-ci * 1.5, max(acf_pos) * 1.1))
  abline(h = 0)
  abline(h =  ci, lty = 2, col = "blue")
  abline(h = -ci, lty = 2, col = "blue")
}

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
plot_ccf_pos(baltimore$Num_Customers,    lancaster$Num_Customers,
             lag.max = 14, main = "CCF: Baltimore vs Lancaster")
plot_ccf_pos(baltimore$Num_Customers,    philadelphia$Num_Customers,
             lag.max = 14, main = "CCF: Baltimore vs Philadelphia")
plot_ccf_pos(baltimore$Num_Customers,    richmond$Num_Customers,
             lag.max = 14, main = "CCF: Baltimore vs Richmond")
plot_ccf_pos(lancaster$Num_Customers,    philadelphia$Num_Customers,
             lag.max = 14, main = "CCF: Lancaster vs Philadelphia")
plot_ccf_pos(lancaster$Num_Customers,    richmond$Num_Customers,
             lag.max = 14, main = "CCF: Lancaster vs Richmond")
plot_ccf_pos(philadelphia$Num_Customers, richmond$Num_Customers,
             lag.max = 14, main = "CCF: Philadelphia vs Richmond")
par(mfrow = c(1, 1))

cat("\nSe o pico do CCF estiver em lag=0 e o gráfico for simétrico\n")
cat("→ as lojas movem-se juntas mas nenhuma antecipa a outra\n")
cat("→ VAR pode ainda ser útil para capturar padrões partilhados\n\n") 


# Decomposição STL
ts_baltimore    <- ts(baltimore$Num_Customers,    frequency = 7) # frequência semanal = 7
ts_lancaster    <- ts(lancaster$Num_Customers,    frequency = 7)
ts_philadelphia <- ts(philadelphia$Num_Customers, frequency = 7)
ts_richmond     <- ts(richmond$Num_Customers,     frequency = 7)

plot(stl(ts_baltimore,    s.window = "periodic"), main = "STL Decomposition - Baltimore")
plot(stl(ts_lancaster,    s.window = "periodic"), main = "STL Decomposition - Lancaster")
plot(stl(ts_philadelphia, s.window = "periodic"), main = "STL Decomposition - Philadelphia")
plot(stl(ts_richmond,     s.window = "periodic"), main = "STL Decomposition - Richmond")



# Teste de Sazonalidade de Friedman
friedman_test_k <- function(ts_data, k, store_name) {
  
  n <- length(ts_data)
  n_complete <- floor(n / k) * k
  
  mat <- matrix(ts_data[1:n_complete], ncol = k, byrow = TRUE)
  
  result <- friedman.test(mat)
  
  cat(store_name, "- Sazonalidade", k, ":\n")
  cat("  Chi-squared =", round(result$statistic, 3),
      "| p-value =", round(result$p.value, 5))
  
  if (result$p.value < 0.05) {
    cat("  => Sazonalidade SIGNIFICATIVA\n\n")
  } else {
    cat("  => Não significativa\n\n")
  }
}

cat("TESTES DE SAZONALIDADE\n\n")

# Baltimore
friedman_test_k(baltimore$Num_Customers, 7,  "Baltimore")
friedman_test_k(baltimore$Num_Customers, 28, "Baltimore")

# Lancaster
friedman_test_k(lancaster$Num_Customers, 7,  "Lancaster")
friedman_test_k(lancaster$Num_Customers, 28, "Lancaster")

# Philadelphia
friedman_test_k(philadelphia$Num_Customers, 7,  "Philadelphia")
friedman_test_k(philadelphia$Num_Customers, 28, "Philadelphia")

# Richmond
friedman_test_k(richmond$Num_Customers, 7,  "Richmond")
friedman_test_k(richmond$Num_Customers, 28, "Richmond")