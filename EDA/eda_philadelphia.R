library(tidyverse)
library(lubridate)
library(corrplot)
library(DataExplorer)

# Configuração inicial dos ficheiros
store_file <- "philadelphia.csv"
store_name <- "Philadelphia"

# Carregar o CSV 
df <- read.csv(store_file)
df$Date <- as.Date(df$Date)
df <- df %>% arrange(Date)

# Feature Engineering
day_levels <- c("segunda-feira", "terça-feira", "quarta-feira",
                "quinta-feira", "sexta-feira", "sábado", "domingo")

df$Month   <- factor(month(df$Date), levels = 1:12,
                     labels = c("Jan","Fev","Mar","Abr","Mai","Jun",
                                "Jul","Ago","Set","Out","Nov","Dez"))
df$DayOfWeek <- factor(weekdays(df$Date), levels = day_levels)
df$Weekend   <- df$DayOfWeek %in% c("sábado", "domingo")
df$Year      <- year(df$Date)

# Visão geral do dataset
cat("Dimensão\n");              print(dim(df))
cat("Estrutura\n");             str(df)
cat("Estatísticas Descritivas\n"); print(summary(df))

# Missing values
cat("Missing Values\n")
print(colSums(is.na(df)))

# Visualização dos missing values
plot_missing(df, title = paste("Missing Values -", store_name))

# Imputação por zero
df$Pct_On_Sale[is.na(df$Pct_On_Sale)] <- 0
cat("Após imputação:", sum(is.na(df$Pct_On_Sale)), "\n")

# Dias com 0 clientes e/ou coom 0 vendas
cat("Dias com 0 clientes:", sum(df$Num_Customers == 0), "\n")
cat("Dias com 0 vendas:",   sum(df$Sales == 0), "\n")


# Distribuição e outliers de Num_Customers
par(mfrow = c(1, 2))
boxplot(df$Num_Customers,
        main = paste("Outliers -", store_name),
        col  = "steelblue", ylab = "Num_Customers")
hist(df$Num_Customers,
     main = paste("Distribuição -", store_name),
     col  = "steelblue", xlab = "Num_Customers")
par(mfrow = c(1, 1)) # Outliers não são removidos porque correspondem a eventos reais



# Correlação entre as variáveis da loja de Philadelphia
numeric_vars <- df %>%
  select(Num_Employees, Num_Customers, Pct_On_Sale, Sales)

cor_matrix <- cor(numeric_vars, use = "complete.obs")

corrplot(cor_matrix,
         method      = "color",
         addCoef.col = "white",
         tl.col      = "black",
         tl.srt      = 45,
         col         = colorRampPalette(c("#d0e8f1", "#2171b5"))(200),
         mar         = c(0, 0, 2, 0),
         title       = paste("Correlação entre Variáveis -", store_name))


# Gráfico Sales vs Num_Costumers
ggplot(df, aes(Sales, Num_Customers)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = paste("Sales vs Customers -", store_name),
       x = "Sales",
       y = "Num_Customers") +
  theme_minimal()


# Relação entre número de empregados e número de clientes
ggplot(df, aes(Num_Employees, Num_Customers)) +
  geom_point(alpha = 0.4, color = "darkgreen") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = paste("Employees vs Customers -", store_name),
       x = "Num_Employees",
       y = "Num_Customers") +
  theme_minimal()


# Série Temporal
# Visualizar tendência, picos e padrões de Philadelphia
ggplot(df, aes(Date, Num_Customers)) +
  geom_line(color = "steelblue") +
  labs(title = paste("Num_Customers Over Time -", store_name),
       x = "Date", y = "Customers") +
  theme_minimal()


# Padrão Semanal
ggplot(df, aes(DayOfWeek, Num_Customers)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  labs(title = paste("Customers by Day of Week -", store_name),
       x = "Day of Week", y = "Num_Customers") +
  theme_minimal() # Detalhe por loja - confirma que DayOfWeek é preditor relevante


# Weekend vs Weekday
ggplot(df, aes(x = factor(Weekend, labels = c("Weekday","Weekend")),
               y = Num_Customers)) +
  geom_boxplot(fill = c("steelblue","tomato"), alpha = 0.7) +
  labs(title = paste("Weekend vs Weekday -", store_name),
       x = "", y = "Num_Customers") +
  theme_minimal()


# Consistência do padrão semanal ao longo dos anos
df %>%
  group_by(Year, DayOfWeek) %>%
  summarise(avg = mean(Num_Customers), .groups = "drop") %>%
  ggplot(aes(x = DayOfWeek, y = avg, fill = factor(Year))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(avg, 0)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  labs(title = paste("Avg Customers by Day of Week and Year -", store_name),
       x = "Day of Week", y = "Average Customers", fill = "Year") +
  theme_minimal()


# Sazonalidade Mensal
df %>%
  group_by(Year, Month) %>%
  summarise(avg = mean(Num_Customers, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Month, y = avg, fill = factor(Year))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(avg, 1)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  labs(title = paste("Avg Customers by Month and Year -", store_name),
       x = "Mês", y = "Média de Clientes", fill = "Ano") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()


# Eventos turísticos e o seu impacto
ggplot(df, aes(TouristEvent, Num_Customers)) +
  geom_boxplot(fill = "purple", alpha = 0.7) +
  labs(title = paste("Impact of Tourist Events -", store_name),
       x = "Tourist Event", y = "Num_Customers") +
  theme_minimal()


# Relação entre quantidade de Promoções e número de Clientes
ggplot(df, aes(Pct_On_Sale, Num_Customers)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = paste("Promotions vs Customers -", store_name),
       x = "Pct_On_Sale", y = "Num_Customers") +
  theme_minimal()

# Capacidade de atendimento
df_cap <- df %>%
  mutate(Cust_per_Emp = Num_Customers / Num_Employees) %>%
  filter(!is.infinite(Cust_per_Emp), !is.na(Cust_per_Emp))

ggplot(df_cap, aes(x = DayOfWeek, y = Cust_per_Emp)) +
  geom_boxplot(fill = "indianred", alpha = 0.7) +
  geom_hline(aes(yintercept = median(Cust_per_Emp, na.rm = TRUE)),
             linetype = "dashed", color = "black") +
  labs(title = paste("Clientes por Empregado -", store_name),
       subtitle = "Linha tracejada = mediana geral",
       x = "Dia da Semana", y = "Rácio (Clientes / Empregado)") +
  theme_minimal()

# Guardar dataset
write.csv(df, paste0(tolower(store_name), "_clean.csv"), row.names = FALSE)
cat("Dataset limpo guardado:", paste0(tolower(store_name), "_clean.csv"), "\n")

