# utils.R
# Funções auxiliares partilhadas por todos os scripts de forecasting

# CARREGAR DADOS
# Lê o CSV de uma loja, prepara os dados e devolve uma lista com
# os objetos necessários para todos os scripts de forecasting
load_store <- function(store_file, store_name) {
  
  df <- read.csv(store_file)          # lê o ficheiro CSV para um dataframe
  df$Date <- as.Date(df$Date)         # converte a coluna Date de texto para formato data
  df <- df[order(df$Date), ]          # ordena por data crescente para garantir índices baixos = passado
  df$Pct_On_Sale[is.na(df$Pct_On_Sale)] <- 0  # imputa o único valor omisso por zero
  
  S  <- df$Num_Customers   # extrai a série temporal alvo - vetor com 714 valores diários
  L  <- length(S)          # tamanho da série - 714
  YR <- diff(range(S))     # range global = max - min, calculado UMA vez aqui para o NMAE ser comparável em todas as iterações
  
  # imprime no console um resumo da loja carregada
  cat("Loja:", store_name,
      "| Dias:", L,
      "| Range:", YR,
      "| Período:", as.character(df$Date[1]),
      "a", as.character(df$Date[L]), "\n\n")
  
  # devolve uma lista com quatro objetos acessíveis via $
  list(df = df, S = S, L = L, YR = YR)
}

# DIVISÃO TREINO/TESTE
# Usado na Fase I para a avaliação simples com uma única semana de teste
train_test_split <- function(S, L, H, K) {
  
  LTR   <- L - H                        
  TR_ts <- ts(S[1:LTR], frequency = K)  
  Y     <- S[(LTR + 1):L]              
  
  list(LTR = LTR, TR_ts = TR_ts, Y = Y) 
}

# MOSTRAR RESULTADO (reshow)
# Calcula as métricas de erro, imprime no console e mostra gráfico
reshow <- function(label, Y, Pred, store_name, YR, plot = TRUE, show_legend = TRUE) {
  
  nmae <- mmetric(Y, Pred, metric = "NMAE", val = YR)
  rmse <- mmetric(Y, Pred, metric = "RMSE")
  r2   <- mmetric(Y, Pred, metric = "R22")
  
  cat(sprintf("  %-20s NMAE: %6.2f%%  RMSE: %8.2f  R2: %.4f\n",
              label, nmae, rmse, r2))
  
  if (plot) {
    
    if (show_legend) {
      mgraph(Y, Pred,
             graph = "REG", Grid = 10,
             main  = paste(store_name, "-", label,
                           "| NMAE =", round(nmae, 1), "%"),
             col   = c("black", "blue"),
             leg   = list(pos = "topleft",
                          leg = c("target", paste(label, "pred."))))
    } else {
      mgraph(Y, Pred,
             graph = "REG", Grid = 10,
             main  = paste(store_name, "-", label,
                           "| NMAE =", round(nmae, 1), "%"),
             col   = c("black", "blue"))
    }
  }
  
  invisible(c(NMAE = nmae, RMSE = rmse, R2 = r2))
}

# Gráfico com legenda exterior (Fase I)
plot_with_legend <- function(label, Y, Pred, store_name, YR, legend_labels = NULL) {
  
  old_par <- par(no.readonly = TRUE)
  par(oma = c(4, 0, 0, 0))
  
  r <- reshow(label, Y, Pred, store_name, YR, show_legend = FALSE)
  
  par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
  plot.new()
  
  if (is.null(legend_labels)) {
    legend_labels <- c("target", paste(label, "pred."))
  }
  
  legend("bottom",
         legend = legend_labels,
         col    = c("black", "blue"),
         lty    = 1,
         lwd    = 2,
         horiz  = TRUE,
         bty    = "n",
         cex    = 0.9)
  
  par(old_par)
  
  invisible(r)
}

# Gráfico GW com legenda (Fase II)
plot_gw_legend <- function(all_Y, all_preds, pred_labels, colors, store_name, model_name) {
  
  old_par <- par(no.readonly = TRUE)
  par(oma = c(3, 0, 0, 0), mar = c(4, 4, 4, 2))
  
  n          <- length(all_Y)
  ylim_range <- range(c(all_Y, unlist(all_preds)))
  
  plot(1:n, all_Y,
       type = "l", col = "black", lwd = 2,
       ylim = ylim_range,
       xlab = "Observações (períodos de teste acumulados)",
       ylab = "Num_Customers",
       main = paste("Growing Window -", model_name, "|", store_name))
  
  for (i in seq_along(all_preds)) {
    lines(all_preds[[i]], col = colors[i], lty = 2)
  }
  
  par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
  plot.new()
  
  legend("bottom",
         legend = c("Real", pred_labels),
         col    = c("black", colors),
         lty    = c(1, rep(2, length(pred_labels))),
         lwd    = 2, horiz = TRUE, bty = "n", cex = 0.9)
  
  par(old_par)
}

# Tabela Resumo
summary_table <- function(results_list, model_names) {
  
  df_res <- data.frame(
    Modelo   = model_names,
    NMAE_pct = round(sapply(results_list, function(r) r["NMAE"]), 4),  # extrai NMAE de cada modelo e arredonda
    RMSE     = round(sapply(results_list, function(r) r["RMSE"]), 2),  # extrai RMSE de cada modelo e arredonda
    R2       = round(sapply(results_list, function(r) r["R2"]),   4),  # extrai R² de cada modelo e arredonda
    row.names = NULL
  )
  
  df_res[order(df_res$NMAE_pct), ]  # ordena por NMAE crescente
}

# Gráfico Growing Window
plot_window <- function(all_Y, pred_list, labels, colors, store_name, type = "Growing") {
  
  n          <- length(all_Y)                       
  ylim_range <- range(c(all_Y, unlist(pred_list)))   
  
  # desenha a linha preta com os valores reais acumulados
  plot(1:n, all_Y,
       type = "l", col = "black", lwd = 2,
       ylim = ylim_range,
       xlab = "Observações (períodos de teste acumulados)",
       ylab = "Num_Customers",
       main = paste(type, "Window — Previsões acumuladas |", store_name))
  
  for (i in seq_along(pred_list)) {
    lines(pred_list[[i]], col = colors[i], lty = 2)  
  }
  
  legend("topright", bty = "n", cex = 0.8,
         legend = c("Real", labels),
         col    = c("black", colors),
         lty    = c(1, rep(2, length(labels)))) 
}

# Guardar previsões finais
save_forecasts <- function(Pred_final, df, store_name, H) {
  
  Pred_final   <- round(pmax(Pred_final, 0))  
  last_date    <- df$Date[length(df$Date)]    
  future_dates <- seq(last_date + 1, by = "day", length.out = H)  
  
  out <- data.frame(
    Store               = store_name,    
    Date                = future_dates,  
    Weekday             = weekdays(future_dates),  
    Predicted_Customers = Pred_final  
  )
  
  output_csv <- paste0("previsoes_", tolower(store_name), ".csv") 
  write.csv(out, output_csv, row.names = FALSE)  
  cat("Previsões guardadas em:", output_csv, "\n")
  print(out, row.names = FALSE)
  
  invisible(out)
}