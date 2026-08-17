library(shiny)
library(ggplot2)
library(dplyr)
library(shinythemes)
library(corrplot)


# 1. Configuração inicial dos ficheiros
lojas <- c("Baltimore", "Lancaster", "Philadelphia", "Richmond")

cost_fixed_vals <- c(
  Baltimore    = 700,
  Lancaster    = 730,
  Philadelphia = 760,
  Richmond     = 800
)

graficos_eda <- c(
  "Série Temporal"            = "serie_temporal",
  "Média por Dia da Semana"   = "media_dia",
  "Distribuição (Boxplot)"    = "boxplot",
  "Correlação entre Lojas"    = "correlacao",
  "CCF entre Lojas"           = "ccf",
  "ACF"                       = "acf",
  "PACF"                      = "pacf",
  "Decomposição STL"          = "stl",
  "Testes de Friedman"        = "friedman"
)

# 2. Carregar dados EDA
eda <- readRDS("eda_data.rds")

# 3. Carregar dados de forecast
carregar_dados <- function() {
  lista_dados <- list()
  for (l in lojas) {
    nome_ficheiro <- paste0("previsoes_gw_", tolower(l), ".csv")
    if (file.exists(nome_ficheiro)) {
      df <- read.csv(nome_ficheiro)
      df$Date <- as.Date(df$Date)
      lista_dados[[l]] <- df
    }
  }
  return(lista_dados)
}

dados_previsoes <- carregar_dados()


# 4. Carregar resultados de otimização

res_O1   <- readRDS("resultados_O1_GA.rds")
res_O2   <- readRDS("resultados_O2_SANN.rds")
res_O3   <- readRDS("resultados_O3_NSGA2.rds")
res_O3_s <- readRDS("resultados_O3_SANN.rds")


# 5. UI

ui <- navbarPage(
  title  = "DSS: Gestão de Retalho",
  theme  = shinythemes::shinytheme("flatly"),
  
  # EDA (PÁGINA INICIAL)
  tabPanel("Análise Exploratória (EDA)",
           sidebarLayout(
             sidebarPanel(
               width = 3,
               h4("Configurações"),
               hr(),
               selectInput("eda_loja", "Selecione a Loja:",
                           choices = lojas,
                           selected = "Baltimore"),
               helpText("Usado nos gráficos ACF, PACF e STL."),
               br(),
               selectInput("eda_grafico", "Selecione o Gráfico:",
                           choices = graficos_eda,
                           selected = "serie_temporal"),
               br(),
               conditionalPanel(
                 condition = "input.eda_grafico == 'acf' ||
                       input.eda_grafico == 'pacf' ||
                       input.eda_grafico == 'stl'",
                 wellPanel(
                   style = "background-color: #eaf4fb; border-left: 3px solid #3498db;",
                   helpText(
                     icon("info-circle"),
                     "Este gráfico usa a loja selecionada acima."
                   )
                 )
               ),
               conditionalPanel(
                 condition = "input.eda_grafico == 'serie_temporal' ||
                       input.eda_grafico == 'media_dia'     ||
                       input.eda_grafico == 'boxplot'       ||
                       input.eda_grafico == 'correlacao'    ||
                       input.eda_grafico == 'ccf'           ||
                       input.eda_grafico == 'friedman'",
                 wellPanel(
                   style = "background-color: #eafaf1; border-left: 3px solid #27ae60;",
                   helpText(
                     icon("info-circle"),
                     "Este gráfico mostra todas as lojas."
                   )
                 )
               )
             ),
             mainPanel(
               width = 9,
               fluidRow(
                 column(12,
                        wellPanel(
                          style = "background-color: #f8f9fa;",
                          htmlOutput("eda_titulo")
                        )
                 )
               ),
               br(),
               # Gráficos ggplot
               conditionalPanel(
                 condition = "input.eda_grafico != 'correlacao' &&
                       input.eda_grafico != 'acf'        &&
                       input.eda_grafico != 'pacf'       &&
                       input.eda_grafico != 'ccf'        &&
                       input.eda_grafico != 'stl'        &&
                       input.eda_grafico != 'friedman'",
                 plotOutput("eda_plot_gg", height = "500px")
               ),
               # Gráficos base R
               conditionalPanel(
                 condition = "input.eda_grafico == 'correlacao' ||
                       input.eda_grafico == 'acf'        ||
                       input.eda_grafico == 'pacf'       ||
                       input.eda_grafico == 'ccf'        ||
                       input.eda_grafico == 'stl'",
                 plotOutput("eda_plot_base", height = "500px")
               ),
               # Tabela Friedman
               conditionalPanel(
                 condition = "input.eda_grafico == 'friedman'",
                 h4("Resultados dos Testes de Sazonalidade de Friedman"),
                 br(),
                 tableOutput("eda_friedman_tabela"),
                 br(),
                 wellPanel(
                   style = "background-color: #fef9e7; border-left: 3px solid #f39c12;",
                   helpText(
                     icon("lightbulb"),
                     strong("Interpretação:"),
                     " p-value < 0.05 indica sazonalidade estatisticamente significativa.",
                     " k=7 testa sazonalidade semanal; k=28 testa sazonalidade mensal."
                   )
                 )
               )
             )
           )
  ),
  

  # FORECASTING
  tabPanel("Dashboard de Previsão",
           sidebarLayout(
             sidebarPanel(
               width = 3,
               h4("Configurações de Visualização"),
               hr(),
               selectInput("loja_sel", "Selecione a Loja:",
                           choices = lojas),
               selectInput("run_sel", "Selecione a Semana (Run):",
                           choices = 1:20),
               br(),
               helpText("Visualização das previsões.")
             ),
             mainPanel(
               width = 9,
               fluidRow(
                 column(12,
                        wellPanel(
                          style = "background-color: #f8f9fa;",
                          htmlOutput("info_semana")
                        )
                 )
               ),
               br(),
               plotOutput("grafico_previsao", height = "400px"),
               br(),
               h4("Detalhes Diários"),
               tableOutput("tabela_previsao"),
               htmlOutput("metricas_previsao")
             )
           )
  ),
  

  # OTIMIZAÇÃO O1/O2
  tabPanel("Otimização (O1/O2)",
           sidebarLayout(
             sidebarPanel(
               width = 3,
               h4("Configurações"),
               hr(),
               selectInput("opt_model", "Modelo:",
                           choices = c("O1 - Lucro Máximo"        = "O1",
                                       "O2 - Restrição de Unidades" = "O2")),
               selectInput("opt_loja", "Selecione a Loja:",
                           choices = lojas),
               selectInput("opt_run", "Semana:",
                           choices = 1:20),
               br(),
               helpText("Plano otimizado e KPIs.")
             ),
             mainPanel(
               width = 9,
               fluidRow(
                 column(12,
                        wellPanel(
                          style = "background-color: #f8f9fa;",
                          htmlOutput("opt_info")
                        )
                 )
               ),
               br(),
               plotOutput("opt_convergencia", height = "300px"),
               br(),
               h4("Detalhes Diários"),
               tableOutput("opt_tabela")
             )
           )
  ),
  

  # DECISÃO MULTI-OBJETIVO O3
  tabPanel("Decisão Multi-Objetivo (O3)",
           sidebarLayout(
             sidebarPanel(
               width = 3,
               h4("Configurações"),
               hr(),
               selectInput("o3_run", "Semana:",
                           choices = 1:20),
               selectInput("o3_loja", "Selecione a Loja:",
                           choices = lojas),
               br(),
               radioButtons("o3_algo_click", "Algoritmo para consultar plano:",
                            choices  = c("NSGA-II" = "nsga2", "SANN" = "sann"),
                            selected = "nsga2"),
               br(),
               helpText("Selecione o algoritmo e clique num ponto da
                  fronteira de Pareto para ver o plano e os KPIs.")
             ),
             mainPanel(
               width = 9,
               fluidRow(
                 column(12,
                        wellPanel(
                          style = "background-color: #f8f9fa;",
                          htmlOutput("o3_info_semana")
                        )
                 )
               ),
               br(),
               h4("Fronteira de Pareto"),
               plotOutput("o3_pareto", height = "400px", click = "o3_click"),
               br(),
               conditionalPanel(
                 condition = "output.o3_ponto_selecionado",
                 wellPanel(
                   style = "background-color: #f8f9fa;",
                   htmlOutput("o3_kpis")
                 ),
                 br(),
                 h4("Detalhes Diários — Ponto Selecionado"),
                 tableOutput("o3_tabela")
               ),
               br(),
               uiOutput("o3_metricas_algo")
             )
           )
  )
)


# 6. Server

server <- function(input, output, session) {
  

  # EDA
  # Título dinâmico
  output$eda_titulo <- renderUI({
    label <- names(graficos_eda)[graficos_eda == input$eda_grafico]
    depende_loja <- input$eda_grafico %in% c("acf", "pacf", "stl")
    loja_txt <- if (depende_loja) paste0(" — ", input$eda_loja) else " — Todas as Lojas"
    HTML(paste0(
      "<div style='font-size:18px;'>",
      "<b>Gráfico:</b> ", label, loja_txt,
      "</div>"
    ))
  })
  
  # Gráficos ggplot2
  output$eda_plot_gg <- renderPlot({
    g <- input$eda_grafico
    day_levels <- eda$day_levels
    
    if (g == "serie_temporal") {
      ggplot(eda$df_all, aes(x = Date, y = Num_Customers, color = Store)) +
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
      
    } else if (g == "media_dia") {
      eda$df_all %>%
        mutate(DayOfWeek = factor(weekdays(Date), levels = day_levels)) %>%
        group_by(Store, DayOfWeek) %>%
        summarise(avg = mean(Num_Customers), .groups = "drop") %>%
        ggplot(aes(DayOfWeek, avg, fill = Store)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(
          title = "Average Customers by Day of Week and Store",
          x     = "Day of Week",
          y     = "Average Customers"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          plot.title  = element_text(face = "bold"),
          axis.text.x = element_text(angle = 20, hjust = 1)
        )
      
    } else if (g == "boxplot") {
      ggplot(eda$df_all, aes(Store, Num_Customers, fill = Store)) +
        geom_boxplot() +
        labs(
          title = "Distribution of Customers by Store",
          x     = "Store",
          y     = "Num_Customers"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          legend.position = "none",
          plot.title      = element_text(face = "bold")
        )
    }
  })
  
  # Gráficos base R
  output$eda_plot_base <- renderPlot({
    g    <- input$eda_grafico
    loja <- input$eda_loja
    
    if (g == "correlacao") {
      corrplot(
        eda$cor_matrix,
        method      = "color",
        addCoef.col = "black",
        tl.col      = "black",
        title       = "Correlation Between Stores - Num_Customers",
        mar         = c(0, 0, 2, 0)
      )
      
    } else if (g == "acf") {
      dados_loja <- eda[[tolower(loja)]]$Num_Customers
      acf(dados_loja, lag.max = 35,
          main = paste("ACF -", loja))
      
    } else if (g == "pacf") {
      dados_loja <- eda[[tolower(loja)]]$Num_Customers
      pacf(dados_loja, lag.max = 35,
           main = paste("PACF -", loja))
      
    } else if (g == "ccf") {
      plot_ccf_pos <- function(x, y, lag.max = 14, main = "") {
        cc       <- ccf(x, y, lag.max = lag.max, plot = FALSE)
        idx      <- cc$lag[,,1] >= 0
        lags_pos <- cc$lag[,,1][idx]
        acf_pos  <- cc$acf[,,1][idx]
        ci       <- qnorm(0.975) / sqrt(length(x))
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
      plot_ccf_pos(eda$baltimore$Num_Customers,    eda$lancaster$Num_Customers,
                   lag.max = 14, main = "CCF: Baltimore vs Lancaster")
      plot_ccf_pos(eda$baltimore$Num_Customers,    eda$philadelphia$Num_Customers,
                   lag.max = 14, main = "CCF: Baltimore vs Philadelphia")
      plot_ccf_pos(eda$baltimore$Num_Customers,    eda$richmond$Num_Customers,
                   lag.max = 14, main = "CCF: Baltimore vs Richmond")
      plot_ccf_pos(eda$lancaster$Num_Customers,    eda$philadelphia$Num_Customers,
                   lag.max = 14, main = "CCF: Lancaster vs Philadelphia")
      plot_ccf_pos(eda$lancaster$Num_Customers,    eda$richmond$Num_Customers,
                   lag.max = 14, main = "CCF: Lancaster vs Richmond")
      plot_ccf_pos(eda$philadelphia$Num_Customers, eda$richmond$Num_Customers,
                   lag.max = 14, main = "CCF: Philadelphia vs Richmond")
      par(mfrow = c(1, 1))
      
    } else if (g == "stl") {
      obj <- switch(loja,
                    Baltimore    = eda$stl_baltimore,
                    Lancaster    = eda$stl_lancaster,
                    Philadelphia = eda$stl_philadelphia,
                    Richmond     = eda$stl_richmond
      )
      plot(obj, main = paste("STL Decomposition -", loja))
    }
  })
  
  # Tabela Friedman
  output$eda_friedman_tabela <- renderTable({
    df <- eda$friedman_df
    df$Sazonalidade <- ifelse(df$k == 7, "Semanal (k=7)", "Mensal (k=28)")
    df$Significativa <- ifelse(df$significant, "✔ Sim", "✘ Não")
    df %>%
      select(
        Loja         = store,
        Sazonalidade,
        `Chi-squared` = chi_squared,
        `p-value`     = p_value,
        Significativa
      )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # FORECASTING
  dados_filtrados <- reactive({
    req(input$loja_sel, input$run_sel)
    df <- dados_previsoes[[input$loja_sel]]
    df %>% filter(Run == as.numeric(input$run_sel))
  })
  
  output$info_semana <- renderUI({
    df <- dados_filtrados()
    total_previsto <- sum(df$Predicted_Customers)
    total_real_html <- ""
    if ("Actual_Customers" %in% names(df)) {
      total_real_html <- paste0(
        " | <b>Total Real:</b> <span style='color:#2c3e50; font-weight:bold;'>",
        sum(df$Actual_Customers), "</span>"
      )
    }
    HTML(paste0(
      "<div style='font-size: 18px;'>",
      "<b>Loja:</b> ", input$loja_sel, " | ",
      "<b>Total Previsto:</b> <span style='color:#3498db; font-weight:bold;'>",
      total_previsto, "</span>",
      total_real_html,
      "</div>"
    ))
  })
  
  output$grafico_previsao <- renderPlot({
    df <- dados_filtrados()
    p <- ggplot(df, aes(x = Date)) +
      theme_minimal() +
      labs(title = paste("Análise de Previsão — Semana", input$run_sel),
           x = "Data", y = "Clientes") +
      theme(legend.position = "bottom")
    
    if ("Actual_Customers" %in% names(df)) {
      p <- p +
        geom_line(aes(y = Predicted_Customers, color = "Previsto"), linewidth = 1.2) +
        geom_point(aes(y = Predicted_Customers, color = "Previsto"), size = 3) +
        geom_line(aes(y = Actual_Customers, color = "Real"), linewidth = 1.2, linetype = "dashed") +
        geom_point(aes(y = Actual_Customers, color = "Real"), size = 3) +
        scale_color_manual(values = c("Previsto" = "#3498db", "Real" = "#e74c3c")) +
        labs(color = "Série:")
    } else {
      p <- p +
        geom_line(aes(y = Predicted_Customers), color = "#3498db", linewidth = 1.2) +
        geom_point(aes(y = Predicted_Customers), color = "#2c3e50", size = 3) +
        geom_area(aes(y = Predicted_Customers), fill = "#3498db", alpha = 0.2)
    }
    p
  })
  
  output$tabela_previsao <- renderTable({
    df <- dados_filtrados()
    df$Date <- as.character(df$Date)
    df
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$metricas_previsao <- renderUI({
    df <- dados_filtrados()
    if ("Actual_Customers" %in% names(df)) {
      real   <- df$Actual_Customers
      prev   <- df$Predicted_Customers
      mae    <- mean(abs(real - prev))
      nmae   <- mae / mean(real)
      ss_res <- sum((real - prev)^2)
      ss_tot <- sum((real - mean(real))^2)
      r2     <- if (ss_tot > 0) 1 - (ss_res / ss_tot) else NA
      HTML(paste0(
        "<div style='background-color:#f4f6f7; border-left:5px solid #34495e;",
        "padding:15px; border-radius:4px; margin-top:15px;'>",
        "<h4 style='margin-top:0; color:#2c3e50; font-weight:bold;'>",
        "Métricas de Avaliação (Random Forest)</h4>",
        "<p style='font-size:15px; margin-bottom:0;'>",
        "<b>NMAE:</b> <span style='color:#e67e22; font-weight:bold;'>",
        round(nmae, 4), "</span>",
        " &nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp; ",
        "<b>R²:</b> <span style='color:#27ae60; font-weight:bold;'>",
        if (is.na(r2)) "N/A" else round(r2, 4), "</span>",
        "</p></div>"
      ))
    } else {
      HTML("<p style='color:#7f8c8d; font-style:italic; margin-top:10px;'>
            Nota: Coluna 'Actual_Customers' não encontrada.</p>")
    }
  })
  

  # OTIMIZAÇÃO O1/O2
  dados_opt <- reactive({
    req(input$opt_model, input$opt_run)
    if (input$opt_model == "O1") {
      res_O1[[as.numeric(input$opt_run)]]
    } else {
      res_O2[[as.numeric(input$opt_run)]]
    }
  })
  
  output$opt_info <- renderUI({
    res      <- dados_opt()
    df       <- res$detalhes
    df_loja  <- df[df$Loja == input$opt_loja, ]
    lucro_loja    <- sum(df_loja$Lucro_Dia) - cost_fixed_vals[input$opt_loja]
    vendas_loja   <- sum(df_loja$Vendas_USD)
    custos_loja   <- sum(df_loja$Custo_HR)
    unidades_loja <- sum(df_loja$Unidades)
    hr_loja       <- sum(df_loja$J + df_loja$X)
    lucro_global    <- sum(df$Lucro_Dia) - sum(cost_fixed_vals)
    unidades_global <- sum(df$Unidades)
    hr_global       <- sum(df$J + df$X)
    HTML(paste0(
      "<div style='font-size:18px;'>",
      "<b>Global (todas as lojas)</b><br>",
      "<b>Lucro:</b> <span style='color:green;'>", round(lucro_global), "</span> USD | ",
      "<b>Unidades:</b> ", unidades_global, " | ",
      "<b>HR:</b> ", hr_global,
      "<br><br>",
      "<b>Loja:</b> ", input$opt_loja, "<br>",
      "<b>Lucro:</b> <span style='color:green;'>", round(lucro_loja), "</span> USD | ",
      "<b>Vendas:</b> ", round(vendas_loja), " USD | ",
      "<b>Custos HR:</b> ", round(custos_loja), " USD | ",
      "<b>Unidades:</b> ", unidades_loja, " | ",
      "<b>HR:</b> ", hr_loja,
      "</div>"
    ))
  })
  
  output$opt_convergencia <- renderPlot({
    res    <- dados_opt()
    curvas <- res$curvas
    if (is.null(curvas)) return(NULL)
    df_plot <- data.frame()
    if (input$opt_model == "O1") {
      for (i in 1:length(curvas)) {
        curva <- curvas[[i]]
        df_plot <- rbind(df_plot, data.frame(
          Iter  = 1:length(curva),
          Valor = curva,
          Serie = paste0("Run ", i)
        ))
      }
      titulo  <- "Convergência GA — O1"
      label_x <- "Geração"
    } else {
      temperaturas <- c(10, 250, 500, 750, 1000)
      for (i in 1:length(curvas)) {
        curva <- curvas[[i]]
        if (is.null(curva) || all(is.na(curva))) next
        df_plot <- rbind(df_plot, data.frame(
          Iter  = 1:length(curva),
          Valor = curva,
          Serie = paste0("Temp ", temperaturas[i])
        ))
      }
      titulo  <- "Convergência SANN — O2"
      label_x <- "Avaliações"
    }
    ggplot(df_plot, aes(x = Iter, y = Valor, color = Serie)) +
      geom_line() +
      theme_minimal() +
      labs(title = titulo, x = label_x, y = "Lucro (USD)", color = NULL)
  })
  
  output$opt_tabela <- renderTable({
    res <- dados_opt()
    df  <- res$detalhes
    df[df$Loja == input$opt_loja, ]
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  

  # O3 - MULTI-OBJETIVO
  dados_o3_nsga <- reactive({
    req(input$o3_run)
    res_O3[[as.numeric(input$o3_run)]]
  })
  
  dados_o3_sann <- reactive({
    req(input$o3_run)
    res_O3_s[[as.numeric(input$o3_run)]]
  })
  
  ponto_idx  <- reactiveVal(NULL)
  ponto_algo <- reactiveVal(NULL)
  
  observeEvent(input$o3_run,        { ponto_idx(NULL); ponto_algo(NULL) })
  observeEvent(input$o3_algo_click, { ponto_idx(NULL); ponto_algo(NULL) })
  
  observeEvent(input$o3_click, {
    click <- input$o3_click
    if (is.null(click)) return()
    res <- if (input$o3_algo_click == "nsga2") dados_o3_nsga() else dados_o3_sann()
    if (length(res$pareto_HR) == 0) return()
    dx   <- res$pareto_HR     - click$x
    dy   <- res$pareto_profit - click$y
    dist <- sqrt(dx^2 + dy^2)
    ponto_idx(which.min(dist))
    ponto_algo(input$o3_algo_click)
  })
  
  output$o3_ponto_selecionado <- reactive({ !is.null(ponto_idx()) })
  outputOptions(output, "o3_ponto_selecionado", suspendWhenHidden = FALSE)
  
  output$o3_info_semana <- renderUI({
    res_n <- dados_o3_nsga()
    res_s <- dados_o3_sann()
    HTML(paste0(
      "<div style='font-size:18px;'>",
      "<b>Semana:</b> ", res_n$semana, " | ",
      "<b>Soluções NSGA-II:</b> ", length(res_n$pontos_pareto), " | ",
      "<b>Soluções SANN:</b> ", length(res_s$pontos_pareto), " | ",
      "<span style='color:#7f8c8d;'>Selecione o algoritmo e clique num ponto.</span>",
      "</div>"
    ))
  })
  
  output$o3_pareto <- renderPlot({
    res_n <- dados_o3_nsga()
    res_s <- dados_o3_sann()
    idx   <- ponto_idx()
    algo  <- ponto_algo()
    p <- ggplot() +
      theme_minimal() +
      labs(title = paste("Fronteira de Pareto — Semana", input$o3_run),
           x = "Total HR (Colaboradores)", y = "Lucro Total (USD)", color = "Algoritmo")
    if (length(res_n$pareto_HR) > 0) {
      df_n <- data.frame(HR = res_n$pareto_HR, Lucro = res_n$pareto_profit, Algoritmo = "NSGA-II")
      p <- p +
        geom_line(data = df_n,  aes(x = HR, y = Lucro, color = Algoritmo), linewidth = 0.8) +
        geom_point(data = df_n, aes(x = HR, y = Lucro, color = Algoritmo), size = 3)
    }
    if (length(res_s$pareto_HR) > 0) {
      df_s <- data.frame(HR = res_s$pareto_HR, Lucro = res_s$pareto_profit, Algoritmo = "SANN")
      p <- p +
        geom_line(data = df_s,  aes(x = HR, y = Lucro, color = Algoritmo), linewidth = 0.8, linetype = "dashed") +
        geom_point(data = df_s, aes(x = HR, y = Lucro, color = Algoritmo), size = 3)
    }
    p <- p + scale_color_manual(values = c("NSGA-II" = "steelblue", "SANN" = "darkorange"))
    if (!is.null(idx) && !is.null(algo)) {
      ponto_sel <- NULL
      if (algo == "nsga2" && length(res_n$pareto_HR) >= idx)
        ponto_sel <- data.frame(HR = res_n$pareto_HR[idx], Lucro = res_n$pareto_profit[idx])
      else if (algo == "sann" && length(res_s$pareto_HR) >= idx)
        ponto_sel <- data.frame(HR = res_s$pareto_HR[idx], Lucro = res_s$pareto_profit[idx])
      if (!is.null(ponto_sel))
        p <- p + geom_point(data = ponto_sel, aes(x = HR, y = Lucro),
                            color = "red", size = 5, shape = 8)
    }
    p
  })
  
  output$o3_kpis <- renderUI({
    req(ponto_idx())
    idx   <- ponto_idx()
    algo  <- ponto_algo()
    res   <- if (algo == "nsga2") dados_o3_nsga() else dados_o3_sann()
    ponto <- res$pontos_pareto[[idx]]
    df    <- ponto$detalhes
    df_loja <- df[df$Loja == input$o3_loja, ]
    lucro_loja     <- sum(df_loja$Lucro_Dia) - cost_fixed_vals[input$o3_loja]
    vendas_loja    <- sum(df_loja$Vendas_USD)
    custos_loja    <- sum(df_loja$Custo_HR)
    unidades_loja  <- sum(df_loja$Unidades)
    hr_loja        <- sum(df_loja$J + df_loja$X)
    lucro_total    <- sum(df$Lucro_Dia) - sum(cost_fixed_vals)
    unidades_total <- sum(df$Unidades)
    hr_total       <- ponto$hr
    HTML(paste0(
      "<div style='font-size:18px;'>",
      "<b>Algoritmo:</b> ", ifelse(algo == "nsga2", "NSGA-II", "SANN"), " | ",
      "<b>Lucro global:</b> <span style='color:green;'>", round(lucro_total), "</span> USD | ",
      "<b>HR global:</b> ", round(hr_total), " | ",
      "<b>Unidades globais:</b> ", unidades_total,
      "<br><br>",
      "<b>Loja:</b> ", input$o3_loja, " — ",
      "<b>Lucro:</b> <span style='color:green;'>", round(lucro_loja), "</span> USD | ",
      "<b>Vendas:</b> ", round(vendas_loja), " USD | ",
      "<b>Custos HR:</b> ", round(custos_loja), " USD | ",
      "<b>Unidades:</b> ", unidades_loja, " | ",
      "<b>HR:</b> ", hr_loja,
      "</div>"
    ))
  })
  
  output$o3_tabela <- renderTable({
    req(ponto_idx())
    algo  <- ponto_algo()
    res   <- if (algo == "nsga2") dados_o3_nsga() else dados_o3_sann()
    ponto <- res$pontos_pareto[[ponto_idx()]]
    df    <- ponto$detalhes
    df[df$Loja == input$o3_loja, ]
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$o3_metricas_algo <- renderUI({
    res_n <- dados_o3_nsga()
    res_s <- dados_o3_sann()
    fmt_usd <- function(x) {
      if (is.null(x) || is.na(x)) return("—")
      paste0("$", format(round(x), big.mark = ","))
    }
    fmt_num <- function(x) {
      if (is.null(x) || is.na(x)) return("—")
      format(round(x, 2), big.mark = ",")
    }
    make_card <- function(label, value, color = "var(--color-text-primary)") {
      tags$div(
        style = paste0(
          "background: var(--color-background-secondary);",
          "border-radius: 8px; padding: 1rem;",
          "display: flex; flex-direction: column; gap: 4px;"
        ),
        tags$span(style = "font-size: 13px; color: var(--color-text-secondary);", label),
        tags$span(style = paste0("font-size: 20px; font-weight: 500; color: ", color, ";"), value)
      )
    }
    make_algo_block <- function(label, res, color_accent) {
      tags$div(
        style = paste0(
          "background: var(--color-background-primary);",
          "border: 0.5px solid var(--color-border-tertiary);",
          "border-left: 3px solid ", color_accent, ";",
          "border-radius: 0 8px 8px 0;",
          "padding: 1rem 1.25rem; margin-bottom: 1rem;"
        ),
        tags$p(style = paste0(
          "font-size: 14px; font-weight: 500; margin: 0 0 12px;",
          "color: var(--color-text-primary);"
        ), label),
        tags$div(
          style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 10px;",
          make_card("Média do lucro",   fmt_usd(res$media_lucro),   "var(--color-text-success)"),
          make_card("Mediana do lucro", fmt_usd(res$mediana_lucro), "var(--color-text-success)"),
          make_card("Lucro máximo",     fmt_usd(max(res$pareto_profit, na.rm = TRUE)), "var(--color-text-success)"),
          make_card("Lucro mínimo",     fmt_usd(min(res$pareto_profit, na.rm = TRUE)), "var(--color-text-secondary)"),
          make_card("HR mínimo",        fmt_num(min(res$pareto_HR, na.rm = TRUE))),
          make_card("HR máximo",        fmt_num(max(res$pareto_HR, na.rm = TRUE))),
          make_card("Soluções Pareto",  as.character(length(res$pontos_pareto))),
          make_card("Hipervolume",      fmt_num(res$hv))
        )
      )
    }
    tagList(
      tags$h4(style = "margin: 0 0 0.75rem; font-size: 16px; font-weight: 500;",
              "Métricas por algoritmo — semana selecionada"),
      make_algo_block("NSGA-II", res_n, "#2980b9"),
      make_algo_block("SANN",    res_s, "#e67e22")
    )
  })
}


# 7. Run APP
shinyApp(ui, server)