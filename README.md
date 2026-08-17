# Predictive & Prescriptive Retail Decision Support System

An end-to-end **Intelligent Decision Support System (IDSS)** for retail store management, combining **time-series forecasting**, **machine learning**, **metaheuristic optimization** and an interactive **R Shiny** dashboard.

The system forecasts daily customer demand for four US stores over a 7-day horizon and converts those predictions into weekly operational decisions for **workforce allocation** and **product promotions**.

> Academic project developed for *Sistemas Adaptativos para Análises Preditivas e Prescritivas* (SAAPP), MSc in Data Engineering and Data Science, University of Minho, 2025/2026.  
> **Project grade: 18.2 / 20**

## Project Overview

The project follows a predictive-prescriptive analytics pipeline:

1. **Exploratory Data Analysis** — identify temporal patterns, seasonality, correlations and the effect of external variables.
2. **Demand Forecasting** — predict `Num_Customers` for each store up to 7 days ahead.
3. **Prescriptive Optimization** — use the forecasts to optimize staffing and promotion decisions under operational constraints.
4. **Decision Support System** — integrate EDA, forecasts and optimization results into an interactive Shiny application.

The four stores are located in:

- Baltimore, MD
- Lancaster, PA
- Philadelphia, PA
- Richmond, VA

## Data

Each store contains daily observations with the following main variables:

| Variable | Description |
|---|---|
| `Date` | Observation date |
| `Num_Employees` | Number of employees in the store |
| `Num_Customers` | Daily number of customers — forecasting target |
| `Pct_On_Sale` | Percentage of products offered at a discount |
| `TouristEvent` | Indicator of a local tourist event |
| `Sales` | Daily product sales |

The data were provided as part of the university project.

## 1. Exploratory Data Analysis

EDA was performed both independently for each store and jointly across stores. The analysis included:

- missing-value and outlier inspection;
- temporal evolution of customer demand;
- weekday and seasonal patterns;
- ACF and PACF analysis;
- STL decomposition;
- Friedman tests for seasonality;
- correlations and cross-correlations between stores;
- analysis of promotions and tourist events.

The analysis identified strong weekly seasonality and temporal dependencies, which informed the lag configurations used during forecasting.

## 2. Demand Forecasting

The forecasting task predicts daily customer demand for a maximum horizon of **7 days**.

### Univariate models

The following approaches were evaluated independently for each store:

- Seasonal Naive;
- Holt-Winters;
- ETS;
- NNETAR;
- AutoARIMA;
- MLPE;
- Random Forest.

### Multivariate models

The multivariate experiments incorporated exogenous variables and information shared across stores using:

- ARIMAX;
- VAR;
- MLPE;
- Random Forest.

Models were evaluated using a temporal **growing-window** procedure across **20 forecast weeks**, avoiding random train/test splits and preserving the chronological structure of the problem.

### Final forecasting model

The best overall model was a **multivariate Random Forest** using lag configuration:

```text
[1, 3, 7, 14, 21, 28]
```

It consistently outperformed the Seasonal Naive baseline across all four stores.

| Store | NMAE | R² | NMAE improvement vs. Seasonal Naive |
|---|---:|---:|---:|
| Baltimore | 1.27% | 0.7713 | 34.54% |
| Lancaster | 1.14% | 0.7444 | 19.72% |
| Philadelphia | 2.02% | 0.8299 | 27.34% |
| Richmond | 1.36% | 0.7795 | 39.56% |

These forecasts were then used as the input to the prescriptive optimization stage.

## 3. Prescriptive Optimization

For each store and day, the optimization stage determines three decision variables:

- number of **junior employees**;
- number of **expert employees**;
- percentage of products placed **on promotion**.

Across four stores and seven days, a complete weekly plan contains **84 decision variables**.

The optimization respects operational rules such as employee capacity, salary costs, store fixed costs, promotion limits and customer demand forecasts.

### Optimization objectives

Three management objectives were considered:

**O1 — Profit maximization**  
Maximize total weekly profit across all stores.

**O2 — Profit maximization with sales constraint**  
Maximize profit while enforcing a global maximum of **10,000 sold units**.

**O3 — Multi-objective optimization**  
Maximize the O2 profit objective while simultaneously minimizing total human resources.

### Algorithms evaluated

The optimization experiments included:

- Monte Carlo search;
- Hill Climbing;
- Simulated Annealing;
- Genetic Algorithm;
- NSGA-II;
- MOEA/D.

For the multi-objective problem, a weighted Simulated Annealing approach was also evaluated.

### Main optimization results

Results were evaluated over the same 20-week growing-window period used by the forecasting pipeline.

| Objective | Best approach | Main result |
|---|---|---|
| O1 | Genetic Algorithm | Highest median weekly profit: **$15,525** |
| O2 | Simulated Annealing | Highest median profit: **$1,306**, with average sales of 9,988 units |
| O3 | NSGA-II | Best overall Pareto quality, lowest median HR usage and highest hypervolume |

The experiments show that no single optimization algorithm dominates every management objective: algorithm suitability depends on the structure and constraints of the decision problem.

## 4. Interactive Decision Support System

The final system is implemented in **R Shiny** and integrates the complete analytical workflow into a single interface.

The application provides:

- interactive EDA visualizations;
- store- and week-specific demand forecasts;
- daily prediction tables and forecasting metrics;
- optimized weekly plans for O1 and O2;
- optimization convergence plots;
- Pareto-front exploration for O3;
- interactive comparison of NSGA-II and Simulated Annealing multi-objective solutions;
- daily staffing, promotion and KPI details for selected solutions.

A recorded demonstration is available here: [YouTube demo](https://youtu.be/K3_i2KM84NQ).

## Repository Structure

```text
.
├── EDA/                         # Store-level and global exploratory analysis
├── baltimore_univariados/      # Univariate forecasting experiments
├── lancaster_univariados/
├── philadelphia_univariados/
├── richmond_univariados/
├── multivariados/              # ARIMAX, VAR, MLPE and Random Forest models
├── otimizacao/                 # Optimization algorithms, outputs and Shiny DSS
│   ├── app.R
│   ├── ga.R
│   ├── hill.R
│   ├── sann.R
│   ├── nsga2.R
│   ├── MOAED.R
│   └── ...
├── utils.R                     # Shared forecasting utilities
├── multi-utils.R               # Shared multivariate utilities
├── Report.pdf                  # Full project report
├── project.pdf                 # Original project specification
└── guia-projeto.pdf            # Project guidelines
```

## Technologies

The project was developed entirely in **R**. Main packages include:

- `forecast`
- `rminer`
- `vars`
- `GA`
- `mco`
- `MOEADr`
- `shiny`
- `shinythemes`
- `ggplot2`
- `dplyr`
- `tidyverse`
- `lubridate`
- `corrplot`
- `DataExplorer`

## Running the Shiny Application

### 1. Install R dependencies

From an R session:

```r
install.packages(c(
  "DataExplorer",
  "GA",
  "MOEADr",
  "corrplot",
  "dplyr",
  "forecast",
  "ggplot2",
  "lubridate",
  "mco",
  "rminer",
  "shiny",
  "shinythemes",
  "tidyverse",
  "vars"
))
```

### 2. Launch the dashboard

From the repository root:

```bash
cd otimizacao
R -e 'shiny::runApp(".")'
```

The application uses the forecast and optimization artifacts already stored inside the `otimizacao/` directory.

## My Contributions

My individual contributions to the team project included:

- exploratory analysis of the **Philadelphia** store;
- development and evaluation of the **AutoARIMA** univariate model;
- implementation and analysis of the **ARIMAX** multivariate forecasting scenarios;
- initial implementation of **Hill Climbing** for O1 and O2;
- experimentation with **Simulated Annealing, Genetic Algorithms, NSGA-II and MOEA/D**;
- adaptation of optimization experiments to the 20-week growing-window evaluation;
- integration of forecasting and optimization results into the **Shiny decision support system**.

## Team

Developed by:

- Beatriz Peixoto
- Diogo Miranda
- Martim Félix
- Sandra Cerqueira

University of Minho — MSc in Data Engineering and Data Science, 2025/2026.

## Documentation

For the complete methodology, mathematical formulation, experiments and detailed results, see [`Report.pdf`](./Report.pdf).
