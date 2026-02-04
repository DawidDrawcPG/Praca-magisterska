 
# 0. Ustawienia i pakiety
 
# ŹRÓDŁA DANYCH:
# - Spółki GPW (KGH, PKO, PKN, CDR, LPP, ALE): Yahoo Finance via quantmod
# - WIG20: wig20_d.csv (dane historyczne)
# - Stopa ref NBP: stopy_proc.csv
# - Inflacja: inflacja.csv (GUS)
# - USD/PLN: usdpln.csv (NBP)
# - Ceny miedzi: miedz.csv (LME/COMEX)
# 
# OKRES: 2020-01-01 do dzisiaj
# METODOLOGIA: Markowitz portfolio optimization + rolling backtest

setwd("C:/Users/daza4/OneDrive/Desktop/PRACA MAGISTERSKA")

options(repos = c(CRAN = "https://cran.r-project.org"))

# Uruchom RAZ, potem zakomentuj:
# install.packages("quantmod")
# install.packages("PerformanceAnalytics")
# install.packages("quadprog")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("zoo")
# install.packages("corpcor")

rm(list = ls())

library(quantmod)
library(PerformanceAnalytics)
library(quadprog)
library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)
library(corpcor)
 
# 1. Pobranie danych z Yahoo + WIG20 z CSV
 

tickers_assets <- c("KGH.WA", "PKO.WA", "PKN.WA",
                    "CDR.WA", "LPP.WA", "ALE.WA")

start_date <- as.Date("2020-01-01")   # ~5 lat danych (okres Allegro)
end_date   <- Sys.Date()

# 1.1. Pobieramy spółki z Yahoo
getSymbols(tickers_assets,
           src  = "yahoo",
           from = start_date,
           to   = end_date,
           auto.assign = TRUE)

# 1.2. Ceny zamknięcia spółek
prices_assets <- do.call(
  merge,
  lapply(tickers_assets, function(sym) Cl(get(sym)))
)
colnames(prices_assets) <- gsub(".WA", "", tickers_assets)

# 1.3. Wczytanie WIG20 z CSV (plik wig20_d.csv musi być w katalogu roboczym)
# Wczytywanie z separatorem przecinek:
wig20_data <- read.csv(
  "wig20_d.csv",
  sep = ",",
  dec = ".",
  stringsAsFactors = FALSE
)

if (!all(c("Data", "Zamkniecie") %in% names(wig20_data))) {
  stop("Plik wig20_d.csv musi mieć kolumny 'Data' i 'Zamkniecie'.")
}

# Konwersja dat z automatycznym fallbackiem na format dzien-miesiac-rok
wig20_data$Data <- as.Date(wig20_data$Data)
if (any(is.na(wig20_data$Data))) {
  wig20_data$Data <- as.Date(wig20_data$Data, format = "%d.%m.%Y")
}

# Upewniamy się, że ceny są numeryczne
wig20_data$Zamkniecie <- suppressWarnings(as.numeric(wig20_data$Zamkniecie))
if (any(is.na(wig20_data$Zamkniecie))) {
  stop("Kolumna 'Zamkniecie' w wig20_d.csv zawiera nienumeryczne wartości.")
}

prices_bench <- xts(wig20_data$Zamkniecie, order.by = wig20_data$Data)
colnames(prices_bench) <- "WIG20"

# 1.4. Wspólny zakres dat i usunięcie NA
all_prices <- na.omit(merge(prices_assets, prices_bench))
prices_assets <- all_prices[, colnames(prices_assets)]
prices_bench  <- all_prices[, "WIG20", drop = FALSE]

 
# 2. Stopy zwrotu (logarytmiczne)
 

rets_assets <- na.omit(Return.calculate(prices_assets, method = "log"))
rets_bench  <- na.omit(Return.calculate(prices_bench,  method = "log"))
rets_bench  <- rets_bench[index(rets_assets), , drop = FALSE]

 
# 3. Funkcje do portfela Markowitza
 

# 3.1. Portfel maksymalnego Sharpe (long-short, sum(w)=1, w dowolne)
tangency_long_short <- function(mu, Sigma, rf = 0) {
  mu_excess <- mu - rf
  w <- solve(Sigma, mu_excess)
  w <- as.numeric(w) / sum(w)
  names(w) <- colnames(Sigma)
  return(w)
}

# 3.2. Portfel minimalnej wariancji (long-only, w >= 0, sum(w)=1)
gmv_long_only <- function(Sigma) {
  n    <- ncol(Sigma)
  Dmat <- 2 * as.matrix(Sigma)
  dvec <- rep(0, n)
  Amat <- cbind(rep(1, n), diag(n))   # sum w=1 oraz w >= 0
  bvec <- c(1, rep(0, n))
  meq  <- 1
  
  sol  <- solve.QP(Dmat, dvec, Amat, bvec, meq)
  w    <- sol$solution
  w    <- as.numeric(w) / sum(w)
  names(w) <- colnames(Sigma)
  return(w)
}

 
# 4. Portfele statyczne, miary efektywności + wykresy
 

rf_annual <- 0.02                  # roczna stopa wolna od ryzyka
rf_daily  <- rf_annual / 252

mu_hat <- colMeans(rets_assets)
Sigma  <- cov(rets_assets)

# Estymator macierzy kowariancji typu shrinkage (Ledoit–Wolff)
Sigma_shrink <- cov.shrink(rets_assets)


# 4.1. Wagi portfela long-short (tangency)
w_ls <- tangency_long_short(mu_hat, Sigma, rf = rf_daily)

# 4.2. Wagi portfela long-only (minimum variance)
w_lo <- gmv_long_only(Sigma)

# 4.3. Stopy zwrotu portfeli na całym okresie
ret_port_ls <- Return.portfolio(rets_assets, weights = w_ls, rebalance_on = NA)
ret_port_lo <- Return.portfolio(rets_assets, weights = w_lo, rebalance_on = NA)
colnames(ret_port_ls) <- "Portfel_LS"
colnames(ret_port_lo) <- "Portfel_LO"

# 4.4. Łączymy z benchmarkiem
all_static <- merge(ret_port_ls, ret_port_lo, rets_bench)
colnames(all_static)[3] <- "WIG20"

# 4.5. Tabelki z miarami efektywności
print(table.AnnualizedReturns(all_static, Rf = rf_daily))
print(SharpeRatio.annualized(all_static, Rf = rf_daily))
print(SortinoRatio(all_static, MAR = rf_daily))
print(maxDrawdown(all_static))

# 4.6. Klasyczny wykres PerformanceAnalytics (3 panele)
charts.PerformanceSummary(all_static,
                          main = "Portfel LS, LO vs WIG20 (statyczne wagi)")

# 4.7. Ładny wykres skumulowanych zwrotów w ggplot2
cumrets_static <- cumprod(1 + all_static) - 1

cumrets_static_df <- data.frame(
  Date = index(cumrets_static),
  coredata(cumrets_static)
) %>%
  pivot_longer(-Date, names_to = "Seria", values_to = "Zwrot")

p1 <- ggplot(cumrets_static_df, aes(x = Date, y = Zwrot, colour = Seria)) +
  geom_line(linewidth = 0.8) +
  labs(
    title  = "Skumulowane stopy zwrotu: Portfel LS, Portfel LO i WIG20",
    x      = "Data",
    y      = "Skumulowany zwrot",
    colour = ""
  ) +
  theme_minimal()

print(p1)
ggsave("wykres_statyczny.png", p1, width = 10, height = 6, dpi = 300)

 
# 5. Rolling backtest (okno 1-roczne, rebalancing miesięczny)
 

window_size <- 252                 # 252 dni ≈ 1 rok
ep <- endpoints(rets_assets, on = "months")
ep <- ep[ep > window_size]

if (length(ep) < 2) {
  stop("Za mało miesięcznych punktów do rolling backtestu (potrzeba > 1 roku danych).")
}

port_ls_roll <- NULL
port_lo_roll <- NULL

for (i in 2:length(ep)) {
  # indeksy okna treningowego
  train_idx_start <- ep[i-1] - window_size + 1
  train_idx_end   <- ep[i-1]
  train_window    <- rets_assets[train_idx_start:train_idx_end, ]
  
  mu_train <- colMeans(train_window)
  Sigma_tr <- cov(train_window)
  
  # wagi wyznaczone na podstawie okna treningowego
  w_ls_i <- tangency_long_short(mu_train, Sigma_tr, rf = rf_daily)
  w_lo_i <- gmv_long_only(Sigma_tr)
  
  # okres testowy = kolejny miesiąc
  test_window <- rets_assets[(ep[i-1] + 1):ep[i], ]
  
  ret_ls_i <- Return.portfolio(test_window, weights = w_ls_i, rebalance_on = NA)
  ret_lo_i <- Return.portfolio(test_window, weights = w_lo_i, rebalance_on = NA)
  
  port_ls_roll <- rbind(port_ls_roll, ret_ls_i)
  port_lo_roll <- rbind(port_lo_roll, ret_lo_i)
}

colnames(port_ls_roll) <- "Portfel_LS_roll"
colnames(port_lo_roll) <- "Portfel_LO_roll"

rets_bench_roll <- rets_bench[index(port_ls_roll), , drop = FALSE]

all_roll <- merge(port_ls_roll, port_lo_roll, rets_bench_roll)
colnames(all_roll)[3] <- "WIG20"

 
# 6. Miary efektywności dla rolling-backtestu + wykresy
 

print(table.AnnualizedReturns(all_roll, Rf = rf_daily))
print(SharpeRatio.annualized(all_roll, Rf = rf_daily))
print(SortinoRatio(all_roll, MAR = rf_daily))
print(maxDrawdown(all_roll))

# 6.1. Klasyczny wykres PerformanceAnalytics
charts.PerformanceSummary(all_roll,
                          main = "Rolling backtest: LS & LO vs WIG20")

# 6.2. ggplot2 – skumulowane zwroty rolling (dla porównania dynamiki)
cumrets_roll <- cumprod(1 + all_roll) - 1

cumrets_roll_df <- data.frame(
  Date = index(cumrets_roll),
  coredata(cumrets_roll)
) %>%
  pivot_longer(-Date, names_to = "Seria", values_to = "Zwrot")

p2 <- ggplot(cumrets_roll_df, aes(x = Date, y = Zwrot, colour = Seria)) +
  geom_line(linewidth = 0.8) +
  labs(
    title  = "Rolling backtest: skumulowane zwroty portfeli LS, LO i WIG20",
    x      = "Data",
    y      = "Skumulowany zwrot",
    colour = ""
  ) +
  theme_minimal()

print(p2)
ggsave("wykres_rolling.png", p2, width = 10, height = 6, dpi = 300)

 
# 7. Tabele: wagi i miary efektywności
 

# Kontrolnie sprawdźmy, że wagi istnieją:
print("Wagi portfeli statycznych - wektor LS:")
print(w_ls)
print("Wagi portfeli statycznych - wektor LO:")
print(w_lo)

 
# 7.1. Tabela wag portfeli statycznych
 
 
# Tabela: wagi portfeli statycznych (udzialy spolek)
 

tabela_wagi <- data.frame(
  Spolka          = names(w_ls),
  Waga_Portfel_LS = as.numeric(w_ls),
  Waga_Portfel_LO = as.numeric(w_lo)
)

# zaokrąglamy tylko kolumny liczbowe
tabela_wagi <- tabela_wagi %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

print("=== Tabela wag portfeli statycznych (LS i LO) ===")
print(tabela_wagi)


# Zapis do pliku CSV (do wklejenia w Word/Excel)
write.csv(tabela_wagi,
          "tabela_wagi_portfeli.csv",
          row.names = FALSE)

 
# 7.2. Tabela 3.1 – portfele statyczne (all_static)
 

ann_static    <- table.AnnualizedReturns(all_static, Rf = rf_daily)
sharpe_static <- SharpeRatio.annualized(all_static, Rf = rf_daily)
sortino_static <- SortinoRatio(all_static, MAR = rf_daily)
mdd_static    <- maxDrawdown(all_static)

tabela_3_1 <- data.frame(
  Strategia    = colnames(all_static),
  Ann_Return   = as.numeric(ann_static["Annualized Return", ]),
  Ann_StdDev   = as.numeric(ann_static["Annualized Std Dev", ]),
  Ann_Sharpe   = as.numeric(sharpe_static),
  Sortino      = as.numeric(sortino_static),
  Max_Drawdown = as.numeric(mdd_static)
)

tabela_3_1 <- tabela_3_1 %>% 
  mutate(across(where(is.numeric), ~ round(.x, 4)))
print("=== Tabela 3.1 – miary dla portfeli statycznych ===")
print(tabela_3_1)

write.csv(tabela_3_1,
          "tabela_3_1_statyczne.csv",
          row.names = FALSE)

 
# 7.3. Tabela 3.2 – rolling backtest (all_roll)
 

ann_roll     <- table.AnnualizedReturns(all_roll, Rf = rf_daily)
sharpe_roll  <- SharpeRatio.annualized(all_roll, Rf = rf_daily)
sortino_roll <- SortinoRatio(all_roll, MAR = rf_daily)
mdd_roll     <- maxDrawdown(all_roll)

tabela_3_2 <- data.frame(
  Strategia    = colnames(all_roll),
  Ann_Return   = as.numeric(ann_roll["Annualized Return", ]),
  Ann_StdDev   = as.numeric(ann_roll["Annualized Std Dev", ]),
  Ann_Sharpe   = as.numeric(sharpe_roll),
  Sortino      = as.numeric(sortino_roll),
  Max_Drawdown = as.numeric(mdd_roll)
)

tabela_3_2 <- tabela_3_2 %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
print("=== Tabela 3.2 – miary dla rolling backtestu ===")
print(tabela_3_2)

write.csv(tabela_3_2,
          "tabela_3_2_rolling.csv",
          row.names = FALSE)

 
# 8. Analiza makroekonomiczna portfeli
 

# 8.0. Pomocnicza funkcja do wczytywania „dziwnych” CSV z NBP
# (linie są w całości w cudzysłowach, np. "Data;Stopa")
clean_macro_csv <- function(path, col_names) {
  raw <- readLines(path, warn = FALSE)
  # usuwamy wszystkie cudzysłowy z linii
  raw <- gsub('"', "", raw, fixed = TRUE)
  con  <- textConnection(raw)
  on.exit(close(con))
  
  df <- read.csv(
    con,
    sep = ";",
    dec = ".",
    stringsAsFactors = FALSE
  )
  
  if (!missing(col_names)) {
    colnames(df) <- col_names
  }
  
  return(df)
}

 
# 8.1. Wczytanie danych makro (CSV w katalogu roboczym)
 

macro_rates <- clean_macro_csv("stopy_proc.csv", c("Data", "Stopa"))
macro_infl  <- clean_macro_csv("inflacja.csv",   c("Data", "Inflacja"))
macro_usd   <- clean_macro_csv("usdpln.csv",     c("Data", "Kurs"))
macro_cu    <- clean_macro_csv("miedz.csv",      c("Data", "Cena"))

macro_rates$Data <- as.Date(macro_rates$Data)
macro_infl$Data  <- as.Date(macro_infl$Data)
macro_usd$Data   <- as.Date(macro_usd$Data)
macro_cu$Data    <- as.Date(macro_cu$Data)

# Debug: sprawdź zakresy dat
cat("\n=== Zakres dat rolling: ", as.character(min(index(all_roll))), " do ", as.character(max(index(all_roll))), " ===")
cat("\n=== Zakres dat makro_rates: ", as.character(min(macro_rates$Data)), " do ", as.character(max(macro_rates$Data)), " ===")
cat("=== Zakres dat makro_usd: ", as.character(min(macro_usd$Data)), " do ", as.character(max(macro_usd$Data)), " ===\n")

# Tworzenie obiektów xts z danych makro
rates_xts <- xts(cbind(Stopa_proc = macro_rates$Stopa),
                 order.by = macro_rates$Data)

infl_xts  <- xts(cbind(Inflacja   = macro_infl$Inflacja),
                 order.by = macro_infl$Data)

usd_xts   <- xts(cbind(USDPLN     = macro_usd$Kurs),
                 order.by = macro_usd$Data)

cu_xts    <- xts(cbind(Miedz      = macro_cu$Cena),
                 order.by = macro_cu$Data)

# Sprawdzenie nakładania zakresów
rolling_start <- min(index(all_roll))
rolling_end   <- max(index(all_roll))
macro_start   <- min(macro_rates$Data)
macro_end     <- max(macro_rates$Data)

cat("\n=== DEBUG ZAKRESÓW ===\n")
cat("rolling_start: ", as.character(rolling_start), " (class: ", class(rolling_start), ")\n")
cat("rolling_end: ", as.character(rolling_end), " (class: ", class(rolling_end), ")\n")
cat("macro_start: ", as.character(macro_start), " (class: ", class(macro_start), ")\n")
cat("macro_end: ", as.character(macro_end), " (class: ", class(macro_end), ")\n")
cat("macro_end < rolling_start? ", macro_end < rolling_start, "\n")
cat("macro_start > rolling_end? ", macro_start > rolling_end, "\n")

# Warunek: jeśli makro kończy się PRZED rozpoczęciem rolling - błąd krytyczny
if (macro_end < rolling_start || macro_start > rolling_end) {
  cat("\n*** BŁĄD KRYTYCZNY: Dane makro (", as.character(macro_start), " - ", as.character(macro_end), 
      ") NIE nakładają się z rolling (", as.character(rolling_start), " - ", as.character(rolling_end), ") ***\n")
  cat("Pominięto całą analizę makro (sekcja 8).\n\n")
} else {
  cat("\n=== Zakresy się nakładają - kontynuuję analizę makro ===\n")
  
  # Jeśli makro kończy się wcześniej niż rolling, ostrzeżenie (ale kontynuujemy z ekstrapolacją)
  if (macro_end < rolling_end) {
    cat("\n*** UWAGA: Dane makro kończą się ", as.character(macro_end), 
        ", ale rolling sięga ", as.character(rolling_end), " ***\n")
    cat("Ostatnia wartość makro będzie ekstrapolowana (forward fill).\n\n")
  }
  
   
  # 8.2. Połączenie serii makro z portfelami rolling
   
  
  cat("\n=== Rozpoczynam dopasowanie danych makro do rolling portfolio ===\n")
  cat("Zakres dat all_roll: ", as.character(min(index(all_roll))), " do ", as.character(max(index(all_roll))), "\n")
  cat("Zakres dat rates_xts: ", as.character(min(index(rates_xts))), " do ", as.character(max(index(rates_xts))), "\n")
  
  # Merge wszystkich makro w jedno (miesięczne)
  macro_all <- merge(rates_xts, infl_xts, usd_xts, cu_xts)
  cat("Macro_all po merge: ", nrow(macro_all), " wierszy (miesięczne)\n")
  
  # Tworzymy dzienną siatkę dat z all_roll
  date_grid <- xts(order.by = index(all_roll))
  cat("Date_grid (dni rolling): ", nrow(date_grid), " wierszy\n")
  
  # Merge makro na dzienną siatkę rolling
  macro_on_grid <- merge(date_grid, macro_all)
  cat("Po merge z date_grid: ", nrow(macro_on_grid), " wierszy\n")
  cat("NA przed na.locf: ", sum(is.na(macro_on_grid)), "\n")
  
  # Forward fill: każda miesięczna wartość rozciąga się na kolejne dni
  macro_on_grid <- na.locf(macro_on_grid, na.rm = FALSE)
  cat("NA po forward fill: ", sum(is.na(macro_on_grid)), "\n")
  
  # Backward fill dla początku (jeśli makro zaczyna się później)
  macro_on_grid <- na.locf(macro_on_grid, fromLast = TRUE, na.rm = FALSE)
  cat("NA po backward fill: ", sum(is.na(macro_on_grid)), "\n")
  
  # Teraz wyciągamy tylko kolumny makro (bez pustego date_grid)
  macro_aligned <- macro_on_grid[, colnames(macro_all)]
  
  cat("Kolumny macro_aligned: ", paste(colnames(macro_aligned), collapse = ", "), "\n")
  cat("Kolumny all_roll: ", paste(colnames(all_roll), collapse = ", "), "\n")
  
  # Merge z portfelami rolling
  dane_makro_roll <- merge(all_roll, macro_aligned, all = FALSE)
  
  cat("Po merge all_roll + macro: ", nrow(dane_makro_roll), " wierszy, ", ncol(dane_makro_roll), " kolumn\n")
  cat("Kolumny dane_makro_roll: ", paste(colnames(dane_makro_roll), collapse = ", "), "\n")
  
  # Sprawdź czy są NA
  if (nrow(dane_makro_roll) > 0) {
    na_count <- sum(is.na(dane_makro_roll))
    cat("Liczba NA w dane_makro_roll: ", na_count, "\n")
    
    if (na_count > 0) {
      # Usuwamy wiersze z NA
      dane_makro_roll <- dane_makro_roll[complete.cases(dane_makro_roll), ]
      cat("Po complete.cases: ", nrow(dane_makro_roll), " wierszy\n")
    }
  }
  
  cat("FINALNIE dane_makro_roll: ", nrow(dane_makro_roll), " wierszy\n")
  
  if (nrow(dane_makro_roll) > 0) {
    cat("Zakres dat: ", as.character(min(index(dane_makro_roll))), " do ", as.character(max(index(dane_makro_roll))), "\n")
    cat("Pierwsze 3 wiersze danych makro:\n")
    print(head(dane_makro_roll, 3))
  }
  
   
  # 8.3. Warunek: czy mamy wystarczająco dużo obserwacji?
   
  
  if (nrow(dane_makro_roll) < 30) {
    cat("\n*** UWAGA: Za mało danych (< 30 obserwacji)! Pominięto analizę makro. ***\n")
    cat("Możliwe przyczyny:\n")
    cat("1. Dane makro nie pokrywają okresu rolling (sprawdź daty powyżej)\n")
    cat("2. Pliki CSV makro mają błędny format lub brakujące wartości\n\n")
    cat("Kontynuuję bez analizy makro...\n\n")
    
  } else {
    
     
    # 8.4. Korelacje i proste regresje makro
     
    
    print("=== Korelacje portfeli rolling i WIG20 z czynnikami makro ===")
    print(cor(dane_makro_roll))
    
    # Przygotowanie ramki danych
    dane_df <- data.frame(
      Date = index(dane_makro_roll),
      coredata(dane_makro_roll)
    )
    
    # Modele liniowe
    model_ls  <- lm(Portfel_LS_roll ~ Stopa_proc + Inflacja + USDPLN + Miedz,
                    data = dane_df)
    model_lo  <- lm(Portfel_LO_roll ~ Stopa_proc + Inflacja + USDPLN + Miedz,
                    data = dane_df)
    model_wig <- lm(WIG20           ~ Stopa_proc + Inflacja + USDPLN + Miedz,
                    data = dane_df)
    
    print("=== Regresja: Portfel LS roll vs czynniki makro ===")
    print(summary(model_ls))
    
    print("=== Regresja: Portfel LO roll vs czynniki makro ===")
    print(summary(model_lo))
    
    print("=== Regresja: WIG20 vs czynniki makro ===")
    print(summary(model_wig))
    
    
    # 8.5. Zapis współczynników regresji do CSV
     
    
    coef_ls  <- as.data.frame(summary(model_ls)$coefficients)
    coef_lo  <- as.data.frame(summary(model_lo)$coefficients)
    coef_wig <- as.data.frame(summary(model_wig)$coefficients)
    
    coef_ls$Zmienna  <- rownames(coef_ls)
    coef_lo$Zmienna  <- rownames(coef_lo)
    coef_wig$Zmienna <- rownames(coef_wig)
    
    rownames(coef_ls)  <- NULL
    rownames(coef_lo)  <- NULL
    rownames(coef_wig) <- NULL
    
    write.csv(coef_ls,
              "tabela_makro_portfel_LS_roll.csv",
              row.names = FALSE)
    write.csv(coef_lo,
              "tabela_makro_portfel_LO_roll.csv",
              row.names = FALSE)
    write.csv(coef_wig,
              "tabela_makro_WIG20.csv",
              row.names = FALSE)
    
    print("=== Zapisano tabele z wynikami regresji makro ===")
  }  # koniec else (gdy >= 30 obserwacji)
}  # koniec else (gdy zakresy się nakładają)

cat("\n=== SKRYPT ZAKOŃCZONY POMYŚLNIE ===\n")
