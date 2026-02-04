 
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

setwd("C:/Users/daza4/OneDrive/Desktop/PRACA MAGISTERSKA/ODP_ Szkic Pracy Magisterskiej Dawid Krzysztof Drawc s202286/Praca-magisterska")


# 0. Ustawienia i pakiety



options(repos = c(CRAN = "https://cran.r-project.org"))

# Uruchom RAZ, potem zakomentuj:
# install.packages("quantmod")
# install.packages("PerformanceAnalytics")
# install.packages("quadprog")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("zoo")

rm(list = ls())

library(quantmod)
library(PerformanceAnalytics)
library(quadprog)
library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)

# --- Foldery na wyniki (auto-tworzenie) ---
dir_png <- file.path(getwd(), "wyniki_png")
dir_csv <- file.path(getwd(), "wyniki_csv")

if (!dir.exists(dir_png)) dir.create(dir_png, recursive = TRUE)
if (!dir.exists(dir_csv)) dir.create(dir_csv, recursive = TRUE)


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
wig20_data <- read.csv("wig20_d.csv")
wig20_data$Data <- as.Date(wig20_data$Data)
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

# BONUS: zapis wykresu PerformanceSummary do PNG
png(file.path(dir_png, "performance_summary_statyczne.png"),
    width = 1600, height = 900, res = 150)
charts.PerformanceSummary(all_static,
                          main = "Portfel LS, LO vs WIG20 (statyczne wagi)")
dev.off()

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
ggsave(filename = file.path(dir_png, "wykres_statyczny.png"),
       plot = p1, width = 10, height = 6, dpi = 300)


# 5. Rolling backtest (okno 1-roczne, rebalancing miesięczny)


window_size <- 252                 # 252 dni ≈ 1 rok
ep <- endpoints(rets_assets, on = "months")
ep <- ep[ep > window_size]

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

# BONUS: zapis wykresu PerformanceSummary do PNG
png(file.path(dir_png, "performance_summary_rolling.png"),
    width = 1600, height = 900, res = 150)
charts.PerformanceSummary(all_roll,
                          main = "Rolling backtest: LS & LO vs WIG20")
dev.off()

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
ggsave(filename = file.path(dir_png, "wykres_rolling.png"),
       plot = p2, width = 10, height = 6, dpi = 300)


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
          file.path(dir_csv, "tabela_wagi_portfeli.csv"),
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
          file.path(dir_csv, "tabela_3_1_statyczne.csv"),
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
          file.path(dir_csv, "tabela_3_2_rolling.csv"),
          row.names = FALSE)

cat("\n=== SKRYPT ZAKOŃCZONY POMYŚLNIE ===\n")
