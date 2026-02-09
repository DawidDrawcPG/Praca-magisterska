# Portfele Markowitza (LS i LO) vs WIG20 + rolling backtest (R)

Ten projekt liczy i porównuje:
- portfel **LS** (long-short) = portfel maksymalnego Sharpe (wagi mogą być ujemne),
- portfel **LO** (long-only) = portfel minimalnej wariancji (wagi >= 0),
- benchmark **WIG20** (z pliku CSV).

Skrypt:
1) pobiera dane spółek z Yahoo Finance (pakiet `quantmod`),
2) wczytuje WIG20 z pliku `wig20_d.csv`,
3) liczy dzienne **logarytmiczne stopy zwrotu**,
4) robi portfele statyczne (jedne wagi na cały okres),
5) robi **rolling backtest** (okno 252 dni, rebalancing miesięczny),
6) zapisuje wykresy i tabele do plików.

---

## Wymagania

- R (np. R 4.x)
- Pakiety:
  - `quantmod`
  - `PerformanceAnalytics`
  - `quadprog`
  - `ggplot2`
  - `dplyr`
  - `tidyr`
  - `zoo`

Jeśli nie masz pakietów, w skrypcie są linie `install.packages(...)` (uruchom raz i potem zakomentuj).

---

## Pliki wejściowe

W katalogu roboczym muszą być:
- `wig20_d.csv` – dane WIG20 z kolumnami:
  - `Data` (format daty, np. `YYYY-MM-DD`)
  - `Zamkniecie` (wartości zamknięcia indeksu)

Skrypt sam pobiera dane spółek z Yahoo dla tickerów:
- `KGH.WA`, `PKO.WA`, `PKN.WA`, `CDR.WA`, `LPP.WA`, `ALE.WA`

Zakres dat:
- start: `2020-01-01`
- koniec: `Sys.Date()` (dzisiejsza data)

---

## Jak uruchomić

1) Wgraj `wig20_d.csv` do folderu projektu.
2) Otwórz skrypt w RStudio.
3) Ustaw katalog roboczy (albo zmień w skrypcie `setwd(...)` na swój).
4) Uruchom cały skrypt.

Uwaga: pierwszy raz może trwać dłużej (pobieranie danych z Yahoo).

---

## Co skrypt liczy (w skrócie)

### Stopy zwrotu
- Liczone jako **logarytmiczne**: `Return.calculate(..., method="log")`

### Portfele statyczne
- **LS (long-short)**: wagi z wzoru na portfel styczny (maks. Sharpe), warunek `sum(w)=1`, wagi mogą być ujemne.
- **LO (long-only)**: minimalna wariancja, warunki `sum(w)=1` i `w >= 0` (QP przez `solve.QP`).

### Rolling backtest
- okno treningowe: `252` sesje (około 1 rok),
- rebalancing: co miesiąc (endpointy `months`),
- wagi liczone na oknie treningowym i stosowane na kolejny miesiąc testowy.

---

## Pliki wyjściowe (co zapisuje)

Wykresy:
- `wykres_statyczny.png` – skumulowane zwroty portfeli statycznych (LS, LO, WIG20)
- `wykres_rolling.png` – skumulowane zwroty w rolling backteście

Tabele CSV:
- `tabela_wagi_portfeli.csv` – wagi portfeli statycznych (LS i LO)
- `tabela_3_1_statyczne.csv` – miary efektywności dla portfeli statycznych
- `tabela_3_2_rolling.csv` – miary efektywności dla rolling backtestu

Miary (wypisywane też w konsoli):
- annualized return
- annualized std dev
- Sharpe (annualized)
- Sortino
- max drawdown

---

## Typowe problemy

1) **Błąd wczytania `wig20_d.csv`**
- sprawdź czy plik jest w katalogu roboczym
- sprawdź nazwy kolumn: `Data` i `Zamkniecie`

2) **Braki danych / NA**
- skrypt robi `na.omit(merge(...))`, więc utnie dni bez wspólnych danych
- jeśli dla jakiejś spółki brakuje notowań w części okresu, zakres może się skrócić

3) **Yahoo czasem blokuje / zwraca błędy**
- spróbuj uruchomić ponownie
- ewentualnie zmień zakres dat na krótszy

---

## Autor / opis do pracy
Skrypt jest przygotowany do użycia w pracy magisterskiej: porównanie strategii Markowitza (LS i LO) z benchmarkiem WIG20 na danych dziennych.

## Dawid Krzysztof Drawc
