# Installiere limSolve, falls nicht vorhanden
if (!requireNamespace("limSolve", quietly = TRUE)) {
  install.packages("limSolve")
}
library(limSolve)

# 1. Beispieldaten laden (Hier als Dataframe generiert, kann durch read.csv() ersetzt werden)
futter_daten <- data.frame(
  Zutat = c(
    "trout", "salmon dry", "trout dry", "sweet potato", "potato", "salmon oil", "salmon bouillon", 
    "minerals", "Fennel", "Chickpeas", "Pumpkin", "Parsnips", "Seaweed", 
    "Jerusalem artichoke", "Thyme", "Marjoram", "Oregano", "Parsley", 
    "Sage", "Tomatoes", "Nettle", "Hawthorn", "Dandelion", "Ginseng", 
    "MOS", "FOS", "Spirulina", "Yucca"
  ),
  Prozent = c(
    25, 6.5, 6.5, 32, NA, NA, NA, 
    NA, NA, NA, NA, NA, 
    0.3, NA, NA, NA, NA, NA, 
    0.1, NA, NA, NA, 0.04, 
    NA, NA, NA, NA, NA
  ),
  stringsAsFactors = FALSE
)

# Heuristik: Mineralstoffe zwingend auf 2.0% setzen
futter_daten$Prozent[futter_daten$Zutat == "minerals"] <- 2.0


# 2. Die neue, dynamische Sampling-Funktion
simulate_lca_proxies_split <- function(df, start_monotony_var, n_iter = 10000) {
  n <- nrow(df)
  
  # A. Gleichungen (E * x = F): Feste Werte und Summe = 100
  known_idx <- which(!is.na(df$Prozent))
  n_known <- length(known_idx)
  
  E <- matrix(0, nrow = 1 + n_known, ncol = n)
  F_vec <- numeric(1 + n_known)
  
  E[1, ] <- 1 # Summe = 100
  F_vec[1] <- 100
  
  for(i in seq_along(known_idx)) {
    idx <- known_idx[i]
    E[i + 1, idx] <- 1
    F_vec[i + 1] <- df$Prozent[idx]
  }
  
  # B. Ungleichungen (G * x >= H)
  # WICHTIG: Wir erzwingen die absteigende Reihenfolge erst ab der definierten Zutat!
  start_idx <- which(df$Zutat == start_monotony_var)
  if(length(start_idx) == 0) stop("Start-Zutat für Monotonie nicht gefunden!")
  
  n_monotony <- n - start_idx
  
  # Matrix G muss Positivität (für alle) + Monotonie (ab start_idx) fassen
  G <- matrix(0, nrow = n_monotony + n, ncol = n)
  H_vec <- numeric(n_monotony + n)
  
  row_counter <- 1
  
  # B1: Monotonie (x_i >= x_{i+1}) NUR für den Block ab 'sweet potato'
  for(i in start_idx:(n-1)) {
    G[row_counter, i] <- 1
    G[row_counter, i+1] <- -1
    row_counter <- row_counter + 1
  }
  
  # B2: Positivität (x_i >= 0) für ALLE Zutaten
  for(i in 1:n) {
    G[row_counter, i] <- 1
    # H_vec bleibt an dieser Stelle 0
    row_counter <- row_counter + 1
  }
  
  # C. Monte-Carlo Sampling (Hit-and-Run)
  message("Starte Sampling mit ", n_iter, " Iterationen ab Zutat: ", start_monotony_var)
  mc_samples <- xsample(E = E, F = F_vec, G = G, H = H_vec, iter = n_iter)
  
  samples_df <- as.data.frame(mc_samples$X)
  colnames(samples_df) <- df$Zutat
  
  return(samples_df)
}

# 3. Ausführen (Wir setzen 'sweet potato' als Ankerpunkt für die absteigende Liste)
set.seed(42)
samples_full <- simulate_lca_proxies_split(futter_daten, start_monotony_var = "sweet potato", n_iter = 5000)

# 4. Filtern & Aggregieren (< 1% vernachlässigen)
mean_vals <- colMeans(samples_full)
relevante_zutaten <- names(mean_vals[mean_vals >= 1.0])

samples_relevant <- samples_full[, relevante_zutaten]

lca_summary <- data.frame(
  Zutat = relevante_zutaten,
  Mittelwert = mean_vals[relevante_zutaten],
  Quantil_5 = apply(samples_relevant, 2, quantile, probs = 0.05),
  Quantil_95 = apply(samples_relevant, 2, quantile, probs = 0.95)
)
rownames(lca_summary) <- NULL

# Ergebnisse anzeigen
print("=== LCA PROXY ÜBERSICHT (Erwartungswerte >= 1%) ===")
print(lca_summary)

# Die Variable 'samples_relevant' enthält nun die N-Iterationen 
# für jede Zutat >= 1% und kann für weitere Varianzanalysen genutzt werden.
# head(samples_relevant)