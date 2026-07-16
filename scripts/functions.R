# Installiere limSolve, falls nicht vorhanden
if (!requireNamespace("limSolve", quietly = TRUE)) {
  install.packages("limSolve")
}
library(limSolve)
library(dplyr)
library(stringr)

# functions for data summary statistics

get_co2_value <- function(product, co2_data, dry_factor = 2.5, protein_factor = 3.5, printing = FALSE) {
  # search for co2 value of product in co2_data
  # use criterion if it applies
  
  criterion <- "raw"
  product_col <- "LCI.Name"
  cat_col <- "Sous.groupe.d.aliment"
  
  # scan for "dry" and "protein" first
  is_dry <- str_detect(product, regex("\\bdry\\b", ignore_case = TRUE))
  is_protein <- str_detect(product, regex("\\bprotein\\b", ignore_case = TRUE))
  is_meat <- str_detect(product, regex("\\bmeat\\b", ignore_case = TRUE))
  is_offal <- str_detect(product, regex("\\boffal\\b", ignore_case = TRUE))
  
  # clean up product name
  if (is_dry) {
    product <- str_remove(product, regex("\\bdry\\b", ignore_case = TRUE)) %>% 
      str_trim()
  }
  if (is_protein) {
    product <- str_remove(product, regex("\\bprotein\\b", ignore_case = TRUE)) %>% 
      str_trim()
  }
  if (is_meat) {
    product <- str_remove(product, regex("\\bmeat\\b", ignore_case = TRUE)) %>% 
      str_trim()
  }
  if (is_offal) {
    product <- str_remove(product, regex("\\boffal\\b", ignore_case = TRUE)) %>% 
      str_trim()
  }
  
  # scan for "bouillon"
  is_bouillon <- str_detect(product, regex("\\bbouillon\\b", ignore_case = TRUE))
  if (is_bouillon) {
    product = "bouillon poultry"
  }
  
  key_words <- str_split(product, "\\s+")[[1]] 

  # sequenzielle Suche
  products_selected <- co2_data
  
  for (word in key_words) {
    exact_word <- paste0("\\b", word)
    
    products_selected <- products_selected %>%
      filter(str_detect(products_selected[[product_col]], regex(exact_word, ignore_case = TRUE)))
  }
  
  # For Meat: exclude offals
  # 1. Deine Ausschlussliste
  exclude_offal <- c("liver", "heart", "kidney", "gizzard", "tripe")
  offal_pattern <- paste(exclude_offal, collapse = "|")
  if (is_meat){
    
    # Filtern des Dataframes (Das ! negiert grepl, sodass nur Zeilen OHNE diese Begriffe bleiben)
    products_selected <- products_selected[
      !grepl(offal_pattern, products_selected[[product_col]], ignore.case = TRUE), 
    ]
  } 
  if (is_offal){
    products_selected <- products_selected[
      grepl(offal_pattern, products_selected[[product_col]], ignore.case = TRUE),
    ]
  }

  # Sicherheitsabfrage
  if (nrow(products_selected) == 0) {
    warning(paste("Das Produkt", product, "wurde nicht gefunden."))
    return(NA) # Gibt "Not Available" zurück, falls es keinen Treffer gibt
  }
  
  # check whether there is "raw" in product name and use this
  products_selected$has_criterion <- grepl(criterion, products_selected[[product_col]], ignore.case = TRUE)
  if (sum(products_selected$has_criterion == TRUE) >= 1) {
    ind <- products_selected$has_criterion == TRUE
    products_selected <- products_selected[ind,]
}
  
  # exclude category "charcuteries"
  products_selected$has_criterion <- grepl("charcuteries", products_selected[[cat_col]], ignore.case = TRUE)
  if (sum(products_selected$has_criterion == TRUE) >= 1) {
    ind <- products_selected$has_criterion == FALSE
    products_selected <- products_selected[ind,]
  }
  
  # exclude category "autres produits à base de viande"
  products_selected$has_criterion <- grepl("autres produits à base de viande", products_selected[[cat_col]], ignore.case = TRUE)
  if (sum(products_selected$has_criterion == TRUE) >= 1) {
    ind <- products_selected$has_criterion == FALSE
    products_selected <- products_selected[ind,]
  }
  
  counts = as.data.frame(table(products_selected$Sous.groupe.d.aliment))
  maxcat = which.max(counts[,2])
  maxcat = counts[maxcat,1]
  products_selected <- products_selected %>%
    filter(.data[[cat_col]] == maxcat)
  

  # CO2-Wert
  co2_col <- products_selected[,c(13,14)] # CO2 values from agricultur and transformation!
  
  # if dry product: multiply with fixed factor:
  if (is_dry) {
    co2_col <- co2_col * dry_factor
  }
  if (is_protein) {
    co2_col <- co2_col * dry_factor
  }
  
  co2_values <- rowSums(co2_col)
  co2_value <- median(co2_values)
  co2_min <- min(co2_values)
  co2_max <- max(co2_values)
  if (printing==TRUE) {
    print(product)
    print(products_selected[[product_col]])
  }
  output <- c(co2_value,co2_min,co2_max)
  
  return(output)
}



# Monte-Carlo-Sampling-Funktion
simulate_lca_proxies_split <- function(df, n_iter = 1000) {
  n <- nrow(df)
  n_vars <- n + 1 # Zuweisung der Slack-Variable für Feuchtigkeit/Asche
  
  known_idx <- which(!is.na(df$Prozent))
  n_known <- length(known_idx)
  
  # A. Gleichungen (E * x = F)
  E <- matrix(0, nrow = 1 + n_known, ncol = n_vars)
  F_vec <- numeric(1 + n_known)
  
  E[1, ] <- 1 # Summe (Zutaten + Slack) = 100
  F_vec[1] <- 100
  
  for(i in seq_along(known_idx)) {
    idx <- known_idx[i]
    E[i + 1, idx] <- 1
    F_vec[i + 1] <- df$Prozent[idx]
  }
  
  # Monotonie-Start
  first_na_idx <- which(is.na(df$Prozent))[1]
  start_idx <- ifelse(!is.na(first_na_idx) && first_na_idx > 1, first_na_idx - 1, 1)
  start_ingredient <- df$Zutat[start_idx]
  n_monotony <- n - start_idx
  
  # B. Ungleichungen (G * x >= H)
  G <- matrix(0, nrow = n_monotony + n + 2, ncol = n_vars)
  H_vec <- numeric(n_monotony + n + 2)
  
  row_counter <- 1
  
  # B1: Monotonie (x_i >= x_{i+1}) exklusive Slack
  if(n_monotony > 0) {
    for(i in start_idx:(n-1)) {
      G[row_counter, i] <- 1
      G[row_counter, i+1] <- -1
      row_counter <- row_counter + 1
    }
  }
  
  # B2: Standard-Positivitätsprüfung angewendet
  for(i in 1:n) {
    G[row_counter, i] <- 1
    row_counter <- row_counter + 1
  }
  
  # B3: Slack-Korridor definieren
  G[row_counter, n_vars] <- 1
  H_vec[row_counter] <- 0
  row_counter <- row_counter + 1
  
  G[row_counter, n_vars] <- -1
  H_vec[row_counter] <- -25
  
  # C. Monte-Carlo Sampling
  message("Starte Sampling mit ", n_iter, " Iterationen ab Zutat: ", start_ingredient)
  mc_samples <- limSolve::xsample(E = E, F = F_vec, G = G, H = H_vec, iter = n_iter)
  
  # Slack-Variable vor Datenübergabe isolieren und entfernen
  samples_df <- as.data.frame(mc_samples$X[, 1:n])
  colnames(samples_df) <- df$Zutat
  
  return(samples_df)
}

filter_ingredients <- function(samples_full) {
  
  mean_vals <- colMeans(samples_full)
  
  # Nur Zutaten auswählen, deren Mittelwert >= 1.0 (Prozent) ist
  relevante_zutaten <- names(mean_vals[mean_vals >= 1.0])
  samples_relevant <- samples_full[, relevante_zutaten, drop = FALSE]
  
  
  return(samples_relevant)
}


prepare_food <- function(data, product) {
  data %>%
    filter(Produkt == product) %>%
    select(Zutat, Prozent) # Die Spalte "Produkt" ausblenden
}