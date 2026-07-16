source("scripts/functions.R")
# CSV Datei mit CO2 Werten von IFEU
#datapath <- "data/IFEU_2020_CO2.csv"
datapath <- "data/agribalyse-31-detail-par-etape.csv"
co2_data <- read.csv(datapath, stringsAsFactors = FALSE)

Nsamples = 1000
# check
#head(co2_data)

#------------------------------------------------------------
##Beispiel: Suchbegriffe definieren
#product <- "salmon dry"
#
#co2_values = get_co2_value(product,co2_data)
#
## FOR NOW: use median co2 value of all categories
#co2_value = co2_values[1]



#-------------------------------------------------#
#-------------- LOOP THROUGH PRODUCTS ------------#
#-------------------------------------------------#

# fixed values for some ingredients not in data base
fixed_co2_values <- list(
  # Anorganische Stoffe (Bergbau/Verarbeitung)
  "minerals"              = 0.3, 
  "calcium carbonate"     = 0.3, 
  "monocalcium phosphate" = 0.3, 
  "sodium chloride"       = 0.3,
  
  # Hefe-Produkte (Schätzwert für Fermentation & Trocknung)
  "brewer's yeast dry"        = 1.5,
  "brewer's yeast"        = 1.5,
  "yeast extract" = 1.5,
  "yeast protein hydrolysate" = 1.5,
  "protein hydrolysate" = 1.5,
  "protein autolysate" = 1.5,
  "plant protein hydrolyzed" = 1.5,
  "yeast cell walls"          = 1.5,
  "yeast hydrolyzed"          = 1.5,
  "pea protein" = 2,
  "potato peeled protein" = 2,
  "potato peeled starch" = 2,
  "chicken protein hydrolysate" = 4.0,
  "herring protein hydrolysate" = 3.0,
  "lamb protein hydrolysate" = 4.0,
  "chicken protein hydrolyzed" = 4.0,
  "lignocellulose" = 0.3,
  "meat hydrolysate" = 4.0
)
  
#----------------------- PROXY MAPPING
proxy_mapping <- list(
  "sunflower protein" = "Soy protein dehydrated",
  "sunflower seeds" = "flaxseed",
  # Pferd wird zu Schwein
  "horse offal protein" = "pork offal protein",
  "flaxseed oil" = "rapeseed oil",
  "water buffalo offal protein" = "beef offal protein",
  "ostrich offal protein" = "turkey offal protein",
  "beef liver" = "beef heart",
  "whitefish dry" = "cod dry",
  "carob"          = "peas",
  "poultry meat" = "chicken meat",
  "poultry meat protein" = "chicken meat protein",
  "poultry meat dry" = "chicken meat dry",
  "poultry fat" = "chicken fat",
  "corn grits" = "corn",
  "psyllium husks" = "flaxseed",
  "lamb tripe" = "lamb liver",
  "psyllium" = "flaxseed",
  "parsnip dry" = "carrot dry",
  "cassava dry" = "potato peeled dry"
)

  # Futterdaten laden
  all_food <- read.csv("data/futtermittel_datenbank.csv", stringsAsFactors = FALSE, na.strings = "NA")
  all_names = unique(all_food$Produkt)
  
  set.seed(42)
  co2_results <- data.frame(
    Product = character(),
    CO2_med = numeric(),
    CO2_5 = numeric(),
    CO2_95 = numeric(),
    CO2_min = numeric(),
    CO2_max = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (current_product in all_names) {
    message("Produkt: ",current_product)
    # einzelnes Futter
    futter_daten <- prepare_food(all_food, current_product)
  
    # Heuristik: Mineralstoffe zwingend auf 2.0% setzen
    futter_daten$Prozent[futter_daten$Zutat == "minerals"] <- 3.0
  
    # MONTE-CARLO-SAMPLING
    simulation_required <- any(is.na(futter_daten$Prozent))
    if (simulation_required) {
      samples_full <- simulate_lca_proxies_split(futter_daten, n_iter = Nsamples)
    } else {
      statische_matrix <- matrix(
        rep(futter_daten$Prozent, times = Nsamples),
        nrow = Nsamples,
        ncol = length(futter_daten$Prozent),
        byrow = TRUE
      )
      
      # In Dataframe umwandeln (damit es sich exakt wie der Output deiner Funktion verhält)
      samples_full <- as.data.frame(statische_matrix)
      
      # Spaltennamen aus den Zutaten übernehmen
      colnames(samples_full) <- futter_daten$Zutat
    }
    ingredients_summary = filter_ingredients(samples_full)
  
    
  
  #-----------------------  GET CO2 VALUES FOR EXAMPLE DRY FOOD
  
    ingredients_co2 <- matrix(NA_real_, nrow = Nsamples, ncol = ncol(ingredients_summary))
    ingredients_co2_min <- matrix(NA_real_, nrow = Nsamples, ncol = ncol(ingredients_summary))
    ingredients_co2_max <- matrix(NA_real_, nrow = Nsamples, ncol = ncol(ingredients_summary))
  
  
  pr = TRUE
  for (i in 1:ncol(ingredients_summary)) {
    
    # produktname
    product <- colnames(ingredients_summary)[i]

    if (product %in% names(fixed_co2_values)) {
      
      co2_values <- c(fixed_co2_values[[product]],fixed_co2_values[[product]],fixed_co2_values[[product]])
      
    } else {
      
      if (product %in% names(proxy_mapping)) {
        search_product <- proxy_mapping[[product]]
      } else {
        search_product <- product
      }
      
      co2_values <- get_co2_value(search_product, co2_data, printing=pr) # this is the function using the AGRIBALYS database
    }
    if (!is.na(co2_values[1])) {
        ingredients_co2[,i] <- co2_values[1] * ingredients_summary[[product]]/100 # AVERAGE
        ingredients_co2_min[,i] <- co2_values[2] * ingredients_summary[[product]]/100 # MIN CO2
        ingredients_co2_max[,i] <- co2_values[3] * ingredients_summary[[product]]/100 # MAX CO2
      } else {
        stop("NaN value detected")
      }
    
  if (i>3) {
    pr = FALSE
  }
  }
  
  
  final_co2 <- rowSums(ingredients_co2)
  final_co2_min <- rowSums(ingredients_co2_min)
  final_co2_max <- rowSums(ingredients_co2_max)
    
  # Quantile von dieser finalen Verteilung und av. über mediane und min max
  mean_co2 <- mean(final_co2)
  co2_05  <- quantile(final_co2, probs = 0.05)
  co2_95  <- quantile(final_co2, probs = 0.95)
  min_co2 <- mean(final_co2_min)
  max_co2 <- mean(final_co2_max)
  
  # add the new data
  new_entry <- data.frame(
    Product = current_product,
    CO2_med = mean_co2, 
    CO2_05 = co2_05,
    CO2_95 = co2_95,
    CO2_min = min_co2,
    CO2_max = max_co2,
    stringsAsFactors = FALSE
  )
  
  co2_results <- rbind(co2_results, new_entry)
  }
  print(co2_results)
  
  
  
  #------------------------- PLOT
  
  # label vegan products
  vegane_produkte <- c(
    "Green Petfood Veggie", 
    "Greta Groß", 
    "Vutter Wie Rind", 
    "VegDog Green Crunch",
    "VegDog Farmers Crunch"
  )
  co2_results <- co2_results %>%
    mutate(Proteinquelle = if_else(Product %in% vegane_produkte, "Vegan", "Tierisch"))
  
  install.packages("ggplot2")
  library("ggplot2")
  
  # Generate the main plot
  ggplot(co2_results, aes(x = reorder(Product, -CO2_med), y = CO2_med, color = Proteinquelle)) +
    geom_point(size = 4) +  # Punkte etwas vergrößern, damit die Farbe gut wirkt
    scale_color_manual(values = c("Vegan" = "#2E8D59", 
                                  "Tierisch" = "#2E86C1")) +
    
    # min-max range
    geom_errorbar(aes(ymin = CO2_min, ymax = CO2_max), 
                  width = 0.2, 
                  color = "darkgray", 
                  linewidth = 0.8) +
    
    # names left
    coord_flip() +
    
    labs(
      title = "CO2-Fußabdruck von Hundefutter",
      subtitle = "Median inkl. Min.-Max. Bereichen",
      x = NULL, 
      y = "kg CO2-Äquivalente pro kg Futter-Inhaltsstoffe"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  # Generate the plot for median <5
  # filter the table
    co2_results_subgroup <- co2_results %>% 
    filter(CO2_med < 5)
    
    ggplot(co2_results_subgroup, aes(x = reorder(Product, -CO2_med), y = CO2_med, color = Proteinquelle)) +
      geom_point(size = 4) +  # Punkte etwas vergrößern, damit die Farbe gut wirkt
      scale_color_manual(values = c("Vegan" = "#2E8D59", 
                                    "Tierisch" = "#2E86C1")) +
    
    # min-max range
    geom_errorbar(aes(ymin = CO2_min, ymax = CO2_max), 
                  width = 0.2, 
                  color = "darkgray", 
                  linewidth = 0.8) +
    
    # names left
    coord_flip() +
    
    labs(
      title = "CO2-Fußabdruck von Hundefutter (niedriger Bereich)",
      subtitle = "Median inkl. Min.-Max. Bereichen",
      x = NULL, 
      y = "kg CO2-Äquivalente pro kg Futter-Inhaltsstoffe"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )
  