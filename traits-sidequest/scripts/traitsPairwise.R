# Load required packages
library(cluster)  # for daisy()
library(dplyr)
library(readr)
library(readxl)


# LOADING/FILTERING DATA --------------------------------------------------
# Read your dataset (replace with your actual path)
df <- read_xlsx("traits-sidequest/data/turtle-traits-ours.xlsx")

# select relevant columns
df_select <- df %>%
  select(Species,Order,Suborder,Family,Genus,                  # TAXONOMY
        Microhabitat, `Habitat type`, `Mean Annual Temperature (°C)`,  # ECOLOGY
        `Temperature Seasonality (standard deviation*100)`, 
        Diet, `Active time`, `Foraging mode`, 
        `Maximum Longevity (years)`, `Sex-determining mechanism (GSD or TSD)`,  # LIFE-HISTORY
        `Mean number of offspring per litter or number of eggs per clutch`,
        `Number of litters or clutches produced per year`,
        `Maximum body mass (g)`, # MORPHOLOGY
        `Maximum length ("SVL", mm)/straight carapace length for turtles ("SCL", mm)`) 

library(dplyr)
library(janitor)  # for clean_names()

df_cleaned <- df_select %>%
  janitor::clean_names() %>%  # standardize column names: lowercase, underscores
  mutate(
    # Character identifiers
    species = as.character(species),
    order = as.character(order),
    suborder = as.character(suborder),
    family = as.character(family),
    genus = as.character(genus),
    
    # Categorical ecological and life history traits
    microhabitat = as.factor(microhabitat),
    habitat_type = as.factor(habitat_type),
    diet = as.factor(diet),
    active_time = as.factor(active_time),
    foraging_mode = as.factor(foraging_mode),
    sex_determining_mechanism_gsd_or_tsd = as.factor(sex_determining_mechanism_gsd_or_tsd),
    
    # Numeric traits
    mean_annual_temperature_c = as.numeric(mean_annual_temperature_c),
    temperature_seasonality_standard_deviation_100 = as.numeric(temperature_seasonality_standard_deviation_100),
    maximum_longevity_years = as.numeric(maximum_longevity_years),
    mean_number_of_offspring_per_litter_or_number_of_eggs_per_clutch = as.numeric(mean_number_of_offspring_per_litter_or_number_of_eggs_per_clutch),
    number_of_litters_or_clutches_produced_per_year = as.numeric(number_of_litters_or_clutches_produced_per_year),
    maximum_body_mass_g = as.numeric(maximum_body_mass_g),
    maximum_length_svl_mm_straight_carapace_length_for_turtles_scl_mm = as.numeric(maximum_length_svl_mm_straight_carapace_length_for_turtles_scl_mm)
  )

# DEFINING HELPER FUNCTIONS  ----------------------------------------------
# Step 1: Define helper functions
# Function to impute numeric values with genus-level mean

library(dplyr)
col <- 'maximum_body_mass_g'
df <- df_cleaned
# numeric
impute_numeric_by_genus <- function(df, col, by = "genus") {
  col_sym <- sym(col)
  by_sym <- sym(by)
  
  df %>%
    group_by(!!by_sym) %>%
    mutate(!!col_sym := if_else(is.na(!!col_sym),
                                median(!!col_sym, na.rm = TRUE),
                                !!col_sym)) %>%
    ungroup()
}

# Function to impute missing values in a categorical column by the mode within genus (or any other group)
impute_categorical_by_genus <- function(df, col, by = "genus") {
  col_sym <- sym(col)
  by_sym <- sym(by)
  
  mode_fn <- function(x) {
    ux <- na.omit(x)
    if (length(ux) == 0) return(NA)
    ux[which.max(tabulate(match(ux, ux)))]
  }
  
  df %>%
    group_by(!!by_sym) %>%
    mutate(!!col_sym := if_else(is.na(!!col_sym), mode_fn(!!col_sym), !!col_sym)) %>%
    ungroup()
}


# IMPUTE ------------------------------------------------------------------

# Step 2: Apply to your selected data
df_imputed <- df_cleaned

# Separate numeric and categorical columns (excluding IDs)
num_cols <- df_imputed %>%
  select(where(is.numeric)) %>%
  names()

cat_cols <- df_imputed %>%
  select(where(~is.character(.) || is.factor(.))) %>%
  select(-species, -order, -suborder, -family, -genus) %>%
  names()

# Step 3: Apply imputations
# Impute numeric columns
for (col in num_cols) {
  df_imputed <- impute_numeric_by_genus(df_imputed, col,by = 'genus')
}
# Impute categorical columns
for (col in cat_cols) {
  df_imputed <- impute_categorical_by_genus(df_imputed, col, by = 'genus' )
}
table(is.na(df_imputed))

# Step 3: Apply imputations
# Impute numeric columns
for (col in num_cols) {
  df_imputed <- impute_numeric_by_genus(df_imputed, col,by = 'family')
}
# Impute categorical columns
for (col in cat_cols) {
  df_imputed <- impute_categorical_by_genus(df_imputed, col, by = 'family' )
}
table(is.na(df_imputed))

# Step 3: Apply imputations
# Impute numeric columns
for (col in num_cols) {
  df_imputed <- impute_numeric_by_genus(df_imputed, col,by = 'order')
}
# Impute categorical columns
for (col in cat_cols) {
  df_imputed <- impute_categorical_by_genus(df_imputed, col, by = 'order' )
}
table(is.na(df_imputed))



# Optional: convert character columns to factors
df_imputed <- df_imputed %>%
  mutate(across(all_of(cat_cols), as.factor))


# OLD ---------------------------------------------------------------------
my_order <- c(
  "Malaclemys terrapin",
  "Chrysemys picta",
  "Terrapene mexicana",
  "Platysternon megacephalum",
  "Cuora mccordi",
  "Cuora amboinensis",
  "Gopherus agassizii",
  "Chelonoidis abingdonii",     # Closest available match to "Chelonoidis abingdonii"
  "Dermochelys coriacea",
  "Chelonia mydas",
  "Dermatemys mawii",
  "Chelydra serpentina",       
  "Carettochelys insculpta",
  "Pelodiscus sinensis",
  "Pelusios castaneus",
  "Podocnemis expansa",
  "Emydura subglobosa",
  "Mesoclemmys tuberculata"
)

# change name for chelonoidis
df_imputed <- df_imputed %>%
  mutate(species = ifelse(species == "Chelonoidis carbonarius", 
                          "Chelonoidis abingdonii", 
                          species))

df_cleaned_subset <- df_imputed %>%
  filter(species %in% my_order)


# TAXA --------------------------------------------------------------------
# Set species names as rownames
species_names <- df_cleaned_subset$species
df_imputed_taxon <- df_cleaned_subset %>% select(family,genus)

df_imputed_taxon <- df_imputed_taxon %>%
  mutate(across(where(is.character), as.factor))

# Compute Gower distance (works with mixed types)
gower_dist <- daisy(df_imputed_taxon, metric = "gower")

# Convert to a matrix and set row/col names
gower_matrix <- as.matrix(gower_dist)
rownames(gower_matrix) <- colnames(gower_matrix) <- species_names


gower_matrix_ordered <- gower_matrix[my_order, my_order]

# Plot the heatmap with clustering
pheatmap(gower_matrix_ordered,
         cluster_rows = F,
         cluster_cols = F,
         main = "Pairwise Trait Distance (Gower)",
         fontsize_row = 8,
         fontsize_col = 8,
         border_color = NA)

# Step 3B: Plot a dendrogram
# Hierarchical clustering using the distance matrix
hc <- hclust(gower_dist, method = "average")
plot(as.phylo(hc), type = "fan", tip.color = "steelblue",
     main = "Dendrogram of Trait-Based Clustering", cex = 0.8)


# ECOLOGY -----------------------------------------------------------------
species_names <- df_cleaned_subset$species
df_imputed_taxon <- df_cleaned_subset %>% select(microhabitat,habitat_type,mean_annual_temperature_c,
                                                 temperature_seasonality_standard_deviation_100,diet,
                                                 active_time,foraging_mode)

df_imputed_taxon <- df_imputed_taxon %>%
  mutate(across(where(is.character), as.factor))

# Compute Gower distance (works with mixed types)
gower_dist <- daisy(df_imputed_taxon, metric = "gower")

# Convert to a matrix and set row/col names
gower_matrix <- as.matrix(gower_dist)
rownames(gower_matrix) <- colnames(gower_matrix) <- species_names


gower_matrix_ordered <- gower_matrix[my_order, my_order]

# Plot the heatmap with clustering
pheatmap(gower_matrix_ordered,
         cluster_rows = F,
         cluster_cols = F,
         main = "Pairwise Trait Distance (Gower)",
         fontsize_row = 8,
         fontsize_col = 8,
         border_color = NA)

# Step 3B: Plot a dendrogram
# Hierarchical clustering using the distance matrix
hc <- hclust(gower_dist, method = "average")
plot(as.phylo(hc), type = "fan", tip.color = "steelblue",
     main = "Dendrogram of Trait-Based Clustering", cex = 0.8)


