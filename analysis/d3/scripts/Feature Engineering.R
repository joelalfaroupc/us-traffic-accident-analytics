# =======================
# FEATURE ENGINEERING
# =======================
library(tidyverse)
library(caret)
library(corrplot)
library(XICOR)
library(vcd)
library(randomForest)

data_final <- readRDS("Datasets/US_Accidents_Final_Outliers.rds")


  # 1. Agrupación de variables
dades_fe <- data_final %>%
  mutate(
    # Weather_Condition
    Weather_Group = case_when(
      is.na(Weather_Condition) ~ NA_character_,
      Weather_Condition %in% c("Clear", "Fair", "Fair / Windy") ~ "Clear",
      Weather_Condition %in% c("Cloudy", "Partly Cloudy", "Scattered Clouds",
                               "Mostly Cloudy", "Overcast", "Mostly Cloudy / Windy",
                               "Cloudy / Windy", "Partly Cloudy / Windy") ~ "Cloudy",
      
      Weather_Condition %in% c("Rain", "Light Rain", "Light Drizzle", "Drizzle",
                               "Light Rain / Windy", "Rain / Windy") ~ "Light_Rain",
      
      Weather_Condition %in% c("Heavy Rain", "Heavy Rain / Windy", "Heavy Drizzle",
                               "Thunderstorms and Rain", "Thunder", "T-Storm",
                               "Thunder in the Vicinity", "Thunderstorm",
                               "Heavy T-Storm / Windy", "Heavy T-Storm",
                               "Light Thunderstorms and Rain", "Thunder / Windy",
                               "Heavy Thunderstorms and Rain", "Thunder and Hail",
                               "Light Rain with Thunder") ~ "Heavy_Storms",
      
      Weather_Condition %in% c("Snow", "Light Snow", "Heavy Snow", "Snow / Windy",
                               "Blowing Snow", "Light Snow / Windy") ~ "Snow",
      
      Weather_Condition %in% c("Light Freezing Drizzle", "Wintry Mix", "Snow and Sleet",
                               "Light Snow and Sleet", "Light Freezing Rain", "Ice Pellets",
                               "Light Snow and Sleet / Windy", "Snow and Sleet / Windy",
                               "Light Sleet") ~ "Ice_Sleet_Mix",
      
      Weather_Condition %in% c("Fog", "Haze", "Shallow Fog", "Patches of Fog",
                               "Drizzle and Fog", "Haze / Windy") ~ "Fog_Haze",
      TRUE ~ "Other"
    ),
    Weather_Group = as.factor(Weather_Group),
    
    # Wind_Direction
    Wind_Group = case_when(
      is.na(Wind_Direction) ~ NA_character_,
      Wind_Direction %in% c("Calm", "CALM") ~ "Calm",
      Wind_Direction %in% c("Var", "VAR", "Variable") ~ "Variable",
      str_starts(Wind_Direction, "N") ~ "North",
      str_starts(Wind_Direction, "S") ~ "South",
      str_starts(Wind_Direction, "E") ~ "East",
      str_starts(Wind_Direction, "W") ~ "West",
      TRUE ~ "Other"
    ),
    Wind_Group = as.factor(Wind_Group)
  ) %>%
  
  select(-Weather_Condition, -Wind_Direction) %>%
  
  # Eliminar distancia, ya que representa información conocida DESPUÉS del accidente
  select(-Distance.mi.) %>%
  
  # Conversión de Fahrenheit a Celsius para interpretabilidad
  mutate(
    Temperature.C. = round((Temperature.F. - 32) * 5 / 9, 1),
    Wind_Chill.C.  = round((Wind_Chill.F. - 32) * 5 / 9, 1)
  ) %>%
  select(-Temperature.F., -Wind_Chill.F.)

###########

numeric_cols <- sapply(dades_fe, is.numeric)
dades_num <- dades_fe[, numeric_cols]
dades_cat <- dades_fe[, !numeric_cols]

# ==============================================================================
# 1. FILTRO: VARIANZA NULA
# ==============================================================================
cat("\n--- 1. Variables con varianza casi nula ---\n")
varianza <- nearZeroVar(dades_fe, saveMetrics = TRUE)
print(varianza[varianza$nzv == TRUE, ]) 
# ELIMINAR Visibility.mi., agrupar el resto en una sola variable: Has_Infrastructure

# ==============================================================================
# 2. FILTRO: CORRELACIÓN DE PEARSON (Numéricas)
# ==============================================================================
cat("\n--- 2. Matriz de Correlación ---\n")
matriz_corr <- cor(dades_num)

# Detectar combinaciones lineales altas (redundancia)
highlyCorDescr <- findCorrelation(matriz_corr, cutoff = 0.85)
if(length(highlyCorDescr) > 0) {
  cat("Redundancia detectada. Considera eliminar:\n")
  print(colnames(dades_num)[highlyCorDescr])
} else {
  cat("No hay variables numéricas altamente correlacionadas (>0.85).\n")
}
corrplot(matriz_corr, method = "circle")

# ==============================================================================
# 3. RELACIONES NO LINEALES (XICOR)
# ==============================================================================
cat("\n--- 3. Coeficiente de Correlación Xi ---\n")

var_temp <- "Temperature.C." #
var_hum <- "Humidity..."     #

if(var_temp %in% names(dades_num) & var_hum %in% names(dades_num)) {
  xi_result <- xicor(dades_num[[var_temp]], dades_num[[var_hum]])
  cat("Coeficiente Xi (", var_temp, " vs ", var_hum, "):", round(xi_result, 4), "\n\n")
} else {
  cat("Advertencia: Las variables especificadas para el ejemplo de Xi no existen.\n\n")
}

cat("Calculando matriz completa de correlación no lineal Xi...\n")

n_vars <- ncol(dades_num)
nombres <- colnames(dades_num)

# Inicializamos una matriz vacía
matriz_xi <- matrix(0, nrow = n_vars, ncol = n_vars, dimnames = list(nombres, nombres))

# Llenamos la matriz cruzando todas las variables numéricas
for(i in 1:n_vars) {
  for(j in 1:n_vars) {
    if(i == j) {
      matriz_xi[i, j] <- 1
    } else {
      # tryCatch evita que el bucle entero colapse si una variable tiene problemas (ej. varianza cero residual)
      matriz_xi[i, j] <- tryCatch({
        xicor(dades_num[[i]], dades_num[[j]])
      }, error = function(e) NA)
    }
  }
}

# Mostrar las relaciones funcionales más fuertes (mayores a 0.6)
cat("Relaciones no lineales fuertes detectadas (Xi > 0.6):\n")
indices_altos <- which(matriz_xi > 0.6 & matriz_xi < 0.99, arr.ind = TRUE)

if(nrow(indices_altos) > 0) {
  for(k in 1:nrow(indices_altos)) {
    var_x <- rownames(matriz_xi)[indices_altos[k, 1]]
    var_y <- colnames(matriz_xi)[indices_altos[k, 2]]
    valor_xi <- matriz_xi[indices_altos[k, 1], indices_altos[k, 2]]
    
    cat(sprintf("%s -> %s : %.3f\n", var_x, var_y, valor_xi))
  }
} else {
  cat("No se detectaron relaciones no lineales fuertemente dependientes.\n")
}

# ==============================================================================
# 4. ASOCIACIÓN CATEGÓRICA (V de Cramer)
# ==============================================================================
cat("\n--- 4. Asociación entre Categóricas (V de Cramer) ---\n")

if("Severity" %in% names(dades_cat)) {
  
  # 1. Definimos las variables a evaluar (todas menos Severity)
  vars_eval <- setdiff(names(dades_cat), "Severity")
  
  # 2. Calculamos V de Cramer vectorizado y directo
  cramer_valores <- sapply(vars_eval, function(v) {
    assocstats(table(dades_cat[[v]], dades_cat$Severity))$cramer
  })
  
  # 3. Montamos el dataframe final y lo ordenamos de mayor a menor
  resultados_cramer <- data.frame(Variable = vars_eval, Cramer_V = cramer_valores)
  resultados_cramer <- resultados_cramer[order(-resultados_cramer$Cramer_V), ]
  
  print(resultados_cramer)
  
} else {
  cat("Error: 'Severity' no está en dades_cat.\n")
}


# ==============================================================================
# 5. RANDOM FOREST
# ==============================================================================
cat("\n--- 5. Importancia de Variables mediante Random Forest ---\n")
# Preparamos el dataset eliminando variables identificativas/texto libre
cols_a_excluir <- c("Start_Time", "Description")
dades_rf <- dades_fe %>% select(-any_of(cols_a_excluir))

if("Severity" %in% names(dades_rf)) {
  dades_rf$Severity <- as.factor(dades_rf$Severity)
  
  set.seed(123)
  rf_modelo <- randomForest(Severity ~ ., data = dades_rf, ntree = 100, importance = TRUE)
  
  print(importance(rf_modelo))
  varImpPlot(rf_modelo, main = "Top Variables por Importancia")
} else {
  cat("Variable 'Severity' no encontrada. Imposible generar Random Forest.\n")
}

#### Ajustes finales a la BD ####

valores_true <- c("True", "TRUE", TRUE, "1", 1)

dades_fe_final <- dades_fe %>%
  # Agrupar las variables de infraestructura
  mutate(
    Has_Infrastructure = as.factor(ifelse(
      (Amenity %in% valores_true) | 
        (Crossing %in% valores_true) | 
        (Give_Way %in% valores_true) | 
        (No_Exit %in% valores_true) | 
        (Railway %in% valores_true) | 
        (Station %in% valores_true) | 
        (Stop %in% valores_true) | 
        (Traffic_Calming %in% valores_true), 
      1, 0))
  ) %>%
  select(
    # Agrupadas en Has_Infrastructure
    -Amenity, -Crossing, -Give_Way, -No_Exit, -Railway, -Station, -Stop, -Traffic_Calming,
    # Eliminada por varianza nula (Paso 1)
    -Visibility.mi.,
    # Eliminadas por redundancia y dispersión de importancia detectada en RF (Paso 5)
    -Astronomical_Twilight, -Nautical_Twilight, -Civil_Twilight
  )

write.csv(dades_fe_final, "US_Accidents_Final_FE.csv", row.names = FALSE)
saveRDS(dades_fe_final, "US_Accidents_Final_FE.rds")