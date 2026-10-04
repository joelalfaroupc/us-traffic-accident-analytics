# ==============================================================================
# PARTE 0: PREPARACIÓN DE LA BASE DE DATOS
# ==============================================================================
library(tidyverse)
library(visdat)
library(inspectdf)
library(skimr)
library(DataExplorer)
library(SmartEDA)

dades <- read.csv("Datasets/US_Accidents_Final.csv")

# Pre: Visualización inicial de los datos crudos
ExpData(data = dades,type = 1)

skim(dades)
plot_intro(dades)

# 2. Evaluación de Missings
plot_missing(dades)
vis_miss(dades, warn_large_data = FALSE)

# 3. Análisis de Distribuciones Numéricas
plot_histogram(dades)

# 4. Análisis de Desbalanceo (binarias y categóricas)
inspect_imb(dades) %>% show_plot()
plot_bar(dades, maxcat = 31)

# Limpieza general
dades_clean <- dades %>%
  
  # 1. Eliminar variables con un solo valor (varianza cero)
  select(-State, -Country, -Timezone, -Bump, -Roundabout, -Turning_Loop) %>%
  
  # 2. Eliminar ID, Source, End_Lat/Lng y Street por cardinalidad/no aporta información
  select(-ID, -Source, -End_Lat, -End_Lng, -Street, -City, -Zipcode) %>%
  
  # 3. Convertir entradas de texto vacío "" a NAs
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  
  # 4. Convertir la variable objetivo a factor ordenado
  mutate(Severity = factor(Severity, levels = c(1, 2, 3, 4), ordered = TRUE)) %>%
  # 5. Eliminar variables temporales redundantes
  select(-Weather_Timestamp, -End_Time) %>%
  
  # 6. Transformar a tipos de variables a factores
  mutate(across(c(Amenity, Crossing, Give_Way, Junction, No_Exit, 
                  Railway, Station, Stop, Traffic_Calming, Traffic_Signal), 
                ~ as.factor(.))) %>%
  
  mutate(across(c(where(is.character), -Description, -Start_Time), as.factor))

# ==============================================================================
# PARTE I: EDA INICIAL
# ==============================================================================
# ------------------------------------------------------------------------------
# EDA BASICO UNIVARIANTE:
# ------------------------------------------------------------------------------
# 1. Visualización inicial de los datos limpios
skim(dades_clean)
plot_intro(dades_clean)

# 2. Análisis de Missings
plot_missing(dades_clean)
vis_miss(dades_clean, warn_large_data = FALSE)

# 3. Análisis de Distribuciones Numéricas
plot_histogram(dades_clean)

# 4. Análisis de Desbalanceo (binarias y categóricas)
inspect_imb(dades_clean) %>% show_plot()

# ------------------------------------------------------------------------------
# EDA BASICO BIVARIANTE:
# ------------------------------------------------------------------------------
# 1. Matriz de correlación vars numericas incluyendo Severity
dades_clean %>%
  mutate(Severity_Num = as.numeric(as.character(Severity))) %>%
  select(where(is.numeric)) %>%
  vis_cor()

# 2. Relación categóricas/binarias vs Severity
plot_bar(dades_clean, by = "Severity")

# 3. Relación numéricas vs Severity
plot_boxplot(dades_clean, by = "Severity")

# ------------------------------------------------------------------------------
# SMARTEDA:
# ------------------------------------------------------------------------------

# Estructura
ExpData(data = dades_clean,type = 2)

### VARIABLES NUMERICAS ###

## Summary statistics by – overall
ExpNumStat(dades_clean,by="A",gp=NULL,Qnt=seq(0,1,0.1),MesofShape=2,Outlier=TRUE,round=2)

## Summary statistics by – overall with correlation 
ExpNumStat(dades_clean,by="A",gp="Severity",Qnt=seq(0,1,0.1),MesofShape=1,Outlier=TRUE,round=2)

## Summary statistics by – category
ExpNumStat(dades_clean,by="GA",gp="Severity",Qnt=seq(0,1,0.1),MesofShape=2,Outlier=TRUE,round=2)

# GRAFICOS

## Density plot
ExpNumViz(dades_clean,target=NULL,type=3,nlim=25)

### VARIABLES CATEGORICAS/BIN ###

## Frequency or custom tables for categorical variables
ExpCTable(dades_clean,Target=NULL,margin=1,clim=10,nlim=5,round=2,bin=NULL,per=T)
## Summary statistics of categorical variables
ExpCatStat(dades_clean,Target="Severity",result = "Stat",clim=10,nlim=5)

# GRAFICOS

## column chart
ExpCatViz(dades_clean,target="Severity",fname=NULL,clim=10,col=NULL,margin=2,sample=12)

# Exportar
write.csv(dades_clean, "US_Accidents_Clean.csv", row.names = FALSE)
saveRDS(dades_clean, "US_Accidents_Clean.rds")
# ==============================================================================
# PARTE II: IMPUTACIÓN CON MICE
# ==============================================================================

library(VIM)
library(mice)

vars_a_imputar <- c(
  "Airport_Code",
  "Temperature.F.", "Wind_Chill.F.", "Humidity...", "Pressure.in.", 
  "Visibility.mi.", "Wind_Speed.mph.", "Precipitation.in.",
  "Wind_Direction", "Weather_Condition", "Sunrise_Sunset", 
  "Civil_Twilight", "Nautical_Twilight", "Astronomical_Twilight"
)

par(mfrow = c(1, 3))

datos_mice <- dades_clean[, vars_a_imputar]

imputed_Data <- mice(datos_mice, m = 5, maxit = 50, method = "pmm", seed = 500)
summary(imputed_Data)

dd_mice <- dades_clean
dd_mice[, names(datos_mice)] <- mice::complete(imputed_Data, action = 1)

# Gráficos MICE (numéricos y categórico)
plot_comparativa <- function(df_imputado, df_original, nombre_metodo) {
  try(dev.off(), silent = TRUE) 
  
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1)) 
  
  plot_densidad <- function(var_name, titulo, color_linea) {
    # Extraer los datos válidos
    datos_orig <- as.numeric(na.omit(df_original[[var_name]]))
    datos_imp <- as.numeric(na.omit(df_imputado[[var_name]]))
        if(length(datos_orig) < 2 || length(datos_imp) < 2) {
      warning(paste("Saltando", var_name, "- No hay suficientes datos numéricos."))
      plot.new()
      title(main = paste("Sin datos para", titulo))
      return()
    }
    
    dens_orig <- density(datos_orig)
    dens_imp <- density(datos_imp)
    
    plot(dens_orig, main = paste(nombre_metodo, "-", titulo), col = "black", lwd = 2, 
         ylim = c(0, max(dens_orig$y, dens_imp$y)), xlab = titulo, ylab = "Densidad")
    lines(dens_imp, col = color_linea, lwd = 2, lty = 2)
    legend("topright", legend=c("Original", "Imputado"), col=c("black", color_linea), 
           lwd=2, lty=c(1,2), cex=0.7)
  }
  
  # 1. Temperatura 
  plot_densidad("Temperature.F.", "Temp (F)", "red")
  
  # 2. Humedad
  plot_densidad("Humidity...", "Humedad (%)", "blue")
  
  # 3. Sensación Térmica (Wind Chill) 
  plot_densidad("Wind_Chill.F.", "Wind Chill (F)", "purple")
  
  # 4. Presión Atmosférica
  plot_densidad("Pressure.in.", "Presión (in)", "darkgreen")
  
  # 5. Velocidad del Viento
  plot_densidad("Wind_Speed.mph.", "Vel. Viento (mph)", "brown")
  
  # 6. Sunrise_Sunset (Categórica)
  if(!is.null(df_original$Sunrise_Sunset) && !is.null(df_imputado$Sunrise_Sunset)) {
    prop_orig <- prop.table(table(na.omit(df_original$Sunrise_Sunset)))
    prop_imp <- prop.table(table(na.omit(df_imputado$Sunrise_Sunset)))
    niveles <- union(names(prop_orig), names(prop_imp))
    
    prop_orig_full <- setNames(rep(0, length(niveles)), niveles)
    prop_imp_full <- setNames(rep(0, length(niveles)), niveles)
    prop_orig_full[names(prop_orig)] <- as.numeric(prop_orig)
    prop_imp_full[names(prop_imp)] <- as.numeric(prop_imp)
    
    mat_sun <- rbind(prop_orig_full, prop_imp_full)
    barplot(mat_sun, beside = TRUE, col = c("gray40", "darkorange"),
            main = paste(nombre_metodo, "- Día/Noche"), ylab = "Proporción",
            ylim = c(0, max(mat_sun) * 1.25), las = 1)
    legend("topright", legend = c("Original", "Imputado"), 
           fill = c("gray40", "darkorange"), cex = 0.7)
  } else {
    plot.new()
    title("Sin datos Sunrise_Sunset")
  }
  
  par(mfrow = c(1, 1))
}


plot_comparativa(dd_mice, dades_clean, "MICE")

skim(dd_mice)
# Guardar MICE
write.csv(dd_mice, "US_Accidents_Final_Imputado_Mice.csv", row.names = FALSE)
saveRDS(dd_mice, "US_Accidents_Final_Imputado_Mice.rds")

# ==============================================================================
# PARTE III: GESTIÓN DE OUTLIERS POST IMPUTACIÓN
# ==============================================================================
library(EnvStats)
library(chemometrics)
library(Rlof)
library(solitude)
library(adamethods)

tipus <- sapply(dd_mice, class)
varNum <- names(tipus)[tipus %in% c("integer","numeric")]

sub <- dd_mice[, varNum]

# Mahalanobis clásica

d_maha <- mahalanobis(sub, colMeans(sub), cov(sub))
cutoff <- qchisq(0.99, df=ncol(sub))
outliers_classicmaha <- which(d_maha > cutoff)
length(outliers_classicmaha)

par(mfrow= c(1,1))
hist(d_maha, main="Distancias Mahalanobis", col="skyblue")

# Mahalanobis robusta

vars_mah <- c("Temperature.F.", "Wind_Chill.F.", "Humidity...", "Pressure.in.")
datos_mah <- dd_mice[, vars_mah]

resultados_mah <- chemometrics::Moutlier(datos_mah, quantile = 0.99, plot = TRUE)

grados_libertad <- ncol(datos_mah)
limite_teorico <- sqrt(qchisq(0.99, df = grados_libertad))

outliers_mahalanobis <- which(resultados_mah$rd > limite_teorico)
cat("Outliers multivariantes detectados:", length(outliers_mahalanobis), "\n")

# LOF scores
lof_scores <- Rlof::lof(sub, k = 5)

sub_con_lof <- sub %>%
  mutate(LOF_Score = lof_scores)

outliers_lof <- sub_con_lof %>%
  filter(LOF_Score > 2) %>%
  arrange(desc(LOF_Score)) 

print(nrow(outliers_lof))
#########################
#KNN
datos_knn <- sub %>%
  scale()

outliers_knn <- do_knno(datos_knn, k = 1, top_n = 200)

#########################
# Isolation Forest

datos_if <- dd_mice %>% select(all_of(varNum))

# Inicialización del modelo
isoforest <- isolationForest$new(
  sample_size = as.integer(nrow(datos_if) / 2),
  num_trees   = 500, 
  replace     = TRUE,
  seed        = 123
)

# Entrenamiento del modelo
isoforest$fit(dataset = datos_if)
predicciones <- isoforest$predict(data = datos_if)

# Unificar datos
data_evaluada <- dd_mice %>% bind_cols(predicciones)

# Visualización
ggplot(data = predicciones, aes(x = average_depth)) +
  geom_histogram(color = "gray40", fill = "skyblue", bins = 50) +
  geom_vline(
    xintercept = quantile(predicciones$average_depth, seq(0, 1, 0.1)),
    color      = "red",
    linetype   = "dashed") +
  labs(
    title = "Distancia promedio en el modelo Isolation Forest",
    x = "Average Depth",
    y = "Frecuencia"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 11))


##########################################################
outliers_extremos <- data_evaluada %>%
  filter(average_depth < 11) %>%
  arrange(average_depth) %>%
  select(average_depth, Start_Lat, Start_Lng, Distance.mi., 
         Temperature.F., Wind_Chill.F., Humidity..., Pressure.in., 
         Visibility.mi., Wind_Speed.mph., Precipitation.in.)

cat("Número total de outliers (<11): ", nrow(outliers_extremos), "\n\n")

###########################################
# Identificar Verdaderos Outliers

votos_outliers <- tibble(ID_Fila = 1:nrow(dd_mice))

votos_outliers <- votos_outliers %>%
  mutate(
    Maha_Clasica = if_else(ID_Fila %in% outliers_classicmaha, 1, 0),
    Maha_Robusta = if_else(ID_Fila %in% outliers_mahalanobis, 1, 0),
    KNN          = if_else(ID_Fila %in% outliers_knn, 1, 0),
    LOF          = if_else(lof_scores > 2, 1, 0),
    IsoForest    = if_else(predicciones$average_depth < 11, 1, 0)
  )

votos_outliers <- votos_outliers %>%
  mutate(Total_Votos = Maha_Clasica + Maha_Robusta + LOF + KNN + IsoForest)

table(votos_outliers$Total_Votos)
# Juntar los datos para visualizar
outliers_consenso <- dd_mice %>%
  mutate(Total_Votos = votos_outliers$Total_Votos) %>%
  filter(Total_Votos >= 2) %>%
  arrange(desc(Total_Votos))

cat("\nTotal de outliers reales:", nrow(outliers_consenso), "\n")

# Detectar outliers anómalos

accidentes_anomalos <- outliers_consenso %>%
  filter(
    # Distancias extremas
    Distance.mi. > 5 |
      
      # Sensación mayor a la temp real con viento
      Wind_Chill.F. > Temperature.F. |
      
      # Bajadas extremas de temperatura
      (Temperature.F. - Wind_Chill.F.) > 62 |
      
      # Lluvias masivas
      Precipitation.in. > 5
  ) %>%
  # Ordenar para ver lo más grave primero y seleccionar variables clave
  arrange(desc(Distance.mi.)) %>%
  select(Distance.mi., Temperature.F., Wind_Chill.F., Precipitation.in., Severity, Start_Time, Description)

cat("Num. total de anomalías detectadas:", nrow(accidentes_anomalos), "\n\n")

# Borrar los outliers sin sentido físico, mantener los outliers que sí que tengan sentido

data_final <- anti_join(dd_mice, accidentes_anomalos, by = c("Start_Time", "Description"))

write.csv(data_final, "US_Accidents_Final_Outliers.csv", row.names = FALSE)
saveRDS(data_final, "US_Accidents_Final_Outliers.rds")

# ==============================================================================
# PARTE IV: FEATURE ENGINEERING
# ==============================================================================
library(caret)
library(corrplot)
library(XICOR)
library(vcd)
library(randomForest)

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

# ----------------------------------------------------
# 1. FILTRO: VARIANZA NULA
# ----------------------------------------------------
cat("\n--- 1. Variables con varianza casi nula ---\n")
varianza <- nearZeroVar(dades_fe, saveMetrics = TRUE)
print(varianza[varianza$nzv == TRUE, ]) 
# ELIMINAR Visibility.mi., agrupar el resto en una sola variable: Has_Infrastructure

# ----------------------------------------------------
# 2. FILTRO: CORRELACIÓN DE PEARSON (Numéricas)
# ----------------------------------------------------
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

# ----------------------------------------------------
# 3. RELACIONES NO LINEALES (XICOR)
# ----------------------------------------------------
cat("\n--- 3. Coeficiente de Correlación Xi ---\n")

n_vars <- ncol(dades_num)
nombres <- colnames(dades_num)

# Inicializar una matriz vacía
matriz_xi <- matrix(0, nrow = n_vars, ncol = n_vars, dimnames = list(nombres, nombres))

# Llenar la matriz cruzando todas las vars numericas
for(i in 1:n_vars) {
  for(j in 1:n_vars) {
    if(i == j) {
      matriz_xi[i, j] <- 1
    } else {
      # tryCatch para evitar errores
      matriz_xi[i, j] <- tryCatch({
        xicor(dades_num[[i]], dades_num[[j]])
      }, error = function(e) NA)
    }
  }
}

# Mostrar las relaciones funcionales mas fuertes
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
  
  # 1. Definir las variables a evaluar (todas menos Severity)
  vars_eval <- setdiff(names(dades_cat), "Severity")
  
  # 2. Calcular V de Cramer vectorizado y directo
  cramer_valores <- sapply(vars_eval, function(v) {
    assocstats(table(dades_cat[[v]], dades_cat$Severity))$cramer
  })
  
  # 3. Montar el df ordenado de mayor a menor
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
# Preparamos el dataset eliminando variables identificativas/texto libre que rompen el RF
cols_a_excluir <- c("Start_Time", "Description")
dades_rf <- dades_fe %>% select(-any_of(cols_a_excluir))

if("Severity" %in% names(dades_rf)) {
  # Aseguramos que Severity es un factor para que R entienda que es Clasificación
  dades_rf$Severity <- as.factor(dades_rf$Severity)
  
  # Ejecutamos con toda la base (8.9k observaciones es manejable)
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

# ==============================================================================
# PARTE V: EDA SOBRE LOS DATOS PROCESADOS + COMPARATIVA ANTES VS DESPUÉS
# ==============================================================================

# Base antes y después del procesamiento final
datos_antes <- dades_clean
datos_despues  <- dades_fe_final

# ------------------------------------------------------------------------------
# 1. EDA SOBRE LA BASE FINAL PROCESADA
# ------------------------------------------------------------------------------

ExpData(data = dades_fe_final, type = 1)
skim(dades_fe_final)
plot_intro(dades_fe_final)

plot_missing(dades_fe_final)
vis_miss(dades_fe_final, warn_large_data = FALSE)

plot_histogram(dades_fe_final)

inspect_imb(dades_fe_final) %>% show_plot()
plot_bar(dades_fe_final, maxcat = 30)

dades_fe_final %>%
  mutate(Severity_Num = as.numeric(as.character(Severity))) %>%
  select(where(is.numeric)) %>%
  vis_cor()

plot_bar(dades_fe_final, by = "Severity")
plot_boxplot(dades_fe_final, by = "Severity")

ExpData(data = dades_fe_final, type = 2)

ExpNumStat(dades_fe_final, by = "A", gp = NULL, Qnt = seq(0,1,0.1),
           MesofShape = 2, Outlier = TRUE, round = 2)

ExpNumStat(dades_fe_final, by = "A", gp = "Severity", Qnt = seq(0,1,0.1),
           MesofShape = 1, Outlier = TRUE, round = 2)

ExpNumStat(dades_fe_final, by = "GA", gp = "Severity", Qnt = seq(0,1,0.1),
           MesofShape = 2, Outlier = TRUE, round = 2)

ExpNumViz(dades_fe_final, target = NULL, type = 3, nlim = 25)

ExpCTable(dades_fe_final, Target = NULL, margin = 1, clim = 10,
          nlim = 5, round = 2, bin = NULL, per = TRUE)

ExpCatStat(dades_fe_final, Target = "Severity", result = "Stat",
           clim = 10, nlim = 5)

ExpCatViz(dades_fe_final, target = "Severity", fname = NULL,
          clim = 10, col = NULL, margin = 2, sample = 6)

# ------------------------------------------------------------------------------
# COMPARATIVA ESTRUCTURAL
# ------------------------------------------------------------------------------

comparativa_estructura <- tibble(
  Dataset = c("Antes", "Después"),
  Filas = c(nrow(datos_antes), nrow(datos_despues)),
  Columnas = c(ncol(datos_antes), ncol(datos_despues)),
  Num_Numericas = c(sum(sapply(datos_antes, is.numeric)),
                    sum(sapply(datos_despues, is.numeric))),
  Num_Factores = c(sum(sapply(datos_antes, is.factor)),
                   sum(sapply(datos_despues, is.factor))),
  Missings_Totales = c(sum(is.na(datos_antes)),
                       sum(is.na(datos_despues))),
  Porcentaje_Missings = c(mean(is.na(datos_antes)) * 100,
                          mean(is.na(datos_despues)) * 100)
)

print(comparativa_estructura)

# ==============================================================================
# PARTE V: CLUSTERING DE LA TS
# ==============================================================================
library(lubridate)
library(reshape2)
library(dtwclust)
library(dendextend)

# Convertimos Start_Time a fecha y extraemos el Año-Mes (Periodo)
data_final <- dades_fe_final %>%
  mutate(
    Start_Time = ymd_hms(Start_Time),
    Periodo = format(Start_Time, "%Y-%m")
  ) %>%
  filter(!is.na(Periodo), !is.na(County)) # Limpieza de seguridad

# Agrupamos por Condado y Periodo, contando el número de accidentes
df_agrupado <- data_final %>%
  group_by(County, Periodo) %>%
  summarise(Total = n(), .groups = "drop")

# Transformamos la bbdd a formato ancho (Filas = County, Columnas = Periodo)
datos_wide <- dcast(df_agrupado, County ~ Periodo, value.var = "Total", fill = 0)
rownames(datos_wide) <- datos_wide[, "County"]
datos_wide[, "County"] <- NULL

# Detectamos si hay alguna columna (Mes) que tenga 0 accidentes en todos los condados
columnas_vacias <- colSums(datos_wide) == 0
if(any(columnas_vacias)){
  datos_wide <- datos_wide[, !columnas_vacias]
}

# 3. Escalamos los datos para comparar "formas" de la serie
# Escalamos las filas (Z-score) para agrupar por patrón de estacionalidad
datos <- t(scale(t(datos_wide)))
datos <- na.omit(datos) # Eliminamos si algún condado no tiene varianza (generaría NAs)

dim(datos)
head(datos)

# 4. Calculamos la distancia DTW y clustering jerárquico base
distMatrix <- proxy::dist(datos, method = "DTW")

# Generamos el clustering jerárquico tradicional
hcc <- hclust(distMatrix, method = "ward.D2")
hc1 <- hclust(distMatrix, method = "complete")
hc2 <- hclust(distMatrix, method = "single")
hc3 <- hclust(distMatrix, method = "average")

# Evaluamos qué método preserva mejor las distancias originales
cor(distMatrix, cophenetic(hcc)) # Ward.D2
cor(distMatrix, cophenetic(hc1)) # Complete
cor(distMatrix, cophenetic(hc2)) # Single
cor(distMatrix, cophenetic(hc3)) # Average


# 5. LOOP PARA VISUALIZAR DIFERENTES VALORES DE K 
valores_k <- 2:6 # Prueba valores de K desde 2 hasta 6

for(k_actual in valores_k) {
  
  cat("\n--------------------------------------------------\n")
  cat("Generando visualizaciones para K =", k_actual, "\n")
  cat("--------------------------------------------------\n")
  
  # ..............................................................................
  # 1. Visualización del Dendrograma
  dend <- hcc %>% as.dendrogram()
  plot(dend %>% 
         set("branches_k_color", k = k_actual) %>% 
         set("labels_cex", 0.7),
       main = paste("Clustering DTW de Condados (K =", k_actual, ")"))
  
  # ..............................................................................
  ## Partitional
  pc <- tsclust(datos, type = "partitional", k = k_actual, 
                distance = "dtw_basic", centroid = "pam", 
                seed = 3247L, trace = FALSE, 
                args = tsclust_args(dist = list(window.size = 20L)))
  plot(pc, main = paste("Partitional Clustering (K =", k_actual, ")"))
  
  # ------------------------------------------------------------------------------
  ## Hierarchical
  hc <- tsclust(datos, type = "hierarchical", k = k_actual, 
                distance = "dtw_basic", trace = FALSE,
                control = hierarchical_control(method = "ward.D2"))
  
  plot(hc, type = "sc", main = paste("Hierarchical - Series Solapadas (K =", k_actual, ")")) 
}

# 6. EXPORTAMOS RESULTADOS FINALES
k_final <- 6

hc_final <- tsclust(datos, type = "hierarchical", k = k_final, 
                    distance = "dtw_basic", trace = FALSE,
                    control = hierarchical_control(method = "ward.D2"))

resultados <- data.frame(County = rownames(datos), Cluster = hc_final@cluster)

# Unimos el cluster asignado a las series temporales escaladas
resultados_completos <- cbind(resultados, datos)

# Guardamos los CSV
write.csv(resultados_completos, paste0("resultados_condados_clusters_k", k_final, ".csv"), row.names = FALSE)
resultados_t <- t(resultados_completos)
write.csv(resultados_t, paste0("resultados_condados_clusters_t_k", k_final, ".csv"))
