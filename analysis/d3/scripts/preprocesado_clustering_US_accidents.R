# ==============================================================================
# PMAAD - PREPROCESADO PREVIO AL CLUSTERING (MUESTRA MULTIESTADO)
# ------------------------------------------------------------------------------
# Objetivo:
#   1) Leer la base original de Kaggle US Accidents.
#   2) Seleccionar 1500 registros de forma equilibrada por estado.
#   3) Aplicar preprocesado de missings siguiendo la lógica del script previo.
#   4) Aplicar detección/tratamiento de outliers siguiendo la lógica del script previo.
#   5) Crear una base final lista para clustering (CURE / k-proto).
#
# Salidas:
#   - 01_sample_1500_raw.csv
#   - 02_sample_1500_imputado.csv
#   - 03_sample_1500_clean.csv
#   - 04_base_clustering.csv
#   - resumen_preprocesado.txt
#
# NOTA:
#   - Este script está pensado para trabajar con el CSV original de Kaggle.
#   - Mantiene la filosofía de los scripts previos de imputación y outliers,
#     pero adaptada a una muestra multiestado para clustering mixto.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. CONFIGURACIÓN
# ------------------------------------------------------------------------------
file_in  <- "US_Accidents_March23.csv"
out_dir  <- "Datasets/clustering"
set.seed(123)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. LIBRERÍAS
# ------------------------------------------------------------------------------
packages <- c(
  "tidyverse", "lubridate", "VIM", "mice", "Hmisc", "mi", "missForest",
  "cluster", "chemometrics", "Rlof", "solitude", "adamethods"
)

new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages, dependencies = TRUE)
invisible(lapply(packages, require, character.only = TRUE))

# ------------------------------------------------------------------------------
# 2. LECTURA DE LA BASE ORIGINAL
# ------------------------------------------------------------------------------
# Usamos check.names = TRUE para conservar el mismo estilo de nombres que ya
# aparece en vuestra base previa: Temperature.F., Humidity..., etc.

dd0 <- read.csv(
  file_in,
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", " "),
  check.names = TRUE
)

cat("Dimensiones base original:", dim(dd0), "\n")

# ------------------------------------------------------------------------------
# 3. SELECCIÓN DE VARIABLES BASE
# ------------------------------------------------------------------------------
# Nos quedamos con el bloque de variables coherente con vuestro trabajo previo.

vars_base <- c(
  "ID", "Severity", "Start_Time", "End_Time",
  "Start_Lat", "Start_Lng", "End_Lat", "End_Lng", "Distance.mi.",
  "Description", "Street", "City", "County", "State", "Zipcode",
  "Country", "Timezone", "Airport_Code", "Weather_Timestamp",
  "Temperature.F.", "Wind_Chill.F.", "Humidity...", "Pressure.in.",
  "Visibility.mi.", "Wind_Direction", "Wind_Speed.mph.",
  "Precipitation.in.", "Weather_Condition",
  "Amenity", "Bump", "Crossing", "Give_Way", "Junction", "No_Exit",
  "Railway", "Roundabout", "Station", "Stop", "Traffic_Calming",
  "Traffic_Signal", "Turning_Loop", "Sunrise_Sunset", "Civil_Twilight",
  "Nautical_Twilight", "Astronomical_Twilight"
)

vars_base <- intersect(vars_base, names(dd0))
dd0 <- dd0[, vars_base, drop = FALSE]

# ------------------------------------------------------------------------------
# 4. MUESTRA EQUILIBRADA DE 1500 REGISTROS
# ------------------------------------------------------------------------------
# Estrategia simple y robusta:
#   - Seleccionamos los 15 estados con mayor volumen.
#   - Tomamos 1000 registros por estado.
#   - Así obtenemos 15000 filas, balanceadas territorialmente.
#
# Es una solución sencilla, defendible y mucho más adecuada para clustering que
# un muestreo totalmente aleatorio dominado por pocos estados.

state_counts <- dd0 %>%
  filter(!is.na(State), State != "") %>%
  count(State, sort = TRUE)

states_validos <- state_counts %>%
  filter(n >= 1000) %>%
  slice_head(n = 15) %>%
  pull(State)

if (length(states_validos) < 15) {
  stop("No hay al menos 15 estados con 1000 o más registros. Revisa el fichero de entrada.")
}

sample_raw <- dd0 %>%
  filter(State %in% states_validos) %>%
  group_by(State) %>%
  slice_sample(n = 1000) %>%
  ungroup()

write.csv(sample_raw,
          file.path(out_dir, "01_sample_15000_raw.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. LIMPIEZA INICIAL Y TIPADO
# ------------------------------------------------------------------------------
sample_raw <- sample_raw %>%
  
  # 1. Eliminar variables con varianza cero y de texto libre
  select(-Country, -Timezone, -Bump, -Roundabout, -Turning_Loop, -Description) %>%
  
  # 2. Eliminar ID y Lat/Lng finales
  select(-ID, -End_Lat, -End_Lng) %>%
  
  # 3. Convertir vacíos a NA
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  
  # 4. Convertir fechas y crear variables temporales de negocio
  mutate(
    Start_Time = suppressWarnings(ymd_hms(Start_Time, quiet = TRUE)),
    End_Time   = suppressWarnings(ymd_hms(End_Time, quiet = TRUE)),
    Duration_Hours = as.numeric(difftime(End_Time, Start_Time, units = "hours")),
    Start_Hour = as.factor(hour(Start_Time)),
    Start_Month = as.factor(month(Start_Time)),
    Is_Weekend = as.factor(ifelse(wday(Start_Time, week_start = 1) %in% c(6, 7), "Weekend", "Weekday"))
  ) %>%
  
  # 5. Eliminar variables temporales originales y de alta cardinalidad
  select(-Start_Time, -End_Time, -Weather_Timestamp, -Street, -Zipcode, -City) %>%
  
  # 6. Agrupaciones de Clima y Viento (TU LÓGICA ORIGINAL)
  mutate(
    Weather_Group = case_when(
      is.na(Weather_Condition) ~ NA_character_,
      Weather_Condition %in% c("Clear", "Fair", "Fair / Windy") ~ "Clear",
      str_detect(tolower(Weather_Condition), "cloud|overcast") ~ "Cloudy",
      str_detect(tolower(Weather_Condition), "rain|drizzle") & !str_detect(tolower(Weather_Condition), "heavy|storm|thunder") ~ "Light_Rain",
      str_detect(tolower(Weather_Condition), "heavy|storm|thunder") ~ "Heavy_Storms",
      str_detect(tolower(Weather_Condition), "snow") ~ "Snow",
      str_detect(tolower(Weather_Condition), "ice|sleet|freez|wintry") ~ "Ice_Sleet_Mix",
      str_detect(tolower(Weather_Condition), "fog|haze") ~ "Fog_Haze",
      TRUE ~ "Other" 
    ),
    Weather_Group = as.factor(Weather_Group),
    
    Wind_Group = case_when(
      is.na(Wind_Direction) ~ NA_character_,
      Wind_Direction %in% c("Calm", "CALM") ~ "Calm",
      Wind_Direction %in% c("Var", "VAR", "Variable") ~ "Variable",
      str_starts(toupper(Wind_Direction), "N") ~ "North",
      str_starts(toupper(Wind_Direction), "S") ~ "South",
      str_starts(toupper(Wind_Direction), "E") ~ "East",
      str_starts(toupper(Wind_Direction), "W") ~ "West",
      TRUE ~ "Other"
    ),
    Wind_Group = as.factor(Wind_Group)
  ) %>%
  select(-Weather_Condition, -Wind_Direction) %>%
  
  # 7. Conversión a Celsius (TU LÓGICA ORIGINAL)
  mutate(
    Temperature.C. = round((Temperature.F. - 32) * 5 / 9, 1),
    Wind_Chill.C.  = round((Wind_Chill.F. - 32) * 5 / 9, 1)
  ) %>%
  select(-Temperature.F., -Wind_Chill.F.) %>%
  
  # 8. Creación de Has_Infrastructure y Tipado
  mutate(
    Severity = factor(Severity, levels = c(1, 2, 3, 4), ordered = TRUE),
    Has_Infrastructure = as.factor(ifelse(
      (Amenity %in% c("True", "TRUE", TRUE, "1", 1)) | 
        (Crossing %in% c("True", "TRUE", TRUE, "1", 1)) | 
        (Give_Way %in% c("True", "TRUE", TRUE, "1", 1)) | 
        (No_Exit %in% c("True", "TRUE", TRUE, "1", 1)) | 
        (Railway %in% c("True", "TRUE", TRUE, "1", 1)) | 
        (Station %in% c("True", "TRUE", TRUE, "1", 1)) | 
        (Stop %in% c("True", "TRUE", TRUE, "1", 1)) | 
        (Traffic_Calming %in% c("True", "TRUE", TRUE, "1", 1)), 
      1, 0)),
    Junction = as.factor(as.logical(Junction)),
    Traffic_Signal = as.factor(as.logical(Traffic_Signal))
  ) %>%
  mutate(across(where(is.character), as.factor)) %>%
  
  # 9. Purga final de variables
  select(
    -Amenity, -Crossing, -Give_Way, -No_Exit, -Railway, -Station, -Stop, -Traffic_Calming,
    -Visibility.mi., 
    -Astronomical_Twilight, -Nautical_Twilight, -Civil_Twilight,
    -Distance.mi.
  )

# Limpieza de duraciones anómalas (negativas o cero)
sample_raw$Duration_Hours[sample_raw$Duration_Hours <= 0] <- NA


# ------------------------------------------------------------------------------
# 6. IMPUTACIÓN DE MISSINGS - DIAGNÓSTICO
# ------------------------------------------------------------------------------
# Variables con missings a imputar: mantenemos la misma lógica del script previo.

vars_a_imputar <- intersect(c(
  "End_Lat", "End_Lng", "Street", "Airport_Code",
  "Temperature.C.", "Wind_Chill.C.", "Humidity...", "Pressure.in.",
  "Wind_Speed.mph.", "Precipitation.in.",
  "Wind_Group", "Weather_Group", "Sunrise_Sunset"
), names(sample_raw))

# Para ser fieles al script original, se conservan varios métodos a modo de
# comparación/diagnóstico. La imputación final será híbrida y práctica para el
# clustering mixto:
#   - Numéricas -> MICE (PMM)
#   - Categóricas -> moda
# Esto mantiene el espíritu del script y deja la base completamente usable.

sample_imp <- sample_raw

# -----------------------------
# 6.1 Diagnóstico visual básico
# -----------------------------
if ("Wind_Chill.C." %in% names(sample_imp)) {
  dens_wind_orig <- tryCatch(
    density(as.numeric(na.omit(sample_imp$Wind_Chill.C.)), na.rm = TRUE),
    error = function(e) NULL
  )

  if (!is.null(dens_wind_orig)) {
    png(file.path(out_dir, "diag_imputacion_windchill.png"), width = 1200, height = 500)
    dd_basic <- sample_imp
    dd_basic$Wind_Mean   <- with(dd_basic, Hmisc::impute(Wind_Chill.C., mean))
    dd_basic$Wind_Random <- with(dd_basic, Hmisc::impute(Wind_Chill.C., "random"))

    par(mfrow = c(1, 2))
    dens_wind_mean   <- density(as.numeric(dd_basic$Wind_Mean), na.rm = TRUE)
    dens_wind_random <- density(as.numeric(dd_basic$Wind_Random), na.rm = TRUE)

    plot(dens_wind_orig, main = "Básica: Wind Chill (Media)", col = "black", lwd = 2,
         ylim = c(0, max(dens_wind_orig$y, dens_wind_mean$y, na.rm = TRUE)),
         xlab = "Sensación Térmica (F)")
    lines(dens_wind_mean, col = "blue", lwd = 2)
    legend("topleft", legend = c("Original", "Imputado Media"),
           col = c("black", "blue"), lwd = 2, cex = 0.8)

    plot(dens_wind_orig, main = "Básica: Wind Chill (Random)", col = "black", lwd = 2,
         ylim = c(0, max(dens_wind_orig$y, dens_wind_random$y, na.rm = TRUE)),
         xlab = "Sensación Térmica (F)")
    lines(dens_wind_random, col = "green", lwd = 2, lty = 2)
    legend("topright", legend = c("Original", "Imputado Random"),
           col = c("black", "green"), lwd = 2, lty = c(1, 2), cex = 0.8)
    dev.off()
  }
}

# -----------------------------
# 6.2 Imputación final práctica
# -----------------------------
# NUMÉRICAS -> MICE (PMM)
num_impute_vars <- intersect(c(
  "End_Lat", "End_Lng", "Temperature.C.", "Wind_Chill.C.", "Humidity...",
  "Pressure.in.", "Wind_Speed.mph.", "Precipitation.in."
), names(sample_imp))

if (length(num_impute_vars) > 0) {
  datos_mice <- sample_imp[, num_impute_vars, drop = FALSE]
  if (sum(is.na(datos_mice)) > 0) {
    mice_model <- mice(datos_mice, m = 1, maxit = 5, method = "pmm", seed = 500, printFlag = FALSE)
    sample_imp[, num_impute_vars] <- mice::complete(mice_model, 1)  }
}

# CATEGÓRICAS -> MODA
mode_value <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

cat_impute_vars <- intersect(c(
  "Street", "Airport_Code", "Wind_Group", "Weather_Group", "Sunrise_Sunset"
), names(sample_imp))

for (v in cat_impute_vars) {
  moda <- mode_value(sample_imp[[v]])
  if (!is.na(moda)) {
    sample_imp[[v]][is.na(sample_imp[[v]]) | sample_imp[[v]] == ""] <- moda
  }
  sample_imp[[v]] <- as.factor(sample_imp[[v]])
}

write.csv(sample_imp,
          file.path(out_dir, "02_sample_15000_imputado.csv"),
          row.names = FALSE)

saveRDS(sample_imp, file.path(out_dir, "02_sample_15000_imputado.rds"))


# ------------------------------------------------------------------------------
# 7. DETECCIÓN Y TRATAMIENTO DE OUTLIERS (ADAPTADO AL CONTEXTO MULTIESTADO)
# ------------------------------------------------------------------------------
# Idea clave:
#   - Mantenemos la filosofía del script original.
#   - Pero NO tratamos Start_Lat/Start_Lng como anómalos, porque ahora la muestra
#     es multiestado y la dispersión geográfica es natural.
#   - Se detectan outliers estadísticos y solo se eliminan los que además tengan
#     "poco sentido físico" según reglas sencillas.

sample_feat <- sample_imp

# Variables numéricas a revisar para outliers
vars_outliers_num <- intersect(c(
  "Temperature.C.", "Wind_Chill.C.", "Humidity...", "Pressure.in.", 
  "Wind_Speed.mph.", "Precipitation.in.", "Duration_Hours"
), names(sample_feat))

# Resumen min/max
sink(file.path(out_dir, "resumen_preprocesado.txt"))
cat("==============================\n")
cat("RESUMEN PREPROCESADO\n")
cat("==============================\n\n")
cat("Dimensión base original: ", paste(dim(dd0), collapse = " x "), "\n")
cat("Dimensión muestra raw: ", paste(dim(sample_raw), collapse = " x "), "\n")
cat("Dimensión muestra imputada: ", paste(dim(sample_imp), collapse = " x "), "\n\n")

cat("Estados incluidos en la muestra (100 por estado):\n")
print(table(sample_raw$State))
cat("\n")

cat("Missings antes de imputar:\n")
print(sort(colSums(is.na(sample_raw)), decreasing = TRUE))
cat("\nMissings después de imputar:\n")
print(sort(colSums(is.na(sample_imp)), decreasing = TRUE))
cat("\n")

cat("Resumen de máximos y mínimos (variables numéricas outliers):\n")
for (nm in vars_outliers_num) {
  x <- sample_feat[[nm]]
  cat("Variable:", nm,
      " | Min:", min(x, na.rm = TRUE),
      " | Max:", max(x, na.rm = TRUE), "\n")
}
cat("\n")

# -----------------------------
# 7.1 Detección univariante
# -----------------------------
IQROutlier <- function(variable) {
  IQ <- IQR(variable, na.rm = TRUE)
  Q1 <- quantile(variable, 0.25, na.rm = TRUE)
  Q3 <- quantile(variable, 0.75, na.rm = TRUE)
  inf <- Q1 - 1.5 * IQ
  sup <- Q3 + 1.5 * IQ
  which(variable < inf | variable > sup)
}

outliers_iqr <- lapply(sample_feat[, vars_outliers_num, drop = FALSE], IQROutlier)
outliers_zscore <- lapply(sample_feat[, vars_outliers_num, drop = FALSE], function(x) {
  z <- scale(x)
  which(abs(z) > 3)
})
outliers_hampel <- lapply(sample_feat[, vars_outliers_num, drop = FALSE], function(x) {
  mediana <- median(x, na.rm = TRUE)
  desv_mad <- mad(x, na.rm = TRUE)
  lower <- mediana - 3 * desv_mad
  upper <- mediana + 3 * desv_mad
  which(x < lower | x > upper)
})

cat("Outliers IQR por variable:\n")
print(sapply(outliers_iqr, length))
cat("\nOutliers Z-Score por variable:\n")
print(sapply(outliers_zscore, length))
cat("\nOutliers Hampel por variable:\n")
print(sapply(outliers_hampel, length))
cat("\n")

# -----------------------------
# 7.2 Detección multivariante
# -----------------------------
sub <- sample_feat[, vars_outliers_num, drop = FALSE]
sub <- na.omit(sub)

# Mahalanobis clásica
outliers_classicmaha <- integer(0)
if (nrow(sub) > ncol(sub) + 1) {
  d_maha <- mahalanobis(sub, colMeans(sub), cov(sub))
  cutoff <- qchisq(0.99, df = ncol(sub))
  outliers_classicmaha <- which(d_maha > cutoff)
  cat("Outliers Mahalanobis clásica:", length(outliers_classicmaha), "\n")
}

# Mahalanobis robusta (bloque clima, fiel al guion previo)
vars_mah <- intersect(c("Temperature.C.", "Wind_Chill.C.", "Humidity...", "Pressure.in."), names(sample_feat))
outliers_mahalanobis <- integer(0)
if (length(vars_mah) >= 3) {
  datos_mah <- sample_feat[, vars_mah, drop = FALSE]
  datos_mah <- na.omit(datos_mah)
  if (nrow(datos_mah) > 10) {
    resultados_mah <- chemometrics::Moutlier(datos_mah, quantile = 0.99, plot = FALSE)
    grados_libertad <- ncol(datos_mah)
    limite_teorico <- sqrt(qchisq(0.99, df = grados_libertad))
    outliers_mahalanobis <- which(resultados_mah$rd > limite_teorico)
    cat("Outliers Mahalanobis robusta:", length(outliers_mahalanobis), "\n")
  }
}

# LOF
lof_scores <- Rlof::lof(scale(sub), k = 5)
outliers_lof_idx <- which(lof_scores > 2)
cat("Outliers LOF:", length(outliers_lof_idx), "\n")

# KNN (misma filosofía)
outliers_knn <- tryCatch({
  datos_knn <- scale(sub)
  adamethods::do_knno(datos_knn, k = 1, top_n = min(200, nrow(sub)))
}, error = function(e) integer(0))
cat("Outliers KNN:", length(outliers_knn), "\n")

# Isolation Forest
predicciones <- NULL
outliers_if_idx <- integer(0)
try({
  datos_if <- sub
  isoforest <- isolationForest$new(
    sample_size = as.integer(nrow(datos_if) / 2),
    num_trees   = 500,
    replace     = TRUE,
    seed        = 123
  )
  isoforest$fit(dataset = datos_if)
  predicciones <- isoforest$predict(data = datos_if)

  # Umbral adaptativo: percentil 10 de average_depth
  umbral_if <- quantile(predicciones$average_depth, 0.10, na.rm = TRUE)
  outliers_if_idx <- which(predicciones$average_depth <= umbral_if)
  cat("Outliers Isolation Forest:", length(outliers_if_idx), "\n")
}, silent = TRUE)

cat("\n")

# -----------------------------
# 7.3 Votación por consenso
# -----------------------------
votos_outliers <- tibble(ID_Fila = 1:nrow(sub)) %>%
  mutate(
    Maha_Clasica = if_else(ID_Fila %in% outliers_classicmaha, 1L, 0L),
    Maha_Robusta = if_else(ID_Fila %in% outliers_mahalanobis, 1L, 0L),
    KNN          = if_else(ID_Fila %in% outliers_knn, 1L, 0L),
    LOF          = if_else(ID_Fila %in% outliers_lof_idx, 1L, 0L),
    IsoForest    = if_else(ID_Fila %in% outliers_if_idx, 1L, 0L)
  ) %>%
  mutate(Total_Votos = Maha_Clasica + Maha_Robusta + KNN + LOF + IsoForest)

cat("Distribución de votos outlier:\n")
print(table(votos_outliers$Total_Votos))
cat("\n")

outliers_consenso_sub <- sub %>%
  mutate(ID_Fila = 1:nrow(sub), Total_Votos = votos_outliers$Total_Votos) %>%
  filter(Total_Votos >= 2) %>%
  arrange(desc(Total_Votos))

cat("Total de candidatos por consenso (>=2 votos):", nrow(outliers_consenso_sub), "\n\n")

# -----------------------------
# 7.4 Reglas físicas / de negocio para eliminar solo anomalías claras
# -----------------------------
sample_feat2 <- sample_feat %>% mutate(ID_ROW = row_number())
sub_ids <- sample_feat2 %>%
  select(ID_ROW, all_of(vars_outliers_num)) %>%
  na.omit()

votos_full <- sub_ids %>%
  select(ID_ROW) %>%
  mutate(Total_Votos = votos_outliers$Total_Votos)

outliers_consenso_full <- sample_feat2 %>%
  left_join(votos_full, by = "ID_ROW") %>%
  mutate(Total_Votos = ifelse(is.na(Total_Votos), 0, Total_Votos)) %>%
  filter(Total_Votos >= 2)

# CORRECCIÓN: Eliminado Distance.mi del filtro ya que fue purgado antes
accidentes_anomalos <- outliers_consenso_full %>%
  filter(
    Wind_Chill.C. > Temperature.C. |
      (Temperature.C. - Wind_Chill.C.) > 16.7 |
      Precipitation.in. > 5 |
      Duration_Hours > 48
  )

cat("Anomalías finales eliminadas por criterio físico/negocio:", nrow(accidentes_anomalos), "\n")
cat("\nDimensión final tras limpieza de outliers: ",
    nrow(sample_feat2) - nrow(accidentes_anomalos), " x ", ncol(sample_feat2), "\n", sep = "")

sink()

sample_clean <- anti_join(sample_feat2, accidentes_anomalos, by = "ID_ROW")

write.csv(sample_clean,
          file.path(out_dir, "03_sample_15000_clean.csv"),
          row.names = FALSE)

saveRDS(sample_clean, file.path(out_dir, "03_sample_15000_clean.rds"))

# ------------------------------------------------------------------------------
# 8. SELECCIÓN FINAL DE VARIABLES PARA CLUSTERING
# ------------------------------------------------------------------------------
vars_cluster <- intersect(c(
  # Numéricas
  "Severity", "Temperature.C.", "Wind_Chill.C.",
  "Humidity...", "Pressure.in.", "Wind_Speed.mph.",
  "Precipitation.in.", "Duration_Hours", "Start_Hour", "Start_Month",
  "Start_Lat", "Start_Lng",
  # Categóricas
  "State", "County", "Weather_Group", "Wind_Group",
  "Sunrise_Sunset", "Is_Weekend", 
  "Traffic_Signal", "Junction", "Has_Infrastructure"
), names(sample_clean))

base_clustering <- sample_clean[, vars_cluster, drop = FALSE]

# CORRECCIÓN: Lista limpia, sin variables borradas ni objetos fantasma
cat_cluster <- c(
  "State", "County", "Weather_Group", "Wind_Group",
  "Sunrise_Sunset", "Is_Weekend", "Traffic_Signal", "Junction", "Has_Infrastructure"
)

cat_cluster <- intersect(cat_cluster, names(base_clustering))
for (v in cat_cluster) {
  base_clustering[[v]] <- as.factor(base_clustering[[v]])
}

# Eliminamos columnas que hayan quedado constantes
niveles_utiles <- sapply(base_clustering, function(x) length(unique(na.omit(x))))
base_clustering <- base_clustering[, niveles_utiles > 1, drop = FALSE]

write.csv(base_clustering,
          file.path(out_dir, "04_base_clustering.csv"),
          row.names = FALSE)

saveRDS(base_clustering, file.path(out_dir, "04_base_clustering.rds"))

# ------------------------------------------------------------------------------
# 9. MENSAJES FINALES
# ------------------------------------------------------------------------------
cat("\n============================================\n")
cat("PREPROCESADO COMPLETADO CORRECTAMENTE\n")
cat("============================================\n")
cat("Archivos generados en:", normalizePath(out_dir), "\n")
cat("- 01_sample_15000_raw.csv\n")
cat("- 02_sample_15000_imputado.csv\n")
cat("- 03_sample_15000_clean.csv\n")
cat("- 04_base_clustering.csv\n")
cat("- resumen_preprocesado.txt\n")
cat("\nLa base 04_base_clustering.csv queda lista para aplicar k-proto y CURE.\n")

