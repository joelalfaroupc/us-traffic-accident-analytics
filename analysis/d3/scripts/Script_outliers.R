library(tidyverse)
library(EnvStats)
library(chemometrics)
library(Rlof)
library(solitude)
library(adamethods)


data <- readRDS("Datasets/US_Accidents_Final_imputado_mice.rds")

tipus <- sapply(data, class)
varNum <- names(tipus)[tipus %in% c("integer","numeric")]

# DETECCIÓN UNIVARIANTE

# Máximos y mínimos

invisible(mapply(function(x, name) {
  cat("Variable:", name, "\n  Min:", min(x, na.rm=TRUE),
      "\n  Max:", max(x, na.rm=TRUE), "\n\n")
}, data[, varNum], varNum))

# IQR

IQROutlier <- function(variable) {
  IQ <- iqr(variable, na.rm = TRUE)
  Q1 <- quantile(variable, 0.25, na.rm=TRUE)
  Q3 <- quantile(variable, 0.75, na.rm=TRUE)
  inf <- Q1 - 1.5*IQ
  sup <- Q3 + 1.5*IQ
  which(variable < inf | variable > sup)
}

outliers_iqr <- lapply(data[, varNum, drop = FALSE], IQROutlier)

sapply(outliers_iqr, length)


# Boxplot

par(mfrow = c(2, 5)) 
for (columna in varNum) {
  boxplot(data[[columna]], 
          main = columna, 
          col = "skyblue", 
          border = "black")
}
par(mfrow = c(1, 1))

# Z-Score

outliers_zscore <- lapply(data[, varNum, drop = FALSE], function(x) {
  z <- scale(x)
  which(abs(z) > 3)
})
sapply(outliers_zscore, length)

# Hampel Identifier

outliers_hampel <- lapply(data[, varNum, drop = FALSE], function(x) {
  mediana <- median(x)
  desv_mad <- mad(x)
  lower <- mediana - 3 * desv_mad
  upper <- mediana + 3 * desv_mad
  
  which(x < lower | x > upper)
})
sapply(outliers_hampel, length)


# DETECCIÓN MULTIVARIANTE

sub <- data[, varNum]
sub <- na.omit(sub)

# Mahalanobis clásica

d_maha <- mahalanobis(sub, colMeans(sub), cov(sub))
cutoff <- qchisq(0.99, df=ncol(sub))
outliers_classicmaha <- which(d_maha > cutoff)
length(outliers_classicmaha)

par(mfrow= c(1,1))
hist(d_maha, main="Distancias Mahalanobis", col="skyblue")

# Mahalanobis robusta

vars_mah <- c("Temperature.F.", "Wind_Chill.F.", "Humidity...", "Pressure.in.")
datos_mah <- data[, vars_mah]

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

datos_if <- data %>% select(all_of(varNum))

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
data_evaluada <- data %>% bind_cols(predicciones)

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

votos_outliers <- tibble(ID_Fila = 1:nrow(data))

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
outliers_consenso <- data %>%
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
      (Temperature.F. - Wind_Chill.F.) > 30 |
      
      # Lluvias masivas
      Precipitation.in. > 5
  ) %>%
  # Ordenamos para ver lo más grave primero y seleccionamos las variables clave
  arrange(desc(Distance.mi.)) %>%
  select(Distance.mi., Temperature.F., Wind_Chill.F., Precipitation.in., Severity, Start_Time, Description)

cat("Num. total de anomalías detectadas:", nrow(accidentes_anomalos), "\n\n")

# Borrar los outliers sin sentido físico, mantener los outliers que sí tengan sentido

data_final <- anti_join(data, accidentes_anomalos, by = c("Start_Time", "Description"))

# Exportar
write.csv(data_final, "US_Accidents_Final_Outliers.csv", row.names = FALSE)
saveRDS(data_final, "US_Accidents_Final_Outliers.rds")
