# ==============================================================================
# PARTE V: CLUSTERING DE LA TS Y PROFILING ESTADÍSTICO
# ==============================================================================

# Cargamos las librerías necesarias
list.of.packages <- c(
  "cluster", "dtw", "dtwclust", "dendextend", "ggpubr",
  "lubridate", "proxy", "psych", "reshape2", "tidyverse"
)
invisible(lapply(list.of.packages, require, character.only = TRUE))

# Definimos los parámetros base
input_path <- "US_Accidents_Final_FE.csv" 
k_values <- 2:8
k_final <- 6
k_visual_compare <- 2:6
min_months_with_data <- 6
seed_value <- 3247L

# Funciones auxiliares de clustering y evaluación
safe_scale <- function(x) as.numeric(scale(x))

silhouette_mean <- function(clusters, dist_mat) {
  sil <- cluster::silhouette(as.integer(clusters), dist_mat)
  mean(sil[, "sil_width"])
}

cluster_balance <- function(clusters) {
  sizes <- as.numeric(table(clusters))
  min(sizes) / max(sizes)
}

evaluate_hclust <- function(dist_mat, linkage_method, method_label, ks) {
  hc <- hclust(dist_mat, method = linkage_method)
  metrics <- purrr::map_dfr(ks, function(k) {
    clusters <- cutree(hc, k = k)
    tibble(
      method = method_label, family = "hierarchical", k = k,
      silhouette = silhouette_mean(clusters, dist_mat),
      balance = cluster_balance(clusters),
      cophenetic_corr = cor(as.numeric(dist_mat), as.numeric(cophenetic(hc))),
      min_cluster_size = min(table(clusters)), max_cluster_size = max(table(clusters))
    )
  })
  clusters_by_k <- purrr::map(ks, ~ cutree(hc, k = .x))
  names(clusters_by_k) <- as.character(ks)
  list(model = hc, metrics = metrics, clusters = clusters_by_k)
}

evaluate_partitional <- function(series_mat, dist_mat, ks, seed_value) {
  models <- vector("list", length(ks))
  names(models) <- as.character(ks)
  metrics <- purrr::map_dfr(ks, function(k) {
    pc <- tsclust(series_mat, type = "partitional", k = k, distance = "dtw_basic", 
                  centroid = "pam", seed = seed_value, trace = FALSE, 
                  args = tsclust_args(dist = list(window.size = 20L)))
    models[[as.character(k)]] <<- pc
    clusters <- pc@cluster
    tibble(
      method = "dtw_partitional_pam", family = "partitional", k = k,
      silhouette = silhouette_mean(clusters, dist_mat), balance = cluster_balance(clusters),
      cophenetic_corr = NA_real_, min_cluster_size = min(table(clusters)), 
      max_cluster_size = max(table(clusters))
    )
  })
  clusters_by_k <- purrr::map(models, ~ .x@cluster)
  list(model = models, metrics = metrics, clusters = clusters_by_k)
}

# Funciones auxiliares para el Profiling
safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  shapiro.test(x)$p.value
}

build_series_features <- function(real_matrix) {
  idx <- seq_len(ncol(real_matrix))
  tibble(
    County = rownames(real_matrix),
    total_accidents = rowSums(real_matrix),
    mean_monthly = rowMeans(real_matrix),
    sd_monthly = apply(real_matrix, 1, sd),
    max_monthly = apply(real_matrix, 1, max),
    active_months = rowSums(real_matrix > 0),
    pct_active_months = rowSums(real_matrix > 0) / ncol(real_matrix),
    mean_nonzero = apply(real_matrix, 1, function(x) mean(x[x > 0])),
    peak_to_mean_ratio = apply(real_matrix, 1, function(x) max(x) / mean(x)),
    trend_slope = apply(real_matrix, 1, function(x) coef(lm(x ~ idx))[2])
  )
}

plot_profiles_scaled <- function(long_data, title_text) {
  resumen <- long_data %>% group_by(Cluster, Periodo, Periodo_fecha) %>% summarise(media = mean(Valor), .groups = "drop")
  p <- ggplot(long_data, aes(x = Periodo_fecha, y = Valor, group = County)) +
    geom_line(alpha = 0.20, color = "grey55", linewidth = 0.30) +
    geom_line(data = resumen, aes(x = Periodo_fecha, y = media, group = 1), linewidth = 1.10, color = "#d7301f", inherit.aes = FALSE) +
    facet_wrap(~ Cluster, scales = "free_y") + theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title_text, subtitle = "Lineas grises: series individuales. Linea roja: perfil cluster", x = "Periodo", y = "Z-score")
  print(p)
}

plot_profiles_real <- function(long_data, title_text) {
  resumen <- long_data %>% group_by(Cluster, Periodo, Periodo_fecha) %>% summarise(media = mean(Valor), .groups = "drop")
  p <- ggplot(resumen, aes(x = Periodo_fecha, y = media, color = Cluster, group = Cluster)) +
    geom_line(linewidth = 1) + geom_point(size = 1.2) + facet_wrap(~ Cluster, scales = "free_y") +
    theme_minimal(base_size = 11) + theme(legend.position = "none", strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title_text, subtitle = "Promedio mensual real", x = "Periodo", y = "Accidentes")
  print(p)
}

# 1. Carga y limpieza directa de datos
dad <- read.csv(input_path, stringsAsFactors = FALSE)

dad$Start_Time <- lubridate::parse_date_time(
  dad$Start_Time,
  orders = c("Y-m-d H:M:S", "Y-m-d H:M:SOS", "Y/m/d H:M:S")
)

df_agrupado <- dad %>%
  filter(!is.na(Start_Time), !is.na(County), County != "") %>%
  mutate(Periodo = format(Start_Time, "%Y-%m")) %>%
  count(County, Periodo, name = "Total")

datos_wide <- dcast(df_agrupado, County ~ Periodo, value.var = "Total", fill = 0)
rownames(datos_wide) <- datos_wide[, "County"]
datos_wide[, "County"] <- NULL

columnas_vacias <- colSums(datos_wide) == 0
datos_wide <- datos_wide[, !columnas_vacias, drop = FALSE]

filas_validas <- rowSums(datos_wide > 0) >= min_months_with_data
datos_wide <- datos_wide[filas_validas, order(colnames(datos_wide)), drop = FALSE]
datos_reales <- as.matrix(datos_wide)

# 2. Escalado
datos <- t(apply(datos_reales, 1, safe_scale))
rownames(datos) <- rownames(datos_reales)
colnames(datos) <- colnames(datos_reales)

datos <- datos[complete.cases(datos), , drop = FALSE]
datos_reales <- datos_reales[rownames(datos), , drop = FALSE]

# 3. Distancia y Clustering Base
dist_dtw <- proxy::dist(datos, method = "DTW")

eval_dtw_ward <- evaluate_hclust(dist_dtw, "ward.D2", "dtw_ward", k_values)
eval_partitional <- evaluate_partitional(datos, dist_dtw, k_values, seed_value)

metricas_clustering <- bind_rows(
  eval_dtw_ward$metrics,
  eval_partitional$metrics
) %>%
  arrange(k, method)

print(metricas_clustering)

# Tabla resumida de silhouettes
tabla_silhouette <- metricas_clustering %>%
  select(method, family, k, silhouette, balance, min_cluster_size, max_cluster_size)

print(tabla_silhouette)

# Si quieres ver solo silhouette por método
tabla_silhouette %>%
  select(method, k, silhouette) %>%
  arrange(method, k) %>%
  print()

# Gráfico silhouette vs k
ggplot(tabla_silhouette, aes(x = k, y = silhouette, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Comparación de silhouette por número de clusters",
    x = "Número de clusters (k)",
    y = "Silhouette media",
    color = "Método"
  )

# 4. LOOP DE VISUALIZACIÓN DE K
for(k_actual in k_visual_compare) {
  # Dendrograma
  dend <- eval_dtw_ward$model %>% as.dendrogram()
  plot(dend %>% set("branches_k_color", k = k_actual) %>% set("labels_cex", 0.7),
       main = paste("Clustering DTW de Condados (K =", k_actual, ")"))
  
  # Partitional
  plot(eval_partitional$model[[as.character(k_actual)]], main = paste("Partitional (K =", k_actual, ")"))
  
  # Hierarchical
  hc_loop <- tsclust(datos, type = "hierarchical", k = k_actual, distance = "dtw_basic", trace = FALSE, control = hierarchical_control(method = "ward.D2"))
  plot(hc_loop, type = "sc", main = paste("Hierarchical Solapadas (K =", k_actual, ")")) 
}

# ==============================================================================
# 5. MODELO FINAL (K=6) Y PROFILING VISUAL
# ==============================================================================
clusters_final <- factor(eval_dtw_ward$clusters[[as.character(k_final)]])

datos_long_scaled <- as.data.frame(datos) %>% rownames_to_column(var = "County") %>% pivot_longer(cols = -County, names_to = "Periodo", values_to = "Valor") %>% mutate(Cluster = factor(clusters_final[County], levels = levels(clusters_final)), Periodo_fecha = as.Date(paste0(Periodo, "-01")))
datos_long_real <- as.data.frame(datos_reales) %>% rownames_to_column(var = "County") %>% pivot_longer(cols = -County, names_to = "Periodo", values_to = "Valor") %>% mutate(Cluster = factor(clusters_final[County], levels = levels(clusters_final)), Periodo_fecha = as.Date(paste0(Periodo, "-01")))

plot_profiles_scaled(datos_long_scaled, title_text = paste("Perfiles temporales escalados (k =", k_final, ")"))
plot_profiles_real(datos_long_real, title_text = paste("Evolucion media real por cluster (k =", k_final, ")"))

# ==============================================================================
# 6. PROFILING ESTADÍSTICO DE VARIABLES DERIVADAS DE LA SERIE
# ==============================================================================
cluster_features <- build_series_features(datos_reales) %>%
  mutate(Cluster = factor(clusters_final[County], levels = levels(clusters_final)))

numeric_features <- setdiff(colnames(cluster_features), c("County", "Cluster"))
significant_features <- character(0)
test_used <- list()

for (feature_name in numeric_features) {
  values <- cluster_features[[feature_name]]
  shapiro_by_cluster <- tapply(values, cluster_features$Cluster, safe_shapiro)
  
  # Gráfico QQ por Cluster
  qq_plot <- ggplot(cluster_features, aes(sample = .data[[feature_name]])) +
    geom_qq() + geom_qq_line() + facet_wrap(~ Cluster, scales = "free_y") +
    theme_minimal(base_size = 11) + labs(title = paste("QQ-plot:", feature_name))
  print(qq_plot)
  
  # Selección del test en base a normalidad (Shapiro > 0.05)
  normal_groups <- all(shapiro_by_cluster > 0.05, na.rm = TRUE)
  
  if (normal_groups) {
    model_aov <- aov(stats::reformulate("Cluster", response = feature_name), data = cluster_features)
    p_value <- summary(model_aov)[[1]][["Pr(>F)"]][1]
    if (p_value <= 0.05) {
      significant_features <- c(significant_features, feature_name)
      test_used[[feature_name]] <- "anova"
    }
  } else {
    kw_test <- kruskal.test(stats::reformulate("Cluster", response = feature_name), data = cluster_features)
    if (kw_test$p.value <= 0.05) {
      significant_features <- c(significant_features, feature_name)
      test_used[[feature_name]] <- "kruskal"
    }
  }
  
  # Visualización de la distribución de la variable (Boxplot + Histograma)
  gr_boxplot <- ggpubr::ggboxplot(cluster_features, "Cluster", feature_name, fill = "Cluster") +
    labs(title = paste("Boxplot:", feature_name))
  gr_hist <- ggpubr::gghistogram(cluster_features, x = feature_name, bins = 10, add = "mean", rug = TRUE, color = "Cluster", fill = "Cluster") +
    labs(title = paste("Hist:", feature_name))
  
  print(ggpubr::ggarrange(gr_boxplot, gr_hist, heights = c(2, 0.9), ncol = 2, nrow = 1, align = "v"))
}

# ==============================================================================
# 7. POST-HOC Y RESÚMENES DE VARIABLES SIGNIFICATIVAS
# ==============================================================================
if (length(significant_features) > 0) {
  print(psych::describeBy(cluster_features[, significant_features, drop = FALSE], cluster_features$Cluster))
  
  for (feature_name in significant_features) {
    if (test_used[[feature_name]] == "anova") {
      model_aov <- aov(stats::reformulate("Cluster", response = feature_name), data = cluster_features)
      print(TukeyHSD(model_aov))
    } else {
      print(pairwise.wilcox.test(cluster_features[[feature_name]], cluster_features$Cluster, exact = FALSE, p.adjust.method = "bonferroni"))
    }
  }
}

# ==============================================================================
# 8. EXPORTAMOS RESULTADOS FINALES
# ==============================================================================
resultados <- data.frame(County = rownames(datos), Cluster = as.integer(clusters_final))
resultados_completos <- cbind(resultados, datos)

write.csv(resultados_completos, paste0("resultados_condados_clusters_k", k_final, ".csv"), row.names = FALSE)
write.csv(t(resultados_completos), paste0("resultados_condados_clusters_t_k", k_final, ".csv"))