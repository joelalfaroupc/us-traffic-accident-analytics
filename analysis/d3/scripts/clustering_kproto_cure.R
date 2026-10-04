
# ==============================================================================
# [PMAAD] - Clustering final ROBUSTO: k-proto + CURE adaptado a datos mixtos
# ------------------------------------------------------------------------------
# Qué mejora respecto a la versión anterior:
#   - Validaciones y coerciones más robustas
#   - Selección final menos "ciega": además de silhouette mira balance de clusters
#   - Más trazas en consola
#   - Plots que se GUARDAN y además se MUESTRAN por pantalla
#   - PCA numérico coloreado por clusters
#   - Mapa aproximado Lat/Lng coloreado por clusters
#   - Barras de tamaño de cluster para comparar interpretabilidad
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# 0. Configuración
# ------------------------------------------------------------------------------
file_base  <- "Datasets/clustering/04_base_clustering.csv"
file_clean <- "Datasets/clustering/03_sample_15000_clean.csv"
dir_out    <- "clustering_results/"

k_grid <- 2:6
prefer_interpretable_if_close <- TRUE
sil_close_threshold <- 0.01
min_cluster_prop_threshold <- 0.08
cure_sample_frac <- 0.35
cure_min_state_n <- 10
cure_n_rep <- 5
cure_h_method <- "average"

# ------------------------------------------------------------------------------
# 1. Librerías
# ------------------------------------------------------------------------------
packages <- c("dplyr", "readr", "cluster", "clustMixType", "ggplot2", "gower")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new.packages) > 0) install.packages(new.packages, dependencies = TRUE)
invisible(lapply(packages, require, character.only = TRUE))

# ------------------------------------------------------------------------------
# 2. Validaciones de rutas
# ------------------------------------------------------------------------------
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
if (!file.exists(file_base))  stop("No se encuentra 04_base_clustering.csv. Ejecuta antes el script de refinamiento.")
if (!file.exists(file_clean)) stop("No se encuentra 03_sample_1500_clean.csv")

# ------------------------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------------------------
safe_as_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

avg_silhouette <- function(cluster_vec, diss_matrix) {
  cluster_vec <- as.integer(as.factor(cluster_vec))
  if (length(unique(cluster_vec)) < 2) return(NA_real_)
  sil <- cluster::silhouette(cluster_vec, diss_matrix)
  mean(sil[, 3], na.rm = TRUE)
}

cluster_balance_metrics <- function(cluster_vec) {
  tb <- table(cluster_vec)
  props <- as.numeric(tb) / sum(tb)
  data.frame(
    n_clusters = length(tb),
    min_size = min(as.numeric(tb)),
    max_size = max(as.numeric(tb)),
    min_prop = min(props),
    max_prop = max(props),
    imbalance_ratio = max(as.numeric(tb)) / min(as.numeric(tb)),
    stringsAsFactors = FALSE
  )
}

sample_ids_by_state <- function(states, frac = 0.35, min_n = 10) {
  ids_all <- split(seq_along(states), states)
  ids_sample <- unlist(lapply(ids_all, function(ix) {
    n_take <- max(min_n, ceiling(length(ix) * frac))
    n_take <- min(n_take, length(ix))
    sample(ix, size = n_take, replace = FALSE)
  }))
  sort(as.integer(ids_sample))
}

select_rep_indices <- function(df_cluster, n_rep = 5) {
  n <- nrow(df_cluster)
  if (n <= 1) return(1)
  
  d <- as.matrix(cluster::daisy(df_cluster, metric = "gower"))
  medoid <- which.min(rowMeans(d))
  selected <- medoid
  
  while (length(selected) < min(n_rep, n)) {
    remaining <- setdiff(seq_len(n), selected)
    scores <- sapply(remaining, function(i) min(d[i, selected]))
    next_i <- remaining[which.max(scores)]
    selected <- c(selected, next_i)
  }
  selected
}

run_cure_mixed <- function(data_all, k, ids_sample, n_rep = 5, h_method = "average") {
  data_sample <- data_all[ids_sample, , drop = FALSE]
  d_sample <- cluster::daisy(data_sample, metric = "gower")
  hc <- hclust(d_sample, method = h_method)
  sample_clusters <- cutree(hc, k = k)
  
  rep_global_ids <- c()
  rep_cluster_labels <- c()
  
  for (cl in sort(unique(sample_clusters))) {
    local_ids <- which(sample_clusters == cl)
    reps_local <- select_rep_indices(data_sample[local_ids, , drop = FALSE], n_rep = n_rep)
    reps_global_cluster <- ids_sample[local_ids[reps_local]]
    rep_global_ids <- c(rep_global_ids, reps_global_cluster)
    rep_cluster_labels <- c(rep_cluster_labels, rep(cl, length(reps_global_cluster)))
  }
  
  full_clusters <- integer(nrow(data_all))
  full_clusters[ids_sample] <- sample_clusters
  
  ids_rest <- setdiff(seq_len(nrow(data_all)), ids_sample)
  if (length(ids_rest) > 0) {
    # Usamos gower_topn para buscar el índice del representante más cercano de forma eficiente
    match_res <- gower::gower_topn(
      x = data_all[ids_rest, , drop = FALSE],
      y = data_all[rep_global_ids, , drop = FALSE],
      n = 1
    )
    
    # match_res$index nos da qué fila de 'y' (representantes) está más cerca de cada fila de 'x'
    nearest_rep <- match_res$index
    full_clusters[ids_rest] <- rep_cluster_labels[nearest_rep]
  }
  
  list(
    cluster = full_clusters,
    representatives = rep_global_ids,
    rep_cluster_labels = rep_cluster_labels
  )
}

save_and_print_plot <- function(plot_obj, filename, width = 7, height = 5) {
  print(plot_obj)
  ggsave(filename = filename, plot = plot_obj, width = width, height = height)
}

# ------------------------------------------------------------------------------
# 4. Carga de datos
# ------------------------------------------------------------------------------
base_final <- read.csv(file_base, stringsAsFactors = FALSE)
base_clean <- read.csv(file_clean, stringsAsFactors = FALSE)

cat("\n==============================\n")
cat("CARGA DE DATOS\n")
cat("==============================\n")
cat("Base final clustering:", nrow(base_final), "x", ncol(base_final), "\n")
cat("Base limpia:", nrow(base_clean), "x", ncol(base_clean), "\n")


# ------------------------------------------------------------------------------
# 5. Preparación de tipos
# ------------------------------------------------------------------------------

# Extraer categóricas explícitas y cualquier otra que sea lógicamente un factor
vars_factor <- intersect(c(
  "Severity", "State", "County", "Weather_Group", "Wind_Group", 
  "Sunrise_Sunset", "Is_Weekend", "Traffic_Signal", "Junction", 
  "Has_Infrastructure", "Start_Hour", "Start_Month"
), names(base_final))

vars_num <- setdiff(names(base_final), vars_factor)

base_final[vars_factor] <- lapply(base_final[vars_factor], factor)
base_final[vars_num]    <- lapply(base_final[vars_num], safe_as_numeric)

if (sum(is.na(base_final)) > 0) {
  na_by_col <- sort(colSums(is.na(base_final))[colSums(is.na(base_final)) > 0], decreasing = TRUE)
  print(na_by_col)
  stop("La base final todavía contiene NA. Hay que revisarlo antes del clustering.")
}

zero_var_num <- vars_num[sapply(base_final[vars_num], function(x) stats::sd(x, na.rm = TRUE) == 0)]
if (length(zero_var_num) > 0) {
  cat("\nSe eliminan numéricas de varianza cero:\n")
  print(zero_var_num)
  base_final <- base_final[, setdiff(names(base_final), zero_var_num), drop = FALSE]
  vars_num <- setdiff(vars_num, zero_var_num)
}

# Eliminar factores con un solo nivel
single_level_factors <- vars_factor[sapply(base_final[vars_factor], function(x) length(unique(na.omit(x))) <= 1)]
if (length(single_level_factors) > 0) {
  cat("\nSe eliminan categóricas de un solo nivel:\n")
  print(single_level_factors)
  base_final <- base_final[, setdiff(names(base_final), single_level_factors), drop = FALSE]
  vars_factor <- setdiff(vars_factor, single_level_factors)
}

gower_full <- cluster::daisy(base_final, metric = "gower")

# ------------------------------------------------------------------------------
# 6. K-PROTO ROBUSTO
# ------------------------------------------------------------------------------
cat("\n==============================\n")
cat("K-PROTO\n")
cat("==============================\n")

kproto_list <- list()
kproto_models <- list()
iter <- 1

for (k in k_grid) {
  # Semilla limpia solo dependiente de k
  set.seed(100 + k)
  
  # Ejecutamos directamente con el lambda por defecto
  modelo <- clustMixType::kproto(base_final, k = k, nstart = 20, verbose = FALSE)
  
  clusters <- modelo$cluster
  sil <- avg_silhouette(clusters, gower_full)
  met <- cluster_balance_metrics(clusters)
  
  # Aplicamos el filtro duro
  sil_valid <- if(!is.na(met$min_prop) && met$min_prop >= min_cluster_prop_threshold) sil else NA_real_
  
  kproto_list[[iter]] <- data.frame(
    k = k,
    lambda = NA, # Se mantiene como NA para no romper el código posterior
    silhouette = sil,
    silhouette_valid = sil_valid,
    tot_within = modelo$tot.withinss,
    min_size = met$min_size,
    max_size = met$max_size,
    min_prop = met$min_prop,
    imbalance_ratio = met$imbalance_ratio,
    stringsAsFactors = FALSE
  )
  
  key <- paste0("k", k, "_lam_default")
  kproto_models[[key]] <- modelo
  
  cat("k =", k,
      "| lambda = default",
      "| silhouette =", round(sil, 4),
      "| min_prop =", round(met$min_prop, 3),
      "| imbalance =", round(met$imbalance_ratio, 2), 
      "| Valid:", !is.na(sil_valid), "\n")
  
  iter <- iter + 1
}

# Unimos los resultados
kproto_results <- do.call(rbind, kproto_list)

# Selección del mejor modelo
if (all(is.na(kproto_results$silhouette_valid))) {
  cat("\nADVERTENCIA: Ningún modelo de K-proto superó el umbral mínimo. Revisa tus parámetros.\n")
  best_row_kproto <- kproto_results[which.max(kproto_results$silhouette), ]
} else {
  best_row_kproto <- kproto_results[which.max(kproto_results$silhouette_valid), ]
}

best_key_kproto <- paste0("k", best_row_kproto$k, "_lam_default")
best_kproto <- kproto_models[[best_key_kproto]]
cluster_kproto <- best_kproto$cluster

# ------------------------------------------------------------------------------
# 7. CURE ADAPTADO ROBUSTO
# ------------------------------------------------------------------------------
cat("\n==============================\n")
cat("CURE ADAPTADO\n")
cat("==============================\n")

set.seed(999)
# Muestreo estratificado (Excelente práctica, la mantenemos intacta)
if ("State" %in% names(base_final)) {
  ids_sample <- sample_ids_by_state(base_final$State, frac = cure_sample_frac, min_n = cure_min_state_n)
} else {
  ids_sample <- sample(seq_len(nrow(base_final)),
                       size = ceiling(cure_sample_frac * nrow(base_final)),
                       replace = FALSE)
}

# Usamos listas en lugar de rbind iterativo
cure_list <- list()
cure_models <- list()
iter <- 1

for (k in k_grid) {
  set.seed(200 + k)

  # ATENCIÓN: Esto asume que ya actualizaste run_cure_mixed en la Sección 3 con gower_topn
  modelo_cure <- run_cure_mixed(base_final, k = k, ids_sample = ids_sample,
                                n_rep = cure_n_rep, h_method = cure_h_method)

  sil <- avg_silhouette(modelo_cure$cluster, gower_full)
  met <- cluster_balance_metrics(modelo_cure$cluster)

  # Filtro duro: invalidamos el modelo si el cluster es demasiado pequeño
  sil_valid <- if(!is.na(met$min_prop) && met$min_prop >= min_cluster_prop_threshold) sil else NA_real_

  cure_list[[iter]] <- data.frame(
    k = k,
    silhouette = sil,
    silhouette_valid = sil_valid, # Métrica de decisión
    n_representatives = length(modelo_cure$representatives),
    min_size = met$min_size,
    max_size = met$max_size,
    min_prop = met$min_prop,
    imbalance_ratio = met$imbalance_ratio,
    stringsAsFactors = FALSE
  )

  cure_models[[paste0("k", k)]] <- modelo_cure

  cat("k =", k,
      "| silhouette =", round(sil, 4),
      "| min_prop =", round(met$min_prop, 3),
      "| imbalance =", round(met$imbalance_ratio, 2),
      "| Valid:", !is.na(sil_valid), "\n")

  iter <- iter + 1
}

# Unimos los resultados
cure_results <- do.call(rbind, cure_list)

# Selección del mejor modelo CURE
if (all(is.na(cure_results$silhouette_valid))) {
  cat("\nADVERTENCIA: Ningún modelo de CURE superó el umbral mínimo de tamaño de clúster.\n")
  best_row_cure <- cure_results[which.max(cure_results$silhouette), ]
} else {
  best_row_cure <- cure_results[which.max(cure_results$silhouette_valid), ]
}

best_k_cure <- best_row_cure$k
best_cure <- cure_models[[paste0("k", best_k_cure)]]
cluster_cure <- best_cure$cluster


# ------------------------------------------------------------------------------
# 8. Comparativa y selección final robusta
# ------------------------------------------------------------------------------
best_kproto_sil <- best_row_kproto$silhouette_valid
best_cure_sil   <- best_row_cure$silhouette_valid

# Castigamos fuertemente a los métodos que no superaron la validación de tamaño
# Si fallaron, su silueta válida se vuelve -1 para que nunca ganen la comparativa final.
if (is.na(best_kproto_sil)) best_kproto_sil <- -1
if (is.na(best_cure_sil)) best_cure_sil <- -1

best_kproto_bal <- cluster_balance_metrics(cluster_kproto)
best_cure_bal   <- cluster_balance_metrics(cluster_cure)

comparativa <- data.frame(
  metodo = c("kproto", "cure_adaptado"),
  best_k = c(best_row_kproto$k, best_row_cure$k),
  best_lambda = c(ifelse(is.na(best_row_kproto$lambda), "default", as.character(best_row_kproto$lambda)), NA),
  silhouette_valid = c(best_kproto_sil, best_cure_sil),
  min_prop = c(best_kproto_bal$min_prop, best_cure_bal$min_prop),
  imbalance_ratio = c(best_kproto_bal$imbalance_ratio, best_cure_bal$imbalance_ratio),
  stringsAsFactors = FALSE
)

# El ganador inicial es el que tenga mayor silhouette válido
metodo_final <- comparativa$metodo[which.max(comparativa$silhouette_valid)]

# Regla de desempate técnico (Excelente lógica, la mantenemos)
if (prefer_interpretable_if_close && abs(best_kproto_sil - best_cure_sil) <= sil_close_threshold) {
  if (best_kproto_bal$imbalance_ratio < best_cure_bal$imbalance_ratio) {
    metodo_final <- "kproto"
  } else {
    metodo_final <- "cure_adaptado"
  }
}

cluster_final <- if (metodo_final == "kproto") cluster_kproto else cluster_cure

cat("\n==============================\n")
cat("MÉTODO FINAL SELECCIONADO\n")
cat("==============================\n")
print(comparativa)
cat("\nMétodo final:", metodo_final, "\n")

# ------------------------------------------------------------------------------
# 9. Guardado de resultados
# ------------------------------------------------------------------------------
write.csv(kproto_results,
          file = file.path(dir_out, "06_kproto_resumen_robusto.csv"),
          row.names = FALSE)

write.csv(cure_results,
          file = file.path(dir_out, "06_cure_resumen_robusto.csv"),
          row.names = FALSE)

write.csv(comparativa,
          file = file.path(dir_out, "06_comparativa_metodos_robusta.csv"),
          row.names = FALSE)


base_clustered <- base_final
base_clustered$cluster_kproto <- cluster_kproto
base_clustered$cluster_cure   <- cluster_cure
base_clustered$cluster_final  <- cluster_final

write.csv(base_clustered,
          file = file.path(dir_out, "06_base_clustering_final_con_clusters_robusta.csv"),
          row.names = FALSE)

base_clean_clusters <- base_clean
base_clean_clusters$cluster_kproto <- cluster_kproto
base_clean_clusters$cluster_cure   <- cluster_cure
base_clean_clusters$cluster_final  <- cluster_final

write.csv(base_clean_clusters,
          file = file.path(dir_out, "06_base_limpia_con_clusters_robusta.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 10. Plots
# ------------------------------------------------------------------------------
p_kproto <- ggplot(kproto_results, aes(x = k, y = silhouette,
                                       group = ifelse(is.na(lambda), "default", as.character(lambda)),
                                       color = ifelse(is.na(lambda), "default", as.character(lambda)))) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  labs(title = "k-proto: silhouette media (Gower)",
       subtitle = "Comparación por k y lambda",
       x = "k", y = "Silhouette media", color = "lambda")

save_and_print_plot(p_kproto, file.path(dir_out, "06_kproto_silhouette_robusta.png"), 8, 5)

p_cure <- ggplot(cure_results, aes(x = k, y = silhouette)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  labs(title = "CURE adaptado: silhouette media (Gower)",
       subtitle = "Comparación por k",
       x = "k", y = "Silhouette media")

save_and_print_plot(p_cure, file.path(dir_out, "06_cure_silhouette_robusta.png"), 8, 5)

df_sizes <- rbind(
  data.frame(metodo = "kproto", cluster = factor(names(table(cluster_kproto))), n = as.numeric(table(cluster_kproto))),
  data.frame(metodo = "cure_adaptado", cluster = factor(names(table(cluster_cure))), n = as.numeric(table(cluster_cure))),
  data.frame(metodo = "final", cluster = factor(names(table(cluster_final))), n = as.numeric(table(cluster_final)))
)

p_sizes <- ggplot(df_sizes, aes(x = cluster, y = n, fill = cluster)) +
  geom_col() +
  facet_wrap(~ metodo, scales = "free_x") +
  theme_bw() +
  labs(title = "Tamaño de clusters por método",
       x = "Cluster", y = "Número de observaciones")

save_and_print_plot(p_sizes, file.path(dir_out, "06_cluster_sizes_robusta.png"), 9, 5)

if (length(vars_num) >= 2) {
  pca_obj <- prcomp(scale(base_final[vars_num]), center = TRUE, scale. = FALSE)
  
  df_pca <- data.frame(
    PC1 = pca_obj$x[, 1],
    PC2 = pca_obj$x[, 2],
    cluster_kproto = factor(cluster_kproto),
    cluster_cure = factor(cluster_cure),
    cluster_final = factor(cluster_final)
  )
  
  p_pca_k <- ggplot(df_pca, aes(x = PC1, y = PC2, color = cluster_kproto)) +
    geom_point(alpha = 0.7, size = 1.8) +
    theme_bw() +
    labs(title = "PCA numérica - k-proto", color = "Cluster")
  save_and_print_plot(p_pca_k, file.path(dir_out, "06_pca_kproto_robusta.png"), 7, 5)
  
  p_pca_c <- ggplot(df_pca, aes(x = PC1, y = PC2, color = cluster_cure)) +
    geom_point(alpha = 0.7, size = 1.8) +
    theme_bw() +
    labs(title = "PCA numérica - CURE adaptado", color = "Cluster")
  save_and_print_plot(p_pca_c, file.path(dir_out, "06_pca_cure_robusta.png"), 7, 5)
  
  p_pca_f <- ggplot(df_pca, aes(x = PC1, y = PC2, color = cluster_final)) +
    geom_point(alpha = 0.7, size = 1.8) +
    theme_bw() +
    labs(title = "PCA numérica - método final", color = "Cluster")
  save_and_print_plot(p_pca_f, file.path(dir_out, "06_pca_final_robusta.png"), 7, 5)
}

if (all(c("Start_Lat", "Start_Lng") %in% names(base_final))) {
  df_map <- data.frame(
    Start_Lat = base_final$Start_Lat,
    Start_Lng = base_final$Start_Lng,
    cluster_kproto = factor(cluster_kproto),
    cluster_cure = factor(cluster_cure),
    cluster_final = factor(cluster_final)
  )
  
  p_map_k <- ggplot(df_map, aes(x = Start_Lng, y = Start_Lat, color = cluster_kproto)) +
    geom_point(alpha = 0.7, size = 1.6) +
    theme_bw() +
    labs(title = "Mapa aproximado Lat/Lng - k-proto", color = "Cluster")
  save_and_print_plot(p_map_k, file.path(dir_out, "06_map_kproto_robusta.png"), 8, 5)
  
  p_map_c <- ggplot(df_map, aes(x = Start_Lng, y = Start_Lat, color = cluster_cure)) +
    geom_point(alpha = 0.7, size = 1.6) +
    theme_bw() +
    labs(title = "Mapa aproximado Lat/Lng - CURE adaptado", color = "Cluster")
  save_and_print_plot(p_map_c, file.path(dir_out, "06_map_cure_robusta.png"), 8, 5)
  
  p_map_f <- ggplot(df_map, aes(x = Start_Lng, y = Start_Lat, color = cluster_final)) +
    geom_point(alpha = 0.7, size = 1.6) +
    theme_bw() +
    labs(title = "Mapa aproximado Lat/Lng - método final", color = "Cluster")
  save_and_print_plot(p_map_f, file.path(dir_out, "06_map_final_robusta.png"), 8, 5)
}

# ------------------------------------------------------------------------------
# 11. Resumen texto
# ------------------------------------------------------------------------------
sink(file.path(dir_out, "06_resumen_clustering_robusto.txt"))
cat("==============================\n")
cat("RESUMEN CLUSTERING FINAL ROBUSTO\n")
cat("==============================\n\n")

cat("Dimensión base final clustering: ", nrow(base_final), "x", ncol(base_final), "\n\n")

cat("Variables categóricas:\n")
print(vars_factor)
cat("\n")

cat("Variables numéricas:\n")
print(vars_num)
cat("\n")

cat("Resultados k-proto (todas las combinaciones):\n")
print(kproto_results)
cat("\nMejor fila k-proto:\n")
print(best_row_kproto)
cat("\nTamaño clusters k-proto:\n")
print(table(cluster_kproto))
cat("\n")

cat("Resultados CURE adaptado:\n")
print(cure_results)
cat("\nMejor fila CURE:\n")
print(best_row_cure)
cat("\nTamaño clusters CURE adaptado:\n")
print(table(cluster_cure))
cat("\n")

cat("Comparativa final:\n")
print(comparativa)
cat("\nMétodo final seleccionado:", metodo_final, "\n\n")

cat("Tamaño clusters método final:\n")
print(table(cluster_final))
cat("\n")
sink()

cat("\nScript finalizado correctamente.\n")
cat("Resultados guardados en:", dir_out, "\n")

