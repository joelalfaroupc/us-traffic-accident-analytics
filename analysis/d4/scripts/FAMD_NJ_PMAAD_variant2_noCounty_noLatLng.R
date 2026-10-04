# ==============================================================================
# PMAAD - FAMD + K-means
# ------------------------------------------------------------------------------
# Objetivo:
#   - Eliminar County para evitar dominancia geográfica por alta cardinalidad
#   - Eliminar Start_Lat y Start_Lng para quitar peso geográfico directo
#   - Reducir redundancias meteorológicas fuertes
#   - Transformar variables temporales cíclicas: mes y hora
#   - Imputar valores ausentes en vez de eliminar columnas completas
#   - Seleccionar dimensiones FAMD según varianza acumulada
#   - Probar K = 3, 4, 5, 6 con K-means
#   - Generar gráficos y tablas de perfil para interpretar los clusters
# ------------------------------------------------------------------------------
# Base esperada:
#   US_Accidents_Final_FE.csv
# ==============================================================================

rm(list = ls())
set.seed(123)

# ------------------------------------------------------------------------------
# 0. Configuración
# ------------------------------------------------------------------------------

candidate_files <- c(
  "US_Accidents_Final_FE.csv",
  "US_Accidents_Final_FE.csv"
)

file_in <- candidate_files[file.exists(candidate_files)][1]

if (is.na(file_in)) {
  stop("No se ha encontrado US_Accidents_Final_FE.csv")
}

out_dir <- "FAMD_KMEANS_RESULTS"

plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
base_dir  <- file.path(out_dir, "bases")

dirs <- c(out_dir, plot_dir, table_dir, base_dir)

for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Configuración del clustering
k_grid <- 3:6

# Número máximo de dimensiones a calcular en FAMD
ncp_famd_total <- 20

# Varianza acumulada mínima que se quiere conservar para K-means
target_cumvar <- 40

# Dimensiones mínimas y máximas para clustering
min_dims_cluster <- 5
max_dims_cluster <- 12

# ------------------------------------------------------------------------------
# 1. Librerías
# ------------------------------------------------------------------------------

packages <- c(
  "FactoMineR",
  "factoextra",
  "ggplot2",
  "dplyr",
  "cluster",
  "clusterCrit",
  "readr",
  "tibble",
  "tidyr",
  "forcats"
)

new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]

if (length(new_packages) > 0) {
  install.packages(new_packages, dependencies = TRUE)
}

invisible(lapply(packages, require, character.only = TRUE))

# ------------------------------------------------------------------------------
# 2. Funciones auxiliares
# ------------------------------------------------------------------------------

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  print(plot_obj)
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
}

safe_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

metric_rank_score <- function(x, decreasing = TRUE) {
  if (decreasing) {
    rank(-x, ties.method = "min")
  } else {
    rank(x, ties.method = "min")
  }
}

save_pairs_plot <- function(df_plot, filename, title_text = "") {
  png(filename = filename, width = 1200, height = 1200, res = 120)
  pairs(
    df_plot,
    main = title_text,
    pch = 19,
    cex = 0.45,
    col = as.numeric(df_plot$cluster)
  )
  dev.off()
}

to_binary_factor <- function(x, false_label, true_label) {
  x_chr <- tolower(as.character(x))
  
  out <- dplyr::case_when(
    x_chr %in% c("false", "f", "0", "no", "n") ~ false_label,
    x_chr %in% c("true", "t", "1", "yes", "y") ~ true_label,
    TRUE ~ NA_character_
  )
  
  factor(out, levels = c(false_label, true_label))
}

# ------------------------------------------------------------------------------
# 3. Carga de datos
# ------------------------------------------------------------------------------

df_raw <- read.csv(file_in, stringsAsFactors = FALSE)

cat("Archivo cargado:", file_in, "\n")
cat("Dimensión base:", nrow(df_raw), "x", ncol(df_raw), "\n")

df <- df_raw

# ------------------------------------------------------------------------------
# 4. Preparación temporal
# ------------------------------------------------------------------------------

if (!"Start_Time" %in% names(df)) {
  stop("No existe la columna Start_Time en la base.")
}

dt <- as.POSIXct(df$Start_Time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

if (sum(is.na(dt)) > 0) {
  dt_alt <- as.POSIXct(df$Start_Time, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
  dt[is.na(dt)] <- dt_alt[is.na(dt)]
}

if (sum(is.na(dt)) > 0) {
  warning("Hay fechas Start_Time que no se han podido parsear.")
}

# Variables temporales cíclicas
hour_num <- as.integer(format(dt, "%H"))
month_num <- as.integer(format(dt, "%m"))

df$Hour_sin <- sin(2 * pi * hour_num / 24)
df$Hour_cos <- cos(2 * pi * hour_num / 24)

df$Month_sin <- sin(2 * pi * month_num / 12)
df$Month_cos <- cos(2 * pi * month_num / 12)

# Variables interpretables para perfiles posteriores
df$Start_Hour <- hour_num
df$Start_Month_Num <- month_num

df$Start_Month_Label <- factor(
  month_num,
  levels = 1:12,
  labels = month.abb
)

df$Is_Weekend <- ifelse(
  weekdays(dt) %in% c("Saturday", "Sunday"),
  "Weekend",
  "Weekday"
)

df$Is_Weekend <- factor(df$Is_Weekend)

# ------------------------------------------------------------------------------
# 5. Reetiquetado de variables categóricas
# ------------------------------------------------------------------------------

if ("Severity" %in% names(df)) {
  df$Severity_num_original <- safe_numeric(df$Severity)
  
  df$Severity <- factor(
    df$Severity,
    levels = c(1, 2, 3, 4),
    labels = c("S1", "S2", "S3", "S4")
  )
}

if ("Junction" %in% names(df)) {
  df$Junction <- to_binary_factor(
    df$Junction,
    false_label = "No_Junction",
    true_label = "Junction"
  )
}

if ("Traffic_Signal" %in% names(df)) {
  df$Traffic_Signal <- to_binary_factor(
    df$Traffic_Signal,
    false_label = "No_Signal",
    true_label = "Signal"
  )
}

if ("Has_Infrastructure" %in% names(df)) {
  df$Has_Infrastructure <- to_binary_factor(
    df$Has_Infrastructure,
    false_label = "No_Infrastructure",
    true_label = "With_Infrastructure"
  )
}

if ("Sunrise_Sunset" %in% names(df)) {
  df$Sunrise_Sunset <- factor(df$Sunrise_Sunset)
}

if ("Weather_Group" %in% names(df)) {
  df$Weather_Group <- factor(df$Weather_Group)
}

if ("Wind_Group" %in% names(df)) {
  df$Wind_Group <- factor(df$Wind_Group)
}

# ------------------------------------------------------------------------------
# 6. Selección de variables para FAMD
# ------------------------------------------------------------------------------

# Se eliminan County, Start_Lat y Start_Lng para reducir el peso geográfico.
# Se elimina Wind_Chill.C. para evitar redundancia con Temperature.C.
# Start_Hour y Start_Month se sustituyen por variables cíclicas.

vars_quanti <- c(
  "Humidity...",
  "Pressure.in.",
  "Wind_Speed.mph.",
  "Precipitation.in.",
  "Temperature.C.",
  "Hour_sin",
  "Hour_cos",
  "Month_sin",
  "Month_cos"
)

vars_quali <- c(
  "Severity",
  "Junction",
  "Traffic_Signal",
  "Sunrise_Sunset",
  "Weather_Group",
  "Wind_Group",
  "Has_Infrastructure",
  "Is_Weekend"
)

vars_quanti <- intersect(vars_quanti, names(df))
vars_quali  <- intersect(vars_quali, names(df))

df_famd <- df[, c(vars_quali, vars_quanti), drop = FALSE]

# Tipos
df_famd[vars_quanti] <- lapply(df_famd[vars_quanti], safe_numeric)
df_famd[vars_quali]  <- lapply(df_famd[vars_quali], factor)

# ------------------------------------------------------------------------------
# 7. Imputación de valores ausentes
# ------------------------------------------------------------------------------

# Numéricas: mediana
for (v in vars_quanti) {
  if (any(is.na(df_famd[[v]]))) {
    med_v <- median(df_famd[[v]], na.rm = TRUE)
    
    if (is.na(med_v)) {
      warning(paste("La variable numérica", v, "tiene todos los valores NA. Se eliminará."))
      df_famd[[v]] <- NULL
    } else {
      df_famd[[v]][is.na(df_famd[[v]])] <- med_v
    }
  }
}

vars_quanti <- intersect(vars_quanti, names(df_famd))

# Cualitativas: nivel Unknown
for (v in vars_quali) {
  if (v %in% names(df_famd)) {
    x <- as.character(df_famd[[v]])
    x[is.na(x) | x == ""] <- "Unknown"
    df_famd[[v]] <- factor(x)
  }
}

vars_quali <- intersect(vars_quali, names(df_famd))

# ------------------------------------------------------------------------------
# 8. Limpieza de variables problemáticas
# ------------------------------------------------------------------------------

zero_var_num <- vars_quanti[
  sapply(df_famd[vars_quanti], function(x) sd(x, na.rm = TRUE) == 0)
]

if (length(zero_var_num) > 0) {
  cat("Se eliminan numéricas de varianza cero:\n")
  print(zero_var_num)
  
  df_famd <- df_famd[, setdiff(names(df_famd), zero_var_num), drop = FALSE]
  vars_quanti <- intersect(vars_quanti, names(df_famd))
}

one_level_fac <- vars_quali[
  sapply(df_famd[vars_quali], function(x) length(unique(na.omit(x))) <= 1)
]

if (length(one_level_fac) > 0) {
  cat("Se eliminan cualitativas con un solo nivel:\n")
  print(one_level_fac)
  
  df_famd <- df_famd[, setdiff(names(df_famd), one_level_fac), drop = FALSE]
  vars_quali <- intersect(vars_quali, names(df_famd))
}

cat("Dimensión base FAMD:", nrow(df_famd), "x", ncol(df_famd), "\n")
cat("Variables cualitativas:", paste(vars_quali, collapse = ", "), "\n")
cat("Variables cuantitativas:", paste(vars_quanti, collapse = ", "), "\n")

write.csv(
  df_famd,
  file.path(base_dir, "00_base_famd.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. FAMD
# ------------------------------------------------------------------------------

ncp_final <- min(ncp_famd_total, ncol(df_famd) - 1)

if (ncp_final < 2) {
  stop("No hay suficientes variables para calcular FAMD.")
}

res.famd <- FactoMineR::FAMD(
  df_famd,
  ncp = ncp_final,
  graph = FALSE
)

capture.output(
  print(res.famd),
  file = file.path(table_dir, "01_famd_print.txt")
)

# Extraer eigenvalues de forma robusta
eig.val <- as.data.frame(factoextra::get_eigenvalue(res.famd))

# Normalizar nombres de columnas
names(eig.val) <- make.names(names(eig.val))

cat("Columnas detectadas en eig.val:\n")
print(names(eig.val))

# Buscar columna de porcentaje de varianza
var_col <- grep("variance|percent|percentage", names(eig.val), ignore.case = TRUE, value = TRUE)

if (length(var_col) == 0) {
  stop("No se ha encontrado la columna de porcentaje de varianza en eig.val.")
}

if ("variance.percent" %in% var_col) {
  var_col <- "variance.percent"
} else {
  var_col <- var_col[1]
}

eig.val$cumvar <- cumsum(eig.val[[var_col]])

write.csv(
  eig.val,
  file.path(table_dir, "02_famd_eigenvalues.csv"),
  row.names = TRUE
)

# ------------------------------------------------------------------------------
# 10. Selección automática de dimensiones para clustering
# ------------------------------------------------------------------------------

dim_target <- which(eig.val$cumvar >= target_cumvar)[1]

if (is.na(dim_target)) {
  dim_target <- nrow(eig.val)
}

ncp_cluster <- max(min_dims_cluster, dim_target)
ncp_cluster <- min(ncp_cluster, max_dims_cluster)
ncp_cluster <- min(ncp_cluster, ncol(res.famd$ind$coord))

cat("Dimensiones usadas para clustering:", ncp_cluster, "\n")
cat("Varianza acumulada aproximada:", eig.val$cumvar[ncp_cluster], "%\n")

write.csv(
  data.frame(
    target_cumvar = target_cumvar,
    ncp_cluster = ncp_cluster,
    cumvar_used = eig.val$cumvar[ncp_cluster]
  ),
  file.path(table_dir, "02b_dims_used_for_clustering.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. Visualización FAMD
# ------------------------------------------------------------------------------

p_eig <- factoextra::fviz_screeplot(res.famd, addlabels = TRUE) +
  theme_bw() +
  ggtitle("FAMD - Scree plot")

save_plot(
  p_eig,
  file.path(plot_dir, "01_famd_screeplot.png"),
  8,
  5
)

p_var_all <- factoextra::fviz_famd_var(res.famd, repel = TRUE) +
  theme_bw() +
  ggtitle("FAMD - Variables")

save_plot(
  p_var_all,
  file.path(plot_dir, "02_famd_var_all.png"),
  9,
  7
)

p_contrib_1 <- factoextra::fviz_contrib(
  res.famd,
  choice = "var",
  axes = 1,
  top = 15
) +
  theme_bw() +
  ggtitle("FAMD - Contribución a Dim 1")

save_plot(
  p_contrib_1,
  file.path(plot_dir, "03_famd_contrib_dim1.png"),
  11,
  5
)

p_contrib_2 <- factoextra::fviz_contrib(
  res.famd,
  choice = "var",
  axes = 2,
  top = 15
) +
  theme_bw() +
  ggtitle("FAMD - Contribución a Dim 2")

save_plot(
  p_contrib_2,
  file.path(plot_dir, "04_famd_contrib_dim2.png"),
  11,
  5
)

if ("Severity" %in% vars_quali) {
  p_ind_severity <- factoextra::fviz_famd_ind(
    res.famd,
    habillage = "Severity",
    addEllipses = TRUE,
    ellipse.type = "confidence",
    repel = FALSE,
    label = "none",
    alpha.ind = 0.45,
    pointsize = 1.2
  ) +
    theme_bw() +
    ggtitle("FAMD - Individuos según Severity")
  
  save_plot(
    p_ind_severity,
    file.path(plot_dir, "05_famd_ind_severity.png"),
    8,
    6
  )
}

# ------------------------------------------------------------------------------
# 12. Biplot manual
# ------------------------------------------------------------------------------

ind_coord <- as.data.frame(res.famd$ind$coord[, 1:2, drop = FALSE])

if ("Severity" %in% names(df_famd)) {
  ind_coord$Severity <- df_famd$Severity
} else {
  ind_coord$Severity <- "Individuo"
}

var_coord <- as.data.frame(res.famd$var$coord[, 1:2, drop = FALSE])
var_coord$label <- rownames(var_coord)

range_ind_1 <- max(ind_coord$Dim.1, na.rm = TRUE) - min(ind_coord$Dim.1, na.rm = TRUE)
range_ind_2 <- max(ind_coord$Dim.2, na.rm = TRUE) - min(ind_coord$Dim.2, na.rm = TRUE)

range_var_1 <- max(var_coord$Dim.1, na.rm = TRUE) - min(var_coord$Dim.1, na.rm = TRUE)
range_var_2 <- max(var_coord$Dim.2, na.rm = TRUE) - min(var_coord$Dim.2, na.rm = TRUE)

mult <- min(range_ind_1 / range_var_1, range_ind_2 / range_var_2) * 0.7

var_coord$Dim.1 <- var_coord$Dim.1 * mult
var_coord$Dim.2 <- var_coord$Dim.2 * mult

p_biplot <- ggplot() +
  geom_point(
    data = ind_coord,
    aes(x = Dim.1, y = Dim.2, color = Severity),
    alpha = 0.45,
    size = 1.3
  ) +
  geom_segment(
    data = var_coord,
    aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
    arrow = grid::arrow(length = grid::unit(0.18, "cm")),
    color = "black",
    linewidth = 0.5
  ) +
  geom_text(
    data = var_coord,
    aes(x = Dim.1, y = Dim.2, label = label),
    size = 3,
    vjust = -0.4
  ) +
  theme_bw() +
  labs(
    title = "FAMD - Biplot",
    x = "Dim 1",
    y = "Dim 2",
    color = "Severity"
  )

save_plot(
  p_biplot,
  file.path(plot_dir, "06_famd_biplot.png"),
  10,
  7
)

# ------------------------------------------------------------------------------
# 13. Coordenadas de individuos para clustering
# ------------------------------------------------------------------------------

coords_ind <- as.data.frame(
  res.famd$ind$coord[, 1:ncp_cluster, drop = FALSE]
)

write.csv(
  coords_ind,
  file.path(base_dir, "01_coords_ind_dims_cluster.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 14. K-means para K = 3:6
# ------------------------------------------------------------------------------

dist_coords <- dist(coords_ind)

results_k <- data.frame(
  k = integer(),
  silhouette = numeric(),
  calinski_harabasz = numeric(),
  dunn = numeric(),
  tot_withinss = numeric(),
  betweenss = numeric(),
  ratio_betweenss_totss = numeric(),
  stringsAsFactors = FALSE
)

models_k <- list()

for (k in k_grid) {
  set.seed(100 + k)
  
  km <- kmeans(
    coords_ind,
    centers = k,
    nstart = 50,
    iter.max = 100
  )
  
  cl <- km$cluster
  
  sil <- mean(cluster::silhouette(cl, dist_coords)[, 3])
  
  int_crit <- clusterCrit::intCriteria(
    traj = as.matrix(coords_ind),
    part = as.integer(cl),
    crit = c("Calinski_Harabasz", "Dunn")
  )
  
  results_k <- rbind(
    results_k,
    data.frame(
      k = k,
      silhouette = sil,
      calinski_harabasz = as.numeric(int_crit$calinski_harabasz),
      dunn = as.numeric(int_crit$dunn),
      tot_withinss = km$tot.withinss,
      betweenss = km$betweenss,
      ratio_betweenss_totss = km$betweenss / km$totss
    )
  )
  
  models_k[[as.character(k)]] <- km
}

results_k$rank_sil  <- metric_rank_score(results_k$silhouette, decreasing = TRUE)
results_k$rank_ch   <- metric_rank_score(results_k$calinski_harabasz, decreasing = TRUE)
results_k$rank_dunn <- metric_rank_score(results_k$dunn, decreasing = TRUE)
results_k$rank_wss  <- metric_rank_score(results_k$tot_withinss, decreasing = FALSE)

results_k$rank_total <- results_k$rank_sil +
  results_k$rank_ch +
  results_k$rank_dunn +
  results_k$rank_wss

results_k <- results_k |>
  dplyr::arrange(rank_total)

write.csv(
  results_k,
  file.path(table_dir, "03_kmeans_metrics_k3_k6.csv"),
  row.names = FALSE
)

best_k <- results_k$k[1]

cat("Mejor K según ranking agregado:", best_k, "\n")

# ------------------------------------------------------------------------------
# 15. Plots de métricas
# ------------------------------------------------------------------------------

metrics_long <- results_k |>
  dplyr::select(k, silhouette, calinski_harabasz, dunn, tot_withinss, ratio_betweenss_totss) |>
  tidyr::pivot_longer(
    cols = -k,
    names_to = "metric",
    values_to = "value"
  )

p_metrics <- ggplot(metrics_long, aes(x = k, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Comparación de métricas para K = 3..6",
    x = "Número de clusters",
    y = "Valor"
  )

save_plot(
  p_metrics,
  file.path(plot_dir, "07_kmeans_metrics_comparison.png"),
  10,
  7
)

p_wss <- factoextra::fviz_nbclust(
  coords_ind,
  kmeans,
  method = "wss",
  k.max = max(k_grid)
) +
  theme_bw() +
  ggtitle("WSS - K-means sobre dimensiones FAMD seleccionadas")

save_plot(
  p_wss,
  file.path(plot_dir, "08_nbclust_wss.png"),
  8,
  5
)

p_sil_aux <- factoextra::fviz_nbclust(
  coords_ind,
  kmeans,
  method = "silhouette",
  k.max = max(k_grid)
) +
  theme_bw() +
  ggtitle("Silhouette - K-means sobre dimensiones FAMD seleccionadas")

save_plot(
  p_sil_aux,
  file.path(plot_dir, "09_nbclust_silhouette.png"),
  8,
  5
)

# ------------------------------------------------------------------------------
# 16. Plots y perfiles para cada K
# ------------------------------------------------------------------------------

for (k in k_grid) {
  km <- models_k[[as.character(k)]]
  
  df_k <- coords_ind
  df_k$cluster <- factor(km$cluster)
  
  # ---------------------------------------------------------------------------
  # 16.1 Plots factoriales
  # ---------------------------------------------------------------------------
  
  p12 <- ggplot(df_k, aes(x = Dim.1, y = Dim.2, color = cluster)) +
    geom_point(alpha = 0.60, size = 1.2) +
    stat_ellipse(type = "norm", linewidth = 0.7, alpha = 0.15) +
    theme_bw() +
    labs(
      title = paste0("K = ", k, " | Dim 1 vs Dim 2"),
      x = "Dim 1",
      y = "Dim 2",
      color = "Cluster"
    )
  
  save_plot(
    p12,
    file.path(plot_dir, paste0("10_k", k, "_dim1_dim2.png")),
    8,
    6
  )
  
  if (ncp_cluster >= 3) {
    p13 <- ggplot(df_k, aes(x = Dim.1, y = Dim.3, color = cluster)) +
      geom_point(alpha = 0.60, size = 1.2) +
      stat_ellipse(type = "norm", linewidth = 0.7, alpha = 0.15) +
      theme_bw() +
      labs(
        title = paste0("K = ", k, " | Dim 1 vs Dim 3"),
        x = "Dim 1",
        y = "Dim 3",
        color = "Cluster"
      )
    
    save_plot(
      p13,
      file.path(plot_dir, paste0("11_k", k, "_dim1_dim3.png")),
      8,
      6
    )
    
    p23 <- ggplot(df_k, aes(x = Dim.2, y = Dim.3, color = cluster)) +
      geom_point(alpha = 0.60, size = 1.2) +
      stat_ellipse(type = "norm", linewidth = 0.7, alpha = 0.15) +
      theme_bw() +
      labs(
        title = paste0("K = ", k, " | Dim 2 vs Dim 3"),
        x = "Dim 2",
        y = "Dim 3",
        color = "Cluster"
      )
    
    save_plot(
      p23,
      file.path(plot_dir, paste0("12_k", k, "_dim2_dim3.png")),
      8,
      6
    )
  }
  
  if (ncp_cluster >= 5) {
    p45 <- ggplot(df_k, aes(x = Dim.4, y = Dim.5, color = cluster)) +
      geom_point(alpha = 0.60, size = 1.2) +
      stat_ellipse(type = "norm", linewidth = 0.7, alpha = 0.15) +
      theme_bw() +
      labs(
        title = paste0("K = ", k, " | Dim 4 vs Dim 5"),
        x = "Dim 4",
        y = "Dim 5",
        color = "Cluster"
      )
    
    save_plot(
      p45,
      file.path(plot_dir, paste0("13_k", k, "_dim4_dim5.png")),
      8,
      6
    )
  }
  
  dims_for_facets <- paste0("Dim.", 2:min(5, ncp_cluster))
  
  df_long <- do.call(
    rbind,
    lapply(dims_for_facets, function(d) {
      data.frame(
        x = df_k$Dim.1,
        y = df_k[[d]],
        comp = gsub("Dim\\.", "Dim", d),
        cluster = df_k$cluster
      )
    })
  )
  
  p_facets <- ggplot(df_long, aes(x = x, y = y, color = cluster)) +
    geom_point(alpha = 0.55, size = 1.0) +
    facet_wrap(~ comp, scales = "free") +
    theme_bw() +
    labs(
      title = paste0("K = ", k, " | Dim 1 frente al resto de dimensiones"),
      x = "Dim 1",
      y = "Componente",
      color = "Cluster"
    )
  
  save_plot(
    p_facets,
    file.path(plot_dir, paste0("14_k", k, "_facets_dim1_vs_rest.png")),
    10,
    7
  )
  
  dims_pairs <- paste0("Dim.", 1:min(5, ncp_cluster))
  
  df_pairs <- df_k[, dims_pairs, drop = FALSE]
  df_pairs$cluster <- as.numeric(df_k$cluster)
  
  save_pairs_plot(
    df_pairs,
    file.path(plot_dir, paste0("15_k", k, "_pairs_dims.png")),
    title_text = paste0("Pairs plot - dimensiones FAMD | K = ", k)
  )
  
  sizes_k <- as.data.frame(table(df_k$cluster))
  names(sizes_k) <- c("cluster", "n")
  
  p_sizes <- ggplot(sizes_k, aes(x = cluster, y = n, fill = cluster)) +
    geom_col() +
    theme_bw() +
    labs(
      title = paste0("K = ", k, " | Tamaño de clusters"),
      x = "Cluster",
      y = "Número de observaciones",
      fill = "Cluster"
    )
  
  save_plot(
    p_sizes,
    file.path(plot_dir, paste0("16_k", k, "_cluster_sizes.png")),
    8,
    5
  )
  
  # ---------------------------------------------------------------------------
  # 16.2 Bases con cluster
  # ---------------------------------------------------------------------------
  
  df_out_k <- df_raw
  df_out_k$cluster_famd_kmeans <- km$cluster
  
  write.csv(
    df_out_k,
    file.path(base_dir, paste0("02_base_con_cluster_k", k, ".csv")),
    row.names = FALSE
  )
  
  coords_out_k <- coords_ind
  coords_out_k$cluster_famd_kmeans <- km$cluster
  
  write.csv(
    coords_out_k,
    file.path(base_dir, paste0("03_coords_cluster_k", k, ".csv")),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 16.3 Perfil numérico de clusters
  # ---------------------------------------------------------------------------
  
  df_profile <- df
  df_profile$cluster <- factor(km$cluster)
  
  if ("Severity" %in% names(df_profile)) {
    df_profile$Severity_numeric_profile <- safe_numeric(gsub("S", "", as.character(df_profile$Severity)))
  } else if ("Severity_num_original" %in% names(df_profile)) {
    df_profile$Severity_numeric_profile <- df_profile$Severity_num_original
  } else if ("Severity" %in% names(df_raw)) {
    df_profile$Severity_numeric_profile <- safe_numeric(df_raw$Severity)
  } else {
    df_profile$Severity_numeric_profile <- NA_real_
  }
  
  if ("Has_Infrastructure" %in% names(df_profile)) {
    df_profile$Infrastructure_binary <- ifelse(
      df_profile$Has_Infrastructure == "With_Infrastructure",
      1,
      0
    )
  } else {
    df_profile$Infrastructure_binary <- NA_real_
  }
  
  if ("Junction" %in% names(df_profile)) {
    df_profile$Junction_binary <- ifelse(
      df_profile$Junction == "Junction",
      1,
      0
    )
  } else {
    df_profile$Junction_binary <- NA_real_
  }
  
  if ("Traffic_Signal" %in% names(df_profile)) {
    df_profile$Signal_binary <- ifelse(
      df_profile$Traffic_Signal == "Signal",
      1,
      0
    )
  } else {
    df_profile$Signal_binary <- NA_real_
  }
  
  profile_num <- df_profile |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      n = dplyr::n(),
      pct = 100 * dplyr::n() / nrow(df_profile),
      
      severity_mean = mean(Severity_numeric_profile, na.rm = TRUE),
      severe_rate_s3_s4 = mean(Severity_numeric_profile >= 3, na.rm = TRUE),
      
      temp_mean = if ("Temperature.C." %in% names(df_profile)) mean(Temperature.C., na.rm = TRUE) else NA_real_,
      humidity_mean = if ("Humidity..." %in% names(df_profile)) mean(Humidity..., na.rm = TRUE) else NA_real_,
      pressure_mean = if ("Pressure.in." %in% names(df_profile)) mean(Pressure.in., na.rm = TRUE) else NA_real_,
      wind_speed_mean = if ("Wind_Speed.mph." %in% names(df_profile)) mean(Wind_Speed.mph., na.rm = TRUE) else NA_real_,
      precipitation_mean = if ("Precipitation.in." %in% names(df_profile)) mean(Precipitation.in., na.rm = TRUE) else NA_real_,
      
      infrastructure_rate = mean(Infrastructure_binary, na.rm = TRUE),
      junction_rate = mean(Junction_binary, na.rm = TRUE),
      traffic_signal_rate = mean(Signal_binary, na.rm = TRUE),
      
      hour_mean = mean(Start_Hour, na.rm = TRUE),
      month_mean = mean(Start_Month_Num, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  write.csv(
    profile_num,
    file.path(table_dir, paste0("profile_numeric_k", k, ".csv")),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 16.4 Perfiles categóricos
  # ---------------------------------------------------------------------------
  
  categorical_profile_vars <- intersect(
    c(
      "Severity",
      "Weather_Group",
      "Wind_Group",
      "Sunrise_Sunset",
      "Is_Weekend",
      "Junction",
      "Traffic_Signal",
      "Has_Infrastructure",
      "Start_Month_Label"
    ),
    names(df_profile)
  )
  
  for (cat_var in categorical_profile_vars) {
    profile_cat <- df_profile |>
      dplyr::count(cluster, .data[[cat_var]]) |>
      dplyr::group_by(cluster) |>
      dplyr::mutate(prop = n / sum(n)) |>
      dplyr::ungroup()
    
    names(profile_cat)[2] <- cat_var
    
    write.csv(
      profile_cat,
      file.path(table_dir, paste0("profile_", cat_var, "_k", k, ".csv")),
      row.names = FALSE
    )
  }
}

# ------------------------------------------------------------------------------
# 17. Profiling comparativo K = 3 vs K = 4
# ------------------------------------------------------------------------------

if ("3" %in% names(models_k) && "4" %in% names(models_k)) {
  
  km3 <- models_k[["3"]]
  km4 <- models_k[["4"]]
  
  df_compare <- df
  df_compare$cluster_k3 <- factor(km3$cluster)
  df_compare$cluster_k4 <- factor(km4$cluster)
  
  # ---------------------------------------------------------------------------
  # 17.1 Tabla de cruce entre K = 3 y K = 4
  # ---------------------------------------------------------------------------
  
  cross_k3_k4 <- as.data.frame(table(
    cluster_k3 = df_compare$cluster_k3,
    cluster_k4 = df_compare$cluster_k4
  ))
  
  cross_k3_k4 <- cross_k3_k4 |>
    dplyr::group_by(cluster_k4) |>
    dplyr::mutate(
      prop_within_k4 = Freq / sum(Freq)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(cluster_k3) |>
    dplyr::mutate(
      prop_within_k3 = Freq / sum(Freq)
    ) |>
    dplyr::ungroup()
  
  write.csv(
    cross_k3_k4,
    file.path(table_dir, "profiling_cross_k3_k4.csv"),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 17.2 Distancias entre centroides de K = 4 y centroides de K = 3
  # ---------------------------------------------------------------------------
  
  centers_k3 <- as.data.frame(km3$centers)
  centers_k3$cluster_k3 <- rownames(centers_k3)
  
  centers_k4 <- as.data.frame(km4$centers)
  centers_k4$cluster_k4 <- rownames(centers_k4)
  
  centroid_distances <- expand.grid(
    cluster_k4 = centers_k4$cluster_k4,
    cluster_k3 = centers_k3$cluster_k3,
    stringsAsFactors = FALSE
  )
  
  centroid_distances$distance <- NA_real_
  
  for (i in seq_len(nrow(centroid_distances))) {
    c4 <- centroid_distances$cluster_k4[i]
    c3 <- centroid_distances$cluster_k3[i]
    
    v4 <- as.numeric(centers_k4[centers_k4$cluster_k4 == c4, names(coords_ind)])
    v3 <- as.numeric(centers_k3[centers_k3$cluster_k3 == c3, names(coords_ind)])
    
    centroid_distances$distance[i] <- sqrt(sum((v4 - v3)^2))
  }
  
  centroid_distances <- centroid_distances |>
    dplyr::group_by(cluster_k4) |>
    dplyr::arrange(distance, .by_group = TRUE) |>
    dplyr::mutate(rank_distance = dplyr::row_number()) |>
    dplyr::ungroup()
  
  write.csv(
    centroid_distances,
    file.path(table_dir, "profiling_centroid_distances_k4_to_k3.csv"),
    row.names = FALSE
  )
  
  closest_centroid <- centroid_distances |>
    dplyr::filter(rank_distance == 1) |>
    dplyr::select(cluster_k4, closest_cluster_k3 = cluster_k3, distance)
  
  write.csv(
    closest_centroid,
    file.path(table_dir, "profiling_closest_centroid_k4_to_k3.csv"),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 17.3 Variables auxiliares para profiling
  # ---------------------------------------------------------------------------
  
  if ("Severity" %in% names(df_compare)) {
    df_compare$Severity_numeric_profile <- safe_numeric(gsub("S", "", as.character(df_compare$Severity)))
  } else if ("Severity_num_original" %in% names(df_compare)) {
    df_compare$Severity_numeric_profile <- df_compare$Severity_num_original
  } else if ("Severity" %in% names(df_raw)) {
    df_compare$Severity_numeric_profile <- safe_numeric(df_raw$Severity)
  } else {
    df_compare$Severity_numeric_profile <- NA_real_
  }
  
  if ("Has_Infrastructure" %in% names(df_compare)) {
    df_compare$Infrastructure_binary <- ifelse(
      df_compare$Has_Infrastructure == "With_Infrastructure",
      1,
      0
    )
  } else {
    df_compare$Infrastructure_binary <- NA_real_
  }
  
  if ("Junction" %in% names(df_compare)) {
    df_compare$Junction_binary <- ifelse(
      df_compare$Junction == "Junction",
      1,
      0
    )
  } else {
    df_compare$Junction_binary <- NA_real_
  }
  
  if ("Traffic_Signal" %in% names(df_compare)) {
    df_compare$Signal_binary <- ifelse(
      df_compare$Traffic_Signal == "Signal",
      1,
      0
    )
  } else {
    df_compare$Signal_binary <- NA_real_
  }
  
  # ---------------------------------------------------------------------------
  # 17.4 Perfil numérico de K = 3
  # ---------------------------------------------------------------------------
  
  profile_k3 <- df_compare |>
    dplyr::group_by(cluster_k3) |>
    dplyr::summarise(
      n = dplyr::n(),
      pct = 100 * dplyr::n() / nrow(df_compare),
      
      severity_mean = mean(Severity_numeric_profile, na.rm = TRUE),
      severe_rate_s3_s4 = mean(Severity_numeric_profile >= 3, na.rm = TRUE),
      
      temp_mean = if ("Temperature.C." %in% names(df_compare)) mean(Temperature.C., na.rm = TRUE) else NA_real_,
      humidity_mean = if ("Humidity..." %in% names(df_compare)) mean(Humidity..., na.rm = TRUE) else NA_real_,
      pressure_mean = if ("Pressure.in." %in% names(df_compare)) mean(Pressure.in., na.rm = TRUE) else NA_real_,
      wind_speed_mean = if ("Wind_Speed.mph." %in% names(df_compare)) mean(Wind_Speed.mph., na.rm = TRUE) else NA_real_,
      precipitation_mean = if ("Precipitation.in." %in% names(df_compare)) mean(Precipitation.in., na.rm = TRUE) else NA_real_,
      
      infrastructure_rate = mean(Infrastructure_binary, na.rm = TRUE),
      junction_rate = mean(Junction_binary, na.rm = TRUE),
      traffic_signal_rate = mean(Signal_binary, na.rm = TRUE),
      
      hour_mean = mean(Start_Hour, na.rm = TRUE),
      month_mean = mean(Start_Month_Num, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  write.csv(
    profile_k3,
    file.path(table_dir, "profiling_numeric_k3.csv"),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 17.5 Perfil numérico de K = 4
  # ---------------------------------------------------------------------------
  
  profile_k4 <- df_compare |>
    dplyr::group_by(cluster_k4) |>
    dplyr::summarise(
      n = dplyr::n(),
      pct = 100 * dplyr::n() / nrow(df_compare),
      
      severity_mean = mean(Severity_numeric_profile, na.rm = TRUE),
      severe_rate_s3_s4 = mean(Severity_numeric_profile >= 3, na.rm = TRUE),
      
      temp_mean = if ("Temperature.C." %in% names(df_compare)) mean(Temperature.C., na.rm = TRUE) else NA_real_,
      humidity_mean = if ("Humidity..." %in% names(df_compare)) mean(Humidity..., na.rm = TRUE) else NA_real_,
      pressure_mean = if ("Pressure.in." %in% names(df_compare)) mean(Pressure.in., na.rm = TRUE) else NA_real_,
      wind_speed_mean = if ("Wind_Speed.mph." %in% names(df_compare)) mean(Wind_Speed.mph., na.rm = TRUE) else NA_real_,
      precipitation_mean = if ("Precipitation.in." %in% names(df_compare)) mean(Precipitation.in., na.rm = TRUE) else NA_real_,
      
      infrastructure_rate = mean(Infrastructure_binary, na.rm = TRUE),
      junction_rate = mean(Junction_binary, na.rm = TRUE),
      traffic_signal_rate = mean(Signal_binary, na.rm = TRUE),
      
      hour_mean = mean(Start_Hour, na.rm = TRUE),
      month_mean = mean(Start_Month_Num, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  write.csv(
    profile_k4,
    file.path(table_dir, "profiling_numeric_k4.csv"),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 17.6 Perfil categórico comparativo
  # ---------------------------------------------------------------------------
  
  categorical_profile_vars <- intersect(
    c(
      "Severity",
      "Weather_Group",
      "Wind_Group",
      "Sunrise_Sunset",
      "Is_Weekend",
      "Junction",
      "Traffic_Signal",
      "Has_Infrastructure",
      "Start_Month_Label"
    ),
    names(df_compare)
  )
  
  for (cat_var in categorical_profile_vars) {
    
    profile_cat_k3 <- df_compare |>
      dplyr::count(cluster_k3, .data[[cat_var]]) |>
      dplyr::group_by(cluster_k3) |>
      dplyr::mutate(prop = n / sum(n)) |>
      dplyr::ungroup()
    
    names(profile_cat_k3)[2] <- cat_var
    
    write.csv(
      profile_cat_k3,
      file.path(table_dir, paste0("profiling_", cat_var, "_k3.csv")),
      row.names = FALSE
    )
    
    profile_cat_k4 <- df_compare |>
      dplyr::count(cluster_k4, .data[[cat_var]]) |>
      dplyr::group_by(cluster_k4) |>
      dplyr::mutate(prop = n / sum(n)) |>
      dplyr::ungroup()
    
    names(profile_cat_k4)[2] <- cat_var
    
    write.csv(
      profile_cat_k4,
      file.path(table_dir, paste0("profiling_", cat_var, "_k4.csv")),
      row.names = FALSE
    )
  }
  
  # ------------------------------------------------------------------------------
  # 18. Visualización del profiling K = 3 vs K = 4
  # ------------------------------------------------------------------------------
  
  profiling_plot_dir <- file.path(plot_dir, "profiling_k3_vs_k4")
  
  if (!dir.exists(profiling_plot_dir)) {
    dir.create(profiling_plot_dir, recursive = TRUE)
  }
  
  # ---------------------------------------------------------------------------
  # 18.1 Heatmap de correspondencia K = 3 vs K = 4
  # ---------------------------------------------------------------------------
  
  p_cross_heatmap <- ggplot(
    cross_k3_k4,
    aes(x = cluster_k4, y = cluster_k3, fill = prop_within_k4)
  ) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = paste0(Freq, "\n", round(prop_within_k4 * 100, 1), "%")),
      size = 3
    ) +
    theme_bw() +
    labs(
      title = "Correspondencia entre clusters K = 3 y K = 4",
      subtitle = "Porcentaje calculado dentro de cada cluster de K = 4",
      x = "Cluster en K = 4",
      y = "Cluster en K = 3",
      fill = "Prop. dentro de K4"
    )
  
  save_plot(
    p_cross_heatmap,
    file.path(profiling_plot_dir, "01_heatmap_cross_k3_k4.png"),
    8,
    5
  )
  
  # ---------------------------------------------------------------------------
  # 18.2 Comparación de tamaños de clusters
  # ---------------------------------------------------------------------------
  
  sizes_k3 <- profile_k3 |>
    dplyr::select(cluster = cluster_k3, n, pct) |>
    dplyr::mutate(k = "K = 3")
  
  sizes_k4 <- profile_k4 |>
    dplyr::select(cluster = cluster_k4, n, pct) |>
    dplyr::mutate(k = "K = 4")
  
  sizes_compare <- dplyr::bind_rows(sizes_k3, sizes_k4)
  
  p_sizes_compare <- ggplot(
    sizes_compare,
    aes(x = cluster, y = pct, fill = cluster)
  ) +
    geom_col() +
    facet_wrap(~ k, scales = "free_x") +
    geom_text(
      aes(label = paste0(round(pct, 1), "%")),
      vjust = -0.3,
      size = 3
    ) +
    theme_bw() +
    labs(
      title = "Tamaño relativo de clusters",
      x = "Cluster",
      y = "% de observaciones",
      fill = "Cluster"
    )
  
  save_plot(
    p_sizes_compare,
    file.path(profiling_plot_dir, "02_cluster_sizes_k3_k4.png"),
    9,
    5
  )
  
  # ---------------------------------------------------------------------------
  # 18.3 Perfil numérico comparativo K = 3
  # ---------------------------------------------------------------------------
  
  profile_k3_long <- profile_k3 |>
    tidyr::pivot_longer(
      cols = -cluster_k3,
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::filter(!variable %in% c("n", "pct"))
  
  p_profile_k3 <- ggplot(
    profile_k3_long,
    aes(x = cluster_k3, y = value, fill = cluster_k3)
  ) +
    geom_col() +
    facet_wrap(~ variable, scales = "free_y") +
    theme_bw() +
    labs(
      title = "Perfil numérico de clusters - K = 3",
      x = "Cluster",
      y = "Valor medio / proporción",
      fill = "Cluster"
    )
  
  save_plot(
    p_profile_k3,
    file.path(profiling_plot_dir, "03_numeric_profile_k3.png"),
    12,
    8
  )
  
  # ---------------------------------------------------------------------------
  # 18.4 Perfil numérico comparativo K = 4
  # ---------------------------------------------------------------------------
  
  profile_k4_long <- profile_k4 |>
    tidyr::pivot_longer(
      cols = -cluster_k4,
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::filter(!variable %in% c("n", "pct"))
  
  p_profile_k4 <- ggplot(
    profile_k4_long,
    aes(x = cluster_k4, y = value, fill = cluster_k4)
  ) +
    geom_col() +
    facet_wrap(~ variable, scales = "free_y") +
    theme_bw() +
    labs(
      title = "Perfil numérico de clusters - K = 4",
      x = "Cluster",
      y = "Valor medio / proporción",
      fill = "Cluster"
    )
  
  save_plot(
    p_profile_k4,
    file.path(profiling_plot_dir, "04_numeric_profile_k4.png"),
    12,
    8
  )
  
  # ---------------------------------------------------------------------------
  # 18.5 Perfil comparativo unificado K = 3 vs K = 4
  # ---------------------------------------------------------------------------
  
  profile_k3_unified <- profile_k3 |>
    dplyr::rename(cluster = cluster_k3) |>
    dplyr::mutate(k = "K = 3")
  
  profile_k4_unified <- profile_k4 |>
    dplyr::rename(cluster = cluster_k4) |>
    dplyr::mutate(k = "K = 4")
  
  profile_unified <- dplyr::bind_rows(profile_k3_unified, profile_k4_unified)
  
  profile_unified_long <- profile_unified |>
    tidyr::pivot_longer(
      cols = -c(cluster, k),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::filter(!variable %in% c("n", "pct"))
  
  p_profile_unified <- ggplot(
    profile_unified_long,
    aes(x = cluster, y = value, fill = cluster)
  ) +
    geom_col() +
    facet_grid(variable ~ k, scales = "free_y") +
    theme_bw() +
    labs(
      title = "Comparación visual de perfiles numéricos: K = 3 vs K = 4",
      x = "Cluster",
      y = "Valor medio / proporción",
      fill = "Cluster"
    )
  
  save_plot(
    p_profile_unified,
    file.path(profiling_plot_dir, "05_numeric_profile_k3_vs_k4.png"),
    12,
    14
  )
  
  # ---------------------------------------------------------------------------
  # 18.6 Comparación de severidad por cluster
  # ---------------------------------------------------------------------------
  
  if (file.exists(file.path(table_dir, "profiling_Severity_k3.csv")) &&
      file.exists(file.path(table_dir, "profiling_Severity_k4.csv"))) {
    
    severity_k3 <- read.csv(file.path(table_dir, "profiling_Severity_k3.csv"))
    severity_k4 <- read.csv(file.path(table_dir, "profiling_Severity_k4.csv"))
    
    names(severity_k3)[1] <- "cluster"
    names(severity_k4)[1] <- "cluster"
    
    severity_k3$k <- "K = 3"
    severity_k4$k <- "K = 4"
    
    severity_compare <- dplyr::bind_rows(severity_k3, severity_k4)
    
    p_severity_compare <- ggplot(
      severity_compare,
      aes(x = cluster, y = prop, fill = Severity)
    ) +
      geom_col(position = "fill") +
      facet_wrap(~ k, scales = "free_x") +
      theme_bw() +
      scale_y_continuous(labels = scales::percent) +
      labs(
        title = "Distribución de severidad por cluster",
        x = "Cluster",
        y = "Proporción",
        fill = "Severity"
      )
    
    save_plot(
      p_severity_compare,
      file.path(profiling_plot_dir, "06_severity_distribution_k3_k4.png"),
      9,
      5
    )
  }
  
  # ---------------------------------------------------------------------------
  # 18.7 Comparación de Weather_Group por cluster
  # ---------------------------------------------------------------------------
  
  if (file.exists(file.path(table_dir, "profiling_Weather_Group_k3.csv")) &&
      file.exists(file.path(table_dir, "profiling_Weather_Group_k4.csv"))) {
    
    weather_k3 <- read.csv(file.path(table_dir, "profiling_Weather_Group_k3.csv"))
    weather_k4 <- read.csv(file.path(table_dir, "profiling_Weather_Group_k4.csv"))
    
    names(weather_k3)[1] <- "cluster"
    names(weather_k4)[1] <- "cluster"
    
    weather_k3$k <- "K = 3"
    weather_k4$k <- "K = 4"
    
    weather_compare <- dplyr::bind_rows(weather_k3, weather_k4)
    
    p_weather_compare <- ggplot(
      weather_compare,
      aes(x = cluster, y = prop, fill = Weather_Group)
    ) +
      geom_col(position = "fill") +
      facet_wrap(~ k, scales = "free_x") +
      theme_bw() +
      scale_y_continuous(labels = scales::percent) +
      labs(
        title = "Distribución de Weather_Group por cluster",
        x = "Cluster",
        y = "Proporción",
        fill = "Weather_Group"
      )
    
    save_plot(
      p_weather_compare,
      file.path(profiling_plot_dir, "07_weather_distribution_k3_k4.png"),
      10,
      6
    )
  }
  
  # ---------------------------------------------------------------------------
  # 18.8 Comparación de Sunrise_Sunset por cluster
  # ---------------------------------------------------------------------------
  
  if (file.exists(file.path(table_dir, "profiling_Sunrise_Sunset_k3.csv")) &&
      file.exists(file.path(table_dir, "profiling_Sunrise_Sunset_k4.csv"))) {
    
    light_k3 <- read.csv(file.path(table_dir, "profiling_Sunrise_Sunset_k3.csv"))
    light_k4 <- read.csv(file.path(table_dir, "profiling_Sunrise_Sunset_k4.csv"))
    
    names(light_k3)[1] <- "cluster"
    names(light_k4)[1] <- "cluster"
    
    light_k3$k <- "K = 3"
    light_k4$k <- "K = 4"
    
    light_compare <- dplyr::bind_rows(light_k3, light_k4)
    
    p_light_compare <- ggplot(
      light_compare,
      aes(x = cluster, y = prop, fill = Sunrise_Sunset)
    ) +
      geom_col(position = "fill") +
      facet_wrap(~ k, scales = "free_x") +
      theme_bw() +
      scale_y_continuous(labels = scales::percent) +
      labs(
        title = "Distribución día/noche por cluster",
        x = "Cluster",
        y = "Proporción",
        fill = "Sunrise_Sunset"
      )
    
    save_plot(
      p_light_compare,
      file.path(profiling_plot_dir, "08_light_distribution_k3_k4.png"),
      9,
      5
    )
  }
}
