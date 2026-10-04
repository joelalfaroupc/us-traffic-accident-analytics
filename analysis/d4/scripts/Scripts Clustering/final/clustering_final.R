# ==============================================================================
# PMAAD - Clustering final para profiling
# Método final: FAMD + submuestra jerárquica + expansión
# Solución fijada: 5 clusters
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# 0. Configuración
# ------------------------------------------------------------------------------
file_base <- "Datasets/clustering/04_base_clustering.csv"
out_dir   <- "clustering_final_results"

ncp_famd <- 8
sample_size_target <- 1000
selected_k <- 5

# ------------------------------------------------------------------------------
# 1. Librerías
# ------------------------------------------------------------------------------
packages <- c(
  "dplyr", "ggplot2", "FactoMineR", "factoextra",
  "forcats", "readr", "tibble"
)

new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages, dependencies = TRUE)
invisible(lapply(packages, require, character.only = TRUE))

# ------------------------------------------------------------------------------
# 2. Carpetas
# ------------------------------------------------------------------------------
if (!file.exists(file_base)) stop(paste("No se encuentra la base:", file_base))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
rds_dir   <- file.path(out_dir, "rds_objects")

if (!dir.exists(plot_dir))  dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)
if (!dir.exists(rds_dir))   dir.create(rds_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------------------------
save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  print(plot_obj)
  ggsave(filename = filename, plot = plot_obj, width = width, height = height)
}

safe_as_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

collapse_rare_levels <- function(x, min_prop = 0.01, other_label = "Other") {
  x <- as.factor(x)
  prop_tab <- prop.table(table(x))
  rare_levels <- names(prop_tab[prop_tab < min_prop])
  x <- forcats::fct_other(x, keep = setdiff(levels(x), rare_levels), other_level = other_label)
  droplevels(x)
}

sample_stratified_ids <- function(df, size_target = 1000, strata_col = "State", seed = 123) {
  set.seed(seed)
  n <- nrow(df)
  if (!(strata_col %in% names(df))) {
    return(sort(sample(seq_len(n), size = min(size_target, n), replace = FALSE)))
  }
  
  strata <- df[[strata_col]]
  split_ids <- split(seq_len(n), strata)
  prop <- sapply(split_ids, length) / n
  take <- pmax(10, round(prop * size_target))
  take <- pmin(take, sapply(split_ids, length))
  
  ids <- unlist(
    mapply(
      function(ix, k) sample(ix, size = k, replace = FALSE),
      split_ids, take, SIMPLIFY = FALSE
    )
  )
  
  sort(unique(as.integer(ids)))
}

top_levels_by_cluster <- function(df, cluster_vec, var_name, top_n = 3) {
  tab <- as.data.frame(prop.table(table(cluster_vec, df[[var_name]]), 1))
  names(tab) <- c("cluster_final", "level", "prop")
  tab$variable <- var_name
  tab |>
    dplyr::group_by(variable, cluster_final) |>
    dplyr::slice_max(order_by = prop, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()
}

# ------------------------------------------------------------------------------
# 4. Carga
# ------------------------------------------------------------------------------
base_raw <- read.csv(file_base, stringsAsFactors = FALSE)

cat("\n==============================\n")
cat("CARGA DE DATOS\n")
cat("==============================\n")
cat("Dimensión base:", nrow(base_raw), "x", ncol(base_raw), "\n")

char_cols <- names(base_raw)[sapply(base_raw, is.character)]
if (length(char_cols) > 0) {
  base_raw[char_cols] <- lapply(base_raw[char_cols], function(x) {
    x[x == ""] <- NA
    x
  })
}

# ------------------------------------------------------------------------------
# 5. Variables auxiliares / trazabilidad
# ------------------------------------------------------------------------------
trace_candidates <- intersect(
  c("ID", "ID_ROW", "Start_Time", "Description", "State", "County", "Start_Lat", "Start_Lng"),
  names(base_raw)
)
base_trace <- base_raw[, trace_candidates, drop = FALSE]

infra_cols <- intersect(
  c("Amenity", "Crossing", "Give_Way", "Junction", "No_Exit",
    "Railway", "Station", "Stop", "Traffic_Calming", "Traffic_Signal"),
  names(base_raw)
)

if (!"Has_Infrastructure" %in% names(base_raw) && length(infra_cols) > 0) {
  infra_num <- sapply(base_raw[infra_cols], function(x) {
    if (is.logical(x)) return(as.integer(x))
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
      x <- ifelse(tolower(x) %in% c("true", "1", "yes"), 1L, 0L)
      return(as.integer(x))
    }
    as.integer(x)
  })
  if (is.vector(infra_num)) infra_num <- matrix(infra_num, ncol = 1)
  base_raw$Has_Infrastructure <- ifelse(rowSums(infra_num, na.rm = TRUE) > 0, "Yes", "No")
}

if ("Start_Time" %in% names(base_raw)) {
  dt_try <- suppressWarnings(as.POSIXct(base_raw$Start_Time, tz = "UTC"))
  if (sum(!is.na(dt_try)) > 0) {
    base_raw$Start_Hour <- as.integer(format(dt_try, "%H"))
    base_raw$Start_Month <- as.integer(format(dt_try, "%m"))
    base_raw$Is_Weekend <- ifelse(weekdays(dt_try) %in% c("Saturday", "Sunday"), "Yes", "No")
    
    base_raw$Time_Period <- cut(
      base_raw$Start_Hour,
      breaks = c(-Inf, 5, 11, 17, 21, Inf),
      labels = c("Night", "Morning", "Afternoon", "Evening", "LateNight"),
      right = TRUE
    )
    
    base_raw$Season <- dplyr::case_when(
      base_raw$Start_Month %in% c(12, 1, 2) ~ "Winter",
      base_raw$Start_Month %in% c(3, 4, 5) ~ "Spring",
      base_raw$Start_Month %in% c(6, 7, 8) ~ "Summer",
      base_raw$Start_Month %in% c(9, 10, 11) ~ "Autumn",
      TRUE ~ NA_character_
    )
  }
}

# ------------------------------------------------------------------------------
# 6. Base activa para clustering
# ------------------------------------------------------------------------------
drop_vars <- intersect(
  c(
    "ID", "ID_ROW", "Description", "Start_Time",
    "County", "State",
    "Start_Lat", "Start_Lng",
    "Civil_Twilight", "Nautical_Twilight", "Astronomical_Twilight",
    "Amenity", "Crossing", "Give_Way", "Junction", "No_Exit",
    "Railway", "Station", "Stop", "Traffic_Calming", "Traffic_Signal"
  ),
  names(base_raw)
)

base_active <- base_raw[, setdiff(names(base_raw), drop_vars), drop = FALSE]

duration_raw_name <- intersect(c("Duration_Hours"), names(base_active))
if (length(duration_raw_name) > 0) {
  base_active$Duration_Log <- log1p(base_active[[duration_raw_name[1]]])
  base_active[[duration_raw_name[1]]] <- NULL
}

temp_keep <- intersect(c("Temperature.C.", "Temperature.F."), names(base_active))
wc_keep   <- intersect(c("Wind_Chill.C.", "Wind_Chill.F."), names(base_active))

preferred_keep <- unique(intersect(
  c(
    "Severity",
    temp_keep,
    wc_keep,
    "Humidity...", "Pressure.in.", "Wind_Speed.mph.", "Precipitation.in.",
    "Duration_Log", "Start_Hour", "Start_Month",
    "Weather_Group", "Wind_Group", "Sunrise_Sunset",
    "Is_Weekend", "Has_Infrastructure", "Time_Period", "Season"
  ),
  names(base_active)
))

if (length(preferred_keep) > 0) {
  base_active <- base_active[, preferred_keep, drop = FALSE]
}

# ------------------------------------------------------------------------------
# 7. Tipos y limpieza final
# ------------------------------------------------------------------------------
if ("Severity" %in% names(base_active)) {
  base_active$Severity <- factor(base_active$Severity)
}

cat_candidates <- intersect(
  c("Weather_Group", "Wind_Group", "Sunrise_Sunset",
    "Is_Weekend", "Has_Infrastructure", "Time_Period", "Season"),
  names(base_active)
)

for (nm in names(base_active)) {
  if (nm %in% cat_candidates || nm == "Severity") {
    base_active[[nm]] <- as.factor(base_active[[nm]])
  } else {
    base_active[[nm]] <- safe_as_numeric(base_active[[nm]])
  }
}

for (nm in cat_candidates) {
  if (nm %in% names(base_active)) {
    base_active[[nm]] <- collapse_rare_levels(base_active[[nm]], min_prop = 0.01, other_label = "Other")
  }
}

na_cols <- names(base_active)[colSums(is.na(base_active)) > 0]
if (length(na_cols) > 0) {
  cat("\nColumnas eliminadas por contener NA:\n")
  print(na_cols)
  base_active <- base_active[, setdiff(names(base_active), na_cols), drop = FALSE]
}

num_cols <- names(base_active)[sapply(base_active, is.numeric)]
zero_var_num <- num_cols[sapply(base_active[num_cols], function(x) sd(x, na.rm = TRUE) == 0)]
if (length(zero_var_num) > 0) {
  cat("\nNuméricas eliminadas por varianza cero:\n")
  print(zero_var_num)
  base_active <- base_active[, setdiff(names(base_active), zero_var_num), drop = FALSE]
}

fac_cols <- names(base_active)[sapply(base_active, is.factor)]
one_level_fac <- fac_cols[sapply(base_active[fac_cols], function(x) length(unique(na.omit(x))) <= 1)]
if (length(one_level_fac) > 0) {
  cat("\nCategóricas eliminadas por un solo nivel:\n")
  print(one_level_fac)
  base_active <- base_active[, setdiff(names(base_active), one_level_fac), drop = FALSE]
}

num_cols <- names(base_active)[sapply(base_active, is.numeric)]
fac_cols <- names(base_active)[sapply(base_active, is.factor)]

cat("\n==============================\n")
cat("BASE ACTIVA PARA CLUSTERING\n")
cat("==============================\n")
cat("Dimensión base activa:", nrow(base_active), "x", ncol(base_active), "\n")
cat("Numéricas:", length(num_cols), "\n")
cat("Categóricas:", length(fac_cols), "\n")
cat("Missings totales:", sum(is.na(base_active)), "\n")

write.csv(base_active, file.path(out_dir, "00_base_active_clustering_final.csv"), row.names = FALSE)

# Base para profiling
profile_keep <- intersect(
  c("State", "County", "Start_Lat", "Start_Lng", "Duration_Hours", names(base_active)),
  names(base_raw)
)
base_profile <- base_raw[, unique(profile_keep), drop = FALSE]

# ------------------------------------------------------------------------------
# 8. FAMD
# ------------------------------------------------------------------------------
cat("\n==============================\n")
cat("FAMD\n")
cat("==============================\n")

if (length(num_cols) == 0 || length(fac_cols) == 0) {
  stop("La base activa no tiene mezcla suficiente de variables numéricas y categóricas para FAMD.")
}

famd_res <- FactoMineR::FAMD(base_active, ncp = ncp_famd, graph = FALSE)

p_scree <- factoextra::fviz_eig(famd_res, addlabels = TRUE) +
  theme_bw() +
  ggtitle("FAMD - Varianza explicada")
save_plot(p_scree, file.path(plot_dir, "01_famd_screeplot.png"), 8, 6)

p_var <- factoextra::fviz_famd_var(famd_res, repel = TRUE) +
  theme_bw() +
  ggtitle("FAMD - Variables")
save_plot(p_var, file.path(plot_dir, "02_famd_variables.png"), 9, 7)

coord_full <- as.data.frame(famd_res$ind$coord[, 1:ncp_famd, drop = FALSE])
coord_full$ID_ROW_PLOT <- seq_len(nrow(coord_full))
p_ind <- ggplot(coord_full, aes(x = Dim.1, y = Dim.2)) +
  geom_point(alpha = 0.35, size = 0.9) +
  theme_bw() +
  labs(title = "FAMD - Individuos", x = "Dimensión 1", y = "Dimensión 2")
save_plot(p_ind, file.path(plot_dir, "03_famd_individuos.png"), 8, 6)

# ------------------------------------------------------------------------------
# 9. Submuestra jerárquica
# ------------------------------------------------------------------------------
cat("\n==============================\n")
cat("SUBMUESTRA JERÁRQUICA\n")
cat("==============================\n")

sample_ids <- sample_stratified_ids(base_raw, size_target = sample_size_target, strata_col = "State", seed = 123)
coord_sample <- coord_full[sample_ids, 1:ncp_famd, drop = FALSE]

cat("Tamaño submuestra:", nrow(coord_sample), "\n")

dist_sample <- dist(coord_sample)
hc_sample <- hclust(dist_sample, method = "ward.D2")

# Guardar un dendrograma del corte seleccionado
png(
  filename = file.path(plot_dir, paste0("04_dendrogram_selected_k", selected_k, ".png")),
  width = 1200, height = 700, res = 120
)
plot(
  hc_sample,
  labels = FALSE,
  hang = -1,
  main = paste0("Dendrograma submuestra | k = ", selected_k),
  xlab = "", sub = "", cex = 0.6
)
rect.hclust(hc_sample, k = selected_k, border = 2:(selected_k + 1))
dev.off()

# ------------------------------------------------------------------------------
# 10. Expansión al total con k-means
# ------------------------------------------------------------------------------
cat("\n==============================\n")
cat("EXPANSIÓN AL TOTAL\n")
cat("==============================\n")

X_full <- as.matrix(coord_full[, 1:ncp_famd, drop = FALSE])
X_sample <- as.matrix(coord_sample[, 1:ncp_famd, drop = FALSE])

cl_sample <- cutree(hc_sample, k = selected_k)

init_centers <- do.call(
  rbind,
  lapply(sort(unique(cl_sample)), function(cl) {
    colMeans(X_sample[cl_sample == cl, , drop = FALSE])
  })
)

set.seed(100 + selected_k)
km_res <- kmeans(X_full, centers = init_centers, iter.max = 60, nstart = 1)
cluster_final <- km_res$cluster

cat("Tamaños de cluster finales:\n")
print(table(cluster_final))

# ------------------------------------------------------------------------------
# 11. Outputs finales para profiling
# ------------------------------------------------------------------------------
# Base final clusterizada para el compañero
base_profile_clustered <- base_profile
base_profile_clustered$cluster_final <- factor(cluster_final)

write.csv(
  base_profile_clustered,
  file.path(out_dir, "10_base_final_para_profiling_k5.csv"),
  row.names = FALSE
)

# Base activa clusterizada
base_active_clustered <- base_active
base_active_clustered$cluster_final <- factor(cluster_final)

write.csv(
  base_active_clustered,
  file.path(out_dir, "10_base_activa_clusterizada_k5.csv"),
  row.names = FALSE
)

# Trazabilidad mínima
if (nrow(base_trace) == nrow(base_active_clustered)) {
  base_trace_clustered <- cbind(base_trace, cluster_final = factor(cluster_final))
  write.csv(
    base_trace_clustered,
    file.path(out_dir, "10_base_trace_clusterizada_k5.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 12. Perfiles numéricos y categóricos
# ------------------------------------------------------------------------------
num_prof_vars <- names(base_profile_clustered)[sapply(base_profile_clustered, is.numeric)]
if (length(num_prof_vars) > 0) {
  numeric_profile <- base_profile_clustered |>
    dplyr::group_by(cluster_final) |>
    dplyr::summarise(dplyr::across(
      dplyr::all_of(num_prof_vars),
      ~ mean(.x, na.rm = TRUE)
    ))
  
  write.csv(
    numeric_profile,
    file.path(table_dir, "11_numeric_profile_k5.csv"),
    row.names = FALSE
  )
}

cat_prof_vars <- names(base_profile_clustered)[sapply(base_profile_clustered, function(x) is.factor(x) || is.character(x))]
cat_prof_vars <- setdiff(cat_prof_vars, c("Description", "Start_Time", "ID", "ID_ROW"))

if (length(cat_prof_vars) > 0) {
  top_list <- lapply(cat_prof_vars, function(v) top_levels_by_cluster(base_profile_clustered, base_profile_clustered$cluster_final, v, top_n = 3))
  categorical_top <- do.call(rbind, top_list)
  
  write.csv(
    categorical_top,
    file.path(table_dir, "12_categorical_top_levels_k5.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 13. Gráficos finales útiles
# ------------------------------------------------------------------------------
# Plano FAMD coloreado
df_plot <- data.frame(
  Dim.1 = X_full[, 1],
  Dim.2 = X_full[, 2],
  cluster = factor(cluster_final)
)

p_clusters <- ggplot(df_plot, aes(x = Dim.1, y = Dim.2, color = cluster)) +
  geom_point(alpha = 0.40, size = 0.9) +
  theme_bw() +
  labs(
    title = "FAMD + clustering final | k = 5",
    x = "Dimensión 1",
    y = "Dimensión 2",
    color = "Cluster"
  )

save_plot(p_clusters, file.path(plot_dir, "06_famd_clusters_k5.png"), 8, 6)

# Tamaños de clusters
sizes_k <- as.data.frame(table(cluster_final))
names(sizes_k) <- c("cluster", "n")

p_sizes <- ggplot(sizes_k, aes(x = factor(cluster), y = n, fill = factor(cluster))) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Tamaño de clusters | k = 5",
    x = "Cluster",
    y = "Número de observaciones",
    fill = "Cluster"
  )

save_plot(p_sizes, file.path(plot_dir, "07_cluster_sizes_k5.png"), 8, 5)

# Mapa aproximado
if (all(c("Start_Lat", "Start_Lng") %in% names(base_profile_clustered))) {
  df_map <- data.frame(
    Start_Lat = base_profile_clustered$Start_Lat,
    Start_Lng = base_profile_clustered$Start_Lng,
    cluster = factor(cluster_final)
  )
  
  p_map <- ggplot(df_map, aes(x = Start_Lng, y = Start_Lat, color = cluster)) +
    geom_point(alpha = 0.40, size = 0.8) +
    theme_bw() +
    labs(
      title = "Distribución geográfica aproximada | k = 5",
      x = "Longitud",
      y = "Latitud",
      color = "Cluster"
    )
  
  save_plot(p_map, file.path(plot_dir, "08_map_clusters_k5.png"), 8, 5)
}

# ------------------------------------------------------------------------------
# 14. Guardar objetos
# ------------------------------------------------------------------------------
saveRDS(famd_res, file.path(rds_dir, "famd_res.rds"))
saveRDS(hc_sample, file.path(rds_dir, "hc_sample_subsample.rds"))
saveRDS(sample_ids, file.path(rds_dir, "sample_ids_for_dendrogram.rds"))
saveRDS(cluster_final, file.path(rds_dir, "cluster_final_k5.rds"))

# ------------------------------------------------------------------------------
# 15. Resumen final
# ------------------------------------------------------------------------------
sink(file.path(out_dir, "00_RESUMEN_FINAL_K5.txt"))
cat("====================================\n")
cat("RESULTADOS FINALES CLUSTERING K=5\n")
cat("====================================\n\n")

cat("Base original usada:\n")
cat("-", file_base, "\n\n")

cat("Base activa para clustering:\n")
cat("- Filas:", nrow(base_active), "\n")
cat("- Columnas:", ncol(base_active), "\n")
cat("- Numéricas:", length(num_cols), "\n")
cat("- Categóricas:", length(fac_cols), "\n\n")

cat("Parámetros:\n")
cat("- ncp_famd:", ncp_famd, "\n")
cat("- sample_size_target:", sample_size_target, "\n")
cat("- selected_k:", selected_k, "\n\n")

cat("Criterio del dendrograma:\n")
cat("- Distancia: Euclídea sobre componentes FAMD\n")
cat("- Enlace: Ward.D2\n\n")

cat("Tamaños de cluster finales:\n")
print(table(cluster_final))
cat("\n")

cat("Archivos clave para profiling:\n")
cat("- 10_base_final_para_profiling_k5.csv\n")
cat("- 11_numeric_profile_k5.csv\n")
cat("- 12_categorical_top_levels_k5.csv\n")
cat("\n")

cat("Gráficos clave:\n")
cat("- 06_famd_clusters_k5.png\n")
cat("- 07_cluster_sizes_k5.png\n")
cat("- 08_map_clusters_k5.png\n")
sink()

cat("\nScript ejecutado correctamente.\n")
cat("Resultados guardados en:", out_dir, "\n")