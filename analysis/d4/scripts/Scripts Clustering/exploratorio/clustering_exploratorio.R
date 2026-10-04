# ==============================================================================
# PMAAD - Re-clustering visual v2 (rápido) para base grande (~15k filas)
# ------------------------------------------------------------------------------
# Variante: REGLA DE TAMAÑO MÍNIMO
# ------------------------------------------------------------------------------
# Mantiene la lógica de corte por altura del dendrograma y añade:
#   - una regla de tamaño mínimo de cluster
#   - si un cluster queda por debajo del umbral, sus observaciones se reasignan
#     al centroide más cercano entre los clusters válidos
# ------------------------------------------------------------------------------
# Mejoras incluidas:
#   1) NO colapsa Severity a "Other"
#   2) Controla el peso de Duration_Hours usando log1p(Duration_Hours)
#   3) Saca Start_Lat / Start_Lng del clustering activo para que no dominen
#   4) Hace FAMD sobre toda la base activa
#   5) Hace clustering jerárquico SOLO en una submuestra estratificada
#   6) Expande al total con k-means sobre coordenadas FAMD
#   7) Aplica regla de tamaño mínimo
#   8) Genera perfiles para varias alturas
#   9) Guarda todo en carpetas
#
# Criterio para el dendrograma:
#   - Distancia euclídea sobre componentes FAMD
#   - Método de enlace Ward.D2
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# 0. Configuración
# ------------------------------------------------------------------------------
file_base <- "Datasets/clustering/04_base_clustering.csv"
out_dir   <- "clustering_famd_fast_results"

# Alturas a evaluar visualmente
h_grid <- c(32, 36, 40)

# Altura provisional seleccionada
selected_h <- 36

# Nº de componentes FAMD
ncp_famd <- 8

# Tamaño de submuestra para dendrograma visual
sample_size_target <- 1000

# Regla de tamaño mínimo
# Puedes usar porcentaje o absoluto; aquí mando porcentaje de la base total
min_cluster_prop <- 0.0   # 5%
# Si prefieres fijarlo manualmente, usa por ejemplo: min_cluster_size_manual <- 400
min_cluster_size_manual <- NA

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
# 2. Validaciones y carpetas
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

reassign_small_clusters <- function(X, clusters, min_size) {
  clusters <- as.integer(as.factor(clusters))
  tb <- table(clusters)
  small_clusters <- as.integer(names(tb)[tb < min_size])
  big_clusters   <- as.integer(names(tb)[tb >= min_size])

  if (length(small_clusters) == 0 || length(big_clusters) == 0) {
    return(list(cluster = clusters, changed = FALSE, small_clusters = small_clusters))
  }

  X <- as.matrix(X)
  new_clusters <- clusters

  # centroides de clusters válidos
  centers_big <- do.call(
    rbind,
    lapply(big_clusters, function(cl) colMeans(X[new_clusters == cl, , drop = FALSE]))
  )
  rownames(centers_big) <- big_clusters

  # reasignar observaciones de clusters pequeños al centroide grande más cercano
  ids_small <- which(new_clusters %in% small_clusters)

  if (length(ids_small) > 0) {
    for (i in ids_small) {
      xi <- X[i, , drop = FALSE]
      dists <- apply(centers_big, 1, function(cc) sum((xi - cc)^2))
      best_big <- as.integer(names(which.min(dists)))
      new_clusters[i] <- best_big
    }
  }

  # reindexar clusters a 1..K
  new_clusters <- as.integer(as.factor(new_clusters))

  list(
    cluster = new_clusters,
    changed = TRUE,
    small_clusters = small_clusters
  )
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

write.csv(base_active, file.path(out_dir, "00_base_active_clustering_v3.csv"), row.names = FALSE)

profile_keep <- intersect(
  c("State", "County", "Start_Lat", "Start_Lng", "Duration_Hours", names(base_active)),
  names(base_raw)
)
base_profile <- base_raw[, unique(profile_keep), drop = FALSE]

# ------------------------------------------------------------------------------
# 8. FAMD sobre toda la base
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
# 9. Dendrograma visual sobre submuestra
# ------------------------------------------------------------------------------
cat("\n==============================\n")
cat("SUBMUESTRA PARA DENDROGRAMA\n")
cat("==============================\n")

sample_ids <- sample_stratified_ids(base_raw, size_target = sample_size_target, strata_col = "State", seed = 123)
coord_sample <- coord_full[sample_ids, 1:ncp_famd, drop = FALSE]

cat("Tamaño submuestra dendrograma:", nrow(coord_sample), "\n")

dist_sample <- dist(coord_sample)
hc_sample <- hclust(dist_sample, method = "ward.D2")

height_df <- data.frame(step = seq_along(hc_sample$height), height = hc_sample$height)
height_df$jump <- c(NA, diff(hc_sample$height))
height_df2 <- height_df[!is.na(height_df$jump), ]

p_height <- ggplot(height_df2, aes(x = step, y = jump)) +
  geom_line() +
  geom_point(size = 1) +
  theme_bw() +
  labs(title = "Saltos de altura del dendrograma (submuestra)",
       x = "Paso de fusión", y = "Incremento de altura")
save_plot(p_height, file.path(plot_dir, "04_dendrogram_height_jumps.png"), 8, 5)

tail_n <- min(150, nrow(height_df2))
height_tail <- tail(height_df2, tail_n)

p_height_tail <- ggplot(height_tail, aes(x = step, y = jump)) +
  geom_line() +
  geom_point(size = 1.2) +
  theme_bw() +
  labs(title = "Saltos de altura - últimos pasos",
       x = "Paso de fusión", y = "Incremento de altura")
save_plot(p_height_tail, file.path(plot_dir, "04b_dendrogram_height_jumps_tail.png"), 8, 5)

for (h_cut in h_grid) {
  cluster_h <- cutree(hc_sample, h = h_cut)
  n_clusters_h <- length(unique(cluster_h))

  png(
    filename = file.path(plot_dir, paste0("05_dendrogram_sample_h", gsub("\\.", "_", as.character(h_cut)), ".png")),
    width = 1200, height = 700, res = 120
  )
  plot(
    hc_sample,
    labels = FALSE,
    hang = -1,
    main = paste0("Dendrograma visual (submuestra) | h = ", h_cut,
                  " | clusters = ", n_clusters_h),
    xlab = "", sub = "", cex = 0.6
  )
  rect.hclust(hc_sample, h = h_cut, border = 2:(n_clusters_h + 1))
  abline(h = h_cut, col = "red", lty = 2, lwd = 2)
  dev.off()
}

# ------------------------------------------------------------------------------
# 10. Expansión al total con k-means + regla de tamaño mínimo
# ------------------------------------------------------------------------------
cat("\n==============================\n")
cat("EXPANSIÓN AL TOTAL CON K-MEANS EN FAMD + REGLA TAMAÑO MÍNIMO\n")
cat("==============================\n")

cluster_sizes_all <- list()
all_cluster_assignments <- list()
rules_summary <- list()

X_full <- as.matrix(coord_full[, 1:ncp_famd, drop = FALSE])
X_sample <- as.matrix(coord_sample[, 1:ncp_famd, drop = FALSE])

# umbral mínimo
if (is.na(min_cluster_size_manual)) {
  min_cluster_size <- ceiling(min_cluster_prop * nrow(X_full))
} else {
  min_cluster_size <- as.integer(min_cluster_size_manual)
}

cat("Tamaño mínimo exigido:", min_cluster_size, "\n")

for (h_cut in h_grid) {
  cl_sample <- cutree(hc_sample, h = h_cut)
  n_clusters_h <- length(unique(cl_sample))

  init_centers <- do.call(
    rbind,
    lapply(sort(unique(cl_sample)), function(cl) {
      colMeans(X_sample[cl_sample == cl, , drop = FALSE])
    })
  )

  set.seed(100 + round(h_cut * 10))
  km_res <- kmeans(X_full, centers = init_centers, iter.max = 60, nstart = 1)
  cluster_k <- km_res$cluster

  # aplicar regla de tamaño mínimo
  reassign_res <- reassign_small_clusters(X_full, cluster_k, min_cluster_size)
  cluster_k_final <- reassign_res$cluster

  all_cluster_assignments[[as.character(h_cut)]] <- cluster_k_final

  sizes_before <- as.data.frame(table(cluster_k))
  names(sizes_before) <- c("cluster", "n_before")
  sizes_before$h_cut <- h_cut

  sizes_k <- as.data.frame(table(cluster_k_final))
  names(sizes_k) <- c("cluster", "n")
  sizes_k$h_cut <- h_cut
  sizes_k$n_clusters_initial <- n_clusters_h
  sizes_k$n_clusters_final <- length(unique(cluster_k_final))
  cluster_sizes_all[[as.character(h_cut)]] <- sizes_k

  rules_summary[[as.character(h_cut)]] <- data.frame(
    h_cut = h_cut,
    n_clusters_initial = n_clusters_h,
    n_clusters_final = length(unique(cluster_k_final)),
    min_cluster_size = min_cluster_size,
    small_clusters_detected = paste(reassign_res$small_clusters, collapse = ","),
    reassignment_applied = reassign_res$changed,
    stringsAsFactors = FALSE
  )

  # plano FAMD
  df_plot <- data.frame(
    Dim.1 = X_full[, 1],
    Dim.2 = X_full[, 2],
    cluster = factor(cluster_k_final)
  )

  p_clusters <- ggplot(df_plot, aes(x = Dim.1, y = Dim.2, color = cluster)) +
    geom_point(alpha = 0.40, size = 0.9) +
    theme_bw() +
    labs(title = paste0("FAMD + clustering expandido | h = ", h_cut,
                        " | clusters finales = ", length(unique(cluster_k_final))),
         x = "Dimensión 1", y = "Dimensión 2", color = "Cluster")

  save_plot(p_clusters, file.path(plot_dir, paste0("06_famd_clusters_h", gsub("\\.", "_", as.character(h_cut)), ".png")), 8, 6)

  # tamaños
  p_sizes <- ggplot(sizes_k, aes(x = factor(cluster), y = n, fill = factor(cluster))) +
    geom_col() +
    theme_bw() +
    labs(title = paste0("Tamaño de clusters | h = ", h_cut,
                        " | clusters finales = ", length(unique(cluster_k_final))),
         x = "Cluster", y = "Número de observaciones", fill = "Cluster")

  save_plot(p_sizes, file.path(plot_dir, paste0("07_cluster_sizes_h", gsub("\\.", "_", as.character(h_cut)), ".png")), 8, 5)

  # mapa
  if (all(c("Start_Lat", "Start_Lng") %in% names(base_profile))) {
    df_map <- data.frame(
      Start_Lat = base_profile$Start_Lat,
      Start_Lng = base_profile$Start_Lng,
      cluster = factor(cluster_k_final)
    )

    p_map <- ggplot(df_map, aes(x = Start_Lng, y = Start_Lat, color = cluster)) +
      geom_point(alpha = 0.40, size = 0.8) +
      theme_bw() +
      labs(title = paste0("Distribución geográfica aproximada | h = ", h_cut,
                          " | clusters finales = ", length(unique(cluster_k_final))),
           x = "Longitud", y = "Latitud", color = "Cluster")

    save_plot(p_map, file.path(plot_dir, paste0("08_map_clusters_h", gsub("\\.", "_", as.character(h_cut)), ".png")), 8, 5)
  }

  # exportar base y perfiles
  base_profile_h <- base_profile
  base_profile_h$cluster_final <- factor(cluster_k_final)

  write.csv(base_profile_h,
            file.path(out_dir, paste0("10_base_profile_clustered_h", gsub("\\.", "_", as.character(h_cut)), ".csv")),
            row.names = FALSE)

  num_prof_vars <- names(base_profile_h)[sapply(base_profile_h, is.numeric)]
  if (length(num_prof_vars) > 0) {
    numeric_profile <- base_profile_h |>
      dplyr::group_by(cluster_final) |>
      dplyr::summarise(dplyr::across(
        dplyr::all_of(num_prof_vars),
        ~ mean(.x, na.rm = TRUE)
      ))

    write.csv(numeric_profile,
              file.path(table_dir, paste0("11_numeric_profile_h", gsub("\\.", "_", as.character(h_cut)), ".csv")),
              row.names = FALSE)
  }

  cat_prof_vars <- names(base_profile_h)[sapply(base_profile_h, function(x) is.factor(x) || is.character(x))]
  cat_prof_vars <- setdiff(cat_prof_vars, c("Description", "Start_Time", "ID", "ID_ROW"))

  if (length(cat_prof_vars) > 0) {
    top_list <- lapply(cat_prof_vars, function(v) top_levels_by_cluster(base_profile_h, base_profile_h$cluster_final, v, top_n = 3))
    categorical_top <- do.call(rbind, top_list)

    write.csv(categorical_top,
              file.path(table_dir, paste0("12_categorical_top_levels_h", gsub("\\.", "_", as.character(h_cut)), ".csv")),
              row.names = FALSE)
  }
}

cluster_sizes_df <- do.call(rbind, cluster_sizes_all)
write.csv(cluster_sizes_df, file.path(table_dir, "09_cluster_sizes_all_h.csv"), row.names = FALSE)

rules_summary_df <- do.call(rbind, rules_summary)
write.csv(rules_summary_df, file.path(table_dir, "09b_min_size_rule_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 11. Selección provisional y exportación final
# ------------------------------------------------------------------------------
cluster_final <- all_cluster_assignments[[as.character(selected_h)]]

cat("\n==============================\n")
cat("SELECCIÓN PROVISIONAL POR ALTURA\n")
cat("==============================\n")
cat("h seleccionada provisionalmente:", selected_h, "\n")
cat("Número de clusters resultante:", length(unique(cluster_final)), "\n")
print(table(cluster_final))

base_final_selected <- base_profile
base_final_selected$cluster_final <- factor(cluster_final)

write.csv(base_final_selected,
          file.path(out_dir, paste0("13_base_final_selected_h", gsub("\\.", "_", as.character(selected_h)), ".csv")),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 12. Guardar objetos
# ------------------------------------------------------------------------------
saveRDS(famd_res, file.path(rds_dir, "famd_res.rds"))
saveRDS(hc_sample, file.path(rds_dir, "hc_sample_subsample.rds"))
saveRDS(sample_ids, file.path(rds_dir, "sample_ids_for_dendrogram.rds"))
saveRDS(all_cluster_assignments, file.path(rds_dir, "all_cluster_assignments_by_h.rds"))
saveRDS(cluster_final, file.path(rds_dir, paste0("cluster_final_h", gsub("\\.", "_", as.character(selected_h)), ".rds")))

# ------------------------------------------------------------------------------
# 13. Resumen final
# ------------------------------------------------------------------------------
sink(file.path(out_dir, "00_RESUMEN_LEEME.txt"))
cat("====================================\n")
cat("RESULTADOS FAMD + SUBMUESTRA JERÁRQUICA + EXPANSIÓN (regla tamaño mínimo)\n")
cat("====================================\n\n")

cat("Base original usada:\n")
cat("-", file_base, "\n\n")

cat("Base activa para clustering v3:\n")
cat("- Filas:", nrow(base_active), "\n")
cat("- Columnas:", ncol(base_active), "\n")
cat("- Numéricas:", length(num_cols), "\n")
cat("- Categóricas:", length(fac_cols), "\n\n")

cat("Mejoras metodológicas:\n")
cat("- Severity se mantiene completa\n")
cat("- Duration_Hours se transforma con log1p\n")
cat("- Start_Lat / Start_Lng se sacan del clustering activo\n")
cat("- Dendrograma solo sobre submuestra\n")
cat("- Corte por altura, no por k\n")
cat("- Regla de tamaño mínimo aplicada\n\n")

cat("FAMD:\n")
cat("- Componentes:", ncp_famd, "\n\n")

cat("Submuestra para dendrograma:\n")
cat("- Tamaño:", length(sample_ids), "\n\n")

cat("Criterio del dendrograma:\n")
cat("- Distancia: Euclídea sobre componentes FAMD\n")
cat("- Enlace: Ward.D2\n")
cat("- Evaluación por alturas:\n")
print(h_grid)
cat("\n")

cat("Regla de tamaño mínimo:\n")
cat("- min_cluster_prop:", min_cluster_prop, "\n")
cat("- min_cluster_size efectivo:", min_cluster_size, "\n\n")

cat("Altura provisional seleccionada:\n")
cat(selected_h, "\n")
cat("Número de clusters resultante:\n")
cat(length(unique(cluster_final)), "\n\n")

cat("Tamaños de cluster para la altura seleccionada:\n")
print(table(cluster_final))
cat("\n")

cat("Archivos generados:\n")
cat("- plots/: screeplot, variables, individuos, dendrogramas, planos FAMD, mapas, tamaños\n")
cat("- tables/: tamaños por altura, perfiles numéricos y categóricos por altura, resumen de regla mínima\n")
cat("- rds_objects/: objetos FAMD, dendrograma submuestra y clusters finales\n")
cat("- bases clusterizadas por altura y base final seleccionada\n")
sink()

cat("\nScript ejecutado correctamente.\n")
cat("Resultados guardados en:", out_dir, "\n")
