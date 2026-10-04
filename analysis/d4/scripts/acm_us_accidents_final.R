args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
bootstrap_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

local_libs <- c(
  file.path(dirname(bootstrap_dir), "r_libs"),
  file.path(dirname(dirname(bootstrap_dir)), "r_libs")
)
local_libs <- local_libs[dir.exists(local_libs)]
if (length(local_libs) > 0) {
  .libPaths(unique(c(normalizePath(local_libs, winslash = "/", mustWork = TRUE), .libPaths())))
}

suppressPackageStartupMessages({
  library(FactoMineR)
  library(factoextra)
  library(ggplot2)
})

find_root <- function(start_dir) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "data", "US_Accidents_Final_FE.csv")) &&
        file.exists(file.path(current, "README.md"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) stop("No se encuentra la raiz del proyecto ACM.", call. = FALSE)
    current <- parent
  }
}

script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

root_dir <- find_root(script_dir)
data_path <- file.path(root_dir, "data", "US_Accidents_Final_FE.csv")
tables_dir <- file.path(root_dir, "outputs", "tables")
figures_dir <- file.path(root_dir, "outputs", "figures")
models_dir <- file.path(root_dir, "outputs", "models")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

write_utf8_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

wrap_label <- function(x, width = 30) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

export_matrix <- function(x, path, id_name) {
  x <- as.data.frame(x, check.names = FALSE)
  ids <- rownames(x)
  if (is.null(ids) || length(ids) != nrow(x)) ids <- seq_len(nrow(x))

  out <- data.frame(
    id = ids,
    x,
    row.names = NULL,
    check.names = FALSE
  )
  names(out)[1] <- id_name
  write_utf8_csv(out, path)
}

export_table_matrix <- function(x, path, id_name = "row_category") {
  out <- data.frame(
    row_category = rownames(x),
    as.data.frame.matrix(x, stringsAsFactors = FALSE),
    row.names = NULL,
    check.names = FALSE
  )
  names(out)[1] <- id_name
  write_utf8_csv(out, path)
}

df <- read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)

audit <- data.frame(
  variable = names(df),
  type = vapply(df, function(x) paste(class(x), collapse = ", "), character(1)),
  n_missing = vapply(df, function(x) sum(is.na(x)), integer(1)),
  n_unique = vapply(df, function(x) length(unique(x)), integer(1)),
  example_values = vapply(df, function(x) paste(head(unique(x), 5), collapse = " | "), character(1)),
  row.names = NULL
)
audit <- audit[order(audit$variable), ]
write_utf8_csv(audit, file.path(tables_dir, "00_audit_variables.csv"))

df_fe <- df
start_time <- as.POSIXct(df_fe$Start_Time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
if (any(is.na(start_time))) stop("Start_Time contiene fechas no convertibles.", call. = FALSE)

month_num <- as.integer(format(start_time, "%m"))
weekday_num <- as.integer(format(start_time, "%u"))

df_fe$Month <- factor(month.abb[month_num], levels = month.abb, ordered = TRUE)
df_fe$Hour <- as.integer(format(start_time, "%H"))
df_fe$Weekday <- factor(
  c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")[weekday_num],
  levels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"),
  ordered = TRUE
)
df_fe$Time_Period <- factor(
  ifelse(df_fe$Hour >= 6 & df_fe$Hour < 12, "Morning",
    ifelse(df_fe$Hour >= 12 & df_fe$Hour < 18, "Afternoon",
      ifelse(df_fe$Hour >= 18 & df_fe$Hour < 24, "Evening", "Night")
    )
  ),
  levels = c("Morning", "Afternoon", "Evening", "Night")
)
df_fe$Season <- factor(
  ifelse(month_num %in% c(12, 1, 2), "Winter",
    ifelse(month_num %in% c(3, 4, 5), "Spring",
      ifelse(month_num %in% c(6, 7, 8), "Summer", "Autumn")
    )
  ),
  levels = c("Winter", "Spring", "Summer", "Autumn")
)
df_fe$Severity_Group <- factor(
  ifelse(df_fe$Severity %in% c(1, 2), "Low_Moderate",
    ifelse(df_fe$Severity == 3, "High", "Very_High")
  ),
  levels = c("Low_Moderate", "High", "Very_High")
)
df_fe$Temp_Bin <- factor(
  ifelse(df_fe$Temperature.C. < 0, "Freezing",
    ifelse(df_fe$Temperature.C. < 10, "Cold",
      ifelse(df_fe$Temperature.C. < 25, "Mild", "Warm")
    )
  ),
  levels = c("Freezing", "Cold", "Mild", "Warm"),
  ordered = TRUE
)
df_fe$WindChill_Bin <- factor(
  ifelse(df_fe$Wind_Chill.C. < 0, "Freezing_Chill",
    ifelse(df_fe$Wind_Chill.C. < 10, "Cold_Chill",
      ifelse(df_fe$Wind_Chill.C. < 25, "Mild_Chill", "Warm_Chill")
    )
  ),
  levels = c("Freezing_Chill", "Cold_Chill", "Mild_Chill", "Warm_Chill"),
  ordered = TRUE
)
df_fe$WindSpeed_Bin <- factor(
  ifelse(df_fe$Wind_Speed.mph. == 0, "Calm",
    ifelse(df_fe$Wind_Speed.mph. < 10, "Low_Wind",
      ifelse(df_fe$Wind_Speed.mph. < 20, "Moderate_Wind", "Strong_Wind")
    )
  ),
  levels = c("Calm", "Low_Wind", "Moderate_Wind", "Strong_Wind"),
  ordered = TRUE
)
df_fe$Precip_Bin <- factor(
  ifelse(df_fe$Precipitation.in. == 0, "No_Precip",
    ifelse(df_fe$Precipitation.in. <= 0.05, "Light_Precip", "Relevant_Precip")
  ),
  levels = c("No_Precip", "Light_Precip", "Relevant_Precip"),
  ordered = TRUE
)
df_fe$Humidity_Bin <- factor(
  ifelse(df_fe$Humidity... < 40, "Low_Humidity",
    ifelse(df_fe$Humidity... < 75, "Medium_Humidity", "High_Humidity")
  ),
  levels = c("Low_Humidity", "Medium_Humidity", "High_Humidity"),
  ordered = TRUE
)

pressure_breaks <- unique(quantile(df_fe$Pressure.in., probs = c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE))
if (length(pressure_breaks) < 4) stop("No hay puntos de corte suficientes para Pressure_Bin.", call. = FALSE)
df_fe$Pressure_Bin <- cut(
  df_fe$Pressure.in.,
  breaks = pressure_breaks,
  labels = c("Low_Pressure", "Medium_Pressure", "High_Pressure"),
  include.lowest = TRUE,
  ordered_result = TRUE
)

derived_vars <- c(
  "Month", "Hour", "Weekday", "Time_Period", "Season", "Severity_Group",
  "Temp_Bin", "WindChill_Bin", "WindSpeed_Bin", "Precip_Bin",
  "Humidity_Bin", "Pressure_Bin"
)

freq_derived <- do.call(rbind, lapply(derived_vars, function(v) {
  tab <- as.data.frame(table(df_fe[[v]], useNA = "ifany"), stringsAsFactors = FALSE)
  names(tab) <- c("modality", "n")
  tab$variable <- v
  tab$prop <- tab$n / sum(tab$n)
  tab[, c("variable", "modality", "n", "prop")]
}))
freq_derived <- freq_derived[order(freq_derived$variable, -freq_derived$n, freq_derived$modality), ]
row.names(freq_derived) <- NULL
write_utf8_csv(freq_derived, file.path(tables_dir, "01_categorical_frequencies.csv"))

plot_balloon <- function(tab, title, path) {
  dat <- as.data.frame(as.table(tab), stringsAsFactors = FALSE)
  names(dat) <- c("row_category", "column_category", "n")
  dat$row_category <- factor(dat$row_category, levels = rev(rownames(tab)))
  dat$column_category <- factor(dat$column_category, levels = colnames(tab))

  p <- ggplot(dat, aes(column_category, row_category)) +
    geom_point(aes(size = n), shape = 21, fill = "#4C78A8", color = "grey20", alpha = 0.75) +
    geom_text(aes(label = n), size = 3) +
    scale_size_area(max_size = 22) +
    labs(title = title, x = NULL, y = NULL, size = "n") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(path, plot = p, width = 9, height = 6, dpi = 150)
}

export_chisq <- function(tab, path) {
  test <- suppressWarnings(chisq.test(tab, correct = FALSE))
  valid <- all(test$expected >= 5)
  lines <- c(
    paste("Chi-square test valid by expected frequencies >= 5:", valid),
    "",
    capture.output(test)
  )
  writeLines(lines, path, useBytes = TRUE)
}

run_cross_eda <- function(row_var, col_var, stub) {
  tab <- table(factor(df_fe[[row_var]]), factor(df_fe[[col_var]]))
  names(dimnames(tab)) <- c(row_var, col_var)

  export_table_matrix(tab, file.path(tables_dir, paste0("02_contingency_", stub, ".csv")), row_var)
  export_table_matrix(round(prop.table(tab, 1), 6), file.path(tables_dir, paste0("02_rowprop_", stub, ".csv")), row_var)
  export_table_matrix(round(prop.table(tab, 2), 6), file.path(tables_dir, paste0("02_colprop_", stub, ".csv")), row_var)
  export_chisq(tab, file.path(tables_dir, paste0("02_chisq_", stub, ".txt")))
  plot_balloon(tab, paste(row_var, "x", col_var), file.path(figures_dir, paste0("02_balloonplot_", stub, ".png")))
}

eda_crosses <- data.frame(
  row_var = c("Severity_Group", "Severity_Group", "Severity_Group", "Severity_Group", "Weather_Group"),
  col_var = c("Weather_Group", "Sunrise_Sunset", "Has_Infrastructure", "Season", "Wind_Group"),
  stub = c(
    "severity_group_by_weather_group",
    "severity_group_by_sunrise_sunset",
    "severity_group_by_has_infrastructure",
    "severity_group_by_season",
    "weather_group_by_wind_group"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(eda_crosses))) {
  run_cross_eda(eda_crosses$row_var[i], eda_crosses$col_var[i], eda_crosses$stub[i])
}

write_ca_block <- function(obj, prefix, stub, row_var, col_var) {
  eig <- data.frame(
    dimension = seq_len(nrow(obj$eig)),
    eigenvalue = obj$eig[, 1],
    percentage_inertia = obj$eig[, 2],
    cumulative_percentage = obj$eig[, 3],
    row.names = NULL
  )
  write_utf8_csv(eig, file.path(tables_dir, paste0(prefix, "_eigenvalues_", stub, ".csv")))
  export_matrix(obj$row$contrib, file.path(tables_dir, paste0(prefix, "_row_contrib_", stub, ".csv")), row_var)
  export_matrix(obj$col$contrib, file.path(tables_dir, paste0(prefix, "_col_contrib_", stub, ".csv")), col_var)
  export_matrix(obj$row$cos2, file.path(tables_dir, paste0(prefix, "_row_cos2_", stub, ".csv")), row_var)
  export_matrix(obj$col$cos2, file.path(tables_dir, paste0(prefix, "_col_cos2_", stub, ".csv")), col_var)
}

save_ca_scree <- function(obj, stub) {
  eig <- data.frame(dimension = seq_len(nrow(obj$eig)), inertia = obj$eig[, 2])
  p <- ggplot(eig, aes(factor(dimension), inertia)) +
    geom_col(fill = "#4C78A8", width = 0.7) +
    geom_text(aes(label = round(inertia, 2)), vjust = -0.35, size = 3.2) +
    labs(title = paste("ACS:", stub), x = "Dimension", y = "Inertia (%)") +
    theme_minimal(base_size = 12)
  ggsave(file.path(figures_dir, paste0("03_ca_screeplot_", stub, ".png")), p, width = 8, height = 5, dpi = 150)
}

save_ca_biplot <- function(obj, stub) {
  row_coord <- as.matrix(obj$row$coord)
  col_coord <- as.matrix(obj$col$coord)

  if (ncol(row_coord) >= 2) {
    p <- fviz_ca_biplot(obj, repel = TRUE, ggtheme = theme_minimal()) +
      labs(title = paste("ACS:", stub))
    ggsave(file.path(figures_dir, paste0("03_ca_biplot_", stub, ".png")), p, width = 9, height = 6, dpi = 150)
  } else {
    dat <- rbind(
      data.frame(category = rownames(row_coord), coord = row_coord[, 1], type = "Row", y = 0.08),
      data.frame(category = rownames(col_coord), coord = col_coord[, 1], type = "Column", y = -0.08)
    )
    p <- ggplot(dat, aes(coord, y, color = type, label = category)) +
      geom_hline(yintercept = 0, color = "grey75") +
      geom_vline(xintercept = 0, color = "grey75", linetype = "dashed") +
      geom_point(size = 3) +
      geom_text(vjust = -0.7, size = 3.4, show.legend = FALSE) +
      scale_y_continuous(NULL, breaks = NULL, limits = c(-0.35, 0.35)) +
      labs(title = paste("ACS:", stub), x = "Dimension 1", color = NULL) +
      theme_minimal(base_size = 12)
    ggsave(file.path(figures_dir, paste0("03_ca_biplot_", stub, ".png")), p, width = 9, height = 5, dpi = 150)
  }
}

run_ca <- function(row_var, col_var, stub) {
  tab <- table(factor(df_fe[[row_var]]), factor(df_fe[[col_var]]))
  obj <- FactoMineR::CA(tab, graph = FALSE)
  write_ca_block(obj, "03_ca", stub, row_var, col_var)
  save_ca_scree(obj, stub)
  save_ca_biplot(obj, stub)
  invisible(obj)
}

ca_crosses <- eda_crosses[1:3, ]
for (i in seq_len(nrow(ca_crosses))) {
  run_ca(ca_crosses$row_var[i], ca_crosses$col_var[i], ca_crosses$stub[i])
}

save_residual_heatmap <- function(stdres, stub) {
  dat <- as.data.frame(as.table(stdres), stringsAsFactors = FALSE)
  names(dat) <- c("row_category", "column_category", "standardized_residual")
  dat$relevant <- abs(dat$standardized_residual) > 2
  dat$label <- ifelse(dat$relevant, paste0(round(dat$standardized_residual, 2), "*"), round(dat$standardized_residual, 2))

  p <- ggplot(dat, aes(column_category, row_category, fill = standardized_residual)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 3.4) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, name = "Std. residual") +
    labs(title = paste("Chi-square standardized residuals:", stub), x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

  ggsave(file.path(figures_dir, paste0("04_chisq_residual_heatmap_", stub, ".png")), p, width = 9, height = 5.5, dpi = 150)
}

export_residuals <- function(row_var, col_var, stub) {
  tab <- table(factor(df_fe[[row_var]]), factor(df_fe[[col_var]]))
  test <- suppressWarnings(chisq.test(tab, correct = FALSE))

  export_table_matrix(tab, file.path(tables_dir, paste0("04_chisq_observed_", stub, ".csv")), row_var)
  export_table_matrix(round(test$expected, 6), file.path(tables_dir, paste0("04_chisq_expected_", stub, ".csv")), row_var)
  export_table_matrix(round(test$stdres, 6), file.path(tables_dir, paste0("04_chisq_standardized_residuals_", stub, ".csv")), row_var)

  relevant <- as.data.frame(as.table(test$stdres), stringsAsFactors = FALSE)
  names(relevant) <- c(row_var, col_var, "standardized_residual")
  relevant$observed <- as.vector(tab)
  relevant$expected <- as.vector(test$expected)
  relevant$abs_standardized_residual <- abs(relevant$standardized_residual)
  relevant <- relevant[relevant$abs_standardized_residual > 2, ]
  relevant <- relevant[order(-relevant$abs_standardized_residual), ]
  write_utf8_csv(relevant, file.path(tables_dir, paste0("04_chisq_relevant_cells_", stub, ".csv")))

  save_residual_heatmap(test$stdres, stub)
}

for (i in seq_len(nrow(ca_crosses))) {
  export_residuals(ca_crosses$row_var[i], ca_crosses$col_var[i], ca_crosses$stub[i])
}

active_vars <- c(
  "Severity_Group", "Sunrise_Sunset", "Time_Period", "Season", "Weather_Group",
  "Wind_Group", "Temp_Bin", "WindSpeed_Bin", "Precip_Bin", "Humidity_Bin",
  "Has_Infrastructure"
)
quali_sup_vars <- c("County", "Airport_Code", "Junction", "Traffic_Signal")
quanti_sup_vars <- c(
  "Start_Lat", "Start_Lng", "Temperature.C.", "Wind_Chill.C.", "Pressure.in.",
  "Wind_Speed.mph.", "Precipitation.in."
)

df_mca <- df_fe[df_fe$Weather_Group != "Other", , drop = FALSE]
filter_summary <- data.frame(
  step = c("before_filter", "removed_weather_other", "after_filter"),
  n = c(nrow(df_fe), sum(df_fe$Weather_Group == "Other"), nrow(df_mca))
)
write_utf8_csv(filter_summary, file.path(tables_dir, "05_mca_weather_other_filter_summary.csv"))

mca_data <- df_mca[, c(active_vars, quali_sup_vars, quanti_sup_vars), drop = FALSE]
mca_data[active_vars] <- lapply(mca_data[active_vars], factor)
mca_data[quali_sup_vars] <- lapply(mca_data[quali_sup_vars], factor)

mca_structure <- data.frame(
  variable = names(mca_data),
  role = ifelse(names(mca_data) %in% active_vars, "active",
    ifelse(names(mca_data) %in% quali_sup_vars, "qualitative_supplementary", "quantitative_supplementary")
  ),
  type = vapply(mca_data, function(x) paste(class(x), collapse = ", "), character(1)),
  n_unique = vapply(mca_data, function(x) length(unique(x)), integer(1)),
  row.names = NULL
)
write_utf8_csv(mca_structure, file.path(tables_dir, "05_mca_data_structure.csv"))

active_freq <- do.call(rbind, lapply(active_vars, function(v) {
  tab <- as.data.frame(table(mca_data[[v]], useNA = "ifany"), stringsAsFactors = FALSE)
  names(tab) <- c("modality", "n")
  tab$variable <- v
  tab$prop <- tab$n / sum(tab$n)
  tab[, c("variable", "modality", "n", "prop")]
}))
active_freq <- active_freq[order(active_freq$variable, -active_freq$n, active_freq$modality), ]
row.names(active_freq) <- NULL
write_utf8_csv(active_freq, file.path(tables_dir, "05_mca_active_modalities_frequencies.csv"))

rare_modalities <- active_freq[active_freq$n < 30 | active_freq$prop < 0.005, ]
write_utf8_csv(rare_modalities, file.path(tables_dir, "05_mca_rare_active_modalities.csv"))

for (v in active_vars) {
  dat <- active_freq[active_freq$variable == v, ]
  dat$modality <- factor(dat$modality, levels = dat$modality[order(dat$n)])
  p <- ggplot(dat, aes(modality, n)) +
    geom_col(fill = "#4C78A8") +
    coord_flip() +
    labs(title = paste("Active MCA modalities:", v), x = NULL, y = "n") +
    theme_minimal(base_size = 12)
  ggsave(file.path(figures_dir, paste0("05_mca_active_freq_", v, ".png")), p, width = 8, height = 5, dpi = 150)
}

quali_sup_idx <- which(names(mca_data) %in% quali_sup_vars)
quanti_sup_idx <- which(names(mca_data) %in% quanti_sup_vars)

res_mca <- FactoMineR::MCA(
  mca_data,
  method = "Indicator",
  quali.sup = quali_sup_idx,
  quanti.sup = quanti_sup_idx,
  graph = FALSE,
  ncp = 10
)
saveRDS(res_mca, file.path(models_dir, "res_mca_main.rds"))

eig_mca <- factoextra::get_eigenvalue(res_mca)
eig_table <- data.frame(dimension = seq_len(nrow(eig_mca)), eig_mca, row.names = NULL, check.names = FALSE)
write_utf8_csv(eig_table, file.path(tables_dir, "06_mca_eigenvalues.csv"))

p_active <- length(active_vars)
retention <- data.frame(
  dimension = eig_table$dimension,
  eigenvalue = eig_table$eigenvalue,
  threshold_1_over_p = 1 / p_active,
  eigen_gt_1_over_p = eig_table$eigenvalue > 1 / p_active,
  p_active_variables = p_active,
  row.names = NULL
)
write_utf8_csv(retention, file.path(tables_dir, "06_mca_dimension_retention_criterion.csv"))

p <- factoextra::fviz_screeplot(res_mca, addlabels = TRUE, ncp = 10) +
  geom_hline(yintercept = 100 / p_active, linetype = "dashed", color = "#D55E00") +
  labs(title = "MCA scree plot", x = "Dimension", y = "Explained inertia (%)") +
  theme_minimal(base_size = 12)
ggsave(file.path(figures_dir, "06_mca_screeplot.png"), p, width = 8, height = 5, dpi = 150)

mca_var <- factoextra::get_mca_var(res_mca)
export_matrix(mca_var$coord, file.path(tables_dir, "07_mca_var_coordinates.csv"), "modality")
export_matrix(mca_var$cos2, file.path(tables_dir, "07_mca_var_cos2.csv"), "modality")
export_matrix(mca_var$contrib, file.path(tables_dir, "07_mca_var_contributions.csv"), "modality")

save_mca_var_map <- function(axes, path, color_mode = NULL) {
  args <- list(X = res_mca, axes = axes, repel = TRUE, ggtheme = theme_minimal())
  if (!is.null(color_mode)) {
    args$col.var <- color_mode
    args$gradient.cols <- c("#00AFBB", "#E7B800", "#FC4E07")
  }
  p <- do.call(factoextra::fviz_mca_var, args) +
    labs(title = paste0("MCA modalities map Dim ", axes[1], "-", axes[2]))
  ggsave(path, p, width = 10, height = 7, dpi = 150)
}

save_mca_var_map(c(1, 2), file.path(figures_dir, "07_mca_modalities_map_dim1_2.png"))
save_mca_var_map(c(1, 3), file.path(figures_dir, "07_mca_modalities_map_dim1_3.png"))
save_mca_var_map(c(2, 3), file.path(figures_dir, "07_mca_modalities_map_dim2_3.png"))
save_mca_var_map(c(1, 2), file.path(figures_dir, "07_mca_modalities_contrib_dim1_2.png"), "contrib")
save_mca_var_map(c(1, 2), file.path(figures_dir, "07_mca_modalities_cos2_dim1_2.png"), "cos2")

save_contrib_plot <- function(axis) {
  dim_name <- paste0("Dim ", axis)
  dat <- data.frame(modality = rownames(mca_var$contrib), contribution = mca_var$contrib[, dim_name])
  dat <- head(dat[order(-dat$contribution), ], 20)
  dat$label <- factor(wrap_label(dat$modality), levels = rev(wrap_label(dat$modality)))

  p <- ggplot(dat, aes(label, contribution)) +
    geom_col(fill = "#4C78A8", width = 0.72) +
    geom_text(aes(label = round(contribution, 2)), hjust = -0.12, size = 3.2) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = paste("Top 20 modality contributions to Dim", axis), x = NULL, y = "Contribution") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 9), panel.grid.major.y = element_blank())

  ggsave(
    file.path(figures_dir, paste0("07_mca_var_contrib_dim", axis, "_top20_horizontal.png")),
    p,
    width = 9,
    height = 7,
    dpi = 150
  )
}

for (axis in 1:3) save_contrib_plot(axis)

dimdesc_mca <- FactoMineR::dimdesc(res_mca, axes = 1:3, proba = 0.05)
capture.output(dimdesc_mca, file = file.path(tables_dir, "07_mca_dimdesc.txt"))
for (axis in 1:3) {
  capture.output(dimdesc_mca[[axis]], file = file.path(tables_dir, paste0("07_mca_dimdesc_axis_", axis, ".txt")))
}

mca_ind <- factoextra::get_mca_ind(res_mca)
export_matrix(mca_ind$coord, file.path(tables_dir, "08_mca_ind_coordinates.csv"), "mca_row_id")
export_matrix(mca_ind$cos2, file.path(tables_dir, "08_mca_ind_cos2.csv"), "mca_row_id")
export_matrix(mca_ind$contrib, file.path(tables_dir, "08_mca_ind_contributions.csv"), "mca_row_id")

save_ind_contrib <- function() {
  p <- factoextra::fviz_mca_ind(
    res_mca,
    axes = c(1, 2),
    geom = "point",
    col.ind = "contrib",
    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    ggtheme = theme_minimal()
  ) +
    labs(title = "MCA individuals Dim 1-2 colored by contribution")
  ggsave(file.path(figures_dir, "08_mca_ind_dim1_2_contrib.png"), p, width = 9, height = 6, dpi = 150)
}

save_ind_group <- function(axes, group, stub, ellipses = FALSE) {
  p <- factoextra::fviz_mca_ind(
    res_mca,
    axes = axes,
    geom = "point",
    habillage = group,
    addEllipses = ellipses,
    ellipse.type = "confidence",
    palette = "jco",
    ggtheme = theme_minimal()
  ) +
    labs(title = paste0("MCA individuals Dim ", axes[1], "-", axes[2], " by ", stub))
  suffix <- if (ellipses) "_ellipses" else ""
  ggsave(file.path(figures_dir, paste0("08_mca_ind_dim", axes[1], "_", axes[2], "_", tolower(stub), suffix, ".png")),
    p, width = 9, height = 6, dpi = 150
  )
}

save_ind_contrib()
save_ind_group(c(1, 2), mca_data$Severity_Group, "severity_group", TRUE)
save_ind_group(c(1, 2), mca_data$Sunrise_Sunset, "sunrise_sunset", TRUE)
save_ind_group(c(1, 2), mca_data$Weather_Group, "weather_group", TRUE)
save_ind_group(c(1, 3), mca_data$Weather_Group, "weather_group", FALSE)
save_ind_group(c(2, 3), mca_data$Weather_Group, "weather_group", FALSE)

ind_contrib <- data.frame(
  mca_row_id = seq_len(nrow(mca_ind$contrib)),
  original_row_id = as.integer(rownames(df_mca)),
  as.data.frame(mca_ind$contrib, check.names = FALSE),
  row.names = NULL,
  check.names = FALSE
)
ind_coord <- data.frame(
  mca_row_id = seq_len(nrow(mca_ind$coord)),
  as.data.frame(mca_ind$coord, check.names = FALSE),
  row.names = NULL,
  check.names = FALSE
)

top_ind <- ind_contrib
top_ind$total_contrib_dim1_2 <- top_ind[["Dim 1"]] + top_ind[["Dim 2"]]
top_ind <- head(top_ind[order(-top_ind$total_contrib_dim1_2), ], 30)

derived_top <- df_mca[
  match(top_ind$original_row_id, as.integer(rownames(df_mca))),
  c(active_vars, quali_sup_vars),
  drop = FALSE
]
names(derived_top) <- paste0("mca_", names(derived_top))

coord_top <- ind_coord[match(top_ind$mca_row_id, ind_coord$mca_row_id), c("Dim 1", "Dim 2", "Dim 3")]
names(coord_top) <- c("coord_Dim1", "coord_Dim2", "coord_Dim3")

top_output <- cbind(
  top_ind[, c("mca_row_id", "original_row_id", "Dim 1", "Dim 2", "Dim 3", "total_contrib_dim1_2")],
  coord_top,
  derived_top,
  df[top_ind$original_row_id, , drop = FALSE]
)
names(top_output)[3:5] <- c("contrib_Dim1", "contrib_Dim2", "contrib_Dim3")
write_utf8_csv(top_output, file.path(tables_dir, "08_mca_top30_individuals_dim1_2_original_rows.csv"))

save_biplot <- function(axes, path, group = NULL) {
  args <- list(
    X = res_mca,
    axes = axes,
    label = "var",
    geom.ind = "point",
    geom.var = c("point", "text"),
    repel = TRUE,
    select.var = list(contrib = 25),
    ggtheme = theme_minimal()
  )
  if (!is.null(group)) args$habillage <- group

  p <- do.call(factoextra::fviz_mca_biplot, args) +
    labs(title = paste0("MCA biplot Dim ", axes[1], "-", axes[2]))
  ggsave(path, p, width = 11, height = 8, dpi = 150)
}

save_biplot(c(1, 2), file.path(figures_dir, "09_mca_biplot_dim1_2_top25_modalities.png"))
save_biplot(c(1, 2), file.path(figures_dir, "09_mca_biplot_dim1_2_by_sunrise_sunset.png"), mca_data$Sunrise_Sunset)
save_biplot(c(1, 2), file.path(figures_dir, "09_mca_biplot_dim1_2_by_weather_group.png"), mca_data$Weather_Group)
save_biplot(c(1, 3), file.path(figures_dir, "09_mca_biplot_dim1_3_top25_modalities.png"))
