# ==============================================================================
# PMAAD - Profiling final adaptado al clustering K=5
# Entrada: salida final del script CLUSTERING/CLUSTERING/final/clustering_final.R
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# 0. Configuracion
# ------------------------------------------------------------------------------
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_file <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_file)) {
      return(dirname(normalizePath(active_file, winslash = "/")))
    }
  }

  normalizePath(getwd(), winslash = "/")
}

project_dir <- normalizePath(getwd(), winslash = "/")
input_file <- file.path(project_dir, "Datasets", "10_base_final_para_profiling_k5.csv")
output_dir <- file.path(project_dir, "salidas")

max_levels_plot <- 25
max_levels_catdes <- 30
selected_k <- 5

# ------------------------------------------------------------------------------
# 1. Librerias
# ------------------------------------------------------------------------------
packages <- c("fmsb", "scales", "ggpubr", "psych", "dplyr", "FactoMineR", "ggplot2")
lib_dir <- file.path(project_dir, "r_libs")
if (!dir.exists(lib_dir)) dir.create(lib_dir, recursive = TRUE)
.libPaths(c(normalizePath(lib_dir, winslash = "/", mustWork = TRUE), .libPaths()))

new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  install.packages(
    new_packages,
    lib = lib_dir,
    repos = "https://cloud.r-project.org",
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}
invisible(lapply(packages, require, character.only = TRUE))

cat("\n==============================\n")
cat("INICIANDO PROFILING FINAL K=5\n")
cat("==============================\n")

# ------------------------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------------------------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
text_dir <- file.path(output_dir, "txt")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)
if (!dir.exists(text_dir)) dir.create(text_dir, recursive = TRUE)

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  substr(x, 1, 120)
}

collapse_levels_by_frequency <- function(x, max_levels = 25, other_label = "Other") {
  x_chr <- as.character(x)
  counts <- sort(table(x_chr, useNA = "no"), decreasing = TRUE)
  if (length(counts) <= max_levels) return(factor(x_chr))

  keep <- names(counts)[seq_len(max_levels)]
  x_chr[!(x_chr %in% keep) & !is.na(x_chr)] <- other_label
  factor(x_chr, levels = c(keep, other_label))
}

safe_chisq_test <- function(tab) {
  first_test <- suppressWarnings(chisq.test(tab, correct = FALSE))
  if (any(first_test$expected < 5)) {
    return(suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 2000)))
  }
  first_test
}

safe_plot_save <- function(plot_obj, filename, width = 8, height = 5) {
  tryCatch(
    ggsave(filename = filename, plot = plot_obj, bg = "white", width = width, height = height),
    error = function(e) cat("No se pudo guardar", filename, ":", e$message, "\n")
  )
}

semantic_group <- function(variable_name) {
  dplyr::case_when(
    grepl("Humidity|Pressure|Wind_Speed|Precipitation|Weather_Group|Wind_Group|Temperature|Wind_Chill",
          variable_name) ~ "1. Meteorologia",
    grepl("Start_Lat|Start_Lng|County|State|Has_Infrastructure",
          variable_name) ~ "2. Geografia y Via",
    grepl("Sunrise_Sunset|Duration|Start_Hour|Start_Month|Is_Weekend|Time_Period|Season",
          variable_name) ~ "3. Temporalidad",
    grepl("Severity", variable_name) ~ "4. Gravedad e Impacto",
    TRUE ~ "5. Otras"
  )
}

shorten_label <- function(x, max_chars = 32) {
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

vtest_to_light <- function(v_test) {
  dplyr::case_when(
    v_test <= -4 ~ "Muy bajo",
    v_test < -1.96 ~ "Bajo",
    v_test <= 1.96 ~ "Neutro",
    v_test < 4 ~ "Alto",
    TRUE ~ "Muy alto"
  )
}

tlp_palette <- c(
  "Muy bajo" = "#f7191c",
  "Bajo" = "#f28c00",
  "Neutro" = "#fff200",
  "Alto" = "#16e632",
  "Muy alto" = "#08751f"
)

build_tlp_class_data <- function(vtest_source, cluster_sizes_source) {
  cluster_info <- cluster_sizes_source |>
    dplyr::mutate(
      Cluster = paste("Cluster", cluster),
      Cluster_Label = paste0("C", cluster, "\n(n=", n, ")")
    )

  variable_info <- vtest_source |>
    dplyr::group_by(Variable) |>
    dplyr::summarise(
      Tipo = dplyr::first(Tipo),
      Semantica = dplyr::first(Semantica),
      max_abs_vtest = max(abs(v_test), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(Semantica, dplyr::desc(max_abs_vtest), Variable)

  complete_grid <- expand.grid(
    Cluster = cluster_info$Cluster,
    Variable = variable_info$Variable,
    stringsAsFactors = FALSE
  )

  complete_grid |>
    dplyr::left_join(vtest_source, by = c("Cluster", "Variable")) |>
    dplyr::left_join(variable_info, by = "Variable", suffix = c("", "_info")) |>
    dplyr::left_join(cluster_info[, c("Cluster", "Cluster_Label", "n")], by = "Cluster") |>
    dplyr::mutate(
      Tipo = ifelse(is.na(Tipo), Tipo_info, Tipo),
      Semantica = ifelse(is.na(Semantica), Semantica_info, Semantica),
      v_test = ifelse(is.na(v_test), 0, v_test),
      Nivel = factor(
        vtest_to_light(v_test),
        levels = c("Muy bajo", "Bajo", "Neutro", "Alto", "Muy alto")
      ),
      Variable = factor(Variable, levels = variable_info$Variable),
      Cluster_Label = factor(Cluster_Label, levels = rev(cluster_info$Cluster_Label))
    ) |>
    dplyr::select(Cluster, Cluster_Label, n, Variable, Tipo, Semantica, v_test, Nivel, max_abs_vtest)
}

make_tlp_class_plot <- function(plot_data, title, subtitle, grouped = FALSE, axis_text_size = 7) {
  gr <- ggplot(plot_data, aes(x = Variable, y = Cluster_Label, fill = Nivel)) +
    geom_tile(color = "black", linewidth = 0.35) +
    scale_fill_manual(
      values = tlp_palette,
      drop = FALSE,
      name = "v.test",
      labels = c(
        "Muy bajo" = "< -4",
        "Bajo" = "-4 a -1.96",
        "Neutro" = "-1.96 a 1.96",
        "Alto" = "1.96 a 4",
        "Muy alto" = "> 4"
      )
    ) +
    scale_x_discrete(labels = function(x) shorten_label(x, 32)) +
    labs(title = title, subtitle = subtitle, x = "Variables", y = "Class / nc") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      axis.text.x = element_text(size = axis_text_size, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 9, face = "bold"),
      axis.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.6),
      strip.text = element_text(size = 9, face = "bold")
    )

  if (grouped) {
    gr <- gr + facet_grid(. ~ Semantica, scales = "free_x", space = "free_x")
  }

  gr
}

plot_catdes_safely <- function(res_catdes, filename, show_value, barplot_value, width = 1200, height = 900) {
  png(filename = filename, width = width, height = height, res = 120)
  tryCatch(
    plot(
      res_catdes,
      show = show_value,
      col.upper = "red",
      col.lower = "blue",
      barplot = barplot_value,
      cex.names = 0.8
    ),
    error = function(e) {
      plot.new()
      title(main = paste("No se pudo generar catdes", show_value))
      text(0.5, 0.5, e$message)
    }
  )
  dev.off()
}

# ------------------------------------------------------------------------------
# 3. Preparacion de la base
# ------------------------------------------------------------------------------
if (!file.exists(input_file)) stop(paste("No se encuentra la base final:", input_file))

base_cluster_final <- read.csv(input_file, stringsAsFactors = FALSE)

cat("\nBase cargada:", input_file, "\n")
cat("Dimension:", nrow(base_cluster_final), "x", ncol(base_cluster_final), "\n")

if (!"cluster_final" %in% names(base_cluster_final)) {
  stop("La base no contiene la variable cluster_final.")
}

d <- base_cluster_final
names(d)[names(d) == "cluster_final"] <- "cluster"

# Variables que deben tratarse como cualitativas aunque entren como enteros/character.
forced_cat <- intersect(
  c(
    "cluster", "Severity", "State", "County", "Weather_Group", "Wind_Group",
    "Sunrise_Sunset", "Is_Weekend", "Has_Infrastructure", "Time_Period", "Season"
  ),
  names(d)
)

for (nm in forced_cat) {
  d[[nm]] <- as.factor(d[[nm]])
}

char_cols <- names(d)[sapply(d, is.character)]
for (nm in char_cols) {
  d[[nm]] <- as.factor(d[[nm]])
}

d$cluster <- as.factor(d$cluster)
d <- dplyr::relocate(d, cluster, .after = dplyr::last_col())

tipos <- sapply(d, class)
var_Num <- names(tipos)[tipos %in% c("integer", "numeric")]
var_Cat <- names(tipos)[tipos %in% c("factor", "character", "logical")]
var_Cat <- setdiff(var_Cat, "cluster")

for (nm in var_Cat) {
  d[[nm]] <- as.factor(d[[nm]])
}

cluster_sizes <- as.data.frame(table(d$cluster))
names(cluster_sizes) <- c("cluster", "n")
write.csv(cluster_sizes, file.path(table_dir, "00_cluster_sizes_k5.csv"), row.names = FALSE)

cat("\nClusters finales:\n")
print(cluster_sizes)
cat("Numericas:", length(var_Num), "\n")
cat("Categoricas:", length(var_Cat), "\n")

# ------------------------------------------------------------------------------
# 4. Stage 1: significancia de variables
# ------------------------------------------------------------------------------
columns_validate <- setdiff(colnames(d), "cluster")
significant_vars <- c()
stage1_results <- data.frame(
  variable = character(),
  type = character(),
  test = character(),
  p_value = numeric(),
  significant = logical(),
  stringsAsFactors = FALSE
)

sink(file = file.path(text_dir, "01_Test_stage1.txt"))
cat("STAGE 1 - Tests de significancia por cluster\n")
cat("Base:", input_file, "\n\n")

for (cV in columns_validate) {
  current_var <- d[[cV]]
  valid <- !is.na(current_var) & !is.na(d$cluster)

  if (is.numeric(current_var) || is.integer(current_var)) {
    values <- current_var[valid]
    groups <- d$cluster[valid]

    cat("============ ", cV, " ================\n")

    if (length(unique(values)) < 2 || length(unique(groups)) < 2) {
      cat("Variable omitida: variacion insuficiente.\n\n")
      next
    }

    test_name <- "Kruskal-Wallis"
    p_value <- NA_real_

    shapiro_values <- values
    if (length(shapiro_values) > 5000) {
      set.seed(123)
      shapiro_values <- sample(shapiro_values, 5000)
    }

    if (length(shapiro_values) >= 3 && length(unique(shapiro_values)) >= 2) {
      testSH <- shapiro.test(shapiro_values)
      cat("Shapiro-Wilk p-value:", testSH$p.value, "\n")

      if (testSH$p.value > 0.05) {
        anova <- aov(values ~ groups)
        print(summary(anova))
        p_value <- summary(anova)[[1]][["Pr(>F)"]][1]
        test_name <- "ANOVA"
      } else {
        test <- kruskal.test(values ~ groups)
        print(test)
        p_value <- test$p.value
      }
    } else {
      test <- kruskal.test(values ~ groups)
      print(test)
      p_value <- test$p.value
    }

    is_sig <- !is.na(p_value) && p_value <= 0.05
    if (is_sig) significant_vars <- c(significant_vars, cV)

    stage1_results <- rbind(
      stage1_results,
      data.frame(variable = cV, type = "Numerica", test = test_name,
                 p_value = p_value, significant = is_sig, stringsAsFactors = FALSE)
    )

    gr_boxplot <- ggpubr::ggboxplot(d, x = "cluster", y = cV, fill = "cluster", outlier.shape = NA) +
      ggplot2::labs(title = cV, x = "Cluster", y = cV) +
      ggplot2::theme_bw()

    gr_hist <- ggpubr::gghistogram(
      d, x = cV, add = "mean", rug = FALSE, color = "cluster", fill = "cluster", bins = 30
    ) +
      ggplot2::theme_bw()

    gr <- ggpubr::ggarrange(gr_boxplot, gr_hist, heights = c(1, 1), ncol = 2, nrow = 1, align = "v")
    safe_plot_save(
      gr,
      file.path(plot_dir, paste0("stage1_num_", sanitize_filename(cV), ".png")),
      width = 10,
      height = 4.8
    )

    cat("\n")
  }

  if (is.factor(current_var) || is.character(current_var) || is.logical(current_var)) {
    values <- as.factor(current_var[valid])
    groups <- d$cluster[valid]

    cat("============ ", cV, " ================\n")

    if (length(levels(droplevels(values))) < 2 || length(unique(groups)) < 2) {
      cat("Variable omitida: niveles insuficientes.\n\n")
      next
    }

    tab <- table(values, groups)
    test <- safe_chisq_test(tab)
    print(test)

    p_value <- test$p.value
    is_sig <- !is.na(p_value) && p_value <= 0.05
    if (is_sig) significant_vars <- c(significant_vars, cV)

    stage1_results <- rbind(
      stage1_results,
      data.frame(variable = cV, type = "Categorica", test = "Chi-square",
                 p_value = p_value, significant = is_sig, stringsAsFactors = FALSE)
    )

    plot_values <- collapse_levels_by_frequency(current_var, max_levels = max_levels_plot)
    tabla <- data.frame(table(Variable = plot_values, cluster = d$cluster))

    gr <- ggplot(tabla, aes(x = Variable, y = Freq, fill = cluster)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = cV, x = "", y = "Frecuencia") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    safe_plot_save(
      gr,
      file.path(plot_dir, paste0("stage1_cat_", sanitize_filename(cV), ".png")),
      width = 10,
      height = 5
    )

    cat("\n")
  }
}

cat("\nVariables significativas:\n")
print(unique(significant_vars))
sink()

significant_vars <- unique(significant_vars)
write.csv(stage1_results, file.path(table_dir, "01_stage1_results.csv"), row.names = FALSE)
write.csv(data.frame(variable = significant_vars), file.path(table_dir, "01_significant_vars.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Centroides y modas
# ------------------------------------------------------------------------------
sink(file = file.path(text_dir, "02_Centroides_y_modas.txt"))
cat("CENTROIDES Y MODAS - Profiling K=5\n\n")

sel_num <- significant_vars[significant_vars %in% var_Num]
if (length(sel_num) > 0) {
  cat("============ Numericas significativas ================\n")
  print(psych::describeBy(d[, sel_num, drop = FALSE], d$cluster))

  numeric_summary <- d |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(sel_num),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          median = ~ median(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )

  write.csv(numeric_summary, file.path(table_dir, "02_numeric_centroids_significant.csv"), row.names = FALSE)
} else {
  cat("No hay variables numericas significativas.\n")
}

sel_cat <- significant_vars[significant_vars %in% var_Cat]
mode_rows <- data.frame(variable = character(), cluster = character(), level = character(), freq = integer())

if (length(sel_cat) > 0) {
  cat("\n============ Categoricas significativas ================\n")

  for (vC in sel_cat) {
    tabla <- data.frame(table(level = d[[vC]], cluster = d$cluster))
    cal_moda <- tabla |>
      dplyr::group_by(cluster) |>
      dplyr::filter(Freq == max(Freq)) |>
      dplyr::ungroup()

    cat("\n============ ", vC, " ================\n")
    print(cal_moda)

    mode_rows <- rbind(
      mode_rows,
      data.frame(
        variable = vC,
        cluster = as.character(cal_moda$cluster),
        level = as.character(cal_moda$level),
        freq = cal_moda$Freq,
        stringsAsFactors = FALSE
      )
    )
  }

  write.csv(mode_rows, file.path(table_dir, "02_categorical_modes_significant.csv"), row.names = FALSE)
} else {
  cat("\nNo hay variables categoricas significativas.\n")
}
sink()

# ------------------------------------------------------------------------------
# 6. Stage 2: significancia de modalidades con catdes
# ------------------------------------------------------------------------------
cat("\nEjecutando catdes (Stage 2)...\n")

d_catdes <- d
for (nm in var_Cat) {
  if (nlevels(d_catdes[[nm]]) > max_levels_catdes) {
    d_catdes[[nm]] <- collapse_levels_by_frequency(d_catdes[[nm]], max_levels = max_levels_catdes)
  }
}
d_catdes$cluster <- as.factor(d_catdes$cluster)
d_catdes <- dplyr::relocate(d_catdes, cluster, .after = dplyr::last_col())

sink(file = file.path(text_dir, "03_catdes_stage2.txt"))
cat("STAGE 2 - catdes\n")
cat("Nota: variables categoricas con mas de", max_levels_catdes,
    "niveles se agrupan en top niveles + Other para mejorar legibilidad.\n\n")
res_catdes <- FactoMineR::catdes(d_catdes, num.var = ncol(d_catdes), proba = 0.05)
print(res_catdes)
sink()

plot_catdes_safely(res_catdes, file.path(plot_dir, "03_catdes_quanti_barplot.png"), "quanti", TRUE, 1200, 900)
plot_catdes_safely(res_catdes, file.path(plot_dir, "03_catdes_quali.png"), "quali", FALSE, 1200, 1200)
plot_catdes_safely(res_catdes, file.path(plot_dir, "03_catdes_all.png"), "all", FALSE, 1200, 1400)

# ------------------------------------------------------------------------------
# 7. Stage 3: radar chart de centroides numericos
# ------------------------------------------------------------------------------
cat("\nGenerando radar chart...\n")

radar_vars <- if (length(sel_num) > 0) sel_num else var_Num

if (length(radar_vars) >= 2) {
  centroids_num <- aggregate(d[, radar_vars, drop = FALSE], list(cluster = d$cluster), mean, na.rm = TRUE)
  rownames(centroids_num) <- paste("Cluster", centroids_num$cluster)
  centroids_num <- centroids_num[, setdiff(names(centroids_num), "cluster"), drop = FALSE]

  col_max <- apply(d[, radar_vars, drop = FALSE], 2, max, na.rm = TRUE)
  col_min <- apply(d[, radar_vars, drop = FALSE], 2, min, na.rm = TRUE)
  data_radar <- rbind(col_max, col_min, centroids_num)

  colors_border <- scales::hue_pal()(nrow(centroids_num))
  colors_fill <- scales::alpha(colors_border, 0.18)

  png(filename = file.path(plot_dir, "04_Radar_Chart_Centroids.png"),
      width = 1300, height = 1300, res = 110)
  par(mar = c(1, 1, 3, 1))
  fmsb::radarchart(
    data_radar,
    axistype = 1,
    pcol = colors_border,
    pfcol = colors_fill,
    plwd = 3,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "grey30",
    vlcex = 0.75,
    title = "Centroides por cluster - variables numericas"
  )
  legend(
    x = "topright",
    legend = rownames(centroids_num),
    bty = "n",
    pch = 20,
    col = colors_border,
    text.col = "black",
    cex = 1.0,
    pt.cex = 2
  )
  dev.off()
} else {
  cat("Radar no generado: se necesitan al menos dos variables numericas.\n")
}

# ------------------------------------------------------------------------------
# 8. Stage 4: categorical profiling con balloon plot
# ------------------------------------------------------------------------------
cat("\nGenerando balloon plot categorico...\n")

cat_data_long <- data.frame()

if (length(sel_cat) > 0) {
  for (cat_var in sel_cat) {
    x_plot <- collapse_levels_by_frequency(d[[cat_var]], max_levels = max_levels_plot)
    tbl <- table(Cluster = d$cluster, Level = x_plot)
    prop_tbl <- prop.table(tbl, 1) * 100

    df_temp <- as.data.frame(as.table(prop_tbl))
    colnames(df_temp) <- c("Cluster", "Level", "Percentage")
    df_temp$Variable <- cat_var
    df_temp$Label <- paste0(cat_var, ": ", df_temp$Level)

    cat_data_long <- rbind(cat_data_long, df_temp)
  }

  cat_data_long <- cat_data_long |>
    dplyr::filter(Percentage > 5)

  write.csv(cat_data_long, file.path(table_dir, "04_balloon_categorical_data.csv"), row.names = FALSE)

  if (nrow(cat_data_long) > 0) {
    balloon_height <- min(max(900, 35 * length(unique(cat_data_long$Label)) + 350), 5000)

    png(filename = file.path(plot_dir, "05_BalloonPlot_Categorical.png"),
        width = 1300, height = balloon_height, res = 120)

    gr_balloon <- ggplot(cat_data_long, aes(x = Cluster, y = Label)) +
      geom_point(aes(size = Percentage, color = Percentage)) +
      scale_size_continuous(range = c(2, 12)) +
      scale_color_gradient(low = "gray80", high = "#d73027") +
      theme_bw() +
      labs(
        title = "Profiling categorico - Balloon Plot",
        subtitle = "Tamano y color = porcentaje de la categoria dentro del cluster",
        y = "Nivel",
        x = "Cluster",
        size = "%",
        color = "%"
      ) +
      theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 12, face = "bold"),
        panel.grid.major.x = element_blank(),
        strip.text.y = element_text(angle = 0)
      ) +
      facet_grid(Variable ~ ., scales = "free_y", space = "free_y")

    print(gr_balloon)
    dev.off()
  }
} else {
  cat("Balloon plot no generado: no hay categoricas significativas.\n")
}

# ------------------------------------------------------------------------------
# 9. Stage 5: Traffic Light Panel (TLP)
# ------------------------------------------------------------------------------
cat("\nGenerando Traffic Light Panel...\n")

vtest_data <- data.frame()

if (!is.null(res_catdes$quanti)) {
  for (k in names(res_catdes$quanti)) {
    mat <- res_catdes$quanti[[k]]
    if (!is.null(mat) && nrow(mat) > 0 && "v.test" %in% colnames(mat)) {
      temp <- data.frame(
        Variable = rownames(mat),
        Cluster = paste("Cluster", k),
        v_test = mat[, "v.test"],
        Tipo = "Numerica",
        stringsAsFactors = FALSE
      )
      vtest_data <- rbind(vtest_data, temp)
    }
  }
}

if (!is.null(res_catdes$category)) {
  for (k in names(res_catdes$category)) {
    mat <- res_catdes$category[[k]]
    if (!is.null(mat) && nrow(mat) > 0 && "v.test" %in% colnames(mat)) {
      temp <- data.frame(
        Variable = rownames(mat),
        Cluster = paste("Cluster", k),
        v_test = mat[, "v.test"],
        Tipo = "Categorica",
        stringsAsFactors = FALSE
      )
      vtest_data <- rbind(vtest_data, temp)
    }
  }
}

if (nrow(vtest_data) > 0) {
  vtest_data <- vtest_data |>
    dplyr::mutate(
      Significativo = ifelse(abs(v_test) >= 1.96, "Si", "No"),
      v_test_capped = pmax(pmin(v_test, 10), -10),
      Semantica = semantic_group(Variable)
    ) |>
    dplyr::filter(Significativo == "Si")

  write.csv(vtest_data, file.path(table_dir, "05_vtest_data_catdes.csv"), row.names = FALSE)

  if (nrow(vtest_data) > 0) {
    tlp_class_data <- build_tlp_class_data(vtest_data, cluster_sizes)
    write.csv(tlp_class_data, file.path(table_dir, "06_tlp_class_style_data.csv"), row.names = FALSE)

    n_vars_tlp <- length(unique(tlp_class_data$Variable))
    tlp_width <- min(max(1400, 18 * n_vars_tlp + 450), 9000)
    tlp_height <- 850

    png(filename = file.path(plot_dir, "06_TLP_Semaforo.png"),
        width = tlp_width, height = tlp_height, res = 120)

    gr_tlp <- make_tlp_class_plot(
      tlp_class_data,
      title = "Traffic Light Panel (TLP)",
      subtitle = "Formato de clase: v.test discretizado en cinco niveles",
      grouped = FALSE,
      axis_text_size = ifelse(n_vars_tlp > 80, 5, 7)
    )

    print(gr_tlp)
    dev.off()

    atlp_width <- min(max(1500, 20 * n_vars_tlp + 500), 9500)

    png(filename = file.path(plot_dir, "07_aTLP_Semaforo_Avanzado.png"),
        width = atlp_width, height = tlp_height, res = 120)

    gr_atlp <- make_tlp_class_plot(
      tlp_class_data,
      title = "Advanced Traffic Light Panel (aTLP)",
      subtitle = "Variables agrupadas por semantica, como en las diapositivas",
      grouped = TRUE,
      axis_text_size = ifelse(n_vars_tlp > 80, 5, 7)
    )

    print(gr_atlp)
    dev.off()

    recortado_vars <- vtest_data |>
      dplyr::filter(abs(v_test) >= 4 | grepl("Severity", Variable)) |>
      dplyr::group_by(Variable, Semantica) |>
      dplyr::summarise(max_abs_vtest = max(abs(v_test), na.rm = TRUE), .groups = "drop") |>
      dplyr::group_by(Semantica) |>
      dplyr::slice_max(order_by = max_abs_vtest, n = 12, with_ties = FALSE) |>
      dplyr::ungroup()

    severity_vars <- vtest_data |>
      dplyr::filter(grepl("Severity", Variable)) |>
      dplyr::distinct(Variable, Semantica)

    recortado_vars <- dplyr::bind_rows(
      recortado_vars[, c("Variable", "Semantica")],
      severity_vars
    ) |>
      dplyr::distinct(Variable, Semantica)

    vtest_data_recortado <- vtest_data |>
      dplyr::semi_join(recortado_vars, by = c("Variable", "Semantica"))

    if (nrow(vtest_data_recortado) > 0) {
      tlp_class_data_recortado <- build_tlp_class_data(vtest_data_recortado, cluster_sizes)
      write.csv(
        tlp_class_data_recortado,
        file.path(table_dir, "07_tlp_class_style_data_recortado.csv"),
        row.names = FALSE
      )

      n_vars_rec <- length(unique(tlp_class_data_recortado$Variable))
      rec_width <- min(max(1400, 36 * n_vars_rec + 500), 6500)
      rec_height <- 850

      png(filename = file.path(plot_dir, "08_aTLP_Semaforo_Avanzado_Recortado.png"),
          width = rec_width, height = rec_height, res = 120)

      gr_atlp_recortado <- make_tlp_class_plot(
        tlp_class_data_recortado,
        title = "Advanced Traffic Light Panel (aTLP) - Filtrado",
        subtitle = "Variables mas discriminantes (|v.test| >= 4) y Severity",
        grouped = TRUE,
        axis_text_size = ifelse(n_vars_rec > 50, 6, 8)
      )

      print(gr_atlp_recortado)
      dev.off()
    }
  }
} else {
  cat("No se ha podido extraer v.test desde catdes.\n")
}

# ------------------------------------------------------------------------------
# 10. Graficos extra: tamanos y mapa
# ------------------------------------------------------------------------------
p_sizes <- ggplot(cluster_sizes, aes(x = cluster, y = n, fill = cluster)) +
  geom_col() +
  theme_bw() +
  labs(title = paste0("Tamano de clusters | k = ", selected_k), x = "Cluster", y = "N")
safe_plot_save(p_sizes, file.path(plot_dir, "09_cluster_sizes_k5.png"), width = 8, height = 5)

if (all(c("Start_Lat", "Start_Lng") %in% names(d))) {
  p_map <- ggplot(d, aes(x = Start_Lng, y = Start_Lat, color = cluster)) +
    geom_point(alpha = 0.35, size = 0.8) +
    theme_bw() +
    labs(
      title = "Distribucion geografica aproximada por cluster",
      x = "Longitud",
      y = "Latitud",
      color = "Cluster"
    )
  safe_plot_save(p_map, file.path(plot_dir, "10_map_clusters_k5.png"), width = 8, height = 5)
}

# ------------------------------------------------------------------------------
# 11. Resumen final
# ------------------------------------------------------------------------------
sink(file.path(output_dir, "00_RESUMEN_PROFILING_K5.txt"))
cat("====================================\n")
cat("PROFILING FINAL K=5\n")
cat("====================================\n\n")
cat("Base utilizada:\n")
cat("-", input_file, "\n\n")
cat("Filas:", nrow(d), "\n")
cat("Columnas analizadas:", length(columns_validate), "\n")
cat("Variables numericas:", length(var_Num), "\n")
cat("Variables categoricas:", length(var_Cat), "\n\n")
cat("Tamanos de cluster:\n")
print(cluster_sizes)
cat("\nVariables significativas Stage 1:\n")
print(significant_vars)
cat("\nOutputs principales:\n")
cat("- tables/01_stage1_results.csv\n")
cat("- tables/02_numeric_centroids_significant.csv\n")
cat("- tables/02_categorical_modes_significant.csv\n")
cat("- tables/05_vtest_data_catdes.csv\n")
cat("- plots/06_TLP_Semaforo.png\n")
cat("- plots/07_aTLP_Semaforo_Avanzado.png\n")
cat("- plots/08_aTLP_Semaforo_Avanzado_Recortado.png\n")
sink()

cat("\nScript de profiling finalizado correctamente.\n")
cat("Resultados guardados en:", output_dir, "\n")
