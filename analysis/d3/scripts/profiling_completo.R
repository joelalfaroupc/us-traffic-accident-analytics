# ==============================================================================
# PROFILING
# ==============================================================================
paquetes <- c("fmsb", "scales", "ggpubr", "psych", "dplyr", "FactoMineR", "ggplot2")
nuevos_paquetes <- paquetes[!(paquetes %in% installed.packages()[, "Package"])]
if (length(nuevos_paquetes) > 0) install.packages(nuevos_paquetes, dependencies = TRUE)
invisible(lapply(paquetes, require, character.only = TRUE))

cat("\n==============================\n")
cat("INICIANDO PROFILING\n")
cat("==============================\n")

# ------------------------------------------------------------------------------
# PREPARACIÓN DE LA BASE
# ------------------------------------------------------------------------------

base_cluster_clean <- read.csv("clustering_results/06_base_limpia_con_clusters_robusta.csv")
# Limpiamos variables que rompen el análisis
d <- base_cluster_clean %>%
  select(-cluster_kproto, -cluster_cure, -ID_ROW)

# Renombrar para que coincida con el script del profesor
colnames(d)[which(names(d) == "cluster_final")] <- "cluster"
d$cluster <- as.factor(d$cluster)

# CRÍTICO: Movemos 'cluster' al final para que ncol(d) funcione en el catdes
d <- d %>% relocate(cluster, .after = last_col())

# Tipos de variables (Añadiendo logical para que no se pierdan)
tipos <- sapply(d, class)
var_Num <- names(tipos)[which(tipos %in% c("integer", "numeric"))]
var_Cat <- names(tipos)[which(tipos %in% c("factor", "character", "logical"))]

for(vC in var_Cat){ d[, vC] <- as.factor(d[, vC]) }

# ------------------------------------------------------------------------------
# STAGE 1: VARIABLE SIGNIFICANCE
# ------------------------------------------------------------------------------
pathProfiling <- "Profiling/"
if (!dir.exists(pathProfiling)) dir.create(pathProfiling)

columns_validate <- colnames(d)[colnames(d) != "cluster"]
significant_vars <- c()

sink(file = paste0(pathProfiling, "Test_stage1.txt"))
for (cV in columns_validate) {
  current_var <- d[[cV]]
  
  # --- Numeric Variables ---
  if (any(class(current_var) %in% c("integer", "numeric"))) {
    
    # Parche obligatorio para Shapiro con N > 5000
    if(length(current_var) > 5000) {
      testSH <- shapiro.test(sample(current_var, 5000))
    } else {
      testSH <- shapiro.test(current_var)
    }
    
    if (testSH$p.value > 0.05) {
      # Normal --> ANOVA
      anova <- aov(current_var ~ d$cluster)
      cat("============ ", cV, " ================\n")
      print(summary(anova)); cat("\n")
      p_valor <- summary(anova)[[1]][["Pr(>F)"]][1]
      if (!is.na(p_valor) && p_valor <= 0.05) {significant_vars <- c(significant_vars, cV)}
    } else {
      # Not Normal --> Kruskal-Wallis
      test <- kruskal.test(current_var ~ d$cluster)
      cat("============ ", cV, " ================\n")
      print(test); cat("\n")
      if (!is.na(test$p.value) && test$p.value <= 0.05) {significant_vars <- c(significant_vars, cV)}
    }
    
    # Graphical representation
    gr_Boxplot <- ggboxplot(d, "cluster", cV, fill = "cluster")
    gr_Hist    <- gghistogram(d, x  = cV, add = "mean", rug = TRUE, color = "cluster", fill = "cluster")
    gr         <- ggarrange(gr_Boxplot, gr_Hist, heights = c(2, 0.7), ncol = 2, nrow = 1, align = "v")
    
    ggsave(filename = paste0(pathProfiling, "clustering_var_", cV, ".png"),
           plot = gr, bg = "white", width = 8, height = 4)
  }
  
  # --- Categorical Variables ---
  if (any(class(current_var) %in% c("factor", "character", "logical"))) {
    
    # Chi2 test
    test <- chisq.test(current_var, d$cluster)
    cat("============ ", cV, " ================\n")
    print(test); cat("\n")
    
    tabla <- data.frame(table(Var1 = current_var, cluster = d$cluster))
    
    if (!is.na(test$p.value) && test$p.value <= 0.05) {significant_vars <- c(significant_vars, cV)}
    
    # Graphical representation
    gr <- ggplot(tabla, aes(x = Var1, y = Freq, fill = cluster)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = "", x = "", y = "") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(filename = paste0(pathProfiling, "clustering_var_", cV, ".png"),
           plot = gr, bg = "white", width = 8, height = 4)
  }
}
sink()

# ------------------------------------------------------------------------------
# CENTROIDS & MODES
# ------------------------------------------------------------------------------
sink(file = paste0(pathProfiling, "Centroides.txt"))

sel_num <- significant_vars[which(significant_vars %in% var_Num)]
if(length(sel_num) > 0) {
  print(describeBy(d[, sel_num, drop=FALSE], d$cluster))
}

sel_cat <- significant_vars[which(significant_vars %in% var_Cat)]
listaModa <- list()
if(length(sel_cat) > 0) {
  for (vC in sel_cat){
    tabla   <- data.frame(table(Var1 = d[[vC]], Var2 = d$cluster))
    calModa <- tabla %>%
      group_by(Var2) %>%
      filter(Freq == max(Freq)) %>%
      select(Var1, Var2)  %>% as.data.frame()
    colnames(calModa) <- c(vC, "cluster")
    cat("\n============ ", vC, " ================\n")
    print(calModa); cat("\n")
    listaModa[[vC]][["moda"]] <- calModa
  }
}
sink()

# ------------------------------------------------------------------------------
# STAGE 2: SIGNIFICANCE OF MODALITIES (catdes)
# ------------------------------------------------------------------------------
cat("\nEjecutando catdes (Stage 2)...\n")
res_catdes <- catdes(d, num.var = ncol(d))
print(res_catdes)

##-- Barplot (Se mostrarán en el visor de RStudio)
plot(res_catdes, show = "quanti", 
     col.upper = "red", col.lower = "blue", 
     barplot = TRUE, cex.names = 1)

par(mfrow = c(1, 1))
plot(res_catdes, show = "quali", 
     col.upper = "red", col.lower = "blue", 
     barplot = FALSE, cex.names = 1)

plot(res_catdes, show = "all", 
     col.upper = "red", col.lower = "blue", 
     barplot = FALSE, cex.names = 1)

cat("\nScript de Profiling finalizado. Gráficos y TXTs generados en la carpeta 'Profiling/'.\n")

# -----------------------------------------------------
# STAGE 3: RADAR/BALLOON PLOTS
# -----------------------------------------------------

# RADAR CHART (CENTROIDES)

# 1. Calcular centroides numéricos usando var_Num (heredado del script principal)
centroids_num <- aggregate(d[, var_Num], list(Group = d$cluster), mean, na.rm = TRUE)

# Arreglar rownames
rownames(centroids_num) <- paste("Clúster", centroids_num$Group)
centroids_num <- centroids_num[, -1] 

# 2. Preparar formato fmsb (Fila 1 = Max, Fila 2 = Min)
col_max <- apply(d[, var_Num], 2, max, na.rm = TRUE)
col_min <- apply(d[, var_Num], 2, min, na.rm = TRUE)

data_radar <- rbind(col_max, col_min, centroids_num)

# 3. Definir colores (Usaremos la paleta azul y amarilla que fijamos antes)
colors_border <- c("#0073C2FF", "#EFC000FF")[1:nrow(centroids_num)]
colors_fill   <- scales::alpha(colors_border, 0.2)

# 4. Generar Radar Chart
png(filename = paste0(pathProfiling, "Radar_Chart_Centroids.png"), 
    width = 1200, height = 1200, res = 100)

par(mar = c(1, 1, 2, 1))

radarchart(data_radar,
           axistype = 1,
           pcol = colors_border,
           pfcol = colors_fill,
           plwd = 3,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           vlcex = 0.8,
           title = "Centroides (Variables Numéricas)"
)

legend(x = "topright", 
       legend = rownames(centroids_num), 
       bty = "n", pch = 20, col = colors_border, 
       text.col = "black", cex = 1.2, pt.cex = 2)

dev.off()
cat("Radar Chart generado con éxito.\n")

# ------------------------------------------------------------------------------
# STAGE 4: CATEGORICAL PROFILING (BALLOON PLOT)
# ------------------------------------------------------------------------------
cat_data_long <- data.frame()

for (cat_var in sel_cat) {
  # Tabla de contingencia directa
  tbl <- table(Cluster = d$cluster, Level = d[[cat_var]])
  
  # Porcentaje por fila (dentro del clúster)
  prop_tbl <- prop.table(tbl, 1) * 100 
  
  df_temp <- as.data.frame(as.table(prop_tbl))
  colnames(df_temp) <- c("Cluster", "Level", "Percentage")
  
  df_temp$Variable <- cat_var
  df_temp$Label <- paste0(cat_var, ": ", df_temp$Level)
  
  cat_data_long <- rbind(cat_data_long, df_temp)
}

# Filtro crítico: Solo mostrar lo que represente al menos el 5% del clúster
cat_data_long <- cat_data_long %>% filter(Percentage > 5)

# Generar Balloon Plot (Ampliamos la altura por la cantidad de variables)
png(filename = paste0(pathProfiling, "BalloonPlot_Categorical.png"), 
    width = 1200, height = 2000, res = 100) 

gr_balloon <- ggplot(cat_data_long, aes(x = Cluster, y = Label)) +
  geom_point(aes(size = Percentage, color = Percentage)) +
  scale_size_continuous(range = c(2, 12)) +       
  scale_color_gradient(low = "gray80", high = "#E41A1C") + 
  theme_bw() +
  labs(
    title = "Profiling Categórico: Balloon Plot",
    subtitle = "Tamaño y color = % de la categoría dentro del clúster",
    y = "Niveles de Variable",
    x = "Clúster",
    size = "%",
    color = "%"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_blank() 
  ) +
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  theme(strip.text.y = element_text(angle = 0)) 

print(gr_balloon)

dev.off()
cat("Balloon Plot generado con éxito.\n")

# ==============================================================================
# STAGE 5: TRAFFIC LIGHT PANEL (TLP)
# ==============================================================================
if (!require(tidyr)) install.packages("tidyr", dependencies = TRUE)
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Extraer los valores test (v.test) del objeto res_catdes
vtest_data <- data.frame()

# Extraer numéricas (quanti)
if (!is.null(res_catdes$quanti)) {
  for (k in names(res_catdes$quanti)) {
    mat <- res_catdes$quanti[[k]]
    if (nrow(mat) > 0) {
      temp <- data.frame(
        Variable = rownames(mat),
        Cluster = paste("Clúster", k),
        v_test = mat[, "v.test"],
        Tipo = "Numérica"
      )
      vtest_data <- rbind(vtest_data, temp)
    }
  }
}

# Extraer categóricas (category)
if (!is.null(res_catdes$category)) {
  for (k in names(res_catdes$category)) {
    mat <- res_catdes$category[[k]]
    if (nrow(mat) > 0) {
      temp <- data.frame(
        Variable = rownames(mat),
        Cluster = paste("Clúster", k),
        v_test = mat[, "v.test"],
        Tipo = "Categórica"
      )
      vtest_data <- rbind(vtest_data, temp)
    }
  }
}

# 2. Preparar los datos para el semáforo
# Filtramos para no mostrar variables que no son significativas en absoluto (entre -1.96 y 1.96)
# y limitamos/capeamos los valores extremos (ej. a -10 y 10) para que la escala de colores no se rompa
vtest_data <- vtest_data %>%
  mutate(
    Significativo = ifelse(abs(v_test) >= 1.96, "Si", "No"),
    v_test_capped = ifelse(v_test > 10, 10, ifelse(v_test < -10, -10, v_test))
  ) %>%
  filter(Significativo == "Si") # Solo nos quedamos con lo que define a los clústers

# 3. Dibujar el TLP (Heatmap)
png(filename = paste0(pathProfiling, "TLP_Semaforo.png"), 
    width = 1000, height = 1200, res = 120)

gr_tlp <- ggplot(vtest_data, aes(x = Cluster, y = Variable, fill = v_test_capped)) +
  geom_tile(color = "grey80", size = 0.5) +
  # Escala de colores tipo SEMÁFORO: Rojo (bajo), Blanco (neutro), Verde (alto)
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850", midpoint = 0, 
                       limits = c(-10, 10), name = "v.test") +
  theme_minimal() +
  labs(
    title = "Traffic Light Panel (TLP)",
    subtitle = "Verde (>1.96): Sobrerrepresentado | Rojo (<-1.96): Infrarrepresentado",
    x = "Clústeres",
    y = "Variables y Categorías"
  ) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, face = "bold")
  )

print(gr_tlp)
dev.off()

# ==============================================================================
# STAGE 6: ADVANCED TRAFFIC LIGHT PANEL (aTLP) - DEFINITIVO
# ==============================================================================

library(dplyr)
library(ggplot2)

# 1. Definir la semántica de las variables
vtest_data <- vtest_data %>%
  mutate(
    # Quitamos las categorías en las variables cualitativas para leer solo el nombre base
    Variable_Base = sub("=.*", "", Variable), # por si catdes usa "="
    Variable_Base = sub("\\..*", "", Variable_Base), # limpia algunos sufijos de factor
    Variable_Base = trimws(Variable_Base),
    
    # Mapeo usando los nombres reales de la base_limpia_con_clusters_robusta.csv
    Semantica = case_when(
      grepl("Humidity|Pressure|Wind_Speed|Precipitation|Weather_Group|Wind_Group|Temperature|Wind_Chill", Variable) ~ "1. Meteorología",
      
      grepl("Start_Lat|Start_Lng|County|State|Airport_Code|Junction|Traffic_Signal|Has_Infrastructure", Variable) ~ "2. Geografía y Vía",
      
      grepl("Sunrise_Sunset|Duration_Hours|Start_Hour|Start_Month|Is_Weekend", Variable) ~ "3. Temporalidad",
      
      grepl("Severity", Variable) ~ "4. Gravedad e Impacto",
      
      TRUE ~ "5. Otras"
    )
  )

# 2. Dibujar el aTLP (Panel Semáforo Avanzado)
png(filename = paste0(pathProfiling, "aTLP_Semaforo_Avanzado.png"), 
    width = 1100, height = 3000, res = 120)

gr_atlp <- ggplot(vtest_data, aes(x = Cluster, y = Variable, fill = v_test_capped)) +
  geom_tile(color = "grey80", size = 0.5) +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850", midpoint = 0, 
                       limits = c(-10, 10), name = "v.test") +
  theme_minimal() +
  labs(
    title = "Advanced Traffic Light Panel (aTLP)",
    subtitle = "Perfil de los Clústeres de Accidentes (agrupados por Semántica)",
    x = "Clústeres",
    y = ""
  ) +
  theme(
    axis.text.x = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    strip.text.y = element_text(size = 11, face = "bold", angle = 0), # Texto de los bloques semánticos
    strip.background = element_rect(fill = "grey90", color = "grey50") # Fondo de los bloques
  ) +
  # Cortamos el gráfico en paneles según la dimensión semántica
  facet_grid(Semantica ~ ., scales = "free_y", space = "free_y")

print(gr_atlp)
dev.off()
