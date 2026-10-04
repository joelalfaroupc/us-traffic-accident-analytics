# ==============================================================================
# PARTE II: IMPUTACIÓN Y COMPARATIVA DE MÉTODOS
# ==============================================================================

list.of.packages <- c("VIM", "mice", "Hmisc", "mi", "missForest", "dplyr", "tidyr", "ggplot2", "cluster", "skimr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(list.of.packages, require, character.only = TRUE))

dd <- read.csv("Datasets/US_Accidents_Clean.csv", 
               header = TRUE, sep = ",", stringsAsFactors = TRUE, na.strings = c("", "NA", " "))

vars_a_imputar <- c(
  "Airport_Code",
  "Temperature.F.", "Wind_Chill.F.", "Humidity...", "Pressure.in.", 
  "Visibility.mi.", "Wind_Speed.mph.", "Precipitation.in.",
  "Wind_Direction", "Weather_Condition", "Sunrise_Sunset", 
  "Civil_Twilight", "Nautical_Twilight", "Astronomical_Twilight"
)

# ==============================================================================
# FUNCIÓN DE VISUALIZACIÓN 
# ==============================================================================
plot_comparativa <- function(df_imputado, df_original, nombre_metodo) {
  try(dev.off(), silent = TRUE) 
  
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1)) 
  
  plot_densidad <- function(var_name, titulo, color_linea) {
    datos_orig <- as.numeric(na.omit(df_original[[var_name]]))
    datos_imp <- as.numeric(na.omit(df_imputado[[var_name]]))
    
    if(length(datos_orig) < 2 || length(datos_imp) < 2) {
      plot.new()
      title(main = paste("Sin datos para", titulo))
      return()
    }
    
    dens_orig <- density(datos_orig)
    dens_imp <- density(datos_imp)
    
    plot(dens_orig, main = paste(nombre_metodo, "-", titulo), col = "black", lwd = 2, 
         ylim = c(0, max(dens_orig$y, dens_imp$y)), xlab = titulo, ylab = "Densidad")
    lines(dens_imp, col = color_linea, lwd = 2, lty = 2)
    legend("topright", legend=c("Original", "Imputado"), col=c("black", color_linea), 
           lwd=2, lty=c(1,2), cex=0.7)
  }
  
  plot_densidad("Temperature.F.", "Temp (F)", "red")
  plot_densidad("Humidity...", "Humedad (%)", "blue")
  plot_densidad("Wind_Chill.F.", "Wind Chill (F)", "purple")
  plot_densidad("Pressure.in.", "Presión (in)", "darkgreen")
  plot_densidad("Wind_Speed.mph.", "Vel. Viento (mph)", "brown")
  
  if(!is.null(df_original$Sunrise_Sunset) && !is.null(df_imputado$Sunrise_Sunset)) {
    prop_orig <- prop.table(table(na.omit(df_original$Sunrise_Sunset)))
    prop_imp <- prop.table(table(na.omit(df_imputado$Sunrise_Sunset)))
    niveles <- union(names(prop_orig), names(prop_imp))
    
    prop_orig_full <- setNames(rep(0, length(niveles)), niveles)
    prop_imp_full <- setNames(rep(0, length(niveles)), niveles)
    prop_orig_full[names(prop_orig)] <- as.numeric(prop_orig)
    prop_imp_full[names(prop_imp)] <- as.numeric(prop_imp)
    
    mat_sun <- rbind(prop_orig_full, prop_imp_full)
    barplot(mat_sun, beside = TRUE, col = c("gray40", "darkorange"),
            main = paste(nombre_metodo, "- Día/Noche"), ylab = "Proporción",
            ylim = c(0, max(mat_sun) * 1.25), las = 1)
    legend("topright", legend = c("Original", "Imputado"), 
           fill = c("gray40", "darkorange"), cex = 0.7)
  } else {
    plot.new()
    title("Sin datos Sunrise_Sunset")
  }
  
  par(mfrow = c(1, 1))
}

# ------------------------------------------------------------------------------
# MÉTODO 1: Imputación Básica
# ------------------------------------------------------------------------------
cat("\nEjecutando Método 1: Básica...\n")
dd_basic <- dd
dd_basic$Temperature.F. <- with(dd_basic, impute(Temperature.F., mean))
dd_basic$Humidity... <- with(dd_basic, impute(Humidity..., mean))
moda_sun <- names(sort(table(dd_basic$Sunrise_Sunset), decreasing = TRUE))[1]
dd_basic$Sunrise_Sunset[is.na(dd_basic$Sunrise_Sunset)] <- moda_sun

plot_comparativa(dd_basic, dd, "Básica")
dd_basic <- NULL

# ------------------------------------------------------------------------------
# MÉTODO 2: Imputación con aregImpute
# ------------------------------------------------------------------------------
cat("Ejecutando Método 2: aregImpute...\n")
dd_areg <- dd
impute_arg <- aregImpute(~ Temperature.F. + Humidity... + Sunrise_Sunset + 
                           Pressure.in. + Wind_Speed.mph., 
                         data = dd_areg, n.impute = 5)

dd_areg$Temperature.F.[is.na(dd_areg$Temperature.F.)] <- rowMeans(impute_arg$imputed$Temperature.F.)
dd_areg$Humidity...[is.na(dd_areg$Humidity...)] <- rowMeans(impute_arg$imputed$Humidity...)
dd_areg$Sunrise_Sunset[is.na(dd_areg$Sunrise_Sunset)] <- impute_arg$imputed$Sunrise_Sunset[,1]

plot_comparativa(dd_areg, dd, "aregImpute")
dd_areg <- NULL

# ------------------------------------------------------------------------------
# MÉTODO 3: Multiple Iterative Regression Imputation (MI)
# ------------------------------------------------------------------------------
cat("\nEjecutando Método 3: MI...\n")
vars_mi <- c("Temperature.F.", "Humidity...", "Sunrise_Sunset")
dd_mi_subset <- dd[, vars_mi]

mi_data <- suppressWarnings(mi(dd_mi_subset, seed = 335))
mi_data_extracted <- mi::complete(mi_data, m = 1)

dd_mi_imputed <- dd
dd_mi_imputed$Temperature.F. <- mi_data_extracted$Temperature.F.
dd_mi_imputed$Humidity... <- mi_data_extracted$Humidity...
dd_mi_imputed$Sunrise_Sunset <- mi_data_extracted$Sunrise_Sunset

plot_comparativa(dd_mi_imputed, dd, "MI (Multiple Imputation)")

dd_mi_subset <- NULL
mi_data <- NULL
dd_mi_imputed <- NULL

# ------------------------------------------------------------------------------
# MÉTODO 4: Media con Variable Objetivo (Target: Día/Noche)
# ------------------------------------------------------------------------------
cat("Ejecutando Método 4: Target...\n")
dd_target <- dd
target <- "Sunrise_Sunset"

for (varImp in c("Temperature.F.", "Humidity...")) {
  means <- aggregate(dd_target[, varImp, drop=FALSE], list(dd_target[, target]), mean, na.rm = TRUE)
  for (g in na.omit(unique(dd_target[, target]))) {
    cond <- dd_target[, target] == g & !is.na(dd_target[, target])
    na_index <- is.na(dd_target[, varImp]) & cond
    valor_imputar <- means[means[, "Group.1"] == g, varImp]
    dd_target[na_index, varImp] <- valor_imputar
  }
}
plot_comparativa(dd_target, dd, "Target")
dd_target <- NULL

# ------------------------------------------------------------------------------
# MÉTODO 5: missForest (Random Forest)
# ------------------------------------------------------------------------------
cat("Ejecutando Método 5: missForest...\n")
vars_rf <- c("Temperature.F.", "Humidity...", "Sunrise_Sunset")
dd_rf_subset <- dd[, vars_rf]

rf_model <- missForest(dd_rf_subset, variablewise = T, verbose = T, ntree = 10) 
dd_rf_imputed <- dd
dd_rf_imputed[, vars_rf] <- rf_model$ximp

plot_comparativa(dd_rf_imputed, dd, "missForest")
dd_rf_subset <- NULL
dd_rf_imputed <- NULL

# ------------------------------------------------------------------------------
# MÉTODO 6: MIMMI (Clustering Imputation)
# ------------------------------------------------------------------------------
cat("Ejecutando Método 6: MIMMI...\n")
uncompleteVar <- function(vector){any(is.na(vector))}
Mode_MIMMI <- function(x) {
  x <- as.factor(x)
  maxV <- which.max(table(x))
  return(levels(x)[maxV])
}

MiMMi <- function(data, priork=-1) {
  colsMiss <- which(sapply(data, uncompleteVar))
  if(length(colsMiss) == 0){
    return(list(imputedData = data, imputation = NULL))
  } else {
    K <- dim(data)[2]
    colsNoMiss <- setdiff(c(1:K), as.vector(colsMiss))
    dissimMatrix <- daisy(data[,colsNoMiss, drop=FALSE], metric = "gower", stand=TRUE)
    distMatrix <- dissimMatrix^2
    hcdata <- hclust(distMatrix, method = "ward.D2")
    nk <- ifelse(priork == -1, 5, priork) 
    partition <- cutree(hcdata, nk)
    CompleteData <- data
    newCol <- K + 1
    CompleteData[,newCol] <- partition
    names(CompleteData)[newCol] <- "ClassAux"
    setOfClasses <- as.numeric(levels(as.factor(partition)))
    imputationTable <- data.frame(row.names=setOfClasses)
    p <- 1
    for(k in colsMiss) {
      rowsWithFullValues <- !is.na(CompleteData[,k])
      if(is.numeric(CompleteData[,k])) {
        imputingValues <- aggregate(CompleteData[rowsWithFullValues,k], by=list(partition[rowsWithFullValues]), FUN=mean)
      } else {
        imputingValues <- aggregate(CompleteData[rowsWithFullValues,k], by=list(partition[rowsWithFullValues]), FUN=Mode_MIMMI)
      }
      for(c in setOfClasses) {
        valor_imp <- imputingValues[imputingValues[,1] == c, 2]
        if(length(valor_imp) > 0) {
          CompleteData[is.na(CompleteData[,k]) & partition==c, k] <- valor_imp
        }
      }
      imputationTable[,p] <- imputingValues[,2]
      names(imputationTable)[p] <- names(data)[k]
      p <- p + 1
    }
    rownames(imputationTable) <- paste0("c", 1:nk)
    out <- list(imputedData = CompleteData, imputation = imputationTable)
    return(out)
  }
}

columnas_numericas <- names(dd)[sapply(dd, is.numeric)]
missings_por_col <- colSums(is.na(dd[, columnas_numericas]))
cols_completas <- names(missings_por_col[missings_por_col == 0])

vars_mimmi <- c("Temperature.F.", "Humidity...", "Sunrise_Sunset", cols_completas[1:2])
dd_mimmi_subset <- dd[, vars_mimmi]

mimmi_results <- MiMMi(dd_mimmi_subset, priork = 5)
dd_mimmi_imputed <- dd
dd_mimmi_imputed[, vars_mimmi] <- mimmi_results$imputedData

plot_comparativa(dd_mimmi_imputed, dd, "MIMMI")
dd_mimmi_subset <- NULL
mimmi_results <- NULL
dd_mimmi_imputed <- NULL

# ------------------------------------------------------------------------------
# MÉTODO 7: KNN (k=5)
# ------------------------------------------------------------------------------
cat("\nEjecutando Método 7: kNN...\n")
dd_limpio <- kNN(dd, variable = vars_a_imputar, k = 5)

plot_comparativa(dd_limpio, dd, "kNN")
write.csv(dd_limpio, "US_Accidents_Final_imputado_knn.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# MÉTODO 8: MICE
# ------------------------------------------------------------------------------
cat("\nEjecutando Método 8: MICE...\n")

datos_mice <- dd[, vars_a_imputar]

imputed_Data <- mice(datos_mice, m = 5, maxit = 50, method = "pmm", seed = 500)
summary(imputed_Data)

dd_mice <- dd
dd_mice[, names(datos_mice)] <- mice::complete(imputed_Data, action = 1)

plot_comparativa(dd_mice, dd, "MICE")
write.csv(dd_mice, "US_Accidents_Final_imputado_mice.csv", row.names = FALSE)
saveRDS(dd_mice, "US_Accidents_Final_imputado_mice.rds")
