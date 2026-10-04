library(FactoMineR)
library(factoextra)
library(tm)
library(Matrix)

df <- readRDS("Datasets/US_Accidents_Final_FE.rds")
corpus <- readRDS("Datasets/corpus_light.rds")

# Construir la DTM
minimumFrequency <- 5
dtm <- DocumentTermMatrix(corpus, control = list(bounds = list(global = c(minimumFrequency, Inf))))

# Convertir a matriz
dtm_mat <- as.matrix(dtm)

# Eliminar documentos vacíos
dtm_mat <- dtm_mat[rowSums(dtm_mat) > 0, ]

# Ver dimensión de la matriz (número de documentos y número de términos)
dim(dtm_mat)

# Realizar CA
res.ca <- CA(dtm_mat, graph = FALSE)
res.ca

# Eigenvalues
res.ca$eig
fviz_screeplot(res.ca, addlabels = TRUE, ylim = c(0, 50))


# Documentos (filas)
row <- get_ca_row(res.ca)
row

# Documentos que más contribuyen a la primera dimensión
fviz_contrib(res.ca, choice = "row", axes = 1, top = 20)

# Documentos que más contribuyen a la primera y segunda dimensión
fviz_contrib(res.ca, choice = "row", axes = 1:2, top = 20)

# Bi-plot de documentos (dimensiones 1 y 2)
plot.CA(res.ca, invisible = "col", axes = c(1,2))

# Bi-plot de documentos sin etiquetas (dimensiones 1 y 2)
plot.CA(res.ca, invisible = "col", label = "none", axes = c(1,2))

# Bi-plot de documentos en las dimensiones 3 y 4
plot.CA(res.ca, invisible = "col", axes = c(3,4))


# Términos (columnas)
col <- get_ca_col(res.ca)
col

# Términos que más contribuyen a la primera dimensión
fviz_contrib(res.ca, choice = "col", axes = 1, top = 20)

# Términos que más contribuyen a la primera y segunda dimensión
fviz_contrib(res.ca, choice = "col", axes = 1:2, top = 20)

# Bi-plot de términos (dimensiones 1 y 2)
plot.CA(res.ca, invisible = "row", axes = c(1,2))

# Bi-plot de términos sin etiquetas (dimensiones 1 y 2)
plot.CA(res.ca, invisible = "row", label = "none", axes = c(1,2))

# Bi-plot de términos en las dimensiones 3 y 4
plot.CA(res.ca, invisible = "row", axes = c(3,4))


# Bi-plot con los 20 documentos y términos más contribuyentes (dimensiones 1 y 2)
plot.CA(res.ca, autoLab = "no", selectRow = "contrib 20", selectCol = "contrib 20", axes = c(1,2))


# # Seleccionar variables categóricas
# quali_vars <- df[, c("Severity", "Weather_Group", "Wind_Group")]
# 
# # Comprobar que son factores
# quali_vars <- data.frame(lapply(quali_vars, factor))
# 
# # Comprobar que el número de filas sea igual para realizar CA-GALT
# dim(dtm_mat)
# dim(quali_vars)
# 
# # Mirar cuál es la fila que aparece en quali_vars y no en dtm_mat
# setdiff(rownames(quali_vars), rownames(dtm_mat))
# 
# # Eliminar esta fila
# quali_vars <- quali_vars[-which(rownames(quali_vars) == "8282"), ]
#
# # Realizar CA-GALT
# res.cagalt <- CaGalt(Y = dtm_mat, X = quali_vars, type = "n", graph = FALSE)
# res.cagalt
# 
# # Eigenvalues
# res.cagalt$eig
# 
# # Información básica de los individuos
# names(res.cagalt$ind)
# res.cagalt$ind$coord
# res.cagalt$ind$cos2
# 
# # Coordenadas y calidad de palabras
# names(res.cagalt$freq)
# res.cagalt$freq$coord
# res.cagalt$freq$cos2
# res.cagalt$freq$contr
# 
# # 10 palabras más contribuyentes
# res.cagalt$freq$contr[order(apply(res.cagalt$freq$contr[,1:2], 1, sum), decreasing = TRUE)[1:10], 1:2]
# 
# # Variables categóricas
# names(res.cagalt$quali.var)
# res.cagalt$quali.var$coord
# res.cagalt$quali.var$cos2
# 
# # Resumen
# summary(res.cagalt)
# 
# # CA-GALT plots
# plot.CaGalt(res.cagalt, choix = "freq", axes = c(1,2))
# plot.CaGalt(res.cagalt, choix = "freq", axes = c(1,2), select = "contrib 20")
# plot.CaGalt(res.cagalt, choix = "quali.var", axes = c(1,2), autoLab = "no")
# 
# plot(res.cagalt, choix = "quali.var", conf.ellip = TRUE, axes = c(1,2))
# plot(res.cagalt, choix = "freq", cex = 1.5, col.freq = "darkgreen", select = "contrib 10")
# plot(res.cagalt, choix = "ind", cex = 1.5, col.ind = "darkgreen", select = "contrib 10")
# 
# par(mfrow=c(1,3))
# plot.CaGalt(res.cagalt, select = "cos2 10")
# plot.CaGalt(res.cagalt, choix="freq", axes=c(1,2), select="contrib 10")
# plot(res.cagalt, choix = "quali.var", conf.ellip = FALSE, axes = c(1,2), select="cos2 5")
# par(mfrow=c(1,1))