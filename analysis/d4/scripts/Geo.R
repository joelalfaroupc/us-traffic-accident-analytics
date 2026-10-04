library("raster")
library("spatstat")

library("sf")
library("maps")

db <- readRDS("Datasets/US_Accidents_Final_FE.rds")

# Cargar los datos en un df alternativo
xy = data.frame(x=db$Start_Lng, y=db$Start_Lat)

plot(xy$x,xy$y,xlab="x",ylab="y",main = "Nube de puntos crudos")

df1 = data.frame(xy, Severity=as.factor(db$Severity), Temperature=db$Temperature.C., 
                 Humidity=db$Humidity..., Weather=db$Weather_Group, Sunrise_Sunset=db$Sunrise_Sunset)

# Definir los sistemas de coordenadas
crs_original <- 4326  # Longitud/Latitud (WGS84)
crs_proyecto <- 32618 # UTM Zona 18N (Proyectado en metros)

# Descargar, definir y proyectar el polígono de NJ
nj_sf <- st_as_sf(maps::map("state", "new jersey", plot = FALSE, fill = TRUE))
st_crs(nj_sf) <- crs_original
nj_sf_proj <- st_transform(nj_sf, crs_proyecto)

nj_win <- as.owin(nj_sf_proj)

# 3. PROYECCIÓN DE LOS PUNTOS
puntos_sf <- st_as_sf(xy, coords = c("x", "y"), crs = crs_original)
puntos_proj <- st_transform(puntos_sf, crs_proyecto)
coords_proj <- st_coordinates(puntos_proj)

# Filtro para descartar puntos fuera del estado
dentro_nj <- inside.owin(x = coords_proj[,1], y = coords_proj[,2], w = nj_win)

coords_proj <- coords_proj[dentro_nj, ]
df1 <- df1[dentro_nj, ]
xy <- xy[dentro_nj, ]

# 4. Crear los objetos ppp y reescalar a km
y <- ppp(x = coords_proj[,1], y = coords_proj[,2], window = nj_win)
y1 <- ppp(x = coords_proj[,1], y = coords_proj[,2], window = nj_win, marks = df1[,3:7])

y <- rescale(y, 1000, "km")
y1 <- rescale(y1, 1000, "km")

# Graficos de puntos según cada marca
plot(y1,which.marks="Severity",size=0.8,cols = c("green", "yellow", "orange", "red"),main="Dispersión de accidentes según Severidad")

plot(y1,which.marks="Sunrise_Sunset",size=0.8,cols = c("gold", "darkblue"),main="Dispersión de accidentes según Dia/Noche")

plot(y1, 
     which.marks = "Temperature", 
     symap = symbolmap(pch = 20, 
                       cex = 0.7,
                       range = range(marks(y1)$Temperature, na.rm = TRUE), 
                       col = c("blue", "yellow", "red")), 
     main = "Dispersión de accidentes según Temperatura")

plot(y1, 
     which.marks = "Humidity", 
     symap = symbolmap(pch = 20, 
                       cex = 0.7,
                       range = range(marks(y1)$Humidity, na.rm = TRUE), 
                       col = c("khaki", "deepskyblue", "darkblue")), 
     main = "Dispersión de accidentes según Humedad")

plot(y1,which.marks="Weather",size=0.8,cols = c("gold", "darkgray", "brown", "purple4", "aquamarine", "navy", "black", "cyan"),main="Dispersión de accidentes según Tiempo")

# Análisis de Hotspots

# Sin marcas
plot(density(y), main="Densidad total de accidentes")
contour(density(y), add=TRUE, col="white")

hist(xy$x,xlab = "x",ylab="Frecuencia",main = "Histograma de frecuencia de x",
     col = "cadetblue3")   # Histograma de la longitud.
hist(xy$y,xlab = "y",ylab="Frecuencia",main = "Histograma de frecuencia de y",
     col = "burlywood1")  #Histograma de la latitud.

# Hotspots con marcas

#### Severidad
marks(y) <- df1$Severity

mapas_separados_sev <- split(y)

# Análisis exclusivo de accidentes de Severidad 1
y_sev1 <- mapas_separados_sev$"1"

summary(y_sev1)

# Gráficas de densidad
plot(density(y_sev1), main = "Densidad de Accidentes: Severidad 1")
contour(density(y_sev1), add = TRUE, col = "white")


# Análisis exclusivo de accidentes de Severidad 2
y_sev2 <- mapas_separados_sev$"2"

summary(y_sev2)

# Gráficas de densidad
plot(density(y_sev2), main = "Densidad de Accidentes: Severidad 2")
contour(density(y_sev2), add = TRUE, col = "white")


# Análisis exclusivo de accidentes de Severidad 3
y_sev3 <- mapas_separados_sev$"3"

summary(y_sev3)

# Gráficas de densidad
plot(density(y_sev3), main = "Densidad de Accidentes: Severidad 3")
contour(density(y_sev3), add = TRUE, col = "white")


# Análisis exclusivo de accidentes de Severidad 4
y_sev4 <- mapas_separados_sev$"4"

summary(y_sev4)

# Gráficas de densidad
plot(density(y_sev4), main = "Densidad de Accidentes: Severidad 4")
contour(density(y_sev4), add = TRUE, col = "white")

#### Momento del día

marks(y) <- df1$Sunrise_Sunset

mapas_separados <- split(y)

#Análisis exclusivo de accidentes nocturnos
y_noche <- mapas_separados$Night

summary(y_noche)

# Gráficas de densidad
plot(density(y_noche), main = "Densidad de Accidentes Nocturnos")
contour(density(y_noche), add = TRUE, col = "white")

#Análisis exclusivo de accidentes diurnos
y_day <- mapas_separados$Day

summary(y_day)

# Gráficas de densidad
plot(density(y_day), main = "Densidad de Accidentes Diurnos")
contour(density(y_day), add = TRUE, col = "white")

#### Clima (Top 4 niveles)

marks(y) <- df1$Weather

mapas_separados_wea <- split(y)

# Análisis exclusivo de accidentes con clima Cloudy (Nublado)
y_cloudy <- mapas_separados_wea$Cloudy

summary(y_cloudy)

# Gráficas de densidad
plot(density(y_cloudy), main = "Densidad de Accidentes: Cloudy")
contour(density(y_cloudy), add = TRUE, col = "white")


# Análisis exclusivo de accidentes con clima Clear (Despejado)
y_clear <- mapas_separados_wea$Clear

summary(y_clear)

# Gráficas de densidad
plot(density(y_clear), main = "Densidad de Accidentes: Clear")
contour(density(y_clear), add = TRUE, col = "white")


# Análisis exclusivo de accidentes con clima Light_Rain (Lluvia ligera)
y_lightrain <- mapas_separados_wea$Light_Rain

summary(y_lightrain)

# Gráficas de densidad
plot(density(y_lightrain), main = "Densidad de Accidentes: Light Rain")
contour(density(y_lightrain), add = TRUE, col = "white")


# Análisis exclusivo de accidentes con clima Heavy_Storms (Tormentas fuertes)
y_heavystorms <- mapas_separados_wea$Heavy_Storms

summary(y_heavystorms)

# Gráficas de densidad
plot(density(y_heavystorms), main = "Densidad de Accidentes: Heavy Storms")
contour(density(y_heavystorms), add = TRUE, col = "white")

#### NUMERICAS ####
#### Humedad

# Discretización
humedad_clase <- ifelse(df1$Humidity > 84, "Humedad_Extrema", "Humedad_Normal")

df1$Humedad_Factor <- factor(humedad_clase)

# Marca y División
marks(y) <- df1$Humedad_Factor
mapas_humedad <- split(y)

# Humedad Extrema
y_hum_extrema <- mapas_humedad$Humedad_Extrema

summary(y_hum_extrema)

plot(density(y_hum_extrema), main = "Densidad: Accidentes con Humedad Extrema (>84%)")
contour(density(y_hum_extrema), add = TRUE, col = "white")

# Humedad Normal
y_hum_normal <- mapas_humedad$Humedad_Normal

plot(density(y_hum_normal), main = "Densidad: Accidentes con Humedad Normal (<=84%)")
contour(density(y_hum_normal), add = TRUE, col = "white")

##### Temperatura

# Discretización
clase_termica <- ifelse(df1$Temperature < 5, "Frio_Severo", 
                        ifelse(df1$Temperature > 25, "Calor_Severo", "Normal"))

df1$Temp_Factor <- factor(clase_termica, levels = c("Frio_Severo", "Normal", "Calor_Severo"))

# Marca y división
marks(y) <- df1$Temp_Factor
mapas_temp <- split(y)

# 1. Frío severo
y_frio <- mapas_temp$Frio_Severo
cat("Intensidad en Frío Severo:\n")
intensity(y_frio)

plot(density(y_frio), main = "Densidad: Accidentes con Frío Severo (< 5ºC)")
contour(density(y_frio), add = TRUE, col = "white")

# 2. Temp. normal
y_norm <- mapas_temp$Normal
cat("Intensidad en Temperaturas Normales:\n")
intensity(y_norm)

plot(density(y_norm), main = "Densidad: Accidentes Temp. Media ( 5ºC < T < 25ºC)")
contour(density(y_norm), add = TRUE, col = "white")

# 3. Calor severo
y_calor <- mapas_temp$Calor_Severo
cat("Intensidad en Calor Severo:\n")
intensity(y_calor)

plot(density(y_calor), main = "Densidad: Accidentes con Calor Severo (> 25ºC)")
contour(density(y_calor), add = TRUE, col = "white")

##### Análisis de Cuadrícula

Q <- quadratcount(y, nx = 4, ny = 4)
#Gráfica de cuadrículas
plot(intensity (Q, image=T))
plot(Q, add=T)

# Gráfico 3D de la densidad con las 5 ciudades de NJ más importantes superpuestas

# 1. Definir ciudades en coordenadas originales (grados)
ciudades_base <- data.frame(
  nombre = c("Newark", "Jersey City", "Paterson", "Elizabeth", "Trenton"),
  lng = c(-74.172, -74.078, -74.172, -74.211, -74.760),
  lat = c(40.736, 40.728, 40.917, 40.664, 40.221)
)

# 2. PROYECTAR las ciudades a UTM 18N (metros) y CONVERTIR a kilómetros
ciudades_sf <- st_as_sf(ciudades_base, coords = c("lng", "lat"), crs = crs_original)
ciudades_proj <- st_transform(ciudades_sf, crs_proyecto)

coords_ciudades <- st_coordinates(ciudades_proj) / 1000 

# 3. Calcular la densidad y LIMPIAR LOS NA
densidad_mapa <- density(y)
densidad_mapa$v[is.na(densidad_mapa$v)] <- 0

altura_base <- max(densidad_mapa$v, na.rm = TRUE)

# 4. Reconstruir el dataframe de ciudades con las coordenadas en KILÓMETROS
ciudades <- data.frame(
  nombre = ciudades_base$nombre,
  x = coords_ciudades[,1], 
  y = coords_ciudades[,2], 
  z = c(
    altura_base * 1.2,  # Newark
    altura_base * 0.9, # Jersey City
    altura_base * 0.8,    # Paterson
    altura_base * 0.6,  # Elizabeth
    altura_base * 0.4   # Trenton
  )
)

# 5. Mapa "crudo" 3D
matriz_3d <- persp(densidad_mapa, theta = 290, phi = 20, col = "lightblue", 
                   main = "Densidad de Accidentes según ciudades", border = NA, shade = 0.6)

# 6. Crear los puntos 3D usando las coordenadas correctas (x, y)
puntos_proyectados <- trans3d(x = ciudades$x, 
                              y = ciudades$y, 
                              z = ciudades$z, 
                              pmat = matriz_3d)

# 7. Proyección de los puntos sobre el mapa
altura_superficie <- densidad_mapa[list(x = ciudades$x, y = ciudades$y)]

suelo_proyectado <- trans3d(x = ciudades$x, 
                            y = ciudades$y, 
                            z = altura_superficie, 
                            pmat = matriz_3d)

# 8. Dibujar elementos
points(puntos_proyectados, col = "red", pch = 16, cex = 1.5)
text(puntos_proyectados, labels = ciudades$nombre, col = "black", pos = 3, font = 2)

segments(x0 = puntos_proyectados$x, y0 = puntos_proyectados$y,
         x1 = suelo_proyectado$x, y1 = suelo_proyectado$y,
         col = "red", lty = 3)


