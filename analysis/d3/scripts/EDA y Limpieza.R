library(tidyverse)
library(visdat)
library(inspectdf)
library(skimr)
library(DataExplorer)
library(SmartEDA)

dades <- read.csv("Datasets/US_Accidents_Final.csv")

#####################################################
#            PLOTS SOBRE LA BASE INICIAL            #
#####################################################

# Pre: Visualización inicial de los datos sucios
ExpData(data = dades,type = 1)

skim(dades)
plot_intro(dades)

# 2. Evaluación de Missings
plot_missing(dades)
vis_miss(dades, warn_large_data = FALSE)

# 3. Análisis de Distribuciones Numéricas
plot_histogram(dades)

# 4. Análisis de Desbalanceo (binarias y categóricas)
inspect_imb(dades) %>% show_plot()
plot_bar(dades, maxcat = 31)

###################################################
#            LIMPIEZA DE DATOS GENERAL            #
###################################################

# Limpieza general
dades_clean <- dades %>%
  
  # 1. Eliminar variables con un solo valor (varianza cero)
  select(-State, -Country, -Timezone, -Bump, -Roundabout, -Turning_Loop) %>%
  
  # 2. Eliminar ID, Source, End_Lat/Lng y Street por cardinalidad/no aporta información
  select(-ID, -Source, -End_Lat, -End_Lng, -Street, -City, -Zipcode) %>%
  
  # 3. Convertir entradas de texto vacío "" a NAs
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  
  # 4. Convertir la variable objetivo a factor ordenado
  mutate(Severity = factor(Severity, levels = c(1, 2, 3, 4), ordered = TRUE)) %>%
  # 5. Eliminar variables temporales redundantes
  select(-Weather_Timestamp, -End_Time) %>%
  
  # 6. Transformar a tipos de variables a factores
  mutate(across(c(Amenity, Crossing, Give_Way, Junction, No_Exit, 
                  Railway, Station, Stop, Traffic_Calming, Traffic_Signal), 
                ~ as.factor(.))) %>%
  
  mutate(across(c(where(is.character), -Description, -Start_Time), as.factor))
  

######################
#         EDA        #
######################

### UNIVARIANTE ###

# 1. Visualización inicial de los datos limpios
skim(dades_clean)
plot_intro(dades_clean)

# 2. Evaluación de Missings
plot_missing(dades_clean)
vis_miss(dades_clean, warn_large_data = FALSE)

# 3. Análisis de Distribuciones Numéricas
plot_histogram(dades_clean)

# 4. Análisis de Desbalanceo (binarias y categóricas)
inspect_imb(dades_clean) %>% show_plot()
plot_bar(dades_clean, maxcat = 31)

### BIVARIANTE ###

# 1. Matriz de correlación vars numericas incluyendo Severity
dades_clean %>%
  mutate(Severity_Num = as.numeric(as.character(Severity))) %>%
  select(where(is.numeric)) %>%
  vis_cor()

# 2. Relación categóricas/binarias vs Severity
plot_bar(dades_clean, by = "Severity")

# 3. Relación numéricas vs Severity
plot_boxplot(dades_clean, by = "Severity")

### SMARTEDA ###

# Overview
ExpData(data = dades_clean,type = 1)
# Estructura
ExpData(data = dades_clean,type = 2)

### VARIABLES NUMERICAS ###

## Summary statistics by – overall
ExpNumStat(dades_clean,by="A",gp=NULL,Qnt=seq(0,1,0.1),MesofShape=2,Outlier=TRUE,round=2)

## Summary statistics by – overall with correlation 
ExpNumStat(dades_clean,by="A",gp="Severity",Qnt=seq(0,1,0.1),MesofShape=1,Outlier=TRUE,round=2)

## Summary statistics by – category
ExpNumStat(dades_clean,by="GA",gp="Severity",Qnt=seq(0,1,0.1),MesofShape=2,Outlier=TRUE,round=2)

# GRAFICOS

## Boxplot by category
ExpNumViz(dades_clean, target="Severity",type=2,nlim=25)
## Density plot
ExpNumViz(dades_clean,target=NULL,type=3,nlim=25)

### VARIABLES CATEGORICAS/BIN ###

## Frequency or custom tables for categorical variables
ExpCTable(dades_clean,Target=NULL,margin=1,clim=10,nlim=5,round=2,bin=NULL,per=T)
## Summary statistics of categorical variables
ExpCatStat(dades_clean,Target="Severity",result = "Stat",clim=10,nlim=5)

# GRAFICOS

## column chart
ExpCatViz(dades_clean,target="Severity",fname=NULL,clim=10,col=NULL,margin=2,sample=14)

write.csv(dades_clean, "US_Accidents_Clean.csv", row.names = FALSE)
saveRDS(dades_clean, "US_Accidents_Clean.rds")