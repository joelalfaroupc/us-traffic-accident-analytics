library(tm)

df <- readRDS("Datasets/US_Accidents_Final_FE.rds")

text_vector <- as.character(df$Description)
corpus <- Corpus(VectorSource(text_vector))

english_stopwords <- readLines("https://slcladal.github.io/resources/stopwords_en.txt", encoding = "UTF-8")

# dot, toc y stmc son nombres de agencias de tráfico
domain_stopwords <- c("accident", "road", "street", "st", "ave", "avenue", "traffic",
                      "highway", "hwy", "due", "new", "jersey", "dot", "toc", "stmc")

# Une carreteras como I-95 a i95, US-1 a us1 ANTES de borrar puntuación
protect_highways <- content_transformer(function(x) {
  # Vías federales/interestatales
  x <- gsub("\\bI[- ]?(\\d+)\\b", "i\\1", x, ignore.case = TRUE)
  x <- gsub("\\bUS[- ]?(\\d+)\\b", "us\\1", x, ignore.case = TRUE)
  # Rutas genericas
  x <- gsub("\\bSR[- ]?(\\d+)\\b", "sr\\1", x, ignore.case = TRUE)
  x <- gsub("\\bSH[- ]?(\\d+)\\b", "sh\\1", x, ignore.case = TRUE)
  x <- gsub("\\bRT[- ]?(\\d+)\\b", "rt\\1", x, ignore.case = TRUE)
  # Rutas específicas de nj
  x <- gsub("\\bNJ[- ]?(\\d+)\\b", "nj\\1", x, ignore.case = TRUE)
  x <- gsub("\\bCR[- ]?(\\d+)\\b", "cr\\1", x, ignore.case = TRUE)
  return(x)
})
# Función para unificar las salidas de autopista
protect_exits <- content_transformer(function(x) {
  # Busca la palabra "exit", seguida de posibles espacios o signos, y un número. 
  # Lo fusiona todo en un solo bloque, ej: "exit33"
  x <- gsub("(?i)\\bexit[[:space:][:punct:]]*(\\d+)\\b", "exit\\1", x, perl = TRUE)
  return(x)
})

replace_symbols_with_space <- content_transformer(function(x) {
  gsub("[/-]", " ", x)
})

remove_short_words <- content_transformer(function(x) {
  gsub("\\b[a-z]{1,2}\\b", "", x)
})

# Función para eliminar SOLO números sueltos, protegiendo alfanuméricos
remove_standalone_numbers <- content_transformer(function(x) {
  gsub("\\b\\d+\\b", "", x)
})

# Función de estandarización de cantidad de carriles
unify_lane_counts <- function(x) {
  
  x <- gsub("\\b1\\b(?=\\s*(left|right|center)?\\s*lanes?)", "one", x, perl = TRUE, ignore.case = TRUE)
  x <- gsub("\\b2\\b(?=\\s*(left|right|center)?\\s*lanes?)", "two", x, perl = TRUE, ignore.case = TRUE)
  x <- gsub("\\b3\\b(?=\\s*(left|right|center)?\\s*lanes?)", "three", x, perl = TRUE, ignore.case = TRUE)
  x <- gsub("\\b4\\b(?=\\s*(left|right|center)?\\s*lanes?)", "four", x, perl = TRUE, ignore.case = TRUE)
  
  return(x)
}

# Preprocessing ----------------------------------------------------------------
## 1. Aplicar función de protección de ent. de trafico (I-, US-, SR-)
corpus <- tm_map(corpus, protect_highways)
corpus <- tm_map(corpus, protect_exits)

## 2. Realizar la disminución de las letras
corpus <- tm_map(corpus, content_transformer(tolower))

## 3. Reemplazar guiones y barras por espacios para evitar fusiones
corpus <- tm_map(corpus, replace_symbols_with_space)

## 4. Eliminar el resto de la puntuación
corpus <- tm_map(corpus, removePunctuation)

## 5* Reemplazar numeros sueltos que hagan referencia a cantidad de vías
corpus <- tm_map(corpus, content_transformer(unify_lane_counts))

## 5. Eliminamos las palabras vacías
corpus <- tm_map(corpus, removeWords, english_stopwords)
corpus <- tm_map(corpus, removeWords, stopwords("en"))
corpus <- tm_map(corpus, removeWords, domain_stopwords) 

## 6. Eliminar numeros sueltos
corpus <- tm_map(corpus, remove_standalone_numbers)

## 7. Aplicamos el stemming 
corpus <- tm_map(corpus, stemDocument, language = "english")

## 8. Eliminación de palabras residuales (Solo letras)
corpus <- tm_map(corpus, remove_short_words)

## 9. Eliminar espacios sobrantes
corpus <- tm_map(corpus, stripWhitespace)

saveRDS(corpus, file = "Datasets/corpus_clean.rds")
