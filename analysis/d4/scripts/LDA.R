library(topicmodels)
library(RColorBrewer)
library(wordcloud)
library(remotes)
library(tm)
library(ggplot2)
library(tidytext)
library(dplyr)
library(tidyr)
remotes::install_github("nikita-moor/ldatuning")

library(ldatuning)
library(lda)

textdata <- readRDS("Datasets/US_Accidents_Final_FE.rds")
processedCorpus <- readRDS("Datasets/corpus_clean.rds")

# Creación de la matriz DTM con terms >= minimumFrequency
minimumFrequency <- 5
DTM <- DocumentTermMatrix(processedCorpus, control = list(bounds = list(global = c(minimumFrequency, Inf))))

# Dim. iniciales DTM
dim(DTM)

# Eliminar filas vacías de la matriz DTM y del dataset
sel_idx <- slam::row_sums(DTM) > 0
DTM <- DTM[sel_idx, ]
textdata <- textdata[sel_idx, ]

#Dim. DTM tras eliminar filas vacías 
dim(DTM)


##Función para hallar el K o número de Tópicos
# create models with different number of topics
result <- ldatuning::FindTopicsNumber(
  DTM,
  topics = seq(from = 2, to = 20, by = 1),
  metrics = c("CaoJuan2009",  "Deveaud2014"),
  method = "Gibbs",
  control = list(seed = 77),
  verbose = TRUE
)
FindTopicsNumber_plot(result)

# CaoJuan minimiza en 4 tópicos, mientras que Deveaud maximiza en 2. Como ambas métricas son aceptables en ambos valores de K
# y analizar tan solo 2 tópicos puede ser poco conclusivo, se llevará a cabo el análisis para K = 2 y K = 4. 

### ANÁLISIS CON K = 2 ###

K <- 2
set.seed(67)

# compute the LDA model, inference via 1000 iterations of Gibbs sampling
topicModel <- LDA(DTM, K, method="Gibbs", control=list(iter = 500, verbose = 25,alpha=0.2))
topicModel
attr(topicModel,"alpha")

# Resultados (dists. a posteriori)
tmResult <- posterior(topicModel)

# format of the resulting object
attributes(tmResult)

beta <- tmResult$terms
dim(beta)                # K topics dist. over Terms
theta <- tmResult$topics 
dim(theta)               # Docs distributions over K topics
summary(theta)

# Top de terminos por topico
terms(topicModel, 15)
topicNames <- apply(lda::top.topic.words(beta, 5, by.score = F), 2, paste, collapse = " ")
topicNames

# Nubes de Palabras

for(topicToViz in 1:K){
  
  top50terms <- sort(tmResult$terms[topicToViz, ], decreasing = TRUE)[1:50]
  words <- names(top50terms)
  probabilities <- sort(tmResult$terms[topicToViz, ], decreasing = TRUE)[1:50]
  mycolors <- brewer.pal(8, "Dark2")
  wordcloud(words, probabilities, random.order = FALSE, color = mycolors)

}

# Ver composición de tópicos en documentos
exampleIds <- c(75, 2700, 4500, 6250, 8750)
N <- length(exampleIds) # Number of example documents

# Get topic proportions from example documents
topicProportionExamples <- theta[exampleIds, ]
colnames(topicProportionExamples) <- topicNames

# Reshape data for visualization
vizDataFrame <- reshape2::melt(
  cbind(data.frame(topicProportionExamples),
        document = factor(1:N)
  ),
  variable.name = "topic",
  id.vars = "document")

ggplot(data = vizDataFrame, aes(topic, value, fill = document), ylab = "proportion") +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() +
  facet_wrap(~document, ncol = N)

# Mostrar el texto de cada doc. para validación
textdata$Description[exampleIds]


### ANÁLISIS CON K = 4 ###

K <- 4
set.seed(67)

# compute the LDA model, inference via 1000 iterations of Gibbs sampling
topicModel <- LDA(DTM, K, method="Gibbs", control=list(iter = 500, verbose = 25,alpha=0.2))
topicModel
attr(topicModel,"alpha")

# Resultados (dists. a posteriori)
tmResult <- posterior(topicModel)

# format of the resulting object
attributes(tmResult)

beta <- tmResult$terms
dim(beta)                # K topics dist. over Terms
theta <- tmResult$topics 
dim(theta)               # Docs distributions over K topics
summary(theta)

# Top de terminos por topico
terms(topicModel, 15)
topicNames <- apply(lda::top.topic.words(beta, 5, by.score = F), 2, paste, collapse = " ")
topicNames

# Nubes de Palabras

for(topicToViz in 1:K){
  
  top50terms <- sort(tmResult$terms[topicToViz, ], decreasing = TRUE)[1:50]
  words <- names(top50terms)
  probabilities <- sort(tmResult$terms[topicToViz, ], decreasing = TRUE)[1:50]
  mycolors <- brewer.pal(8, "Dark2")
  wordcloud(words, probabilities, random.order = FALSE, color = mycolors)
  
}

# Ver composición de tópicos en documentos
exampleIds <- c(75, 2700, 4500, 6250, 8750)
N <- length(exampleIds) # Number of example documents

# Get topic proportions from example documents
topicProportionExamples <- theta[exampleIds, ]
colnames(topicProportionExamples) <- topicNames

# Reshape data for visualization
vizDataFrame <- reshape2::melt(
  cbind(data.frame(topicProportionExamples),
        document = factor(1:N)
  ),
  variable.name = "topic",
  id.vars = "document")

ggplot(data = vizDataFrame, aes(topic, value, fill = document), ylab = "proportion") +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() +
  facet_wrap(~document, ncol = N)

# Mostrar el texto de cada doc. para validación
textdata$Description[exampleIds]



