################ MALDI-TOF ANALISIS CLP ########################################
################ 2) ns_m122_d1247  #############################################
#
# ns:    No Supervisado
# m122:  Utiliza las wells (122 muestras)
# d1247: Utiliza todos los días de lectura para los análisis
#
# Autor: Bioing. Facundo Urteaga (IBB-CONICET)
#
#
### CARGA DE LIBRERIAS #########################################################
################################################################################


library("readBrukerFlexData")
library("binda")
library("fs")
library("readxl")
library("MALDIquant")
library("MALDIquantForeign")
library("MALDIrppa")
library("tidyverse")
library("dplyr")
library("clValid")
library(cluster)
library(factoextra)
library(ggplot2)
library(gridExtra)


### CARGA DE ARCHIVOS ##########################################################
################################################################################

data_path <- file.path("Data")

# Load the Rdata files using the relative path
load(file.path(data_path, "matint_122_dico.Rdata"))
load(file.path(data_path, "matint_51_dico.Rdata"))
load(file.path(data_path,"matint_122.Rdata"))
load(file.path(data_path,"matint_51.Rdata"))

df_metadata_prom_rep$dia <- as.integer(gsub("[^0-9]", "", 
                                            df_metadata_prom_rep$dia))


### SELECCIÓN DE PICOS #########################################################
################################################################################


# Selección de picos para binary discriminant analysis (BDA)
factor_tipo <- factor(df_metadata_prom_rep$tipo)
is.binaryMatrix(matint_122_dico) # TRUE
br <- binda.ranking(matint_122_dico, factor_tipo, verbose = FALSE)

# Gráfico de picos vs score 
nueva_columna <- c()
matriz <- matrix(br, nrow = 218, ncol = 4) #218 es la cantidad de picos
for (i in 1:218) {
    nuevo_valor <- colnames(matint_122_dico)[br[i]]
    nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)
plot(df_br$nueva_columna, df_br$V2, 
     xlab = "m/z", ylab = "Score", 
     main = "Ranking de picos de los espectros días 1, 2, 4 y 7")
# Crear un gradiente de colores (por ejemplo, de azul a rojo)
colores <- colorRampPalette(c("green4", "red2"))(218)
# Agregar puntos con colores en forma de gradiente
for (i in 1:218) {
  points(df_br$nueva_columna[i], df_br$V2[i], col = colores[i]) 
}
# Agregar puntos con relleno de colores en forma de gradiente
for (i in 1:218) {
  points(df_br$nueva_columna[i], df_br$V2[i], pch = 19, col = colores[i]) 
}

# Selección de picos mas preponderantes
top.b5 <- br[1:5]  ## primeros 5 picos
top.b10 <- br[1:10]  ## primeros 10 picos
top.b15 <- br[1:15]  ## primeros 15 picos
top.b20 <- br[1:20]  ## primeros 20 picos 
top_actual <- top.b20

# Elección de mejores algoritmos de clustering
comparacion <- clValid(
  obj        = matint_122_dico[, top_actual],
  nClust     = 2:6,
  clMethods  = c("hierarchical", "kmeans", "pam"),
  validation = c("stability", "internal")
)
summary(comparacion)
optimalScores(comparacion) #Se puede ir probando con distintos top picos


### ALGORITMO DE CLUSTERING ####################################################
################################################################################

# k-means clustering con top20 y 2 clusters

top_actual <- top.b20
K.num <- 2 # clusters
var2 = 0.95

km.res20 <- kmeans(matint_122_dico[, top_actual], 
                   K.num, 
                   nstart = 25)

cluster_kmean20 <- fviz_cluster(km.res20, ellipse.type = "convex", 
                                data = matint_122_dico[, top_actual],
                                ellipse.level = var2,
                                show.clust.cent = F, 
                                geom = "point", main = "CLP vs SHAM - días: 1,2,4 y 7 - kmeans - Top 20 - 2 cluster")

# Personalización del ploteo
cluster_kmean20 <- cluster_kmean20 + 
  geom_point(data = cluster_kmean20$data, 
             aes(x = x, y = y, color = factor_tipo,
                 size = df_metadata_prom_rep$dia)) +
  scale_color_manual(values = c("maroon","steelblue4", "maroon", "steelblue4")) +
  scale_size_continuous(range = c(2, 4)) + 
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

# Muestra el gráfico
print(cluster_kmean20)


### TASA DE ACIERTO POR CLUSER #################################################
################################################################################


for(i in 1:length(km.res20[[1]])) {
  km.res20$names[i] = sub("_.*", "", attr(km.res20[[1]][i],"names"))
  km.res20$dia[i] = df_metadata_prom_rep$dia[i]
}

# Calcular la tasa de acierto para todo el conjunto de datos
total_CLP <- sum(km.res20$names == "CLP")  # Total de muestras CLP
total_SHAM <- sum(km.res20$names == "SH")  # Total de muestras SHAM

# Muestras CLP Y SHAM clasificadas correctamente
acierto_CLP <- sum(km.res20$names == "CLP" & km.res20$cluster == 2)
acierto_SHAM <- sum(km.res20$names == "SH" & km.res20$cluster == 1)

total_cluster2 <- sum(km.res20$cluster == 2)
total_cluster1 <- sum(km.res20$cluster == 1)

# Cálculo tasa de acierto por cluster
tasa_acierto_CLP <- acierto_CLP / total_cluster2
tasa_acierto_SHAM <- acierto_SHAM / total_cluster1

# Crear un dataframe con los resultados para el gráfico general
datos_grafico_general <- data.frame(Categoria = c("CLP", "SH"),
                                    Aciertos = c(acierto_CLP, acierto_SHAM),
                                    Total = c(total_cluster2, total_cluster1),
                                    Tasa_Acierto = c(tasa_acierto_CLP, 
                                                     tasa_acierto_SHAM),
                                    Dia = "Total")

# Función que calcula los parámetros por día
datos_grafico_dias <- lapply(c(1,2,4,7), function(d) {
  # Obtener los clusters asignados a cada observación
  clusters <- km.res20$cluster
  names <- km.res20$names
  # Filtrar los clusters según el valor de dia
  clusters_dia <- clusters[km.res20$dia == d]
  names_dia <- names[km.res20$dia == d]
  datos_dia <- list(cluster = clusters_dia, names = names_dia)
  
  total_CLP <- sum(datos_dia$names == "CLP")
  total_SHAM <- sum(datos_dia$names == "SH")
  acierto_CLP <- sum(datos_dia$names == "CLP" & datos_dia$cluster == 2)
  acierto_SHAM <- sum(datos_dia$names == "SH" & datos_dia$cluster == 1)
  total_cluster2 <- sum(datos_dia$cluster == 2)
  total_cluster1 <- sum(datos_dia$cluster == 1)
  tasa_acierto_CLP <- acierto_CLP / total_cluster2
  tasa_acierto_SHAM <- acierto_SHAM / total_cluster1
  return(data.frame(Categoria = c("CLP", "SH"),
                    Aciertos = c(acierto_CLP, acierto_SHAM),
                    Total = c(total_cluster2, total_cluster1),
                    Tasa_Acierto = c(tasa_acierto_CLP, tasa_acierto_SHAM),
                    Dia = as.factor(d)))
})

# Unir los datos de los días en un único dataframe
datos_grafico_dias <- do.call(rbind, datos_grafico_dias)

# Gráfico general
grafico_general <- ggplot(datos_grafico_general, 
                          aes(x = Categoria, y = Tasa_Acierto, 
                              fill = Categoria)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
  labs(title = "Tasa de Acierto Total",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CLP" = "blue4", "SH" = "red4")) +
  theme_minimal()

plot(grafico_general)

# Gráficos por día
graficos_dias <- lapply(unique(datos_grafico_dias$Dia), function(d) {
  datos_dia <- subset(datos_grafico_dias, Dia == d)
  ggplot(datos_dia, aes(x = Categoria, y = Tasa_Acierto, fill = Categoria)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
    labs(title = paste("Tasa de Acierto - Día", d),
         x = "Categoría",
         y = "Tasa de Acierto") +
    scale_fill_manual(values = c("CLP" = "blue4", "SH" = "red4")) +
    theme_minimal()
})

# Mostrar los gráficos
grid.arrange(grafico_general, grobs = graficos_dias, ncol = 2, 
             top = "Tasa de Acierto por Día")

plot

### OBTENCIÓN DE MÉTRICAS  #####################################################
################################################################################

# 1. Matriz de confusión

# Calcular los elementos de la matriz de confusión
TP <- acierto_CLP
TN <- acierto_SHAM
FP <- total_SHAM - acierto_SHAM
FN <- total_CLP - acierto_CLP

# Matriz de confusión
conf_matrix <- matrix(c(TN, FP, FN, TP), nrow = 2, byrow = TRUE,
                      dimnames = list('Referencia' = c('SHAM', 'CLP'),
                                      'Predicho' = c('cluster1', 'cluster2')))
print(conf_matrix)

# Cálculo de métricas
precision <- TP / (TP + FP)
recall <- TP / (TP + FN)
f1_score <- 2 * (precision * recall) / (precision + recall)
accuracy <- (TP + TN) / (total_CLP + total_SHAM)

# Mostrar las métricas
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1-Score:", f1_score, "\n")
cat("Accuracy:", accuracy, "\n")

# Ver parte de MATRIZ CONFUSION
# Prueba de distintas métricas

# 2. Silhouette Score
silhouette_score <- silhouette(km.res20$cluster, dist(matint_122_dico[, top_actual]))
fviz_silhouette(silhouette_score)

# 3. Within-cluster Sum of Squares (WCSS)
wcss <- km.res20$tot.withinss
print(wcss)

# 4. Between-cluster Sum of Squares (BCSS)
bcss <- km.res20$betweenss
print(bcss)
#
#
#
### FIN ########################################################################
################################################################################