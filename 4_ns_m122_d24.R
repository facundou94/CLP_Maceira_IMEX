################ MALDI-TOF ANALISIS CLP ########################################
################ 4) ns_m122_d24  ###############################################
#
# ns:    No Supervisado
# m122:  Utiliza las wells (122 muestras)
# d1247: Utiliza los días 2 y 4
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


### CARGA DE ARCHIVOS ##########################################################
################################################################################


data_path <- file.path("Data")

# Load the Rdata files using the relative path
load(file.path(data_path, "matint_122_dico.Rdata"))
load(file.path(data_path, "matint_122.Rdata"))


df_metadata_prom_rep$dia <- as.integer(gsub("[^0-9]", "", 
                                            df_metadata_prom_rep$dia))

# Filtrado de muestras de día 2 y día 4
filas_filtradas <- rownames(matint_122_dico)[grepl("D2|D4", rownames(matint_122_dico))]
matint_122_dico_d2d4 <- matint_122_dico[filas_filtradas, ]
filas_filtradas <- rownames(matint_122)[grepl("D2|D4", rownames(matint_122))]
matint_122_d2d4 <- matint_122[filas_filtradas, ]

df_prom_rep_d2d4 <- df_metadata_prom_rep %>% dplyr::filter(dia == 2 | dia == 4)


### SELECCIÓN DE PICOS #########################################################
################################################################################


# Selección de picos para binary discriminant analysis (BDA)

factor_tipo <- factor(df_prom_rep_d2d4$factor1)
is.binaryMatrix(matint_122_dico_d2d4) # TRUE
br <- binda.ranking(matint_122_dico_d2d4, factor_tipo, verbose = FALSE)
indices_especificos <- br[1:20]

# Gráfico de picos vs score 
nueva_columna <- c()
matriz <- matrix(br, nrow = 218, ncol = 4)
for (i in 1:218) {
  nuevo_valor <- colnames(matint_122_dico_d2d4)[br[i]]
  nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)

plot(df_br$nueva_columna, df_br$V2, 
     xlab = "m/z", ylab = "Score", 
     main = "Ranking de picos de los espectros, 4 clusters")
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
top.b10 <- br[1:10]  ## primeros 10 picos
top.b15 <- br[1:15]  ## primeros 15 picos
top.b20 <- br[1:20]  ## primeros 20 picos 
top.b30 <- br[1:30]  ## primeros 30 picos
top_actual <- top.b30

# Elección de mejores algoritmos de clustering
comparacion <- clValid(
  obj        = matint_122_dico_d2d4[, top_actual],
  nClust     = 2:6,
  clMethods  = c("hierarchical", "kmeans", "pam"),
  validation = c("stability", "internal")
)

summary(comparacion)
optimalScores(comparacion) # Se puede ir probando con distintos top picos


### ALGORITMOS DE CLUSTERING: Con top30 HKMEANS clustering 2 clusters ##########
################################################################################

top_actual <- top.b30
K.num <- 2 # clusters
var2 = 0.95

hkmeans.top30.k2 <- hkmeans(matint_122_dico_d2d4[, top_actual], 
                      K.num)

cluster.hkmeans.top30.k2 <- fviz_cluster(
                                  hkmeans.top30.k2, 
                                  ellipse.type = "convex", 
                                  data = matint_122_dico_d2d4[, top_actual],
                                  ellipse.level = var2,
                                  show.clust.cent = F, 
                                  geom = "point", 
                                  main = "hkmeans - Top 30 - 2 clusters")

cluster.hkmeans.top30.k2 <- cluster.hkmeans.top30.k2 + 
  geom_point(data = cluster.hkmeans.top30.k2$data, 
             aes(x = x, y = y, color = df_prom_rep_d2d4$factor1, 
                 size = df_prom_rep_d2d4$dia)) +
  scale_color_manual(values = c("blueviolet","aquamarine3","maroon1", "maroon4", 
                                               "blue1","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +  
  theme(legend.position = "right")

print(cluster.hkmeans.top30.k2)

### TASA DE ACIERTO POR CLUSTER #################################################

print("Resultados de hkmeans.top30.k2: ")

for(i in 1:length(hkmeans.top30.k2[[1]])) {
  hkmeans.top30.k2$names[i] = sub("_.*", "", attr(hkmeans.top30.k2[[1]][i],"names"))
  hkmeans.top30.k2$dia[i] = df_metadata_prom_rep$dia[i]
}

# Calcular la tasa de acierto para todo el conjunto de datos
total_CLP <- sum(hkmeans.top30.k2$names == "CLP")  # Total de muestras CLP
total_SHAM <- sum(hkmeans.top30.k2$names == "SH")  # Total de muestras SHAM

# Muestras CLP Y SHAM clasificadas correctamente
acierto_CLP <- sum(hkmeans.top30.k2$names == "CLP" & hkmeans.top30.k2$cluster == 1)
acierto_SHAM <- sum(hkmeans.top30.k2$names == "SH" & hkmeans.top30.k2$cluster == 2)

total_cluster2 <- sum(hkmeans.top30.k2$cluster == 2)
total_cluster1 <- sum(hkmeans.top30.k2$cluster == 1)

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

# Gráfico general
grafico_general <- ggplot(datos_grafico_general, 
                          aes(x = Categoria, y = Tasa_Acierto, 
                              fill = Categoria)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
  labs(title = "Tasa de Acierto Total",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CLP" = "red4", "SH" = "blue4")) +
  theme_minimal()

plot(grafico_general)


### OBTENCIÓN DE MÉTRICAS  #####################################################

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

# 2. Silhouette Score
silhouette_score <- silhouette(hkmeans.top30.k2$cluster, dist(matint_122_dico_d2d4[, top_actual]))
fviz_silhouette(silhouette_score)

# 3. Within-cluster Sum of Squares (WCSS)
wcss <- hkmeans.top30.k2$tot.withinss
print(wcss)

# 4. Between-cluster Sum of Squares (BCSS)
bcss <- hkmeans.top30.k2$betweenss
print(bcss)

### ALGORITMOS DE CLUSTERING: Con top20 PAM clustering 3 clusters ##############
################################################################################

top_actual <- top.b20
K.num <- 3 # clusters
var2 = 0.95

pam.top20.k3 <- pam(x= matint_122_dico_d2d4[, top_actual], K.num,
                    metric = "manhattan")

cluster.pam.top20.k3 <- fviz_cluster(pam.top20.k3, ellipse.type = "convex",
                                data = matint_122_dico_d2d4[, top_actual],
                                ellipse.level = var2,
                                show.clust.cent = F,
                                geom = "point",
                                main = "PAM - Top 20 - 3 cluster")


cluster.pam.top20.k3 <- cluster.pam.top20.k3 +
  geom_point(data = cluster.pam.top20.k3$data,
             aes(x = x, y = y, color = df_prom_rep_d2d4$factor1,
                 size = df_prom_rep_d2d4$dia)) +
  scale_color_manual(values = c("blueviolet","aquamarine3","blue4","maroon1",
                                "maroon4", "blue1","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.pam.top20.k3)

### OBTENCIÓN DE MÉTRICAS  #####################################################

# Solo se calcula el parámetro silhouette en este análisis en particular

print("Resultados de pam.top20.k3: ")

# 1. Silhouette Score
silhouette_score <- silhouette(pam.top20.k3$cluster, dist(matint_122_dico_d2d4[, top_actual]))
fviz_silhouette(silhouette_score)

### ALGORITMOS DE CLUSTERING: Con top20 PAM clustering 2 clusters ##############
################################################################################

# Con top20: PAM clustering 2 clusters
top_actual <- top.b20
K.num <- 2 # clusters
var2 = 0.95

pam.top20.k2 <- pam(x= matint_122_dico_d2d4[, top_actual], K.num,
                    metric = "manhattan")

cluster.pam.top20.k2 <- fviz_cluster(pam.top20.k2, ellipse.type = "convex",
                              data = matint_122_dico_d2d4[, top_actual],
                              ellipse.level = var2,
                              show.clust.cent = F,
                              geom = "point", main = "PAM - Top 20 - 2 cluster")

cluster.pam.top20.k2 <- cluster.pam.top20.k2 +
  geom_point(data = cluster.pam.top20.k2$data,
             aes(x = x, y = y, color = df_prom_rep_d2d4$factor1,
                 size = df_prom_rep_d2d4$dia)) +
  scale_color_manual(values = c("blueviolet","aquamarine3","maroon1", "maroon4",
                                "blue1","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.pam.top20.k2)

### TASA DE ACIERTO POR CLUSTER #################################################

print("Resultados de pam.top20.k2: ")

# Obtener los nombres de las muestras
sample_names <- names(pam.top20.k2$clustering)

# Calcular la tasa de acierto para las muestras que contienen "CLP"
total_CLP <- sum(grepl("CLP", sample_names))

# Calcular la tasa de acierto para las muestras que contienen "SH"
total_SHAM <- sum(grepl("SH", sample_names))

# Muestras CLP_D2, CLP_D4 y SHAM clasificadas correctamente
acierto_CLP <- sum(grepl("CLP", sample_names) & pam.top20.k2$clustering == 2)
acierto_SHAM <- sum(grepl("SH", sample_names) & pam.top20.k2$clustering == 1)

# Totales por cluster
total_cluster1 <- sum(pam.top20.k2$clustering == 1)
total_cluster2 <- sum(pam.top20.k2$clustering == 2)

# Cálculo tasa de acierto por cluster
tasa_acierto_CLP <- acierto_CLP / total_cluster1
tasa_acierto_SHAM <- acierto_SHAM / total_cluster2

# Crear un dataframe con los resultados para el gráfico general
datos_grafico_general <- data.frame(
  Categoria = c("CLP", "SH"),
  Aciertos = c(acierto_CLP, acierto_SHAM),
  Total = c(total_cluster1, total_cluster2),
  Tasa_Acierto = c(tasa_acierto_CLP, tasa_acierto_SHAM),
  Dia = "Total"
)

# Gráfico general
grafico_general <- ggplot(datos_grafico_general, 
                          aes(x = Categoria, y = Tasa_Acierto, 
                              fill = Categoria)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
  labs(title = "Tasa de Acierto Total",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CLP" = "red4", "SH" = "blue4")) +
  theme_minimal()

plot(grafico_general)


### OBTENCIÓN DE MÉTRICAS ######################################################

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

# 2. Silhouette Score
silhouette_score <- silhouette(pam.top20.k2$cluster, dist(matint_122_dico_d2d4[, top_actual]))
fviz_silhouette(silhouette_score)

# 3. Within-cluster Sum of Squares (WCSS)
# Centroides de los clusters
centroides <- pam.top20.k2$medoids

# Inicializar WCSS
wcss <- 0

# Calcular WCSS sumando la distancia de cada punto a su centroide
for (i in 1:K.num) {
  cluster_data <- matint_122_dico_d2d4[pam.top20.k2$clustering == i, top_actual, drop = FALSE]
  distances <- apply(cluster_data, 1, function(x) sum((x - centroides[i, ])^2))
  wcss <- wcss + sum(distances)
}

cat("Within-cluster Sum of Squares (WCSS):", wcss, "\n")

# 4. Between-cluster Sum of Squares (BCSS)
# Calcular el centroide global
centroide_global <- colMeans(matint_122_dico_d2d4[, top_actual])

# Inicializar BCSS
bcss <- 0

# Calcular BCSS sumando las distancias entre centroides y el centroide global, ponderadas por el número de puntos en cada cluster
for (i in 1:K.num) {
  n_i <- sum(pam.top20.k2$clustering == i)
  centroid <- centroides[i, ]
  bcss <- bcss + n_i * sum((centroid - centroide_global)^2)
}

cat("Between-cluster Sum of Squares (BCSS):", bcss, "\n")


### ALGORITMOS DE CLUSTERING: Con top15 KMEAN clustering 2 clusters ############
################################################################################
 
top_actual <- top.b15
K.num <- 2 # clusters
var2 = 0.95

kmeans.top15.k2 <- kmeans(matint_122_dico_d2d4[, top_actual],
                    K.num, nstart = 25)

cluster.kmeans.top15.k2 <- fviz_cluster(kmeans.top15.k2,
                                        ellipse.type = "convex",
                                data = matint_122_dico_d2d4[, top_actual],
                                ellipse.level = var2,
                                show.clust.cent = F,
                                geom = "point",
                                main = "KMEANS - Top 15 - 2 cluster")

cluster.kmeans.top15.k2 <- cluster.kmeans.top15.k2 +
  geom_point(data = cluster.kmeans.top15.k2$data,
             aes(x = x, y = y, color = df_prom_rep_d2d4$factor1,
                 size = df_prom_rep_d2d4$dia)) +
  scale_color_manual(values = c("blueviolet","aquamarine3","maroon1", "maroon4",
                                "blue1","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.kmeans.top15.k2)

### TASA DE ACIERTO POR CLUSTER #################################################

print("Resultados de kmeans.top15.k2: ")

for(i in 1:length(kmeans.top15.k2[[1]])) {
  kmeans.top15.k2$names[i] = sub("_.*", "", attr(kmeans.top15.k2[[1]][i],"names"))
  kmeans.top15.k2$dia[i] = df_metadata_prom_rep$dia[i]
}

# Calcular la tasa de acierto para todo el conjunto de datos
total_CLP <- sum(kmeans.top15.k2$names == "CLP")  # Total de muestras CLP
total_SHAM <- sum(kmeans.top15.k2$names == "SH")  # Total de muestras SHAM

# Muestras CLP Y SHAM clasificadas correctamente
acierto_CLP <- sum(kmeans.top15.k2$names == "CLP" & kmeans.top15.k2$cluster == 1)
acierto_SHAM <- sum(kmeans.top15.k2$names == "SH" & kmeans.top15.k2$cluster == 2)

total_cluster2 <- sum(kmeans.top15.k2$cluster == 2)
total_cluster1 <- sum(kmeans.top15.k2$cluster == 1)

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

# Gráfico general
grafico_general <- ggplot(datos_grafico_general, 
                          aes(x = Categoria, y = Tasa_Acierto, 
                              fill = Categoria)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
  labs(title = "Tasa de Acierto Total",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CLP" = "red4", "SH" = "blue4")) +
  theme_minimal()

plot(grafico_general)


### OBTENCIÓN DE MÉTRICAS  #####################################################

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

# 2. Silhouette Score
silhouette_score <- silhouette(kmeans.top15.k2$cluster, dist(matint_122_dico_d2d4[, top_actual]))
fviz_silhouette(silhouette_score)

# 3. Within-cluster Sum of Squares (WCSS)
wcss <- kmeans.top15.k2$tot.withinss
print(wcss)

# 4. Between-cluster Sum of Squares (BCSS)
bcss <- kmeans.top15.k2$betweenss
print(bcss)


### ALGORITMOS DE CLUSTERING: Con top10 KMEAN clustering 2 clusters ############
################################################################################


top_actual <- top.b10
K.num <- 2 # clusters
var2 = 0.95

kmeans.top10.k2 <- kmeans(matint_122_dico_d2d4[, top_actual],
                    K.num, nstart = 25)

cluster.kmeans.top10.k2 <- fviz_cluster(kmeans.top10.k2,
                                        ellipse.type = "convex",
                                data = matint_122_dico_d2d4[, top_actual],
                                ellipse.level = var2,
                                show.clust.cent = F,
                                geom = "point",
                                main = "KMEANS - Top 10 - 2 cluster")

cluster.kmeans.top10.k2 <- cluster.kmeans.top10.k2 +
  geom_point(data = cluster.kmeans.top10.k2$data,
             aes(x = x, y = y, color = df_prom_rep_d2d4$factor1,
                 size = df_prom_rep_d2d4$dia)) +
  scale_color_manual(values = c("blueviolet","aquamarine3","maroon1", "maroon4",
                                "blue1","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.kmeans.top10.k2)

### TASA DE ACIERTO POR CLUSTER #################################################

print("Resultados de kmeans.top10.k2: ")

for(i in 1:length(kmeans.top10.k2[[1]])) {
  kmeans.top10.k2$names[i] = sub("_.*", "", attr(kmeans.top10.k2[[1]][i],"names"))
  kmeans.top10.k2$dia[i] = df_metadata_prom_rep$dia[i]
}

# Calcular la tasa de acierto para todo el conjunto de datos
total_CLP <- sum(kmeans.top10.k2$names == "CLP")  # Total de muestras CLP
total_SHAM <- sum(kmeans.top10.k2$names == "SH")  # Total de muestras SHAM

# Muestras CLP Y SHAM clasificadas correctamente
acierto_CLP <- sum(kmeans.top10.k2$names == "CLP" & kmeans.top10.k2$cluster == 1)
acierto_SHAM <- sum(kmeans.top10.k2$names == "SH" & kmeans.top10.k2$cluster == 2)

total_cluster2 <- sum(kmeans.top10.k2$cluster == 2)
total_cluster1 <- sum(kmeans.top10.k2$cluster == 1)

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

# Gráfico general
grafico_general <- ggplot(datos_grafico_general, 
                          aes(x = Categoria, y = Tasa_Acierto, 
                              fill = Categoria)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
  labs(title = "Tasa de Acierto Total",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CLP" = "red4", "SH" = "blue4")) +
  theme_minimal()

plot(grafico_general)


### OBTENCIÓN DE MÉTRICAS  #####################################################

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

# 2. Silhouette Score
silhouette_score <- silhouette(kmeans.top10.k2$cluster, dist(matint_122_dico_d2d4[, top_actual]))
fviz_silhouette(silhouette_score)

# 3. Within-cluster Sum of Squares (WCSS)
wcss <- kmeans.top10.k2$tot.withinss
print(wcss)

# 4. Between-cluster Sum of Squares (BCSS)
bcss <- kmeans.top10.k2$betweenss
print(bcss)

### ALGORITMOS DE CLUSTERING: Con top10 PAM clustering 2 clusters ##############
################################################################################

# Con top10: PAM CLUSTERING 2 clusters
top_actual <- top.b10
K.num <- 2 # clusters
var2 = 0.95

pam.top10.k2 <- pam(x = matint_122_dico_d2d4[, top_actual], K.num,
                    metric = "manhattan")

cluster.pam.top10.k2 <- fviz_cluster(pam.top10.k2, ellipse.type = "convex",
                                     data = matint_122_dico_d2d4[, top_actual],
                                     ellipse.level = var2,
                                     show.clust.cent = F,
                                     geom = "point",
                                     main = "PAM - Top 10 - 2 cluster")

cluster.pam.top10.k2 <- cluster.pam.top10.k2 +
  geom_point(data = cluster.pam.top10.k2$data,
             aes(x = x, y = y, color = df_prom_rep_d2d4$factor1,
                 size = df_prom_rep_d2d4$dia)) +
  scale_color_manual(values = c("blueviolet","aquamarine3","maroon1", "maroon4",
                                "blue1","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.pam.top10.k2)

### TASA DE ACIERTO POR CLUSTER #################################################

print("Resultados de pam.top10.k2: ")

for(i in 1:length(pam.top10.k2[[1]])) {
  pam.top10.k2$names[i] = sub("_.*", "", attr(pam.top10.k2[[1]][i],"names"))
  pam.top10.k2$dia[i] = df_metadata_prom_rep$dia[i]
}

# Calcular la tasa de acierto para todo el conjunto de datos
total_CLP <- sum(pam.top10.k2$names == "CLP")  # Total de muestras CLP
total_SHAM <- sum(pam.top10.k2$names == "SH")  # Total de muestras SHAM

# Muestras CLP Y SHAM clasificadas correctamente
acierto_CLP <- sum(pam.top10.k2$names == "CLP" & pam.top10.k2$cluster == 1)
acierto_SHAM <- sum(pam.top10.k2$names == "SH" & pam.top10.k2$cluster == 2)

total_cluster2 <- sum(pam.top10.k2$cluster == 2)
total_cluster1 <- sum(pam.top10.k2$cluster == 1)

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

# Gráfico general
grafico_general <- ggplot(datos_grafico_general, 
                          aes(x = Categoria, y = Tasa_Acierto, 
                              fill = Categoria)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
  labs(title = "Tasa de Acierto Total",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CLP" = "red4", "SH" = "blue4")) +
  theme_minimal()

plot(grafico_general)

### OBTENCIÓN DE MÉTRICAS #####################################################

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

# 2. Silhouette Score
silhouette_score <- silhouette(pam.top10.k2$cluster, dist(matint_122_dico_d2d4[, top_actual]))
fviz_silhouette(silhouette_score)

# 3. Within-cluster Sum of Squares (WCSS)
# Centroides de los clusters
centroides <- pam.top10.k2$medoids

# Inicializar WCSS
wcss <- 0

# Calcular WCSS sumando la distancia de cada punto a su centroide
for (i in 1:K.num) {
  cluster_data <- matint_122_dico_d2d4[pam.top10.k2$clustering == i, top_actual, drop = FALSE]
  distances <- apply(cluster_data, 1, function(x) sum((x - centroides[i, ])^2))
  wcss <- wcss + sum(distances)
}

cat("Within-cluster Sum of Squares (WCSS):", wcss, "\n")

# 4. Between-cluster Sum of Squares (BCSS)
# Calcular el centroide global
centroide_global <- colMeans(matint_122_dico_d2d4[, top_actual])

# Inicializar BCSS
bcss <- 0

# Calcular BCSS sumando las distancias entre centroides y el centroide global, ponderadas por el número de puntos en cada cluster
for (i in 1:K.num) {
  n_i <- sum(pam.top10.k2$clustering == i)
  centroid <- centroides[i, ]
  bcss <- bcss + n_i * sum((centroid - centroide_global)^2)
}

cat("Between-cluster Sum of Squares (BCSS):", bcss, "\n")



#### PRUEBA 3 GRUPOS: SH(D2 Y D4) VS CLP_D4 VS CLP_D2 ##########################
################################################################################

# Agrupo factor SH
df_prom_rep_d2d4 <- df_prom_rep_d2d4 %>%
  mutate(factor1 = ifelse(grepl("SH", factor1), "SH", factor1))

# Selección de picos para binary discriminant analysis (BDA)
factor_tipo <- factor(df_prom_rep_d2d4$factor1)
is.binaryMatrix(matint_122_dico_d2d4) # TRUE
br <- binda.ranking(matint_122_dico_d2d4, factor_tipo, verbose = FALSE)
indices_especificos <- br[1:20]

# Gráfico de picos vs score
nueva_columna <- c()
matriz <- matrix(br, nrow = 218, ncol = 4)
for (i in 1:218) {
  nuevo_valor <- colnames(matint_122_dico_d2d4)[br[i]]
  nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)

plot(df_br$nueva_columna, df_br$V2,
     xlab = "m/z", ylab = "Score",
     main = "Ranking de picos de los espectros, 3 clusters")

# Crear un gradiente de colores (por ejemplo, de azul a rojo)
colores <- colorRampPalette(c("green4", "red4"))(218)

# Agregar puntos con colores en forma de gradiente
for (i in 1:218) {
  points(df_br$nueva_columna[i], df_br$V2[i], col = colores[i])
}

# Agregar puntos con relleno de colores en forma de gradiente
for (i in 1:218) {
  points(df_br$nueva_columna[i], df_br$V2[i], pch = 19, col = colores[i])
}

# Selección de picos
top.b10 <- br[1:10]  ## primeros 10 picos
top.b15 <- br[1:15]  ## primeros 15 picos
top.b20 <- br[1:20]  ## primeros 20 picos
top.b30 <- br[1:30]  ## primeros 30 picos
top.b40 <- br[1:40]  ## primeros 30 picos
top_actual <- top.b10 # Ir probando

# Elección del mejor algoritmo de clustering
comparacion <- clValid(
  obj        = matint_122_dico_d2d4[, top_actual],
  nClust     = 2:6,
  clMethods  = c("hierarchical", "kmeans", "pam"),
  validation = c("stability", "internal")
)

summary(comparacion)
optimalScores(comparacion)


# Con top30: PAM CLUSTERING 3 clusters
top_actual <- top.b30
K.num <- 3 # clusters
var2 = 0.95

pam.top30.k4.g3 <- pam(matint_122_dico_d2d4[, top_actual], metric = "euclidean",
                          K.num)

cluster.pam.top30.k4.g3 <- fviz_cluster(pam.top30.k4.g3,
                                        ellipse.type = "convex",
                                        data = matint_122_dico_d2d4[, top_actual],
                                        ellipse.level = var2,
                                        show.clust.cent = F,
                                        geom = "point",
                                        main = "PAM - Top 30 - 3 cluster")

# Ajustar el tamaño de los puntos según los valores de la columna "dia"
cluster.pam.top30.k4.g3 <- cluster.pam.top30.k4.g3 +
  geom_point(data = cluster.pam.top30.k4.g3$data,
             aes(x = x, y = y, color = df_prom_rep_d2d4$factor1,
                 size = df_prom_rep_d2d4$dia)) +
  scale_color_manual(values = c("maroon1","aquamarine3","blueviolet","blue1",
                                "blue4", "maroon4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

# Muestra el gráfico
print(cluster.pam.top30.k4.g3)

### TASA DE ACIERTO POR CLUSTER #################################################

print("Resultados de pam.top30.k4.g3: ")

for(i in 1:length(pam.top30.k4.g3[[1]])) {
  pam.top30.k4.g3$names[i] = sub("_.*", "", attr(pam.top30.k4.g3[[1]][i],"names"))
  pam.top30.k4.g3$dia[i] = df_metadata_prom_rep$dia[i]
}

# Obtener los nombres de las muestras
sample_names <- names(pam.top30.k4.g3$clustering)

# Calcular la tasa de acierto para las muestras que contienen "CLP_D2"
total_CLP_D2 <- sum(grepl("CLP_D2", sample_names))

# Calcular la tasa de acierto para las muestras que contienen "CLP_D4"
total_CLP_D4 <- sum(grepl("CLP_D4", sample_names))

# Calcular la tasa de acierto para las muestras que contienen "SH"
total_SHAM <- sum(grepl("SH", sample_names))

# Muestras CLP_D2, CLP_D4 y SHAM clasificadas correctamente
acierto_CLP_D2 <- sum(grepl("CLP_D2", sample_names) & pam.top30.k4.g3$clustering == 1)
acierto_CLP_D4 <- sum(grepl("CLP_D4", sample_names) & pam.top30.k4.g3$clustering == 2)
acierto_SHAM <- sum(grepl("SH", sample_names) & pam.top30.k4.g3$clustering == 3)

# Totales por cluster
total_cluster1 <- sum(pam.top30.k4.g3$clustering == 1)
total_cluster2 <- sum(pam.top30.k4.g3$clustering == 2)
total_cluster3 <- sum(pam.top30.k4.g3$clustering == 3)

# Cálculo tasa de acierto por cluster
tasa_acierto_CLP_D2 <- acierto_CLP_D2 / total_cluster1
tasa_acierto_CLP_D4 <- acierto_CLP_D4 / total_cluster2
tasa_acierto_SHAM <- acierto_SHAM / total_cluster3

# Crear un dataframe con los resultados para el gráfico general
datos_grafico_general <- data.frame(
  Categoria = c("CLP_D2", "CLP_D4", "SH"),
  Aciertos = c(acierto_CLP_D2, acierto_CLP_D4, acierto_SHAM),
  Total = c(total_cluster1, total_cluster2, total_cluster3),
  Tasa_Acierto = c(tasa_acierto_CLP_D2, tasa_acierto_CLP_D4, tasa_acierto_SHAM),
  Dia = "Total"
)

# Gráfico general
grafico_general <- ggplot(datos_grafico_general, 
                          aes(x = Categoria, y = Tasa_Acierto, 
                              fill = Categoria)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
  labs(title = "Tasa de Acierto Total",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CLP_D2" = "red4", "CLP_D4" = "orange4", "SH" = "blue4")) +
  theme_minimal()

plot(grafico_general)

### OBTENCIÓN DE MÉTRICAS  #####################################################

# 1. Matriz de confusión

# Calcular los elementos de la matriz de confusión
TP_CLP_D2 <- acierto_CLP_D2
TP_CLP_D4 <- acierto_CLP_D4
TP_SHAM <- acierto_SHAM

FP_CLP_D2 <- total_cluster1 - acierto_CLP_D2
FP_CLP_D4 <- total_cluster2 - acierto_CLP_D4
FP_SHAM <- total_cluster3 - acierto_SHAM

FN_CLP_D2 <- total_CLP_D2 - acierto_CLP_D2
FN_CLP_D4 <- total_CLP_D4 - acierto_CLP_D4
FN_SHAM <- total_SHAM - acierto_SHAM

# Matriz de confusión (para CLP_D2 y CLP_D4 se consideran clusters por separado)
conf_matrix_CLP_D2 <- matrix(c(TP_CLP_D2, FP_CLP_D2, FN_CLP_D2, TP_CLP_D2), nrow = 2, byrow = TRUE,
                             dimnames = list('Referencia' = c('No_CLP_D2', 'CLP_D2'),
                                             'Predicho' = c('No_cluster1', 'cluster1')))
conf_matrix_CLP_D4 <- matrix(c(TP_CLP_D4, FP_CLP_D4, FN_CLP_D4, TP_CLP_D4), nrow = 2, byrow = TRUE,
                             dimnames = list('Referencia' = c('No_CLP_D4', 'CLP_D4'),
                                             'Predicho' = c('No_cluster2', 'cluster2')))
conf_matrix_SHAM <- matrix(c(TP_SHAM, FP_SHAM, FN_SHAM, TP_SHAM), nrow = 2, byrow = TRUE,
                           dimnames = list('Referencia' = c('No_SHAM', 'SHAM'),
                                           'Predicho' = c('No_cluster3', 'cluster3')))

# Mostrar matrices de confusión
print("Matriz de Confusión para CLP_D2:")
print(conf_matrix_CLP_D2)
print("Matriz de Confusión para CLP_D4:")
print(conf_matrix_CLP_D4)
print("Matriz de Confusión para SHAM:")
print(conf_matrix_SHAM)

# Cálculo de métricas para cada grupo

# Precision, recall, f1-score y accuracy para CLP_D2
precision_CLP_D2 <- TP_CLP_D2 / (TP_CLP_D2 + FP_CLP_D2)
recall_CLP_D2 <- TP_CLP_D2 / (TP_CLP_D2 + FN_CLP_D2)
f1_score_CLP_D2 <- 2 * (precision_CLP_D2 * recall_CLP_D2) / (precision_CLP_D2 + recall_CLP_D2)
accuracy_CLP_D2 <- (TP_CLP_D2 + (total_SHAM - FP_SHAM - FN_SHAM)) / (total_cluster1 + total_SHAM)

# Precision, recall, f1-score y accuracy para CLP_D4
precision_CLP_D4 <- TP_CLP_D4 / (TP_CLP_D4 + FP_CLP_D4)
recall_CLP_D4 <- TP_CLP_D4 / (TP_CLP_D4 + FN_CLP_D4)
f1_score_CLP_D4 <- 2 * (precision_CLP_D4 * recall_CLP_D4) / (precision_CLP_D4 + recall_CLP_D4)
accuracy_CLP_D4 <- (TP_CLP_D4 + (total_SHAM - FP_SHAM - FN_SHAM)) / (total_cluster2 + total_SHAM)

# Precision, recall, f1-score y accuracy para SHAM
precision_SHAM <- TP_SHAM / (TP_SHAM + FP_SHAM)
recall_SHAM <- TP_SHAM / (TP_SHAM + FN_SHAM)
f1_score_SHAM <- 2 * (precision_SHAM * recall_SHAM) / (precision_SHAM + recall_SHAM)
accuracy_SHAM <- (TP_SHAM + (total_CLP_D2 + total_CLP_D4 - FP_CLP_D2 - FN_CLP_D2 - FP_CLP_D4 - FN_CLP_D4)) / (total_cluster3 + total_CLP_D2 + total_CLP_D4)

# Mostrar las métricas
cat("Métricas para CLP_D2:\n")
cat("Precision:", precision_CLP_D2, "\n")
cat("Recall:", recall_CLP_D2, "\n")
cat("F1-Score:", f1_score_CLP_D2, "\n")
cat("Accuracy:", accuracy_CLP_D2, "\n\n")

cat("Métricas para CLP_D4:\n")
cat("Precision:", precision_CLP_D4, "\n")
cat("Recall:", recall_CLP_D4, "\n")
cat("F1-Score:", f1_score_CLP_D4, "\n")
cat("Accuracy:", accuracy_CLP_D4, "\n\n")

cat("Métricas para SHAM:\n")
cat("Precision:", precision_SHAM, "\n")
cat("Recall:", recall_SHAM, "\n")
cat("F1-Score:", f1_score_SHAM, "\n")
cat("Accuracy:", accuracy_SHAM, "\n\n")

# 2. Silhouette Score
silhouette_score <- silhouette(pam.top30.k4.g3$clustering, dist(matint_122_dico_d2d4[, top_actual]))
fviz_silhouette(silhouette_score)

# 3. Within-cluster Sum of Squares (WCSS)
# Centroides de los clusters
centroides <- pam.top30.k4.g3$medoids

# Inicializar WCSS
wcss <- 0

# Calcular WCSS sumando la distancia de cada punto a su centroide
for (i in 1:K.num) {
  cluster_data <- matint_122_dico_d2d4[pam.top30.k4.g3$clustering == i, top_actual, drop = FALSE]
  distances <- apply(cluster_data, 1, function(x) sum((x - centroides[i, ])^2))
  wcss <- wcss + sum(distances)
}

cat("Within-cluster Sum of Squares (WCSS):", wcss, "\n")

# 4. Between-cluster Sum of Squares (BCSS)
# Calcular el centroide global
centroide_global <- colMeans(matint_122_dico_d2d4[, top_actual])

# Inicializar BCSS
bcss <- 0

# Calcular BCSS sumando las distancias entre centroides y el centroide global, ponderadas por el número de puntos en cada cluster
for (i in 1:K.num) {
  n_i <- sum(pam.top30.k4.g3$clustering == i)
  centroid <- centroides[i, ]
  bcss <- bcss + n_i * sum((centroid - centroide_global)^2)
}

cat("Between-cluster Sum of Squares (BCSS):", bcss, "\n")

# ### EXPORTAR MATRICES CON 30 PICOS DE INTERÉS ##################################
# ################################################################################
# 
# 
# 
# write.csv(matint_122_dico_d2d4[, top_actual], "top30_107_dico_d24.csv",
#           row.names = TRUE)
# write.csv(matint_122_d2d4[, top_actual], "top30_107_d24.csv", row.names = TRUE)
# write.csv(df_prom_rep_d2d4, "metadata_107.csv", row.names = TRUE)
# #
# #
# #
# ### FIN ########################################################################
# ################################################################################
