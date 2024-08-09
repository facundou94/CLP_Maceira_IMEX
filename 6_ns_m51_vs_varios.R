################ MALDI-TOF ANALISIS CLP ########################################
################ 6) ns_m51_vs_varios  ##########################################
#
# ns:    No Supervisado
# m51:  Utiliza las muestras (51)
# vs_varios: Utiliza los días 2 y 4 y realiza análisis vs entre CLP y SHAM
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
load(file.path(data_path, "matint_51_dico.Rdata"))
df_metadata_unicas$dia <- as.integer(gsub("[^0-9]", "", df_metadata_unicas$dia))


### PRIMER ANÁLISIS: CLP_D2 vs CLP_D4 ##########################################
################################################################################


filas_filtradas <- rownames(matint_51_dico)[grepl("CLP_D2|CLP_D4", rownames(matint_51_dico))]
matint_51_dico_clp_d2d4 <- matint_51_dico[filas_filtradas, ]
df_unicas_clp_d2d4 <- df_metadata_unicas %>%
  dplyr::filter(factor1 == "CLP_D2" | factor1 == "CLP_D4")

# Selección de picos para binary discriminant analysis (BDA)
factor_tipo <- factor(df_unicas_clp_d2d4$factor1)
is.binaryMatrix(matint_51_dico_clp_d2d4) # TRUE
br <- binda.ranking(matint_51_dico_clp_d2d4, factor_tipo, verbose = FALSE)
indices_especificos <- br[1:20]

# Gráfico de picos vs score
nueva_columna <- c()
matriz <- matrix(br, nrow = 218, ncol = 4)
for (i in 1:218) {
  nuevo_valor <- colnames(matint_51_dico_clp_d2d4)[br[i]]
  nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)

plot(df_br$nueva_columna, df_br$V2, 
     xlab = "m/z", ylab = "Score", 
     main = "Ranking de picos de los espectros CLP_D2 vs CLP_D4")

# Crear un gradiente de colores (por ejemplo, de azul a rojo)
colores <- colorRampPalette(c("blue", "red"))(218)

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
top_actual <- top.b30 # Ir probando

# Algoritmo PAM 2 clusters con TOP15
top_actual <- top.b15
K.num <- 2 # clusters
var2 = 0.95

pam.top15.k2 <- pam(matint_51_dico_clp_d2d4[, top_actual], 
                        metric = "manhattan",
                        K.num)

cluster.pam.top15.k2 <- fviz_cluster(pam.top15.k2, ellipse.type = "convex", 
                                         data = matint_51_dico_clp_d2d4[, top_actual],
                                         ellipse.level = var2,
                                         show.clust.cent = F, 
                                         geom = "point", 
                        main = "PAM - Top 15 - 2 clusters - CLP_D2 vs CLP_D4")

# Ajustar el tamaño de los puntos según los valores de la columna "dia"
cluster.pam.top15.k2 <- cluster.pam.top15.k2 + 
  geom_point(data = cluster.pam.top15.k2$data, 
             aes(x = x, y = y, color = df_unicas_clp_d2d4$factor1, 
                 size = df_unicas_clp_d2d4$dia)) +
  scale_color_manual(values = c("maroon1","aquamarine3","maroon4","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.pam.top15.k2)

### TASA DE ACIERTO POR CLUSTER #################################################

print("Resultados de pam.top20.k2: ")

# Obtener los nombres de las muestras
sample_names <- names(pam.top15.k2$clustering)

# Calcular la tasa de acierto para las muestras que contienen "CLP"
total_CLP_D2 <- sum(grepl("D2", sample_names))

# Calcular la tasa de acierto para las muestras que contienen "SH"
total_CLP_D4 <- sum(grepl("D4", sample_names))

# Muestras CLP_D2, CLP_D4 y SHAM clasificadas correctamente
acierto_CLP_D2 <- sum(grepl("D2", sample_names) & pam.top15.k2$clustering == 1)
acierto_CLP_D4 <- sum(grepl("D4", sample_names) & pam.top15.k2$clustering == 2)

# Totales por cluster
total_cluster1 <- sum(pam.top15.k2$clustering == 1)
total_cluster2 <- sum(pam.top15.k2$clustering == 2)

# Cálculo tasa de acierto por cluster
tasa_acierto_CLP_D2 <- acierto_CLP_D2 / total_CLP_D2
tasa_acierto_CLP_D4 <- acierto_CLP_D4 / total_CLP_D4

# Crear un dataframe con los resultados para el gráfico general
datos_grafico_general <- data.frame(
  Categoria = c("CLP_D2", "CLP_D4"),
  Aciertos = c(acierto_CLP_D2, acierto_CLP_D4),
  Total = c(total_cluster1, total_cluster2),
  Tasa_Acierto = c(tasa_acierto_CLP_D2, tasa_acierto_CLP_D4),
  Dia = "Total"
)

# Gráfico general
grafico_general <- ggplot(datos_grafico_general, 
                          aes(x = Categoria, y = Tasa_Acierto, 
                              fill = Categoria)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
  labs(title = "Tasa de Acierto top15 PAM 2 clusters",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CLP_D2" = "red4", "CLP_D4" = "blue4")) +
  theme_minimal()

plot(grafico_general)


### OBTENCIÓN DE MÉTRICAS ######################################################

# 1. Silhouette Score
silhouette_score <- silhouette(pam.top15.k2$cluster, dist(matint_51_dico_clp_d2d4[, top_actual]))
fviz_silhouette(silhouette_score)

# 2. Within-cluster Sum of Squares (WCSS)
# Centroides de los clusters
centroides <- pam.top15.k2$medoids

# Inicializar WCSS
wcss <- 0

# Calcular WCSS sumando la distancia de cada punto a su centroide
for (i in 1:K.num) {
  cluster_data <- matint_51_dico_clp_d2d4[pam.top15.k2$clustering == i, top_actual, drop = FALSE]
  distances <- apply(cluster_data, 1, function(x) sum((x - centroides[i, ])^2))
  wcss <- wcss + sum(distances)
}

cat("Within-cluster Sum of Squares (WCSS):", wcss, "\n")

# 3. Between-cluster Sum of Squares (BCSS)
# Calcular el centroide global
centroide_global <- colMeans(matint_51_dico_clp_d2d4[, top_actual])

# Inicializar BCSS
bcss <- 0

# Calcular BCSS sumando las distancias entre centroides y el centroide global, ponderadas por el número de puntos en cada cluster
for (i in 1:K.num) {
  n_i <- sum(pam.top15.k2$clustering == i)
  centroid <- centroides[i, ]
  bcss <- bcss + n_i * sum((centroid - centroide_global)^2)
}

cat("Between-cluster Sum of Squares (BCSS):", bcss, "\n")




### SEGUNDO ANÁLISIS: CLP_D2 vs SHAM_D2 ########################################
################################################################################


filas_filtradas <- rownames(matint_51_dico)[grepl("CLP_D2|SH_D2", rownames(matint_51_dico))]
matint_51_dico_clpd2_shamd2 <- matint_51_dico[filas_filtradas, ]
df_unicas_clpd2_shamd2 <- df_metadata_unicas %>%
  dplyr::filter(factor1 == "CLP_D2" | factor1 == "SH_D2")

# Selección de picos para binary discriminant analysis (BDA)
factor_tipo <- factor(df_unicas_clpd2_shamd2$factor1)
is.binaryMatrix(matint_51_dico_clpd2_shamd2) # TRUE
br <- binda.ranking(matint_51_dico_clpd2_shamd2, factor_tipo, verbose = FALSE)
indices_especificos <- br[1:20]

# Gráfico de picos vs score
nueva_columna <- c()
matriz <- matrix(br, nrow = 218, ncol = 4)
for (i in 1:218) {
  nuevo_valor <- colnames(matint_51_dico_clpd2_shamd2)[br[i]]
  nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)

plot(df_br$nueva_columna, df_br$V2, 
     xlab = "m/z", ylab = "Score", 
     main = "Ranking de picos de los espectros CLP_D2 vs SHAM_D2")

# Crear un gradiente de colores (por ejemplo, de azul a rojo)
colores <- colorRampPalette(c("blue", "red"))(218)

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
top_actual <- top.b30 # Ir probando

# Algoritmo PAM 2 clusters con TOP20
top_actual <- top.b20
K.num <- 2 # clusters
var2 = 0.95

pam.top20.k2 <- pam(matint_51_dico_clpd2_shamd2[, top_actual], 
                    metric = "manhattan",
                    K.num)

cluster.pam.top20.k2 <- fviz_cluster(pam.top20.k2, ellipse.type = "convex", 
                                     data = matint_51_dico_clpd2_shamd2[, top_actual],
                                     ellipse.level = var2,
                                     show.clust.cent = F, 
                                     geom = "point", 
                        main = "PAM - Top 20 - 2 clusters - CLP_D2 vs SHAM_D2 ")

# Ajustar el tamaño de los puntos según los valores de la columna "dia"
cluster.pam.top20.k2 <- cluster.pam.top20.k2 + 
  geom_point(data = cluster.pam.top20.k2$data, 
             aes(x = x, y = y, color = df_unicas_clpd2_shamd2$factor1, 
                 size = df_unicas_clpd2_shamd2$dia)) +
  scale_color_manual(values = c("maroon1","aquamarine3","maroon4","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.pam.top20.k2)


### TERCER ANÁLISIS: CLP_D4 vs SHAM_D4 #########################################
################################################################################


filas_filtradas <- rownames(matint_51_dico)[grepl("CLP_D4|SH_D4", 
                                                  rownames(matint_51_dico))]
matint_51_dico_clpd4_shamd4 <- matint_51_dico[filas_filtradas, ]
df_unicas_clpd4_shamd4 <- df_metadata_unicas %>%
  dplyr::filter(factor1 == "CLP_D4" | factor1 == "SH_D4")

# Selección de picos para binary discriminant analysis (BDA)
factor_tipo <- factor(df_unicas_clpd4_shamd4$factor1)
is.binaryMatrix(matint_51_dico_clpd4_shamd4) # TRUE
br <- binda.ranking(matint_51_dico_clpd4_shamd4, factor_tipo, verbose = FALSE)
indices_especificos <- br[1:20]

# Gráfico de picos vs score
nueva_columna <- c()
matriz <- matrix(br, nrow = 218, ncol = 4)
for (i in 1:218) {
  nuevo_valor <- colnames(matint_51_dico_clpd4_shamd4)[br[i]]
  nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)

plot(df_br$nueva_columna, df_br$V2, 
     xlab = "m/z", ylab = "Score", 
     main = "Ranking de picos de los espectros CLP_D4 vs SHAM_D4")

# Crear un gradiente de colores (por ejemplo, de azul a rojo)
colores <- colorRampPalette(c("blue", "red"))(218)

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
top_actual <- top.b30 # Ir probando

# Algoritmo PAM 2 clusters con TOP20
top_actual <- top.b20
K.num <- 2 # clusters
var2 = 0.95

pam.top15.k2 <- pam(matint_51_dico_clpd4_shamd4[, top_actual],
                        metric = "manhattan",
                        K.num)

cluster.pam.top15.k2 <- fviz_cluster(pam.top15.k2, 
                                         ellipse.type = "convex", 
                                         data = matint_51_dico_clpd4_shamd4[, top_actual],
                                         ellipse.level = var2,
                                         show.clust.cent = F, 
                                         geom = "point",
                            main = "PAM - Top 15 - 2 clusters - CLP_D4 vs SHAM_D4 ")

cluster.pam.top15.k2 <- cluster.pam.top15.k2 + 
  geom_point(data = cluster.pam.top15.k2$data, 
             aes(x = x, y = y, color = df_unicas_clpd4_shamd4$factor1,
                 size = df_unicas_clpd4_shamd4$dia)) +
  scale_color_manual(values = c("maroon1","aquamarine3","maroon4","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.pam.top15.k2)


### TASA DE ACIERTO POR CLUSTER #################################################

print("Resultados de pam.top15.k2: ")

# Obtener los nombres de las muestras
sample_names <- names(pam.top15.k2$clustering)

# Calcular la tasa de acierto para las muestras que contienen "CLP"
total_CLP <- sum(grepl("CLP", sample_names))

# Calcular la tasa de acierto para las muestras que contienen "SH"
total_SHAM <- sum(grepl("SH", sample_names))

# Muestras CLP_D2, CLP_D4 y SHAM clasificadas correctamente
acierto_CLP <- sum(grepl("CLP", sample_names) & pam.top15.k2$clustering == 1)
acierto_SHAM <- sum(grepl("SH", sample_names) & pam.top15.k2$clustering == 2)

# Totales por cluster
total_cluster1 <- sum(pam.top15.k2$clustering == 1)
total_cluster2 <- sum(pam.top15.k2$clustering == 2)

# Cálculo tasa de acierto por cluster
tasa_acierto_CLP <- acierto_CLP / total_CLP
tasa_acierto_SHAM <- acierto_SHAM / total_SHAM

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
  labs(title = "Tasa de Acierto top20 PAM 2 clusters",
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
silhouette_score <- silhouette(pam.top15.k2$cluster, dist(matint_51_dico_clpd4_shamd4[, top_actual]))
fviz_silhouette(silhouette_score)

# 3. Within-cluster Sum of Squares (WCSS)
# Centroides de los clusters
centroides <- pam.top15.k2$medoids

# Inicializar WCSS
wcss <- 0

# Calcular WCSS sumando la distancia de cada punto a su centroide
for (i in 1:K.num) {
  cluster_data <- matint_51_dico_clpd4_shamd4[pam.top15.k2$clustering == i, top_actual, drop = FALSE]
  distances <- apply(cluster_data, 1, function(x) sum((x - centroides[i, ])^2))
  wcss <- wcss + sum(distances)
}

cat("Within-cluster Sum of Squares (WCSS):", wcss, "\n")

# 4. Between-cluster Sum of Squares (BCSS)
# Calcular el centroide global
centroide_global <- colMeans(matint_51_dico_clpd4_shamd4[, top_actual])

# Inicializar BCSS
bcss <- 0

# Calcular BCSS sumando las distancias entre centroides y el centroide global, ponderadas por el número de puntos en cada cluster
for (i in 1:K.num) {
  n_i <- sum(pam.top15.k2$clustering == i)
  centroid <- centroides[i, ]
  bcss <- bcss + n_i * sum((centroid - centroide_global)^2)
}

cat("Between-cluster Sum of Squares (BCSS):", bcss, "\n")

#
#
#
### FIN ########################################################################
################################################################################