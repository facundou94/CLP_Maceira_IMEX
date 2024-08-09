################ MALDI-TOF ANALISIS CLP ########################################
################ 5) ns_m51_d24  ################################################
#
# ns:  No Supervisado
# m51: Utiliza las muestras (51)
# d24: Utiliza los días 2 y 4
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
load(file.path(data_path, "matint_51.Rdata"))

df_metadata_unicas$dia <- as.integer(gsub("[^0-9]", "", 
                                          df_metadata_unicas$dia))

# Filtrado de muestras de día 2 y día 4 y agrupo los SH
filas_filtradas <- rownames(matint_51_dico)[grepl("D2|D4", rownames(matint_51_dico))]
matint_51_dico_d2d4 <- matint_51_dico[filas_filtradas, ]
df_unicas_d2d4 <- df_metadata_unicas %>%
  dplyr::filter(dia == 2 | dia == 4)
df_unicas_d2d4 <- df_unicas_d2d4 %>%
  mutate(factor1 = ifelse(grepl("SH", factor1), "SH", factor1))

filas_filtradas <- rownames(matint_51)[grepl("D2|D4", rownames(matint_51))]
matint_51_d2d4 <- matint_51[filas_filtradas, ]



### SELECCIÓN DE PICOS #########################################################
################################################################################


factor_tipo <- factor(df_unicas_d2d4$factor1)
is.binaryMatrix(matint_51_dico_d2d4) # TRUE
br <- binda.ranking(matint_51_dico_d2d4, factor_tipo, verbose = FALSE)
indices_especificos <- br[1:20]

# Gráfico de picos vs score 
nueva_columna <- c()
matriz <- matrix(br, nrow = 218, ncol = 4)
for (i in 1:218) {
  nuevo_valor <- colnames(matint_51_dico_d2d4)[br[i]]
  nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)

plot(df_br$nueva_columna, df_br$V2, 
     xlab = "m/z", ylab = "Score", 
     main = "Ranking de picos de los espectros CLP_D2 vs CLP_D4 vs SH")

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
top_actual <- top.b10 # Ir probando


# Elección de mejores algoritmos de clustering
comparacion <- clValid(
  obj        = matint_51_dico_d2d4[, top_actual],
  nClust     = 2:6,
  clMethods  = c("hierarchical", "kmeans", "pam"),
  validation = c("stability", "internal")
)

summary(comparacion)
optimalScores(comparacion)


### ALGORITMOS DE CLUSTERING ###################################################
################################################################################


# Con top15: HKMEANS clustering 3 clusters

top_actual <- top.b15
K.num <- 3 # clusters
var2 = 0.95

hkmeans.top15.k3 <- hkmeans(matint_51_dico_d2d4[, top_actual], 
                            K.num)

cluster.hkmeans.top15.k3 <- fviz_cluster(hkmeans.top15.k3, ellipse.type = "convex", 
                                         data = matint_51_dico_d2d4[, top_actual],
                                         ellipse.level = var2,
                                         show.clust.cent = F, 
                                         geom = "point", 
                                         main = "hkmeans - Top 15 - 3 clusters")

# Ajustar el tamaño de los puntos según los valores de la columna "dia"
cluster.hkmeans.top15.k3 <- cluster.hkmeans.top15.k3 + 
  geom_point(data = cluster.hkmeans.top15.k3$data, 
             aes(x = x, y = y, color = df_unicas_d2d4$factor1,
                 size = df_unicas_d2d4$dia)) +
  scale_color_manual(values = c("maroon1","blueviolet","aquamarine3", "blue1",
                                "blue4", "maroon4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.hkmeans.top15.k3)



print("Resultados de hkmeans.top15.k3: ")

# Obtener los nombres de las muestras
sample_names <- names(hkmeans.top15.k3$cluster)

# Calcular la tasa de acierto para las muestras que contienen "CLP_D2"
total_CLP_D2 <- sum(grepl("CLP_D2", sample_names))

# Calcular la tasa de acierto para las muestras que contienen "CLP_D4"
total_CLP_D4 <- sum(grepl("CLP_D4", sample_names))

# Calcular la tasa de acierto para las muestras que contienen "SH"
total_SHAM <- sum(grepl("SH", sample_names))

# Muestras CLP_D2, CLP_D4 y SHAM clasificadas correctamente
acierto_CLP_D2 <- sum(grepl("CLP_D2", sample_names) & hkmeans.top15.k3$cluster == 1)
acierto_CLP_D4 <- sum(grepl("CLP_D4", sample_names) & hkmeans.top15.k3$cluster == 3)
acierto_SHAM <- sum(grepl("SH", sample_names) & hkmeans.top15.k3$cluster == 2)

# Totales por cluster
total_cluster1 <- sum(hkmeans.top15.k3$cluster == 1)
total_cluster2 <- sum(hkmeans.top15.k3$cluster == 2)
total_cluster3 <- sum(hkmeans.top15.k3$cluster == 3)

# Cálculo tasa de acierto por cluster
tasa_acierto_CLP_D2 <- acierto_CLP_D2 / total_CLP_D2
tasa_acierto_CLP_D4 <- acierto_CLP_D4 / total_CLP_D4
tasa_acierto_SHAM <- acierto_SHAM / total_SHAM

# Crear un dataframe con los resultados para el gráfico general
datos_grafico_general <- data.frame(
  Categoria = c("CLP_D2", "CLP_D4", "SHAM"),
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
  labs(title = "Tasa de Acierto top15: HKMEANS clustering 3 clusters",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CLP_D2" = "red4", "CLP_D4" = "orange4", "SHAM" = "blue4")) +
  theme_minimal()

plot(grafico_general)


# 1. Silhouette Score
silhouette_score <- silhouette(hkmeans.top15.k3$cluster, dist(matint_51_dico_d2d4[, top_actual]))
fviz_silhouette(silhouette_score)


# NOTA: Este mismo esquema con algoritmo PAM da resultados similares

# # Con top15: PAM clustering 3 clusters
# 
# top_actual <- top.b15
# K.num <- 3 # clusters
# var2 = 0.95
# 
# pam.top15.k3 <- pam(matint_51_dico_d2d4[, top_actual], metric = "manhattan",
#                             K.num)
# 
# cluster.pam.top15.k3 <- fviz_cluster(pam.top15.k3, ellipse.type = "convex", 
#                                      data = matint_51_dico_d2d4[, top_actual],
#                                      ellipse.level = var2,
#                                      show.clust.cent = F, 
#                                      geom = "point",
#                                      main = "PAM - Top 15 - 3 clusters")
# 
# # Ajustar el tamaño de los puntos según los valores de la columna "dia"
# cluster.pam.top15.k3 <- cluster.pam.top15.k3 + 
#   geom_point(data = cluster.pam.top15.k3$data, 
#              aes(x = x, y = y, color = df_unicas_d2d4$factor1,
#                  size = df_unicas_d2d4$dia)) +
#   scale_color_manual(values = c("maroon1","aquamarine3","blueviolet", "blue1",
#                                 "blue4", "maroon4")) +
#   scale_size_continuous(range = c(2, 3)) +
#   labs(color = "Cluster", size = "Día") +
#   theme(legend.position = "right")
# 
# print(cluster.pam.top15.k3)

#
#
#
### FIN ########################################################################
################################################################################