################ MALDI-TOF ANALISIS CLP #####################################################################################
################ 2) NO SUPERVISADO    #######################################################################################
##
## Autor: Bioing. Facundo Urteaga (IBB-CONICET)
##
##
##################### CARGA DE LIBRERIAS ####################################################################################
#
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
#
#
#################### CARGA DE ARCHIVOS ######################################################################################
#
load("matint_122_dico.Rdata")
load("matint_51_dico.Rdata")
load("matint_122.Rdata")
load("matint_51.Rdata")
#
df_metadata_prom_rep$dia <- as.integer(gsub("[^0-9]", "", df_metadata_prom_rep$dia))

### Analisis no supervisado
#
### Selección de picos para binary discriminant analysis (BDA)
#
#
factor_tipo <- factor(df_metadata_prom_rep$tipo)
is.binaryMatrix(matint_122_dico) # TRUE
br <- binda.ranking(matint_122_dico, factor_tipo, verbose = FALSE)
indices_especificos <- br[1:20]
#
#
# GRAFICO DISPERSIÓN
#
nueva_columna <- c()
matriz <- matrix(br, nrow = 218, ncol = 4)
for (i in 1:218) {
  nuevo_valor <- colnames(matint_122_dico)[br[i]]
  nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)

plot(df_br$nueva_columna, df_br$V2, 
     xlab = "m/z", ylab = "Score", 
     main = "Ranking de picos de los espectros")

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
#
#
#
# SELECCIÓN DE PICOS MAS PREPONDERANTES
#
top.b5 <- br[1:5]  ## primeros 5 picos
top.b10 <- br[1:10]  ## primeros 10 picos
top.b15 <- br[1:15]  ## primeros 15 picos
top.b20 <- br[1:20]  ## primeros 20 picos 
#
top_actual <- top.b20
#
# ELECCIÓN MEJOR ALGORITMO DE CLUSTERING
#
#
comparacion <- clValid(
  obj        = matint_122_dico[, top_actual],
  nClust     = 2:6,
  clMethods  = c("hierarchical", "kmeans", "pam"),
  validation = c("stability", "internal")
)
#
summary(comparacion)
#
optimalScores(comparacion)
#
# Obtengo k means 2 grupos
#
#
### k-means clustering

K.num <- 2 # clusters
var2 = 0.95

km.res20 <- kmeans(matint_122_dico[, top_actual], 
                   K.num, 
                   nstart = 25)

cluster_kmean20 <- fviz_cluster(km.res20, ellipse.type = "convex", 
                                data = matint_122_dico[, top_actual],
                                ellipse.level = var2,
                                show.clust.cent = F, 
                                geom = "point", main = "Top 20 - 2 cluster")

# Ajustar el tamaño de los puntos según los valores de la columna "dia"
cluster_kmean20 <- cluster_kmean20 + 
  geom_point(data = cluster_kmean20$data, 
             aes(x = x, y = y, color = factor_tipo, size = df_metadata_prom_rep$dia)) +
  scale_color_manual(values = c("red", "blue", "steelblue4", "maroon","green")) +
  scale_size_continuous(range = c(1, 6)) +  # Ajusta el rango de tamaño de los puntos
  labs(color = "Cluster", size = "Día") +  # Etiquetas de las leyendas
  theme(legend.position = "right")  # Posición de las leyendas

# Muestra el gráfico
print(cluster_kmean20)

for(i in 1:length(km.res20[[1]])) {
  km.res20$names[i] = sub("_.*", "", attr(km.res20[[1]][i],"names"))
  km.res20$dia[i] = df_metadata_prom_rep$dia[i]
}


#
#
#

# Calcular la tasa de acierto para todo el conjunto de datos
total_CLP <- sum(km.res20$names == "CLP")  # Total de muestras CLP
total_SHAM <- sum(km.res20$names == "SH")  # Total de muestras SHAM

acierto_CLP <- sum(km.res20$names == "CLP" & km.res20$cluster == 2)  # Muestras CLP clasificadas correctamente
acierto_SHAM <- sum(km.res20$names == "SH" & km.res20$cluster == 1)  # Muestras SHAM clasificadas correctamente

total_cluster2 <- sum(km.res20$cluster == 2)
total_cluster1 <- sum(km.res20$cluster == 1)

tasa_acierto_CLP <- acierto_CLP / total_cluster2  # Calcular tasa de acierto total
tasa_acierto_SHAM <- acierto_SHAM / total_cluster1 

# Crear un dataframe con los resultados para el gráfico general
datos_grafico_general <- data.frame(Categoria = c("CLP", "SH"),
                                    Aciertos = c(acierto_CLP, acierto_SHAM),
                                    Total = c(total_cluster2, total_cluster1),
                                    Tasa_Acierto = c(tasa_acierto_CLP, tasa_acierto_SHAM),
                                    Dia = "Total")

# Calcular la tasa de acierto por día
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

# Graficar
library(ggplot2)
library(gridExtra)

# Gráfico general
grafico_general <- ggplot(datos_grafico_general, aes(x = Categoria, y = Tasa_Acierto, fill = Categoria)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
  labs(title = "Tasa de Acierto Total",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CLP" = "blue", "SH" = "red")) +
  theme_minimal()

# Gráficos por día
graficos_dias <- lapply(unique(datos_grafico_dias$Dia), function(d) {
  datos_dia <- subset(datos_grafico_dias, Dia == d)
  ggplot(datos_dia, aes(x = Categoria, y = Tasa_Acierto, fill = Categoria)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
    labs(title = paste("Tasa de Acierto - Día", d),
         x = "Categoría",
         y = "Tasa de Acierto") +
    scale_fill_manual(values = c("CLP" = "blue", "SH" = "red")) +
    theme_minimal()
})

# Mostrar los gráficos
grid.arrange(grafico_general, grobs = graficos_dias, ncol = 2, top = "Tasa de Acierto por Día y Total")