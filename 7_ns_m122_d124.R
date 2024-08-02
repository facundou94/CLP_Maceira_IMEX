################ MALDI-TOF ANALISIS CLP ########################################
################ 4) ns_m122_d24  ###############################################
#
# ns:    No Supervisado
# m122:  Utiliza las wells (122 muestras)
# d1247: Utiliza los días 1, 2 y 4 (junta 1 con 2)
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


#script_path <- "C:/Users/Facundo/Documents/Proyectos/Data"
#setwd(script_path) session -> set working directory
# Define the relative path to the Data folder

data_path <- file.path("Data")

# Load the Rdata files using the relative path
load(file.path(data_path, "matint_122_dico.Rdata"))
load(file.path(data_path, "matint_122.Rdata"))


df_metadata_prom_rep$dia <- as.integer(gsub("[^0-9]", "", 
                                            df_metadata_prom_rep$dia))

# Filtrado de muestras de día 1 2 y 4
filas_filtradas <- rownames(matint_122_dico)[grepl("D1|D2|D4", rownames(matint_122_dico))]
matint_122_dico_d1d2d4 <- matint_122_dico[filas_filtradas, ]
filas_filtradas <- rownames(matint_122)[grepl("D1|D2|D4", rownames(matint_122))]
matint_122_d1d2d4 <- matint_122[filas_filtradas, ]

df_prom_rep_d1d2d4 <- df_metadata_prom_rep %>% dplyr::filter(dia == 1 | dia == 2 | dia == 4)

# Agrego columna "factor4" que agrupe los dias 1 y 2 y los distinga del dia 4
df_prom_rep_d1d2d4 <- df_metadata_prom_rep %>% 
  dplyr::filter(dia == 1 | dia == 2 | dia == 4) %>%
  dplyr::mutate(factor4 = ifelse(dia == 4, 4, 2))

# Agregar la columna "factor5" basada en las condiciones de "factor1"
df_prom_rep_d1d2d4 <- df_prom_rep_d1d2d4 %>%
  dplyr::mutate(factor5 = case_when(
    grepl("SH", factor1) ~ "CTRL",
    factor1 %in% c("CLP_D1", "CLP_D2") ~ "D1D2",
    factor1 == "CLP_D4" ~ "D4",
    TRUE ~ NA_character_  # Valor por defecto si ninguna condición se cumple
  ))


### SELECCIÓN DE PICOS #########################################################
################################################################################


# Selección de picos para binary discriminant analysis (BDA)

factor_tipo <- factor(df_prom_rep_d1d2d4$factor4)
is.binaryMatrix(matint_122_dico_d1d2d4) # TRUE
br <- binda.ranking(matint_122_dico_d1d2d4, factor_tipo, verbose = FALSE)
indices_especificos <- br[1:20]

# Gráfico de picos vs score 
nueva_columna <- c()
matriz <- matrix(br, nrow = 218, ncol = 4)
for (i in 1:218) {
  nuevo_valor <- colnames(matint_122_dico_d1d2d4)[br[i]]
  nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)

plot(df_br$nueva_columna, df_br$V2, 
     xlab = "m/z", ylab = "Score", 
     main = "Ranking de picos de los espectros, 4 clusters")
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
top_actual <- top.b15

# Elección de mejores algoritmos de clustering
comparacion <- clValid(
  obj        = matint_122_dico_d1d2d4[, top_actual],
  nClust     = 2:6,
  clMethods  = c("hierarchical", "kmeans", "pam"),
  validation = c("stability", "internal")
)

summary(comparacion)
optimalScores(comparacion) # Se puede ir probando con distintos top picos

## Obtenemos que con TOP 10: PAM 2 Clusters, TOP15: hkmeans 2, PAM 6, kmeans 3
## TOP 20: PAM 2, kmeans 3


### ALGORITMOS DE CLUSTERING ###################################################
################################################################################


# Con top15: KMEANS clustering 3 clusters
top_actual <- top.b15
K.num <- 3 # clusters
var2 = 0.95

kmeans.top15.k3 <- kmeans(matint_122_dico_d1d2d4[, top_actual], 
                            K.num)

cluster.kmeans.top15.k3 <- fviz_cluster(
  kmeans.top15.k3, 
  ellipse.type = "convex", 
  data = matint_122_dico_d1d2d4[, top_actual],
  ellipse.level = var2,
  show.clust.cent = F, 
  geom = "point", 
  main = "kmeans - Top 15 - 3 clusters")

cluster.kmeans.top15.k3 <- cluster.kmeans.top15.k3 + 
  geom_point(data = cluster.kmeans.top15.k3$data, 
             aes(x = x, y = y, color = df_prom_rep_d1d2d4$factor5, 
                 size = df_prom_rep_d1d2d4$dia)) +
  scale_color_manual(values = c("blueviolet","aquamarine3","maroon1", "maroon4", 
                                "blue1","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +  
  theme(legend.position = "right")

print(cluster.kmeans.top15.k3)

# Con top15: HKMEANS clustering 2 clusters
top_actual <- top.b15
K.num <- 2 # clusters
var2 = 0.95

hkmeans.top15.k2 <- hkmeans(matint_122_dico_d1d2d4[, top_actual], 
                           K.num)

cluster.hkmeans.top15.k2 <- fviz_cluster(
  hkmeans.top15.k2, 
  ellipse.type = "convex", 
  data = matint_122_dico_d1d2d4[, top_actual],
  ellipse.level = var2,
  show.clust.cent = F, 
  geom = "point", 
  main = "kmeans - Top 15 - 3 clusters")

cluster.hkmeans.top15.k2 <- cluster.hkmeans.top15.k2 + 
  geom_point(data = cluster.hkmeans.top15.k2$data, 
             aes(x = x, y = y, color = df_prom_rep_d1d2d4$factor5, 
                 size = df_prom_rep_d1d2d4$dia)) +
  scale_color_manual(values = c("blueviolet","aquamarine3","maroon1", "maroon4", 
                                "blue1","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +  
  theme(legend.position = "right")

print(cluster.hkmeans.top15.k2)



# Con top20: PAM 2 clusters
top_actual <- top.b20
K.num <- 2 # clusters
var2 = 0.95

pam.top20.k2 <- pam(matint_122_dico_d1d2d4[, top_actual], 
                           K.num)

cluster.pam.top20.k2 <- fviz_cluster(
  pam.top20.k2, 
  ellipse.type = "convex", 
  data = matint_122_dico_d1d2d4[, top_actual],
  ellipse.level = var2,
  show.clust.cent = F, 
  geom = "point", 
  main = "PAM - Top 20 - 2 clusters")

cluster.pam.top20.k2 <- cluster.pam.top20.k2 + 
  geom_point(data = cluster.pam.top20.k2$data, 
             aes(x = x, y = y, color = df_prom_rep_d1d2d4$factor5, 
                 size = df_prom_rep_d1d2d4$dia)) +
  scale_color_manual(values = c("blueviolet","maroon1", "maroon4", 
                                "blue1","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +  
  theme(legend.position = "right")

print(cluster.pam.top20.k2)


# Con top20: kmeans 3 clusters
top_actual <- top.b20
K.num <- 3 # clusters
var2 = 0.95

kmeans.top20.k3 <- kmeans(matint_122_dico_d1d2d4[, top_actual], 
                       K.num)

cluster.kmeans.top20.k3 <- fviz_cluster(
  kmeans.top20.k3, 
  ellipse.type = "convex", 
  data = matint_122_dico_d1d2d4[, top_actual],
  ellipse.level = var2,
  show.clust.cent = F, 
  geom = "point", 
  main = "kmeans - Top 20 - 3 clusters")

cluster.kmeans.top20.k3 <- cluster.kmeans.top20.k3 + 
  geom_point(data = cluster.kmeans.top20.k3$data, 
             aes(x = x, y = y, color = df_prom_rep_d1d2d4$factor5, 
                 size = df_prom_rep_d1d2d4$dia)) +
  scale_color_manual(values = c("blueviolet","aquamarine3","maroon1", "maroon4", 
                                "blue1","blue4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +  
  theme(legend.position = "right")

print(cluster.kmeans.top20.k3)