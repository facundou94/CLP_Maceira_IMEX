################ MALDI-TOF ANALISIS CLP ########################################
################ 5) ns_m51_d24  ################################################
#
# ns:    No Supervisado
# m51:  Utiliza las muestras (51)
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


#script_path <- "C:/Users/Facundo/Documents/Proyectos/Data"
#setwd(script_path) session -> set working directory
load("matint_51_dico.Rdata")
load("matint_51.Rdata")

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
  scale_color_manual(values = c("maroon1","aquamarine3","blueviolet", "blue1",
                                "blue4", "maroon4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.hkmeans.top15.k3)

# Con top15: PAM clustering 3 clusters

top_actual <- top.b15
K.num <- 3 # clusters
var2 = 0.95

pam.top15.k3 <- pam(matint_51_dico_d2d4[, top_actual], metric = "manhattan",
                            K.num)

cluster.pam.top15.k3 <- fviz_cluster(pam.top15.k3, ellipse.type = "convex", 
                                     data = matint_51_dico_d2d4[, top_actual],
                                     ellipse.level = var2,
                                     show.clust.cent = F, 
                                     geom = "point",
                                     main = "PAM - Top 15 - 3 clusters")

# Ajustar el tamaño de los puntos según los valores de la columna "dia"
cluster.pam.top15.k3 <- cluster.pam.top15.k3 + 
  geom_point(data = cluster.pam.top15.k3$data, 
             aes(x = x, y = y, color = df_unicas_d2d4$factor1,
                 size = df_unicas_d2d4$dia)) +
  scale_color_manual(values = c("maroon1","aquamarine3","blueviolet", "blue1",
                                "blue4", "maroon4")) +
  scale_size_continuous(range = c(2, 3)) +
  labs(color = "Cluster", size = "Día") +
  theme(legend.position = "right")

print(cluster.pam.top15.k3)

#
#
#
### FIN ########################################################################
################################################################################