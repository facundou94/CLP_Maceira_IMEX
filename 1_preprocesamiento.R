################ MALDI-TOF ANALISIS CLP ########################################
################ 1) PREPROCESAMIENTO  ##########################################
#
# Autor: Bioing. Facundo Urteaga (IBB-CONICET)
#
#
### CARGA DE LIBRERIAS #########################################################
################################################################################


library(binda)
library(here)
library(dplyr)
library(readBrukerFlexData)
library(MALDIquant)
library(MALDIquantForeign)
library(MALDIrppa)
library(stringr)


### CARGA DE ESPECTROS #########################################################
################################################################################


# Creación de la ruta relativa de los archivos
ruta_proyecto <- here("Datos CLP-SHAM D1-2-4-7")
ruta_datos <- file.path(ruta_proyecto)

# Importar espectros
Spectra_list <- importBrukerFlex(file.path(ruta_datos), verbose=FALSE)


### OBTENCIÓN DE METADATA DE ESPECTROS #########################################
################################################################################


# Creación de columnas vacías
col_dia <- c()
col_tipo <- c()
col_numero <- c()
col_well <- c()
col_rep <- c()

# Patrones auxiliares para buscar el día de la muestra
patron <- "D1-2-4-7"
patron2 <- "D(1|2|4|7)"

# Patrón para extraer el nombre (todo antes de los números)
patron_nombre <- "^[A-Za-z]+"

# Patrón para extraer el número (todo después de las letras)
patron_numero <- "\\d+$"

# Ciclo que extrae dia, tipo, numero, well y réplica de cada muestra
for(i in 1:length(Spectra_list)) {
    nombre <- Spectra_list[[i]]@metaData$file
    
    # Encuentra la posición del patrón en el texto
    posicion <- str_locate(nombre, patron)[1, 2]
    
    # Extrae la parte de la ruta después del patrón
    resultado <- substr(nombre, (posicion+2), nchar(nombre))
    resultado <- str_extract(resultado, patron2)
    
    # Almaceno en columnas
    col_dia <- c(col_dia, resultado)
    partes <- unlist(strsplit(Spectra_list[[i]]@metaData$sampleName, "_+"))
    tipo <- str_extract(partes[1], patron_nombre)
    numero <- str_extract(partes[1], patron_numero)
    col_tipo <- c(col_tipo, tipo)
    col_numero <- c(col_numero, numero)
    col_well <- c(col_well, partes[2])
    col_rep <- c(col_rep, partes[3])
}

# Data Frame con los datos limpios
df_metadata <- data.frame(dia = col_dia, tipo = col_tipo, numero = col_numero, 
                          well = col_well, replica = col_rep)

# Creación de factores de agrupamiento para su uso posterior
df_metadata$factor1 <- paste0(df_metadata$tipo, "_", df_metadata$dia)
df_metadata$factor2 <- paste0(df_metadata$tipo, "_", df_metadata$dia, "_", 
                              df_metadata$numero)
df_metadata$factor3 <- paste0(df_metadata$tipo, "_", df_metadata$dia, "_", 
                              df_metadata$numero, "_", df_metadata$well)

# Asigno metaData a espectros
for(i in 1:length(Spectra_list)) {
    Spectra_list[[i]]@metaData$tipo = df_metadata$tipo[i]
    Spectra_list[[i]]@metaData$dia = df_metadata$dia[i]
    Spectra_list[[i]]@metaData$numero = df_metadata$numero[i]
    Spectra_list[[i]]@metaData$well = df_metadata$well[i]
    Spectra_list[[i]]@metaData$replica = df_metadata$replica[i]
    Spectra_list[[i]]@metaData$factor1 = df_metadata$factor1[i]
    Spectra_list[[i]]@metaData$factor2 = df_metadata$factor2[i]
    Spectra_list[[i]]@metaData$factor3 = df_metadata$factor3[i]
}


### CONTROL DE CALIDAD Y LIMPIEZA DE ESPECTROS #################################
################################################################################


# Screening inicial: Detección de espectros de baja calidad
sc.results <- screenSpectra(Spectra_list, meta = df_metadata)
summary(sc.results)
plot(sc.results, labels = TRUE)

# plot(Spectra_list[[253]]) # Ploteo de espectros ruidosos

# Descartamos espectros defectuosos 
Spectra_list_f1 <- sc.results$fspectra # Filtramos espectros
df_metadata_f1 <- sc.results$fmeta # Filtramos metadatos


### FILTRADO Y TRANSFORMACIÓN DE ESPECTROS #####################################
################################################################################


# Parámetros de procesamiento de espectros
thScale <- 10 # Smoothing
ite <- 105 # Baseline correction
SigNoi <- 2.5 # Peak extraction
hws <- 20 # Peak extraction
tol <- 0.03 # Peak binning

# Transformación/filtrado/corrección de espectros con parámetros definidos
# 1) Transformación de intensidad por medio de función sqrt
Spectra_list_f1 <- transfIntensity(Spectra_list_f1, fun = sqrt)
plot(Spectra_list_f1[[30]])
# 2) Suavizado del espectro mediante el método Wavelet
Spectra_list_f1 <- wavSmoothing(Spectra_list_f1, method = "Wavelet", n.levels = 4)
plot(Spectra_list_f1[[30]])
# Detección de la linea de base
baseline <- estimateBaseline(Spectra_list_f1[[30]], method = "SNIP",
                             iterations = ite)
plot(Spectra_list_f1[[30]])
lines(baseline, col="red", lwd=2)
# 3) Remoción de linea de base mediante método SNIP
Spectra_list_f2 <- removeBaseline(Spectra_list_f1, method = "SNIP", 
                                  iterations = ite)
plot(Spectra_list_f2[[30]])
# 4) Calibración de intensidad mediante método PQN
Spectra_list_f2 <- calibrateIntensity(Spectra_list_f2, method = "PQN")
plot(Spectra_list_f2[[30]])
# 5) Alineación de espectros
Spectra_list_f3 <- alignSpectra(Spectra_list_f2,
                                halfWindowSize=20,
                                SNR=2,
                                tolerance=0.02, # Parámetro sensible
                                warpingMethod="lowess")
plot(Spectra_list_f3[[30]])

# 6) Por control descarto espectros que detectaron una cantidad anómala de picos
indices_a_preservar <- c(1:296, 299)
Spectra_list_f3 <- Spectra_list_f3[indices_a_preservar]
df_metadata_f1 <- df_metadata_f1 %>% 
  slice(indices_a_preservar)


### PROMEDIO DE LECTURAS DE UNA MISMA WELL Y UNA MISMA MUESTRA #################
################################################################################


# Promedio de lecturas de una misma well
Spectra_list_prom_rep <- averageMassSpectra(Spectra_list_f3,
                                      labels = factor(df_metadata_f1$factor3), 
                                      method = "mean")

# Creo la nueva metadata de los espectros promediados
df_metadata_prom_rep <- df_metadata_f1 %>% 
  distinct(factor3, .keep_all = TRUE)

# Promedio de wells de una misma muestra
Spectra_list_prom_muestra <- averageMassSpectra(Spectra_list_prom_rep,
                                            labels = factor(df_metadata_prom_rep$factor2), 
                                            method = "mean")


### EXTRACCIÓN DE PICOS Y ALINEACIÓN ###########################################
################################################################################


# Análisis de la SNR en espectros para chequear que utilizamos el valor correcto
noise <- estimateNoise(Spectra_list_prom_rep[[122]])
plot(Spectra_list_prom_rep[[122]], xlim=c(4000, 20000), ylim=c(0, 0.0005))
lines(noise, col="red")
lines(noise[,1], noise[, 2]*2, col="blue") # Se ve que es correcto el 2

# Detección de picos a partir del umbral definido de SNR
peaks <- detectPeaks(Spectra_list_prom_rep, 
                     SNR = SigNoi, 
                     halfWindowSize = 40)

# Ploteo de picos detectados en un espectro de ejemplo
plot(Spectra_list_prom_rep[[122]], xlim=c(4000, 20000), ylim=c(0, 0.0005))
points(peaks[[122]], col="red", pch=4)

# Alineado de picos
peaks <- alignPeaks(peaks, 
                    minFreq = 0.8, 
                    tolerance = tol)

#summaryPeaks(peaks[1:10])  # resumen estadistico de picos (primeros 10)

# Conteo de picos por perfil
cP <- countPeaks(peaks)

# Gráfico de picos
plot(cP, type = "n")
text(cP, label = 1:length(cP))


# Patrones de picos
peakPatterns(peaks)

# Filtrado de picos de baja frecuencia de aparición
picos_filtrados <- filterPeaks(peaks, 
                               minFreq = 0.25, 
                               labels = df_metadata_prom_rep$factor1 ) #labels

# Metadata de picos filtrados
df_metadata_unicas <- df_metadata_prom_rep %>% 
  distinct(factor2, .keep_all = TRUE)

# Patrones de picos
peakPatterns(picos_filtrados)

# Conteo de picos por perfil
cP2 <- countPeaks(picos_filtrados)

# Gráfico
plot(cP2, type = "n")
text(cP2, label = 1:length(cP2))

# Fusión de picos de la misma muestra
picos_fusion_muestra <- mergeMassPeaks(picos_filtrados, 
                                       labels = df_metadata_prom_rep$factor2, 
                                       method = "median")

# Patrones de picos
peakPatterns(picos_fusion_muestra)


### CREACIÓN DE MATRIZ DE INTENSIDADES Y DICOTÓMICA ############################
################################################################################


# Matriz de intensidades de 122 wells
# matint_na_122 <- intensityMatrix(picos_filtrados) # con valores NA
matint_122 <- intensityMatrix(picos_filtrados, 
                              Spectra_list_prom_rep) # sin valores NA

# Matriz de intensidades de 51 muestras
# matint_na_51 <- intensityMatrix(picos_fusion_muestra) # con valores NA
matint_51 <- intensityMatrix(picos_fusion_muestra, 
                             Spectra_list_prom_muestra) # sin valores NA

# Definición de umbrales
thr1 <- optimizeThreshold(matint_122,
                         df_metadata_prom_rep$tipo, 
                         verbose = T)
thr2 <- optimizeThreshold(matint_51,
                          df_metadata_unicas$tipo, 
                          verbose = T)

# Dicotomización
matint_122_dico <- dichotomize(matint_122, thr1)
matint_51_dico <- dichotomize(matint_51, thr2)

# Agrego nombres a las filas de cada df
rownames(matint_122_dico) <- df_metadata_prom_rep$factor3
rownames(matint_51_dico) <- df_metadata_unicas$factor2
rownames(matint_122) <- df_metadata_prom_rep$factor3
rownames(matint_51) <- df_metadata_unicas$factor2


### GUARDAR DATOS ##############################################################
################################################################################


# Guardar matrices y metadata asociada como archivo .Rdata
save(matint_122_dico,df_metadata_prom_rep, file = "matint_122_dico.Rdata")
save(matint_51_dico, df_metadata_unicas, file = "matint_51_dico.Rdata")
save(matint_122,df_metadata_prom_rep, file = "matint_122.Rdata")
save(matint_51, df_metadata_unicas, file = "matint_51.Rdata")

# Guardar matrices y metadata asociada como archivo .csv
write.csv(matint_122_dico, "matint_122_dico.csv", row.names = TRUE)
write.csv(matint_122, "matint_122.csv", row.names = TRUE)
write.csv(matint_51_dico, "matint_51_dico.csv", row.names = TRUE)
write.csv(matint_51, "matint_51.csv", row.names = TRUE)
write.csv(df_metadata_prom_rep, "df_122.csv", row.names = TRUE)
write.csv(df_metadata_unicas, "df_51.csv", row.names = TRUE)

#
#
#
### FIN ########################################################################
################################################################################