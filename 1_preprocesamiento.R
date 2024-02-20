################ MALDI-TOF ANALISIS CLP #####################################################################################
##
##
## Autor: Bioing. Facundo Urteaga (IBB-CONICET)
##
##
##################### CARGA DE LIBRERIAS ####################################################################################
#
#
library(here)
library(readBrukerFlexData)
library(MALDIquant)
library(MALDIquantForeign)
library(MALDIrppa)
library(stringr)
#
#
##################### CARGA DE ESPECTROS ####################################################################################
#
#
# Creación de la ruta relativa de los archivos
#
ruta_proyecto <- here("OneDrive","Documents","Proyectos RStudio","Maceira_CLP","Datos CLP-SHAM D1-2-4-7")
ruta_datos <- file.path(ruta_proyecto)
#
# Importar espectros
Spectra_list <- importBrukerFlex(file.path(ruta_datos), verbose=FALSE)
#
# Lectura de nombre de espectros, características y ploteo (Opcional/Chequeo)
#
#for(i in 1:length(Spectra_list)) {
#  print(Spectra_list[[i]]@metaData$sampleName)
#}
#summarySpectra(Spectra_list[1:10])
#plot(Spectra_list[[1]])
#
#
###### OBTENCIÓN DE METADATA DE ESPECTROS ###################################################################################
#
# Creación de columnas vacías
#
col_dia <- c()
col_tipo <- c()
col_numero <- c()
col_well <- c()
col_rep <- c()
#
# Patrones auxiliares para buscar el día de la muestra
patron <- "D1-2-4-7"
patron2 <- "D(1|2|4|7)"
# Patrón para extraer el nombre (todo antes de los números)
patron_nombre <- "^[A-Za-z]+"
# Patrón para extraer el número (todo después de las letras)
patron_numero <- "\\d+$"
#
# Ciclo que se encarga de extraer dia, tipo, numero, well y réplica de cada muestra
#
for(i in 1:length(Spectra_list)) {
  nombre <- Spectra_list[[i]]@metaData$file
  
  # Encuentra la posición del patrón en el texto
  posicion <- str_locate(nombre, patron)[1, 2]
  
  # Extrae la parte de la ruta después del patrón
  resultado <- substr(nombre, (posicion+2), nchar(nombre))
  resultado <- str_extract(resultado, patron2)
  
  col_dia <- c(col_dia, resultado)
  
  partes <- unlist(strsplit(Spectra_list[[i]]@metaData$sampleName, "_+"))
  tipo <- str_extract(partes[1], patron_nombre)
  numero <- str_extract(partes[1], patron_numero)
  col_tipo <- c(col_tipo, tipo)
  col_numero <- c(col_numero, numero)
  col_well <- c(col_well, partes[2])
  col_rep <- c(col_rep, partes[3])
  
}
#
#
# Data Frame con los datos limpios
#
df_metadata <- data.frame(dia = col_dia, tipo = col_tipo, numero = col_numero, well = col_well,
                        replica = col_rep)
#
# Creación de factores de agrupamiento para su uso posterior
#
df_metadata$factor1 <- paste0(df_metadata$tipo, "_", df_metadata$dia)
df_metadata$factor2 <- paste0(df_metadata$tipo, "_", df_metadata$dia, "_", df_metadata$numero)
df_metadata$factor3 <- paste0(df_metadata$tipo, "_", df_metadata$dia, "_", df_metadata$numero, "_", df_metadata$well)
#
#
######### ASIGNAR A METADATA DE ESPECTROS ##############################################################################
#
#
# Creo nueva metaData en espectro
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
#
#
############### LIMPIEZA  ##############################################################################################
############### Y TRANSFORMACIÓN #######################################################################################
############### DE ESPECTROS ###########################################################################################
#
#
#
#### SCREENING INICIAL #################################################################################################
# detección de espectros de baja calidad
#
sc.results <- screenSpectra(Spectra_list, meta = df_metadata)
summary(sc.results)
plot(sc.results, labels = TRUE)
#
# Ploteo que no pasan control de calidad y uno correcto (Aprox 6%, bien)
# plot(Spectra_list[[253]]) # Espectro muy ruidoso
#
## Descartamos espectros defectuosos 
#
Spectra_list_f1 <- sc.results$fspectra # Filtramos espectros
df_metadata_f1 <- sc.results$fmeta # Filtramos metadatos
#
#
####### PROCESAMIENTO Y FILTRADO DE ESPECTROS ##########################################################################
# Parámetros de procesamiento de espectros
#
thScale <- 10 # Smoothing
#smoothLevel <- 4 # Wavelet smoothing level
ite <- 105 # Baseline correction
SigNoi <- 2 # Peak extraction
hws <- 20 # Peak extraction
tol <- 0.03 # Peak binning
#
# proceso correccion de espectros
# utiliza los parámetros definidos anteriormente
#
Spectra_list_f1 <- transfIntensity(Spectra_list_f1, fun = sqrt)
plot(Spectra_list_f1[[20]])
Spectra_list_f1 <- wavSmoothing(Spectra_list_f1, method = "Wavelet", n.levels = 4)
plot(Spectra_list_f1[[20]])
#spectra <- wavSmoothing(spectra, method = "Wavelet", n.levels = smoothLevel)
#
# Detección de la linea de base
#
baseline <- estimateBaseline(Spectra_list_f1[[20]], method = "SNIP",
                             iterations = ite)
plot(Spectra_list_f1[[20]])
lines(baseline, col="red", lwd=2)
#
#
Spectra_list_f2 <- removeBaseline(Spectra_list_f1, method = "SNIP", iterations = ite)
plot(Spectra_list_f2[[20]])
Spectra_list_f2 <- calibrateIntensity(Spectra_list_f2, method = "PQN")
plot(Spectra_list_f2[[20]])

Spectra_list_f3 <- alignSpectra(Spectra_list_f2,
                                halfWindowSize=20,
                                SNR=2,
                                tolerance=0.02, #Ojo, ver parámetros
                                warpingMethod="lowess")
#
plot(Spectra_list_f3[[20]])
#
### OJO !!!!
# for (i in 1:length(Spectra_list_f3)) {
#   valor1 <- Spectra_list_f3[[i]]@metaData$sampleName
#   valor2 <- Spectra_list_f3[[i]]@metaData$factor3
#   print(valor1) #Ambos n de muestra deben coincidir
#   print(valor2)
# }
# 
# #Buscador de indices
# #for (i in 1:length(SL_transformado)) {
# #  if (SL_transformado[[i]]@metaData$muestra == 101)
# #  print(i)
# #}








