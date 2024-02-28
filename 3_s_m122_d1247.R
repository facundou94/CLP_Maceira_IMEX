################ MALDI-TOF ANALISIS CLP ########################################
################ 3) s_m122_d1247  #############################################
#
# s:     Supervisado
# m122:  Utiliza las wells (122 muestras)
# d1247: Utiliza todos los días de lectura para los análisis
#
# Autor: Bioing. Facundo Urteaga (IBB-CONICET)
#
#
### CARGA DE LIBRERIAS #########################################################
################################################################################


library(tidyverse)
library(crossval)
library(binda)
library(caret)
library(pROC)
library(ranger)
library(klaR)


### CARGA DE ARCHIVOS ##########################################################
################################################################################


#setwd(script_path) o session -> set working directory -> carpeta data
load("matint_122.Rdata")
load("matint_122_dico.Rdata")

### ANÁLISIS DE MACHINE LEARNING ###############################################
################################################################################


# Train/test splits
etiqueta <- df_metadata_prom_rep %>% pull(tipo)
etiqueta <- factor(etiqueta)
rownames(matint_122_dico) <- etiqueta

set.seed(123) # semilla para reproducibilidad

ind <- sample(2, nrow(matint_122_dico), replace = TRUE, 
              prob = c(0.6, 0.4))

# Separo en datos de entrenamiento y de testeo
Train <- matint_122_dico[ind == 1, ] 
Test <- matint_122_dico[ind == 2, ]
Ytrain1 <- rownames(Train)
Ytest1 <- rownames(Test)
Ytrain1 <- as.factor(Ytrain1)
Ytest1 <- as.factor(Ytest1)


### BINDA ######################################################################


# Binda ranking
tipo <- df_metadata_prom_rep %>% pull(tipo)
br <- binda.ranking(matint_122_dico, tipo, verbose = FALSE)
top20 <- br[1:20] #Uso el top20

tipo <- data.frame(tipo = tipo)
datos <- cbind(data.frame(matint_122_dico), tipo)

set.seed(1)
ind <- sample(2, nrow(datos), replace = TRUE, 
              prob = c(0.6, 0.4))

Train <- datos[ind == 1, ] 
Test <- datos[ind == 2, ]

Train %>% count(tipo)
Test %>% count(tipo)

Train <- Train[, c(top20, 219)]
Test <- Test[, c(top20, 219)]

particiones  <- 10
repeticiones <- 5

hiperparametros <- data.frame(lambda.freqs = c(0, 0.25, 0.5, 0.75, 1))

set.seed(123)

seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)

for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros)) 
}

seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)

control_train <- trainControl(method = "repeatedcv", 
                              number = particiones,
                              repeats = repeticiones, 
                              seeds = seeds,
                              returnResamp = "final", 
                              classProbs = T,
                              verboseIter = FALSE,
                              summaryFunction = twoClassSummary
)

set.seed(342)
Train_x = data.frame(Train[, -21])
Train_y = as.factor(Train[, 21])

modelo_binda <- train(Train_x, Train_y,
                      method = "binda",
                      metric = "Accuracy",
                      trControl = control_train)

Test_x = data.frame(Test[, -21])
Test_y = as.factor(Test[, 21])

pred.cv <- predict(modelo_binda, 
                   newdata = Test_x,
                   type = "prob")


r_bda <- roc(response = factor(Test_y), 
             predictor = pred.cv$CLP,
             levels = c("CLP", "SH")) 

auc(r_bda)

ci.auc(r_bda)

plot.roc(r_bda, print.auc=T,
         col="blue",xlab="1-Especificidad",ylab="Sensibilidad")


### RANDOM FOREST ##############################################################


particiones  <- 10
repeticiones <- 5

hiperparametros <- expand.grid(mtry = c(3, 4, 5, 7),
                               min.node.size = c(2, 3, 4, 5, 10, 15, 20, 30),
                               splitrule = "gini")

set.seed(123)
seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)

control_train <- trainControl(method = "repeatedcv", 
                              number = particiones,
                              repeats = repeticiones, 
                              seeds = seeds,
                              returnResamp = "final", 
                              classProbs = T,
                              verboseIter = FALSE,
                              allowParallel = TRUE,
                              summaryFunction = twoClassSummary)

modelo_rf <- train(Train_x, Train_y,
                   method = "ranger",
                   tuneGrid = hiperparametros,
                   metric = "ROC",
                   trControl = control_train,
                   # Número de árboles ajustados
                   num.trees = 500)

pred.cv <- predict(modelo_rf, 
                   newdata = Test_x,
                   type = "prob")

# Cálculo de la curva
curva_roc <- roc(response = factor(Test_y), 
                 predictor = pred.cv$CLP,
                 levels = c("CLP", "SH")) 

r_rf <- roc(response = factor(Test_y), 
            predictor = pred.cv$CLP,
            levels = c("CLP", "SH"), auc=T,ci=T) 

auc(r_rf)

ci.auc(r_rf)

plot.roc(r_rf, print.auc=T,
         col="blue",xlab="1-Especificidad",ylab="Sensibilidad")


### kNN ########################################################################


particiones  <- 10
repeticiones <- 5

hiperparametros <- data.frame(k = c(1, 2, 5, 10, 15, 20, 30, 50))

set.seed(1)
seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros)) 
}
seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)

control_train <- trainControl(method = "repeatedcv", 
                              number = particiones,
                              repeats = repeticiones, 
                              seeds = seeds,
                              returnResamp = "final", 
                              classProbs = T,
                              summaryFunction = twoClassSummary,
                              verboseIter = FALSE,
                              allowParallel = TRUE)

set.seed(342)
modelo_knn <- train(Train_x, Train_y,
                    method = "knn",
                    tuneGrid = hiperparametros,
                    metric = "ROC",
                    trControl = control_train)
modelo_knn

pred.cv <- predict(modelo_knn, 
                   newdata = Test_x,
                   type = "prob")

r_knn <- roc(response = factor(Test_y), 
             predictor = pred.cv$CLP,
             levels = c("CLP", "SH"))

auc(r_knn)

ci.auc(r_knn)

plot.roc(r_knn, print.auc=T,
         col="blue",xlab="1-Especificidad",ylab="Sensibilidad")


### SVM RADIAL #################################################################


particiones  <- 10
repeticiones <- 5

hiperparametros <- expand.grid(sigma = c(0.001, 0.01, 0.1, 0.5, 1),
                               C = c(1 , 20, 50, 100, 200, 500, 700))

set.seed(123)
seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)

control_train <- trainControl(method = "repeatedcv", 
                              number = particiones,
                              repeats = repeticiones, 
                              seeds = seeds,
                              returnResamp = "final", 
                              classProbs = T,
                              summaryFunction = twoClassSummary,
                              verboseIter = FALSE,
                              allowParallel = TRUE)

set.seed(342)
modelo_svmrad <- train(Train_x, Train_y,
                       method = "svmRadial",
                       tuneGrid = hiperparametros,
                       metric = "ROC",
                       trControl = control_train)

modelo_svmrad


pred.cv <- predict(modelo_svmrad, 
                   newdata = Test_x,
                   type = "prob")

r_svm <- roc(response = factor(Test_y), 
             predictor = pred.cv$CLP,
             levels = c("CLP", "SH"))

auc(r_svm)

ci.auc(r_svm)

plot.roc(r_svm, print.auc=T,
         col="blue",xlab="1-Especificidad",ylab="Sensibilidad")


### GRÁFICO DE TODAS LAS CURVAS ################################################


# Dibuja las curvas ROC para cada clasificador
plot.roc(r_svm, col = "blue4", xlab = "1-Especificidad", ylab = "Sensibilidad")
plot.roc(r_knn, col = "red4", add = TRUE)
plot.roc(r_bda, col = "green4", add = TRUE)
plot.roc(r_rf, col = "yellow4", add = TRUE)

# Muestra el valor del AUC para cada clasificador en una ubicación deseada
text(0.0, 0.7, labels = paste("AUC (SVM) =", round(auc(r_svm), 2)), 
     pos = 1, offset = 0.5, col = "blue4")
text(0.1, 0.6, labels = paste("AUC (KNN) =", round(auc(r_knn), 2)), 
     pos = 1, offset = 0.5, col = "red4")
text(0.2, 0.5, labels = paste("AUC (BDA) =", round(auc(r_bda), 2)), 
     pos = 1, offset = 0.5, col = "green4")
text(0.3, 0.4, labels = paste("AUC (RF) =", round(auc(r_rf), 2)), 
     pos = 1, offset = 0.5, col = "yellow4")

# Añade la leyenda
legend("bottomright", legend = c("ROC SVM", "ROC KNN", "ROC BDA", "ROC RF"), 
       col = c("blue4", "red4", "green4", "yellow4"), lwd = 2)
title("Curvas ROC para Clasificadores", line = 2.5)


#
#
#
### FIN ########################################################################
################################################################################