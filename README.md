<h1 align="center"> Análisis de muestras de MALDI-TOF de ratones.  </h1>
<h3 align="center"> Origen de las muestras: Instituto Malbrán </h3>
<h3 align="center"> Responsables: Lautaro Maceira, Bárbara Rearte  </h3>

<p align="left">
   <img src="https://img.shields.io/badge/ESTADO-EN%20DESAROLLO-green">
   </p>

## Índice

*[Descripción del proyecto](#descripción-del-proyecto)

*[Características de las muestras](#características-de-las-muestras)

*[Características y nomenclatura de  los archivos](#características-y-nomenclatura-de-los-archivos)

*[Pre-procesamiento](#pre-procesamiento)

*[Seguimiento de cantidad de muestras](#seguimiento-de-cantidad-de-muestras)

*[Algoritmos no supervisados](#algoritmos-no-supervisados)

*[Algoritmos supervisados](#algoritmos-supervisados)

*[Conclusión](#conclusión)

## Descripción del proyecto

<p align="justify">
El proyecto consiste en el procesamiento y análisis de espectros MALDI-TOF obtenidos a partir de plasma de ratones a los que se le ha realizado una ligadura y punción cecal (CLP) como modelo de sepsis, y ratones impostores (SHAM). El procesamiento consiste en el filtrado, acondicionamiento y transformación de los espectros. El análisis consiste en la utilización de algoritmos de aprendizaje maquinal que permitan detectar características distintivas de estos espectros y clasificar tanto las muestras de CLP como sus diferentes estadíos dentro del modelo de la patología.
</p>

## Características de las muestras

* Tipo de muestras: CLP y SHAM
* Días de adquisición: 1, 2, 4 y 7
* Cantidad de muestras iniciales: 303
</p>
* Metodología de adquisición: Cada muestra es sub-dividida en "réplicas biológicas" que son depositadas cada una en un pozo o "well" del equipo. Cada una de estas wells puede ser adquirida o leída mas de una vez, obteniéndose así "réplicas técnicas". La cantidad de réplicas biológicas y técnicas por muestra es variable, llegando a un máximo de 3 de cada una. Es decir, en el caso mas extremo, una muestra podría ser replicada tres veces biológicamente, y cada una adquirida otras tres veces, llegando así a un número de nueve adquisiciones correspondientes a una misma muestra.
</p>

## Características y nomenclatura de  los archivos

</p>
Los lenguajes utilizados para el procesamiento y el análisis de los espectros fueron pyhton y R. Los archivos están enumerados por orden de procesamiento. A excepción del primer archivo de pre-procesamiento, los archivos están nomenclados de la siguiente manera:
</p>

*x_alg_muestras_dias*

Donde

* _x: numeración_
* _alg: algoritmos utilizados, supervisados (s) o no supervisados (ns)_
* _muestras: cantidad de muestras, indicando si se usan las réplicas biológicas (m122) o las muestras ya promediadas (m51)_
* _dias: días utilizados para el análisis, todos (d1247) o solo días 2 y 4 (d24)_

## Pre-procesamiento

### Archivo: 1_preprocesamiento.R

El pre-procesamiento está compuesto por las siguientes etapas:

* Carga de los espectros y su metadata correspondiente

  <p align="center">
     <img src="Imagenes/1_pre_crudo.jpeg" width="400">
   </p>
   <p align="center">
     <em>Figura 1: Espectro cargado sin transformar</em>
   </p>

* <p align="justify"> Control de calidad de los espectros mediante el uso de un estimador robusto Q. A fines prácticos, este control de calidad filtra espectros ruidosos o con el espectro "planchado" debido al fenómeno de supresión iónica en la etapa de adquisición. </p>

* <p align="justify"> Transformación de los espectros. Transformación de intensidad mediante la función raíz cuadrada (sqrt), suavizado del espectro mediante la función "wavelet", detección y remoción de la linea de base y alineamiento de los espectros, en ese orden. </p>

   <p align="center">
   <img src='Imagenes/1_pre_baseline.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 2: Detección de la linea de base</em>
   </p>
   
   <p align="center">
   <img src='Imagenes/1_pre_baseline_removed.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 3: Espectro con la línea de base removida</em>
   </p>

* Se realiza un promediado de las réplicas técnicas y biológicas.

* <p align="justify">Extracción de picos preponderantes de cada espectro. Esto se logra definiendo un umbral a partir del cual se comienzan a detectar los picos. Este umbral se define a partir de dos veces la relación señal a ruido del espectro (SNR).</p>

   <p align="center">
   <img src='Imagenes/1_pre_snr.jpeg' width='400'>
   </p>
   <p align="center">
     <em>Figura 4: Espectro con la detección del nivel de ruido (en rojo) y la definición del umbral (en azul)</em>
   </p>

   <p align="center">
   <img src='Imagenes/1_pre_picos.jpeg' width='400'>
   </p>
   <p align="center">
     <em>Figura 5: Detección de picos por encima del umbral establecido</em>
   </p>

* <p align="justify">A partir de la detección de los picos en cada espectro, se crea la matriz de intensidad, las cual contiene en sus filas las muestras y en las columnas los picos detectados. Esta matriz también es sujeta a transformaciones para preservar los picos con mayor frecuencia de aparición en los espectros y eliminar los picos "extraños", ya que lo que buscamos en esta instancia es que las variables (en este caso los picos) aporten información al sistema para su posterior análisis. Se crea también la matriz dicotómica, la cual se origina a partir de la definición de un umbral en la matriz de intensidades que transforma los valores de los picos en 1 y 0 segun la presencia o ausencia de cada pico en cada muestra.</p>

   <p align="center">
   <img src='Imagenes/1_pre_matriz_grafica.jpeg' width='400'>
   </p>
   <p align="center">
     <em>Figura 6: Representación gráfica de la matriz de intensidades dicotómica. Las filas corresponden a las muestras y las columnas a los picos. El color celeste indica presencia del pico en esa muestra</em>
   </p>

   
* <p align="justify">Por último, guardar las matrices y los metadatos. Se obtienen matrices de intensidades y dicotómicas tanto para las muestras individuales (matriz de 51 filas x 218 columnas) como para las réplicas biológicas (matriz de 122 filas x 218 columnas)</p>

## Seguimiento de cantidad de muestras

* Muestras iniciales: 303

## Algoritmos No Supervisados

El procedimiento para la realización de los algoritmos No Supervisados fue el siguiente:
   1) Elección de conjunto de muestras (Réplicas biológicas o muestras independientes)
   2) Elección de los tiempos de muestreo (Todos los días o solo los días 2 y 4)
   3) Por medio de la matriz dicotómica, se aplica la función *bindaranking* con la cual se obtienen los picos que mejor variabilidad aportan a partir de un factor que se ingresa como variable de entrada. Este factor puede ser CLP vs SHAM, o días por ejemplo.
   4) Se realizan simulaciones de los modelos variando la cantidad de X primeros picos del análisis realizado en (3) y los algoritmos de clustering. Se realizaron pruebas con kmeans, HKmeans y PAM.
   5) Una vez realizada la clasificación No Supervisada, se comparan los puntos clasificados con su etiqueta de interés real (CLP, SHAM).

Resultados:

   ### Archivo: 2_ns_m122_d1247

   <p align="center">
   <img src='Imagenes/2_ranking_picos.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 7: Picos mas preponderantes seleccionados por el algoritmo bindaranking a partir del factor CLP vs SHAM</em>
   </p>

   <p align="center">
   <img src='Imagenes/2_CLP_vs_SHAM_kmeans.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 8: Clustering - CLP vs SHAM - Días 1, 2, 4 y 7 - TOP 20 picos - Algoritmo: kmeans</em>
   </p>

   <p align="center">
   <img src='Imagenes/2_tasa_acierto_total.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 9: Tasa de acierto total</em>
   </p>

   <p align="center">
   <img src='Imagenes/2_tasa_acierto_por_dia.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 10: Tasa de acierto por día</em>
   </p>
   
   ### Archivo: 4_ns_m122_d24

   <p align="center">
   <img src='Imagenes/4_picos_4clusters.jpeg' width='400'>
   </p>
   <p align="center">
     <em>Figura 11: Picos mas preponderantes seleccionados por el algoritmo bindaranking a partir del factor CLP_D2 vs CLP_D4 vs SHAM_D4 vs SHAM_D2</em>
   </p>

   <p align="center">
   <img src='Imagenes/4_hkmeans_4_grupos.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 12: Clustering - CLP vs SHAM - Días 2 y 4 - TOP 30 picos - Algoritmo: Hkmeans</em>
   </p>

   <p align="center">
   <img src='Imagenes/4_kmeans_2_clusters_4_grupos.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 12: Clustering - CLP vs SHAM - Días 2 y 4 - TOP 15 picos - Algoritmo: kmeans</em>
   </p>

   <p align="center">
   <img src='Imagenes/4_kmeans_top10.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 13: Clustering - CLP vs SHAM - Días 2 y 4 - TOP 10 picos - Algoritmo: kmeans</em>
   </p>

   <p align="center">
   <img src='Imagenes/4_pam_2_clusters_4_grupos.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 14: Clustering - CLP vs SHAM - Días 2 y 4 - TOP 20 picos - Algoritmo: PAM</em>
   </p>

   <p align="center">
   <img src='Imagenes/4_pam_3_clusters_4_grupos.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 15: Clustering - CLP vs SHAM - Días 2 y 4 - 3 clusters - TOP 20 picos - Algoritmo: PAM</em>
   </p>

   <p align="center">
   <img src='Imagenes/4_pam_3_grupos_3_clusters.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 16: Clustering - CLP vs SHAM - Días 2 y 4 - 3 clusters - TOP 30 picos - Algoritmo: PAM</em>
   </p>

   <p align="center">
   <img src='Imagenes/4_pam_top10.jpg' width='400'>
   </p>
   <p align="center">
     <em>Figura 17: Clustering - CLP vs SHAM - Días 2 y 4 - TOP 10 picos - Algoritmo: PAM</em>
   </p>
   

   ### Archivo: 5_ns_m51_d24

   
   ### Archivo: 6_ns_m51_vs_varios

   
## Algoritmos Supervisados

## Conclusión
